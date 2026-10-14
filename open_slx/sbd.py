import struct
import time
from collections.abc import Iterator
from pathlib import Path

import h5py
import numpy as np
import zstandard

from open_slx.errors import UnsupportedFormat

TRAILER = 17
SUPPORTED_VERSIONS = {36}


def is_url(source: str) -> bool:
    return str(source).startswith(("http://", "https://"))


def decode_mz(blob: bytes, n_peaks: int) -> np.ndarray:
    d = zstandard.ZstdDecompressor().decompressobj()
    planes = d.decompress(blob)
    trailer = d.unused_data
    if len(trailer) != TRAILER:
        raise ValueError(f"m/z trailer has {len(trailer)} bytes, expected {TRAILER}")
    k, scale, mz0 = struct.unpack("<Bdd", trailer)
    if len(planes) != k * n_peaks:
        raise ValueError(f"m/z blob holds {len(planes)} bytes, expected {k * n_peaks}")
    arr = np.frombuffer(planes, np.uint8).reshape(k, n_peaks)
    gaps = np.zeros(n_peaks, np.uint64)
    for p in range(k):
        gaps |= arr[p].astype(np.uint64) << np.uint64(8 * p)
    return mz0 + np.cumsum(gaps, dtype=np.float64) / scale


def encode_mz(mz: np.ndarray, scale: float = 500000.0) -> bytes:
    mz = np.asarray(mz, np.float64)
    q = np.round((mz - mz[0]) * scale).astype(np.int64)
    gaps = np.diff(q, prepend=q[0]).astype(np.uint64)
    k = max(1, (int(gaps.max()).bit_length() + 7) // 8)
    planes = np.stack([((gaps >> np.uint64(8 * p)) & np.uint64(255)).astype(np.uint8) for p in range(k)])
    return zstandard.ZstdCompressor().compress(planes.tobytes()) + struct.pack("<Bdd", k, scale, mz[0])


class _RetryFile:
    def __init__(self, raw, tries: int = 6):
        self.raw = raw
        self.tries = tries
        self._pos = 0

    def _retry(self, name: str, *args):
        for attempt in range(self.tries):
            try:
                return getattr(self.raw, name)(*args)
            except Exception:
                if attempt == self.tries - 1:
                    raise
                time.sleep(2**attempt)
                self.raw.seek(self._pos)

    def seek(self, offset, whence=0):
        out = self.raw.seek(offset, whence)
        self._pos = self.raw.tell()
        return out

    def tell(self):
        return self.raw.tell()

    def read(self, n=-1):
        self._pos = self.raw.tell()
        return self._retry("read", n)

    def readinto(self, buf):
        self._pos = self.raw.tell()
        return self._retry("readinto", buf)

    def seekable(self):
        return True

    def readable(self):
        return True

    def close(self):
        self.raw.close()


class Sbd:
    def __init__(self, source: str | Path):
        self.source = str(source)
        self._handle = None
        if is_url(self.source):
            import fsspec

            self._handle = _RetryFile(fsspec.open(self.source, block_size=2**20, cache_type="blockcache").open())
            self.file = h5py.File(self._handle, "r")
        else:
            self.file = h5py.File(self.source, "r")
        try:
            self._check_version()
        except UnsupportedFormat:
            self.close()
            raise
        sg = self.file["SpectraGroup"]
        self._mz_blob = sg["compressedMz"]
        self._intensity = sg["intensities"]
        self._mz_offsets = sg["mzOffsets"]
        self._int_offsets = sg["intensityOffsets"]
        reg = self.file["Registrations/0"]
        self._world = reg["Coordinates"]
        self._original = reg["OriginalCoordinates"]
        self._original_index = reg["OriginalSpectrumIndices"]
        self._raster = self.file["MeasurementRasterIndices"]
        self._rasters = reg["Rasters"]
        self.n_spots = int(self._mz_offsets.shape[0])
        self._transforms: dict[int, np.ndarray] = {}
        self._axis = None

    def _check_version(self) -> None:
        version = int(np.asarray(self.file["Version"][()]).ravel()[0]) if "Version" in self.file else None
        if version not in SUPPORTED_VERSIONS:
            raise UnsupportedFormat(
                f".sbd Version {version} is not supported (supported: {', '.join(map(str, sorted(SUPPORTED_VERSIONS)))})"
            )
        spectra_type = self.attrs.get("SpectraType")
        if spectra_type != "Peaks":
            raise UnsupportedFormat(f".sbd SpectraType {spectra_type} is not supported (supported: Peaks)")

    def close(self) -> None:
        self.file.close()
        if self._handle is not None:
            self._handle.close()

    def __enter__(self) -> "Sbd":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    @property
    def attrs(self) -> dict[str, str]:
        out = {}
        for k, v in self.file.attrs.items():
            v = v[0] if getattr(v, "ndim", 0) else v
            out[k] = v.decode() if isinstance(v, bytes) else str(v)
        return out

    @property
    def mz_axis(self) -> np.ndarray:
        if self._axis is None:
            self._axis = self.file["SamplePositions/GlobalMassAxis/SamplePositions"][:]
        return self._axis

    def transform(self, raster: int) -> np.ndarray:
        if raster not in self._transforms:
            self._transforms[raster] = self._rasters[str(raster)]["PixelToWorld"][:]
        return self._transforms[raster]

    def pixel_size(self, raster: int) -> tuple[float, float]:
        m = self.transform(raster)
        return abs(float(m[0, 0])), abs(float(m[1, 1]))

    def raster_indices(self, start: int = 0, stop: int | None = None) -> np.ndarray:
        return self._raster[0, start:stop].astype(np.int64)

    def world_coordinates(self, start: int = 0, stop: int | None = None) -> np.ndarray:
        return self._world[:2, start:stop].astype(np.float64)

    def original_spectrum_indices(self, start: int = 0, stop: int | None = None) -> np.ndarray:
        return self._original_index[start:stop].astype(np.int64)

    def pixel_coordinates(self, start: int = 0, stop: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raster = self.raster_indices(start, stop)
        world = self.world_coordinates(start, stop)
        x = np.zeros(len(raster), np.int64)
        y = np.zeros(len(raster), np.int64)
        for r in np.unique(raster):
            m = self.transform(int(r))
            sel = raster == r
            x[sel] = np.round((world[0, sel] - m[0, 3]) / m[0, 0]).astype(np.int64) + 1
            y[sel] = np.round((world[1, sel] - m[1, 3]) / m[1, 1]).astype(np.int64) + 1
        return raster, x, y

    def _bounds(self, offsets, total: int, start: int, stop: int) -> np.ndarray:
        got = offsets[start : min(stop + 1, self.n_spots)].astype(np.int64)
        if stop >= self.n_spots:
            got = np.append(got, total)
        return got

    def spectra(self, start: int = 0, stop: int | None = None, batch: int = 256) -> Iterator[tuple[int, int, np.ndarray, np.ndarray]]:
        stop = self.n_spots if stop is None else min(stop, self.n_spots)
        for a in range(start, stop, batch):
            b = min(a + batch, stop)
            _, xs, ys = self.pixel_coordinates(a, b)
            mzb = self._bounds(self._mz_offsets, self._mz_blob.shape[0], a, b)
            ib = self._bounds(self._int_offsets, self._intensity.shape[0], a, b)
            blob = self._mz_blob[mzb[0] : mzb[-1]].tobytes()
            inten = self._intensity[ib[0] : ib[-1]]
            for i in range(b - a):
                n = int(ib[i + 1] - ib[i])
                mz = decode_mz(blob[mzb[i] - mzb[0] : mzb[i + 1] - mzb[0]], n)
                yield int(xs[i]), int(ys[i]), mz, inten[ib[i] - ib[0] : ib[i + 1] - ib[0]]
