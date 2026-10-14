import sqlite3

import h5py
import numpy as np
import pytest

from open_slx.sbd import encode_mz

RASTERS = np.array([0, 0, 0, 1, 1, 1, 1])
PIXELS = [(0, 0), (1, 0), (2, 1), (3, 2), (4, 2), (3, 3), (4, 3)]
TRANSFORMS = {0: (50.0, 100.0, 1000.0), 1: (25.0, 5000.0, 2000.0)}


@pytest.fixture
def synthetic(tmp_path):
    return _build_sbd(tmp_path / "synthetic.sbd")


def _build_sbd(path, version=36, spectra_type=b"Peaks"):
    rng = np.random.default_rng(1)
    spectra = []
    for _ in RASTERS:
        n = int(rng.integers(5, 40))
        mz = np.sort(rng.uniform(250, 1100, n))
        spectra.append((mz, rng.uniform(1e3, 1e6, n).astype(np.float32)))
    blobs = [encode_mz(mz) for mz, _ in spectra]
    axis = 245.69 * (1100.01 / 245.69) ** (np.arange(1000) / 999)
    world = np.zeros((3, len(RASTERS)), np.float32)
    for i, (r, (px, py)) in enumerate(zip(RASTERS, PIXELS)):
        s, tx, ty = TRANSFORMS[int(r)]
        world[0, i] = tx + s * px
        world[1, i] = ty - s * py
    with h5py.File(path, "w") as f:
        f["Version"] = np.array([version], np.uint64)
        f.attrs["nSpots"] = np.array([str(len(RASTERS)).encode()])
        f.attrs["AxisType"] = np.array([b"MRMS"])
        f.attrs["SpectraType"] = np.array([spectra_type])
        sg = f.create_group("SpectraGroup")
        sg["compressedMz"] = np.frombuffer(b"".join(blobs), np.int8)
        sg["mzOffsets"] = np.cumsum([0] + [len(b) for b in blobs[:-1]]).astype(np.uint64)
        sg["intensities"] = np.concatenate([i for _, i in spectra])
        sg["intensityOffsets"] = np.cumsum([0] + [len(m) for m, _ in spectra][:-1]).astype(np.uint64)
        f["SamplePositions/GlobalMassAxis/SamplePositions"] = axis
        f["MeasurementRasterIndices"] = RASTERS.astype(np.int32)[None, :]
        reg = f.create_group("Registrations/0")
        reg["Coordinates"] = world
        reg["OriginalCoordinates"] = np.zeros((2, len(RASTERS)), np.uint64)
        reg["OriginalSpectrumIndices"] = np.arange(len(RASTERS), dtype=np.uint64)
        for r, (s, tx, ty) in TRANSFORMS.items():
            m = np.diag([s, -s, 1.0, 1.0])
            m[0, 3], m[1, 3] = tx, ty
            reg[f"Rasters/{r}/PixelToWorld"] = m
    return path, spectra


@pytest.fixture
def make_sbd(tmp_path):
    def build(name="bad.sbd", **kwargs):
        path = tmp_path / name
        _build_sbd(path, **kwargs)
        return path

    return build


def _meta(**pairs):
    return b"".join(f"{len(k)}:{k}{len(v)}:{v}".encode() for k, v in pairs.items())


@pytest.fixture
def make_slx(tmp_path):
    def build(name="demo.slx", n_spots=len(RASTERS), version=12):
        path = tmp_path / name
        con = sqlite3.connect(path)
        con.executescript(
            """
            create table Master (KeyColumn text, ValueColumn text);
            create table MetaInformation (ID text, MetaInformationString blob);
            create table ObjectRelations (ID text, type text, parent text, own int, relations text);
            create table DataSources (ID text, SourcePath text, OriginalImportPath text, Type text, Polarity text, SpotNameFormat text);
            create table Rasters (ID text, ScilsRasterIdx int, OriginalRasterIdx int, RasterXOffset int, RasterYOffset int);
            create table Labels (SpotList blob, MergeList text);
            create table SamplePositions (ID text, MzValues blob);
            create table FeatureTable (ID text, ColumnType int, ColumnData blob, ColumnIndex int);
            """
        )
        con.execute("insert into Master values ('Version', ?)", (str(version),))
        con.executemany(
            "insert into MetaInformation values (?, ?)",
            [
                ("file", _meta(nSpots=str(n_spots), Comment="synthetic")),
                ("c0", _meta(ContentId="Basic:Id")),
                ("c1", _meta(ContentId="Basic:Mz")),
            ],
        )
        con.execute("insert into DataSources values ('ds', 'src', 'src', 'TSFSource', 'negative', 'fmt')")
        con.executemany("insert into Rasters values (?, ?, ?, ?, ?)", [("r0", 0, 0, 0, 0), ("r1", 1, 1, 0, 0)])
        con.execute("insert into Labels values (?, null)", (np.arange(n_spots, dtype="<u8").tobytes(),))
        con.execute("insert into SamplePositions values ('sp', ?)", (np.linspace(250, 1100, 10).astype("<f8").tobytes(),))
        con.execute("insert into FeatureTable values ('c0', 3, ?, 0)", (np.array([0, 1], "<i8").tobytes(),))
        con.execute("insert into FeatureTable values ('c1', 1, ?, 1)", (np.array([300, 301, 500, 501], "<f8").tobytes(),))
        con.commit()
        con.close()
        return path

    return build
