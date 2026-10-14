import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from open_slx.errors import UnsupportedFormat
from open_slx.imzml import ImzMLWriter
from open_slx.sbd import Sbd, is_url
from open_slx.slx import Slx, read_slx


def _find_sbd(slx_path: Path, slx: Slx | None) -> Path:
    same_stem = slx_path.with_suffix(".sbd")
    if same_stem.exists() or slx is None:
        return same_stem
    named = slx.project.get("storage", {}).get("file")
    if named:
        candidate = slx_path.parent / Path(str(named).replace("\\", "/")).name
        if candidate.exists():
            return candidate
    return same_stem


def _write_features(slx: Slx, path: Path) -> None:
    f = slx.features
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        tags = sorted(f.tags)
        w.writerow(["row", "id", "mz_low", "mz_high", "mz_center", "mz_width", "intensity", "color", "name", *tags])
        for i in range(len(f)):
            w.writerow(
                [
                    i,
                    int(f.ids[i]),
                    f"{f.mz_low[i]:.8f}",
                    f"{f.mz_high[i]:.8f}",
                    f"{f.mz_center[i]:.8f}",
                    f"{f.mz_width[i]:.8f}",
                    "" if f.intensity is None else f"{f.intensity[i]:.6g}",
                    f.colors[i] if f.colors else "",
                    f.names[i] if f.names else "",
                    *[int(f.tags[t][i]) for t in tags],
                ]
            )


def _slx_metadata(slx: Slx) -> dict:
    mz = slx.mz_axis
    return {
        "slx": slx.path,
        "schema_version": slx.schema_version,
        "n_spots": slx.n_spots,
        "mz_axis": {"points": len(mz), "min": float(mz[0]), "max": float(mz[-1])},
        "data_source": {
            "type": slx.data_source.type,
            "polarity": slx.data_source.polarity,
            "source_path": slx.data_source.source_path,
            "spot_name_format": slx.data_source.spot_name_format,
        },
        "rasters": [
            {"index": r.scils_index, "original_index": r.original_index, "x_offset": r.x_offset, "y_offset": r.y_offset}
            for r in slx.rasters
        ],
        "storage": slx.project.get("storage", {}),
        "mass_axis": slx.project.get("massaxis", {}),
        "n_features": len(slx.features),
        "import_report": slx.import_report,
    }


def _convert(sbd: Sbd, out: Path, mode: str, max_spots: int | None) -> list[tuple[int, int, tuple[float, float], Path]]:
    stop = sbd.n_spots if max_spots is None else min(max_spots, sbd.n_spots)
    raster_of = sbd.raster_indices(0, stop)
    rasters = np.unique(raster_of)
    axis = sbd.mz_axis if mode == "continuous" else None
    writers: dict[int, ImzMLWriter] = {}
    try:
        for r in rasters:
            stem = out if len(rasters) == 1 else out.with_name(f"{out.name}_r{int(r):02d}")
            writers[int(r)] = ImzMLWriter(stem, mode=mode, mz=axis, pixel_size=sbd.pixel_size(int(r)))
        for i, (x, y, mz, inten) in enumerate(sbd.spectra(0, stop)):
            w = writers[int(raster_of[i])]
            if mode == "continuous":
                dense = np.zeros(len(axis), np.float32)
                idx = np.clip(np.searchsorted(axis, mz), 0, len(axis) - 1)
                left = np.clip(idx - 1, 0, len(axis) - 1)
                idx = np.where(np.abs(axis[left] - mz) < np.abs(axis[idx] - mz), left, idx)
                np.add.at(dense, idx, inten)
                w.add_spectrum(x, y, None, dense)
            else:
                w.add_spectrum(x, y, mz, inten)
        for w in writers.values():
            w.close()
    except BaseException:
        for w in writers.values():
            w.discard()
        raise
    return [(r, w.n_spectra, w.pixel_size, w.imzml_path) for r, w in writers.items()]


def _positive_int(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid integer: {text!r}")
    if value < 1:
        raise argparse.ArgumentTypeError("must be 1 or more")
    return value


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="open-slx")
    p.add_argument("slx", nargs="?", help="SCiLS .slx project (path)")
    p.add_argument("--sbd", help="path or http(s) URL of the .sbd file")
    p.add_argument("-o", "--output", help="output path without extension")
    p.add_argument("--mode", choices=["continuous", "processed"], default="processed")
    p.add_argument("--max-spots", type=_positive_int)
    args = p.parse_args(argv)
    if not args.slx and not args.sbd:
        p.error("give a .slx file, --sbd, or both")
    try:
        return _run(args)
    except UnsupportedFormat as e:
        print(f"error: {e}", file=sys.stderr)
        return 4
    except Exception as e:
        print(f"error: {type(e).__name__}: {e}", file=sys.stderr)
        return 1


def _run(args: argparse.Namespace) -> int:
    slx_path = Path(args.slx) if args.slx and not is_url(args.slx) else None
    slx = None
    if slx_path is not None and slx_path.exists():
        slx = read_slx(slx_path)
    elif not args.sbd:
        print(f"error: .slx file not found: {args.slx} (pass --sbd to convert without it)", file=sys.stderr)
        return 2

    if args.sbd:
        sbd_source = args.sbd
    else:
        sbd_source = str(_find_sbd(slx_path, slx))

    if args.output:
        out = Path(args.output)
        if out.suffix.lower() in (".imzml", ".ibd"):
            out = out.with_name(out.stem)
    else:
        base = slx_path if slx_path is not None else Path(sbd_source.split("?")[0])
        out = Path(base.stem)

    if slx is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        _write_features(slx, out.with_name(out.name + "_features.csv"))
        out.with_name(out.name + "_metadata.json").write_text(json.dumps(_slx_metadata(slx), indent=2), encoding="utf-8")

    if not is_url(sbd_source) and not Path(sbd_source).exists():
        print(f"error: .sbd file not found: {sbd_source}", file=sys.stderr)
        return 2

    out.parent.mkdir(parents=True, exist_ok=True)
    with Sbd(sbd_source) as sbd:
        if slx is not None and sbd.n_spots != slx.n_spots:
            print(f"error: .sbd has {sbd.n_spots} spots, .slx has {slx.n_spots}", file=sys.stderr)
            return 3
        if slx is None:
            out.with_name(out.name + "_metadata.json").write_text(
                json.dumps({"sbd": sbd.source, "n_spots": sbd.n_spots, "attributes": sbd.attrs}, indent=2), encoding="utf-8"
            )
        for r, n, (px, py), path in _convert(sbd, out, args.mode, args.max_spots):
            print(f"raster {r}: {n} spots, pixel size {px:g} x {py:g} um, {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
