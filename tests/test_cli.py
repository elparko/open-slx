import itertools
import json
import os
import shutil
import socket
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser

from open_slx.cli import main
from open_slx.sbd import Sbd

PRIVATE = Path(os.environ.get("SLX_TEST_FILE", ""))
ZENODO = "https://zenodo.org/records/21251439/files/Romics_paper_SCils_file_case_2.sbd?download=1"


def _online() -> bool:
    try:
        socket.create_connection(("zenodo.org", 443), timeout=3).close()
        return True
    except OSError:
        return False


def test_sbd_only_conversion_splits_rasters(synthetic, tmp_path):
    sbd, spectra = synthetic
    out = tmp_path / "out" / "demo"
    assert main(["--sbd", str(sbd), "-o", str(out)]) == 0
    with ImzMLParser(str(tmp_path / "out" / "demo_r00.imzML")) as p:
        assert len(p.coordinates) == 3
        assert p.coordinates[2][:2] == (3, 2)
        mz, inten = p.getspectrum(1)
        assert np.allclose(mz, spectra[1][0], atol=1e-6)
    with ImzMLParser(str(tmp_path / "out" / "demo_r01.imzML")) as p:
        assert len(p.coordinates) == 4
        assert p.imzmldict["pixel size x"] == 25.0
    assert json.loads((tmp_path / "out" / "demo_metadata.json").read_text())["n_spots"] == 7


def test_max_spots_and_continuous(synthetic, tmp_path):
    sbd, _ = synthetic
    out = tmp_path / "c"
    assert main(["--sbd", str(sbd), "-o", str(out), "--mode", "continuous", "--max-spots", "2"]) == 0
    with ImzMLParser(str(tmp_path / "c.imzML")) as p:
        assert len(p.coordinates) == 2
        mz, inten = p.getspectrum(0)
        assert len(mz) == 1000
        assert inten.sum() > 0


@pytest.mark.skipif(
    not os.environ.get("SLX_TEST_FILE") or not PRIVATE.is_file(), reason="SLX_TEST_FILE unset or missing"
)
def test_private_project_without_sbd_exits_2_and_writes_slx_outputs(tmp_path, capsys):
    slx = tmp_path / "project.slx"
    shutil.copy(PRIVATE, slx)
    out = tmp_path / "res"
    assert main([str(slx), "-o", str(out)]) == 2
    assert "project.sbd" in capsys.readouterr().err
    meta = json.loads((tmp_path / "res_metadata.json").read_text())
    assert meta["n_spots"] == 94713
    assert meta["mz_axis"]["points"] == 460519
    rows = (tmp_path / "res_features.csv").read_text().splitlines()
    assert len(rows) == 323


def test_missing_sbd_exits_2_with_synthetic_slx(make_slx, tmp_path, capsys):
    slx = make_slx()
    out = tmp_path / "res"
    assert main([str(slx), "-o", str(out)]) == 2
    assert "demo.sbd" in capsys.readouterr().err
    assert json.loads((tmp_path / "res_metadata.json").read_text())["n_spots"] == 7
    assert len((tmp_path / "res_features.csv").read_text().splitlines()) == 3
    assert not list(tmp_path.glob("*.imzML"))


def test_spot_count_mismatch_exits_3(make_slx, synthetic, tmp_path, capsys):
    sbd, _ = synthetic
    slx = make_slx("synthetic.slx", n_spots=9)
    assert main([str(slx), "-o", str(tmp_path / "res")]) == 3
    assert ".sbd has 7 spots, .slx has 9" in capsys.readouterr().err
    assert not list(tmp_path.glob("*.imzML"))


def test_slx_and_sbd_together_convert(make_slx, synthetic, tmp_path):
    sbd, _ = synthetic
    slx = make_slx("synthetic.slx")
    assert main([str(slx), "-o", str(tmp_path / "res")]) == 0
    assert (tmp_path / "res_r00.imzML").exists() and (tmp_path / "res_r01.imzML").exists()


def test_unsupported_sbd_version_exits_4(make_sbd, tmp_path, capsys):
    sbd = make_sbd(version=99)
    assert main(["--sbd", str(sbd), "-o", str(tmp_path / "o")]) == 4
    assert "Version 99" in capsys.readouterr().err
    assert not list(tmp_path.glob("o*"))


def test_unsupported_slx_version_exits_4(make_slx, tmp_path, capsys):
    assert main([str(make_slx(version=13)), "-o", str(tmp_path / "o")]) == 4
    assert "schema version 13" in capsys.readouterr().err


def test_dotted_output_name_with_two_rasters(synthetic, tmp_path):
    sbd, _ = synthetic
    assert main(["--sbd", str(sbd), "-o", str(tmp_path / "2023.05.12_kidney")]) == 0
    assert sorted(p.name for p in tmp_path.glob("2023*")) == [
        "2023.05.12_kidney_metadata.json",
        "2023.05.12_kidney_r00.ibd",
        "2023.05.12_kidney_r00.imzML",
        "2023.05.12_kidney_r01.ibd",
        "2023.05.12_kidney_r01.imzML",
    ]


def test_typed_imzml_suffix_is_stripped(synthetic, tmp_path):
    sbd, _ = synthetic
    assert main(["--sbd", str(sbd), "-o", str(tmp_path / "typed.imzML"), "--max-spots", "3"]) == 0
    assert (tmp_path / "typed.imzML").exists()
    assert (tmp_path / "typed.ibd").exists()


def test_summary_line_per_output(synthetic, tmp_path, capsys):
    sbd, _ = synthetic
    out = tmp_path / "demo"
    assert main(["--sbd", str(sbd), "-o", str(out)]) == 0
    assert capsys.readouterr().out.splitlines() == [
        f"raster 0: 3 spots, pixel size 50 x 50 um, {out}_r00.imzML",
        f"raster 1: 4 spots, pixel size 25 x 25 um, {out}_r01.imzML",
    ]


def test_error_during_conversion_deletes_partial_outputs(synthetic, tmp_path, monkeypatch, capsys):
    sbd, _ = synthetic

    def broken(self, start=0, stop=None, batch=256):
        yield from itertools.islice(real(self, start, stop, batch), 2)
        raise RuntimeError("boom")

    real = Sbd.spectra
    monkeypatch.setattr(Sbd, "spectra", broken)
    assert main(["--sbd", str(sbd), "-o", str(tmp_path / "o")]) == 1
    assert "boom" in capsys.readouterr().err
    assert not list(tmp_path.glob("o*.imzML*")) and not list(tmp_path.glob("o*.ibd")) and not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("value", ["0", "-1", "x"])
def test_max_spots_must_be_positive(synthetic, tmp_path, value, capsys):
    sbd, _ = synthetic
    with pytest.raises(SystemExit) as e:
        main(["--sbd", str(sbd), "-o", str(tmp_path / "o"), "--max-spots", value])
    assert e.value.code == 2


def test_missing_slx_without_sbd_exits_2(tmp_path, capsys):
    assert main([str(tmp_path / "none.slx")]) == 2


@pytest.mark.network
@pytest.mark.skipif(os.environ.get("SLX_NETWORK") != "1" or not _online(), reason="needs SLX_NETWORK=1 and network")
def test_zenodo_first_200_spots(tmp_path):
    out = tmp_path / "zen"
    assert main(["--sbd", ZENODO, "-o", str(out), "--max-spots", "200"]) == 0
    with ImzMLParser(str(tmp_path / "zen.imzML")) as p:
        assert len(p.coordinates) == 200
        xs = [c[0] for c in p.coordinates]
        ys = [c[1] for c in p.coordinates]
        assert min(xs) >= 1 and max(xs) <= 60 and min(ys) >= 1 and max(ys) <= 36
        lo, hi = np.inf, -np.inf
        for i in range(200):
            mz, inten = p.getspectrum(i)
            lo, hi = min(lo, mz.min()), max(hi, mz.max())
            assert (inten > 0).any()
        assert 245 <= lo and hi <= 1101
        tic = np.array([p.getspectrum(i)[1].sum() for i in range(200)], np.float64)
    with Sbd(ZENODO) as sbd:
        for g in sbd.file["SpectraGroup/Normalizations"].values():
            if g.attrs["name"] == b"Total Ion Count":
                expected = g["Values"][:200] * g["TargetValue"][0]
    assert np.allclose(tic, expected, rtol=1e-4)


def test_no_arguments_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as e:
        main([])
    assert e.value.code == 2
    assert "give a .slx file, --sbd, or both" in capsys.readouterr().err


def test_unreadable_slx_prints_error_and_exits_1(tmp_path, capsys):
    bad = tmp_path / "bad.slx"
    bad.write_text("not a database")
    assert main([str(bad), "-o", str(tmp_path / "o")]) == 1
    assert capsys.readouterr().err.startswith("error: ")


def test_failed_csv_write_prints_error_and_exits_1(make_slx, tmp_path, capsys):
    slx = make_slx()
    (tmp_path / "res_features.csv").mkdir()
    assert main([str(slx), "-o", str(tmp_path / "res")]) == 1
    assert capsys.readouterr().err.startswith("error: ")


def test_failed_json_write_prints_error_and_exits_1(make_slx, tmp_path, capsys):
    slx = make_slx()
    (tmp_path / "res_metadata.json").mkdir()
    assert main([str(slx), "-o", str(tmp_path / "res")]) == 1
    assert capsys.readouterr().err.startswith("error: ")
