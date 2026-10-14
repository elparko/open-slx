import hashlib
import uuid

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser

from open_slx.imzml import ImzMLWriter


def _spectra():
    rng = np.random.default_rng(0)
    out = []
    for i, (x, y) in enumerate([(1, 1), (2, 1), (1, 2), (3, 2)]):
        n = 10 + i
        out.append((x, y, np.sort(rng.uniform(300, 900, n)), rng.uniform(1, 1e5, n).astype(np.float32)))
    return out


def test_processed_round_trip(tmp_path):
    data = _spectra()
    with ImzMLWriter(tmp_path / "a", "processed", pixel_size=(25.0, 25.0)) as w:
        for s in data:
            w.add_spectrum(*s)
    with ImzMLParser(str(tmp_path / "a.imzML")) as p:
        assert len(p.coordinates) == 4
        assert p.mzPrecision == "d"
        for i, (x, y, mz, inten) in enumerate(data):
            assert p.coordinates[i][:2] == (x, y)
            pm, pi = p.getspectrum(i)
            assert np.allclose(pm, mz)
            assert np.array_equal(pi, inten)


def test_continuous_round_trip(tmp_path):
    axis = np.linspace(300, 900, 50)
    rng = np.random.default_rng(0)
    rows = [(1, 1, rng.uniform(0, 9, 50)), (2, 1, rng.uniform(0, 9, 50))]
    with ImzMLWriter(tmp_path / "c", "continuous", mz=axis, pixel_size=(50.0, 50.0)) as w:
        for x, y, v in rows:
            w.add_spectrum(x, y, None, v)
    assert (tmp_path / "c.ibd").stat().st_size == 16 + 50 * 8 + 2 * 50 * 4
    with ImzMLParser(str(tmp_path / "c.imzML")) as p:
        for i, (_, _, v) in enumerate(rows):
            pm, pi = p.getspectrum(i)
            assert np.allclose(pm, axis)
            assert np.allclose(pi, v.astype(np.float32))


def test_uuid_and_sha1_and_size(tmp_path):
    with ImzMLWriter(tmp_path / "u", "processed", pixel_size=(10.0, 20.0)) as w:
        for s in _spectra():
            w.add_spectrum(*s)
        uid = w.uuid
    ibd = (tmp_path / "u.ibd").read_bytes()
    assert ibd[:16] == uid.bytes
    xml = (tmp_path / "u.imzML").read_text()
    assert "{" + str(uid).upper() + "}" in xml
    assert hashlib.sha1(ibd).hexdigest().upper() in xml
    assert uuid.UUID(bytes=ibd[:16]) == uid
    with ImzMLParser(str(tmp_path / "u.imzML")) as p:
        assert p.imzmldict["max count of pixels x"] == 3
        assert p.imzmldict["max count of pixels y"] == 2
    assert 'name="pixel size (x)" value="10"' in xml
    assert 'name="pixel size y" value="20"' in xml


def test_continuous_needs_axis(tmp_path):
    with pytest.raises(ValueError):
        ImzMLWriter(tmp_path / "x", "continuous")


def test_dotted_names_and_two_rasters_do_not_collide(tmp_path):
    for suffix in ("_r00", "_r01"):
        with ImzMLWriter(tmp_path / f"2023.05.12_kidney{suffix}", "processed") as w:
            w.add_spectrum(*_spectra()[0])
    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == [
        "2023.05.12_kidney_r00.ibd",
        "2023.05.12_kidney_r00.imzML",
        "2023.05.12_kidney_r01.ibd",
        "2023.05.12_kidney_r01.imzML",
    ]


def test_streamed_xml_counts_and_no_temp_file(tmp_path):
    with ImzMLWriter(tmp_path / "s", "processed", pixel_size=(5.0, 5.0)) as w:
        for s in _spectra():
            w.add_spectrum(*s)
        assert w.n_spectra == 4
    assert sorted(p.name for p in tmp_path.iterdir()) == ["s.ibd", "s.imzML"]
    assert '<spectrumList count="4"' in (tmp_path / "s.imzML").read_text()


def test_discard_removes_all_files(tmp_path):
    w = ImzMLWriter(tmp_path / "d", "processed")
    w.add_spectrum(*_spectra()[0])
    w.discard()
    assert list(tmp_path.iterdir()) == []


def test_schema_location_matches_root_element(tmp_path):
    with ImzMLWriter(tmp_path / "x", "processed") as w:
        w.add_spectrum(*_spectra()[0])
    xml = (tmp_path / "x.imzML").read_text()
    assert "<mzML " in xml
    assert "mzML1.1.0.xsd" in xml and "_idx" not in xml


def test_failed_second_open_removes_the_first_file(tmp_path):
    (tmp_path / "o.imzML.body.tmp").mkdir()
    with pytest.raises(OSError):
        ImzMLWriter(tmp_path / "o")
    assert not (tmp_path / "o.ibd").exists()
