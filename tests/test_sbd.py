import numpy as np
import pytest

from open_slx.errors import UnsupportedFormat
from open_slx.sbd import Sbd, decode_mz, encode_mz
from tests.conftest import PIXELS, RASTERS


def test_mz_codec_round_trip():
    mz = np.array([249.65, 250.1, 301.123456, 1099.9])
    out = decode_mz(encode_mz(mz), len(mz))
    assert np.allclose(out, mz, atol=1e-6)


def test_mz_codec_uses_four_planes_for_large_gap():
    blob = encode_mz(np.array([250.0, 1000.0]))
    assert blob[-17] == 4


def test_spectra_match_input(synthetic):
    path, spectra = synthetic
    with Sbd(path) as sbd:
        got = list(sbd.spectra(batch=3))
    assert len(got) == len(spectra)
    for (x, y, mz, inten), (mz0, in0) in zip(got, spectra):
        assert np.allclose(mz, mz0, atol=1e-6)
        assert np.array_equal(inten, in0)


def test_pixel_coordinates_and_sizes(synthetic):
    path, _ = synthetic
    with Sbd(path) as sbd:
        raster, x, y = sbd.pixel_coordinates()
        assert np.array_equal(raster, RASTERS)
        assert list(zip(x - 1, y - 1)) == PIXELS
        assert sbd.pixel_size(0) == (50.0, 50.0)
        assert sbd.pixel_size(1) == (25.0, 25.0)
        assert sbd.n_spots == 7
        assert sbd.mz_axis.shape == (1000,)


def test_spectra_slice(synthetic):
    path, spectra = synthetic
    with Sbd(path) as sbd:
        got = list(sbd.spectra(2, 5))
    assert len(got) == 3
    assert np.array_equal(got[0][3], spectra[2][1])


def test_unknown_version_is_rejected(make_sbd):
    with pytest.raises(UnsupportedFormat, match="Version 37"):
        Sbd(make_sbd(version=37))


def test_profile_spectra_are_rejected(make_sbd):
    with pytest.raises(UnsupportedFormat, match="SpectraType Profile"):
        Sbd(make_sbd(spectra_type=b"Profile"))
