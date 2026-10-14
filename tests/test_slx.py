import os
from pathlib import Path

import numpy as np
import pytest

from open_slx import parse_meta, read_slx
from open_slx.errors import UnsupportedFormat

SLX_PATH = Path(os.environ.get("SLX_TEST_FILE", ""))

needs_file = pytest.mark.skipif(
    not os.environ.get("SLX_TEST_FILE") or not SLX_PATH.is_file(), reason="SLX_TEST_FILE unset or missing"
)


@pytest.fixture(scope="module")
def slx():
    return read_slx(SLX_PATH)


def test_parse_meta_counts_utf8_bytes():
    raw = "4:name15:1/K0 width [±]4:type13:MzTableColumn".encode()
    assert parse_meta(raw) == {"name": "1/K0 width [±]", "type": "MzTableColumn"}


@needs_file
def test_schema_and_data_source(slx):
    assert slx.schema_version == 12
    assert slx.data_source.type == "TSFSource"
    assert slx.data_source.polarity == "negative"
    assert slx.data_source.spot_name_format == "0_R%02dX%03dY%03d"
    assert slx.data_source.meta["Source/Number of spectra in dataset"] == "94738"


@needs_file
def test_rasters(slx):
    assert [(r.scils_index, r.x_offset, r.y_offset) for r in slx.rasters] == [
        (0, 1131, 1243),
        (1, 2504, 1417),
        (2, 4093, 2059),
        (3, 6124, 1544),
    ]


@needs_file
def test_spots(slx):
    assert slx.n_spots == 94713
    assert slx.spot_ids.dtype == np.uint64
    assert np.array_equal(slx.spot_ids, np.arange(94713))


@needs_file
def test_spot_tree_spans_all_spots(slx):
    tree = slx.spot_tree
    assert tree.edges.shape == (94712, 2)
    assert tree.edges.max() == 94712
    assert tree.distances.min() > 0
    assert tree.distances.max() == 2.0
    assert len(np.unique(tree.edges)) == 94713


@needs_file
def test_mz_axis(slx):
    axis = slx.mz_axis
    assert axis.shape == (460519,)
    assert axis[0] == 200.0
    assert axis[-1] == pytest.approx(2000.0)
    expected = 200.0 * 10.0 ** (np.arange(axis.size) / (axis.size - 1))
    assert np.allclose(axis, expected, rtol=1e-12)
    assert len(slx.sample_positions) == 2
    first, second = slx.sample_positions.values()
    assert np.array_equal(first, second)


@needs_file
def test_features(slx):
    f = slx.features
    assert len(f) == 322
    assert (f.ids >= 0).sum() == 36
    assert np.array_equal(f.ids[:36], np.arange(36))
    assert (f.ids[36:] == -1).all()
    assert f.mz_low[0] == pytest.approx(255.23107067)
    assert f.mz_high[0] == pytest.approx(255.23492596)
    assert f.mz_center[0] == pytest.approx(255.23299832, abs=1e-8)
    assert f.mz_width[0] == pytest.approx(0.00385529, abs=1e-8)
    assert (f.mz_width > 0).all()
    assert f.mz_low.min() > 200 and f.mz_high.max() < 2000
    assert f.intensity.shape == (322,)
    assert len(f.colors) == 322 and f.colors[0] == "#ff0000"
    assert f.names is None
    assert f.tags["manual"].sum() == 32
    assert f.tags["test"].sum() == 4
    assert f.tags["test"][32:36].all()


@needs_file
def test_region_tree(slx):
    root = slx.region_root
    assert root.id == "Regions"
    assert root.children == ["747deb05-20a9-4b56-8ffd-0eaef5a1f935"]
    group = slx.regions[root.children[0]]
    assert len(group.children) == 4
    assert all(slx.regions[c].children == [] for c in group.children)


@needs_file
def test_project_json_and_report(slx):
    assert slx.project["scils-lab-0"]["massaxis"]["datapoints"] == "460519"
    (source,) = slx.project["scils-lab-0"]["sources"].values()
    names = [r["name"] for r in source["Regions"]]
    assert names == ["01", "02", "03", "04"]
    assert "25 spectra were skipped" in slx.import_report


@needs_file
def test_file_is_not_modified():
    before = SLX_PATH.stat().st_mtime_ns
    read_slx(SLX_PATH)
    assert SLX_PATH.stat().st_mtime_ns == before


def test_synthetic_slx_reads(make_slx):
    slx = read_slx(make_slx())
    assert slx.schema_version == 12
    assert slx.n_spots == 7
    assert len(slx.features) == 2


def test_unknown_schema_version_is_rejected(make_slx):
    with pytest.raises(UnsupportedFormat, match="schema version 13"):
        read_slx(make_slx(version=13))
