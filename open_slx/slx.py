import json
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from open_slx.errors import UnsupportedFormat

SUPPORTED_SCHEMA_VERSIONS = {12}
STRING, DOUBLE, INTERVAL, INT64, BOOL, UNKNOWN = range(6)


def parse_meta(raw: bytes) -> dict[str, str]:
    parts = []
    i = 0
    while i < len(raw):
        j = raw.index(b":", i)
        n = int(raw[i:j])
        parts.append(raw[j + 1 : j + 1 + n].decode("utf-8"))
        i = j + 1 + n
    return dict(zip(parts[::2], parts[1::2]))


def parse_relations(text: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for group in filter(None, text.split(";")):
        *ids, kind = group.split(",")
        out[kind] = ids
    return out


@dataclass
class DataSource:
    id: str
    source_path: str
    original_import_path: str
    type: str
    polarity: str
    spot_name_format: str
    meta: dict[str, str]


@dataclass
class Raster:
    id: str
    scils_index: int
    original_index: int
    x_offset: int
    y_offset: int


@dataclass
class ObjectRelation:
    id: str
    type: str
    parent: str
    own: bool
    relations: dict[str, list[str]]


@dataclass
class RegionNode:
    id: str
    parent: str
    children: list[str]


@dataclass
class FeatureTable:
    ids: np.ndarray
    mz_low: np.ndarray
    mz_high: np.ndarray
    intensity: np.ndarray | None
    colors: list[str]
    names: list[str] | None
    tags: dict[str, np.ndarray]

    @property
    def mz_center(self) -> np.ndarray:
        return (self.mz_low + self.mz_high) / 2

    @property
    def mz_width(self) -> np.ndarray:
        return self.mz_high - self.mz_low

    def __len__(self) -> int:
        return len(self.ids)


@dataclass
class SpanningTree:
    edges: np.ndarray
    distances: np.ndarray


@dataclass
class Slx:
    path: str
    schema_version: int
    meta: dict[str, dict[str, str]]
    objects: dict[str, ObjectRelation]
    data_source: DataSource
    rasters: list[Raster]
    n_spots: int
    spot_ids: np.ndarray
    spot_tree: SpanningTree | None
    features: FeatureTable
    sample_positions: dict[str, np.ndarray]
    regions: dict[str, RegionNode]
    project: dict = field(default_factory=dict)
    import_report: str = ""

    @property
    def mz_axis(self) -> np.ndarray:
        return next(iter(self.sample_positions.values()))

    @property
    def region_root(self) -> RegionNode:
        return next(r for r in self.regions.values() if r.parent not in self.regions)


def _strings(blob: bytes) -> list[str]:
    n = int(np.frombuffer(blob, "<u8", 1)[0])
    out = []
    pos = 8
    for _ in range(n):
        size = int(np.frombuffer(blob, "<u8", 1, pos)[0])
        out.append(blob[pos + 8 : pos + 8 + size].decode("utf-8"))
        pos += 8 + size
    return out


def _read_features(con: sqlite3.Connection, meta: dict[str, dict[str, str]]) -> FeatureTable:
    by_content = {}
    tags_raw = {}
    for col_id, ctype, data in con.execute(
        "select ID, ColumnType, ColumnData from FeatureTable order by ColumnIndex"
    ):
        m = meta.get(col_id, {})
        content = m.get("ContentId")
        if content:
            by_content[content] = data
        elif int(ctype) == BOOL:
            tags_raw[m.get("name", col_id)] = data

    ids = np.frombuffer(by_content["Basic:Id"], "<i8")
    n = len(ids)
    pairs = np.frombuffer(by_content["Basic:Mz"], "<f8").reshape(-1, 2)
    intensity_key = next((k for k in by_content if k.startswith("Intensity:")), None)
    intensity = None
    if intensity_key and by_content[intensity_key]:
        intensity = np.frombuffer(by_content[intensity_key], "<f8")
    colors = _strings(by_content["Basic:Color"]) if by_content.get("Basic:Color") else []
    names = _strings(by_content["Basic:Name"]) if by_content.get("Basic:Name") else None
    tags = {}
    for name, data in tags_raw.items():
        flags = np.zeros(n, bool)
        values = np.frombuffer(data or b"", "u1").astype(bool)
        flags[: len(values)] = values
        tags[name] = flags
    return FeatureTable(ids, pairs[:, 0], pairs[:, 1], intensity, colors, names, tags)


def _read_spot_tree(text: str | None) -> SpanningTree | None:
    if not text:
        return None
    _, rest = text.split(" ", 1)
    arr = np.array(rest.split(), dtype=float).reshape(-1, 3)
    return SpanningTree(arr[:, :2].astype(np.int64), arr[:, 2])


def read_slx(path: str | Path) -> Slx:
    path = Path(path).resolve()
    con = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    try:
        schema_version = int(con.execute("select ValueColumn from Master where KeyColumn='Version'").fetchone()[0])

        if schema_version not in SUPPORTED_SCHEMA_VERSIONS:
            raise UnsupportedFormat(
                f".slx schema version {schema_version} is not supported "
                f"(supported: {', '.join(map(str, sorted(SUPPORTED_SCHEMA_VERSIONS)))})"
            )

        con.text_factory = bytes
        meta = {
            k.decode(): parse_meta(v)
            for k, v in con.execute("select ID, MetaInformationString from MetaInformation")
            if k != b"secure" and v
        }
        con.text_factory = str

        objects = {
            r[0]: ObjectRelation(r[0], r[1], r[2], bool(r[3]), parse_relations(r[4]))
            for r in con.execute("select ID, type, parent, own, relations from ObjectRelations")
        }

        ds = con.execute(
            "select ID, SourcePath, OriginalImportPath, Type, Polarity, SpotNameFormat from DataSources"
        ).fetchone()
        data_source = DataSource(*ds, meta=meta.get(ds[0], {}))

        rasters = [
            Raster(*r)
            for r in con.execute(
                "select ID, ScilsRasterIdx, OriginalRasterIdx, RasterXOffset, RasterYOffset "
                "from Rasters order by ScilsRasterIdx"
            )
        ]

        spot_blob, merge_text = con.execute("select SpotList, MergeList from Labels").fetchone()
        spot_ids = np.frombuffer(spot_blob, "<u8")

        sample_positions = {
            meta.get(i, {}).get("name", i): np.frombuffer(b, "<f8")
            for i, b in con.execute("select ID, MzValues from SamplePositions order by ID")
        }

        regions = {
            o.id: RegionNode(o.id, o.parent, o.relations.get("Region", []))
            for o in objects.values()
            if o.type == "Region"
        }

        file_meta = meta.get("file", {})
        return Slx(
            path=str(path),
            schema_version=schema_version,
            meta=meta,
            objects=objects,
            data_source=data_source,
            rasters=rasters,
            n_spots=int(file_meta["nSpots"]),
            spot_ids=spot_ids,
            spot_tree=_read_spot_tree(merge_text),
            features=_read_features(con, meta),
            sample_positions=sample_positions,
            regions=regions,
            project=json.loads(file_meta["Import/Project"]) if "Import/Project" in file_meta else {},
            import_report=file_meta.get("Comment", ""),
        )
    finally:
        con.close()
