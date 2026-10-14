# SCiLS Lab .slx file format

Source file: a private SCiLS Lab 14 project (schema 12, TSF source, 94713 spots), 21,790,720 bytes. Opened read-only.
Reader: `open_slx/slx.py`. Tests: `tests/test_slx.py`.

## Answer to the main question

The .slx file holds no per-spot intensity data.
It also holds no per-spot x/y coordinates and no region spot lists.
These live in the `.sbd` file (HDF5) and in the Bruker `.d` folder. Neither is present.
The `.slx` holds the m/z axis, the feature table, the spot count, raster offsets, and the region tree.

Evidence:

- Tables `Registrations`, `SpotImages`, `Regions`, `HeatMap`, `MeanSpectra`, `PeakLists`, `OpticalImages` have 0 rows.
- `DataSources.Rasters` is empty.
- `Labels.SpotList` holds only the numbers 0 to 94712.
- The project JSON names the storage file `<project>.sbd` (type `HDF5_PEAKS`).
- The free pages hold only deleted 0.0/1.0 masks (see "Free pages").

## Container

| Item | Value |
| --- | --- |
| Format | SQLite 3, journal mode delete |
| Page size / pages | 4096 / 5320 |
| Free pages | 2698 (about 11 MB, half the file) |
| `Master` | `Version` = 12, `open` = false |
| `ObjectSchemaVersion` | `MetaInformation` = 3 |
| Blob compression | none found. All blobs are raw little-endian arrays |
| Row `MetaInformation.secure` | 692-character base64 text. After a second base64 decode it is 368 opaque bytes. Not decoded |

## Non-empty tables

| Table | Rows | Content |
| --- | --- | --- |
| `Master` | 2 | Version 12, open flag |
| `ObjectSchemaVersion` | 1 | schema version of `MetaInformation` |
| `MetaInformation` | 37 | key/value strings per object (see below) |
| `ObjectRelations` | 60 | object types, parents, relations |
| `DataSources` | 1 | source of the raw data |
| `Rasters` | 4 | raster offsets |
| `FeatureTable` | 15 | one row per column, data in a blob |
| `Labels` | 1 | spot list and a spanning tree of the spots |
| `SamplePositions` | 2 | m/z axis, stored twice |

All other tables have 0 rows.

## MetaInformation string encoding

Each string is a list of `<length>:<text>` items.
The length counts UTF-8 bytes, not characters. `±` counts as 2.
Items alternate between key and value.
Example: `4:name15:1/K0 width [±]4:type13:MzTableColumn` gives `{name: "1/K0 width [±]", type: "MzTableColumn"}`.
Parser: `parse_meta`.

| Row ID | Content |
| --- | --- |
| `file` | project settings, import report, project JSON |
| DataSource ID | `Source/Number of spectra in dataset` = 94738, `Source/Spectrum type` = Centroid, `Source/Method` (path of the Bruker method), `Source/TimsLibraryVersion` = 2.24.0.0 |
| Raster IDs | `name` = 0 to 3 |
| Feature column IDs | `ContentId`, `name`, `flags` |
| `<id>_BookmarkInfo` | view settings (colors, zoom). 8 rows. One lists the region colors |
| Normalization ID | `Total Ion Count`, active |
| InteractiveLabels ID | segmentation `ParkerSegTest1`: Bisecting k-Means, Correlation Distance, Peak Area, Total Ion Count |
| SamplePositions IDs | names `SamplePositions of 1v3roc` and `SamplePositions of 2v4roc` |

Key values in `file`:

| Key | Value |
| --- | --- |
| `nSpots` | 94713 |
| `AxisType` | OTOF |
| `SpectraType` | Peaks |
| `IntervalWidth` | 20 (mode 2) |
| `OptimalPeakWidth` | 15.105 |
| `IntervalProcessingMode` | 3 |
| `StorageVersion` | 16 |
| `Produced by` | SCiLS Lab Version 14.00.17781 |
| `Import/Project` | JSON text (see below) |
| `Comment` | import report text |

The import report says 25 of 94738 spectra were empty and skipped. This gives 94738 − 25 = 94713 spots.
It lists the first 20 skipped spectrum numbers: 925, 1328, 1332, 1433, 1434, 1542, ...

## DataSources (1 row)

| Column | Value |
| --- | --- |
| `Type` | `TSFSource` |
| `Polarity` | negative |
| `SpotNameFormat` | `0_R%02dX%03dY%03d` (source 0, raster R, x, y) |
| `SourcePath` | Windows path of the `.d` folder. The folder is not present |
| `OriginalImportPath` | Windows path of the `.mis` file |
| `Rasters`, `RecalibrationObjId` | empty |

## Rasters (4 rows)

| ScilsRasterIdx | OriginalRasterIdx | RasterXOffset | RasterYOffset |
| --- | --- | --- | --- |
| 0 | 0 | 1131 | 1243 |
| 1 | 1 | 2504 | 1417 |
| 2 | 2 | 4093 | 2059 |
| 3 | 3 | 6124 | 1544 |

Offsets are integers. The unit is the scan grid step. The meaning is inferred from the column names. It is not verified.

## Project JSON (`file` / `Import/Project`)

| Path | Value |
| --- | --- |
| `sources.<name>.transform-matrix` | 4x4 matrix. Unit scale. Translation x = 26340.5838, y = −7120.7200 |
| `sources.<name>.Regions` | four raster regions: names `01` to `04`, ids 0 to 3 |
| `massaxis` | type OTOF, min 200, max 2000, `resmz` 1000, `datapoints` 460519 |
| `storage` | `HDF5_PEAKS`, compression ON, file `<project>.sbd` |
| `baselineremoval` | None |
| `featurefinding` | disabled |

No pixel size field exists. The method file name contains `10um`. The reader does not use this.

## Spots (`Labels.SpotList`)

| Property | Value |
| --- | --- |
| dtype | uint64, little-endian |
| Length | 94713 (757704 bytes) |
| Values | 0, 1, 2, ..., 94712. Consecutive |

The spot number is the import order of the non-empty spectra.
The file does not map a spot number to a raster, x, or y.

Open point: the skipped spectrum numbers and the row widths suggest row-major scan order.
This is not enough to rebuild coordinates.

## Spot spanning tree (`Labels.MergeList`)

| Property | Value |
| --- | --- |
| SQL type | text (2,202,685 characters). `Labels.Labels` is NULL |
| Layout | `NoRemove:94712 a b d a b d ...` |
| Items | 94712 triples: spot a, spot b, distance d |
| Example | `70960 71308 1.37676e-05` |
| Range of d | 9.66e-06 to 2.0 |

The edges join all 94713 spots in one component. This is a spanning tree.
Evidence: 94712 edges, 94712 successful union-find merges, 1 component.
The distance range matches a correlation distance (1 − r). The segmentation used this metric.
The file does not store segment labels.
The name `NoRemove` and the purpose of the tree are inferred.

## m/z axis (`SamplePositions`, 2 rows)

| Property | Value |
| --- | --- |
| dtype | float64, little-endian |
| Length | 460519 (3,684,152 bytes) |
| Min / max | 200.0 / 2000.0 |
| First values | 200.0, 200.001, 200.00200001, 200.00300002 |
| Spacing | geometric. Ratio of neighbors 1.000005 |
| Formula | `mz[k] = 200 * 10 ** (k / 460518)`, error below 1e-12 relative |
| Step | 0.001 at 200 and 0.01 at 2000 |
| The two rows | byte-identical |

`ObjectRelations` also lists an object `GlobalMassAxis` of type `SamplePositions`.
It has no table row. The two stored rows are copies of the same axis.

## FeatureTable (15 columns, 322 features)

Each row is one column. `ColumnData` is a blob. `ColumnType` is stored as text.
`MetaInformation` gives the meaning through `ContentId`.
A NULL blob means the column holds only default values.

| ColumnType | Encoding |
| --- | --- |
| 0 | strings: uint64 count, then per string a uint64 length and the UTF-8 bytes |
| 1 | float64 array |
| 2 | float64 pairs (low, high) |
| 3 | int64 array |
| 4 | uint8 flags. Trailing false values are cut off |
| 5 | not seen with data |

| Idx | ContentId / name | Type | Bytes | Decoded |
| --- | --- | --- | --- | --- |
| 0 | `Basic:Id` / ID | 3 | 2576 | 322 int64. Rows 0 to 35 = 0 to 35. Rows 36 to 321 = −1 |
| 1 | `Basic:Mz` / m/z | 2 | 5152 | 322 pairs (low, high). Row 0 = 255.23107067, 255.23492596 |
| 2 | `IntervalWidths:Mz` / m/z width | 1 | 2576 | 322 float64. Equal to high − low (max difference 0). Row 0 = 0.00385529. Range 0.00386 to 0.02900 |
| 3 | `Basic:OneOverK0` | 2 | NULL | no data |
| 4 | `IntervalWidths:OneOverK0` | 1 | NULL | no data |
| 5 | `Basic:CCS` | 2 | NULL | no data |
| 6 | `IntervalWidths:CCS` | 1 | NULL | no data |
| 7 | `Basic:Color` | 0 | 4838 | 322 strings, 24 distinct. Row 0 = `#ff0000` |
| 8 | `Basic:Name` | 0 | NULL | no names |
| 9 | `Intensity:RootRegionIntensity` / Intensity [Regions] | 1 | 2576 | 322 float64: average intensity over the root region, Peak Area mode. Row 0 = 1106.018. Range 14.29 to 81723.2 |
| 10 | `Basic:RelatedSpotImageID` | 0 | NULL | no data |
| 11 | `Grouping:Parent` | 3 | NULL | no data |
| 12 | `Grouping:Type` | 5 | NULL | no data |
| 13 | tag column `manual` | 4 | 32 | rows 0 to 31 true, the rest false |
| 14 | tag column `test` | 4 | 36 | rows 32 to 35 true, the rest false |

Feature m/z values: low 255.2311, high 1449.9896.
The centers of the first three rows are 255.23300, 279.23303, 281.24794.
The reader returns `mz_center = (low + high) / 2` and `mz_width = high − low`.
The first 36 rows have an ID. They are the features that the user tagged. The other 286 rows have ID −1.
This reading of the ID and the tag columns is inferred.
The rows are not sorted by m/z.

## ObjectRelations (60 rows)

Columns: `ID`, `type`, `parent`, `own`, `relations`.
`relations` is a list of groups separated by `;`. Each group is `id,id,...,<type>`.
Example: `c2cf86c8-...,Normalization;292194ef-...,Region;GlobalMassAxis,SamplePositions`.
Parser: `parse_relations`.

| Type | Count | Note |
| --- | --- | --- |
| `MeanSpectrum` | 18 | links a region, a normalization, and the axis. The spectra are not stored |
| `Region` | 6 | the region tree |
| `AttributeRegion`, `AttributeValueRegion` | 1, 2 | region attribute sets. No table rows |
| `Normalization` | 2 | one is `Total Ion Count` |
| `OpticalImage` | 2 | the image data is not stored (`OpticalImages` is empty) |
| `ThresholdMzRange` | 1 | no table row |
| `Registration` | 1 | ID `0`. `Registrations` is empty |
| `SpectraGroup` | 1 | |
| `SamplePositions` | 3 | `GlobalMassAxis` (`own` = 0) and the two stored rows (`own` = 1) |
| `DataSource`, `Raster`, `MzTable`, `MzTableColumn`, `InteractiveLabels`, `Bookmark` | 1, 4, 1, 15, 1, 1 | rows with `own` = 1 |

## Region tree

The tree has 6 nodes. The file stores no region name and no spot list.
The `Regions` table is empty and no `MetaInformation` row names a region.

| Node ID | Parent | Children |
| --- | --- | --- |
| `Regions` | none | `747deb05-20a9-4b56-8ffd-0eaef5a1f935` |
| `747deb05-20a9-4b56-8ffd-0eaef5a1f935` | `Regions` | the four nodes below |
| `0dff06cf-3c4c-4728-aec0-17e18ebfa0f0` | `747deb05` | none |
| `195aa1ae-ee05-49db-abd2-877ae085210d` | `747deb05` | none |
| `292194ef-75e9-4d44-a3b0-821cef340cd6` | `747deb05` | none |
| `cfd4c1a9-568f-4e7b-8a1f-5913640e8f65` | `747deb05` | none |

The four leaves probably match the four rasters (`01` to `04` in the project JSON).
The file does not state which leaf matches which raster.
The segmentation took `747deb05` as input.

## Free pages

SQLite free pages: 2698. They form 6 chains of deleted overflow pages.
Each page starts with a 4-byte big-endian pointer to the next page. The rest is payload.

| First page | Pages | Content |
| --- | --- | --- |
| 2622 | 118 | float64 values 0.0 and 1.0 (60357 values) |
| 4422 | 235 | float64 values 0.0 and 1.0 (120202 values) |
| 4659 | 663 | float64 values 0.0 and 1.0. The pages are not contiguous |
| 2741 | 781 | float64 values 0.0 and 1.0. The pages are not contiguous |
| 3523 | 897 | float64 values 0.0, 0.5, 1.0. Last non-zero index 458811 (close to the 460519-point axis) |
| 3522 | 1 | 5 non-zero bytes |

The values are masks. They are deleted data. They are not spectra and not intensities.
What each mask was is not known. A mask along the m/z axis is a possible reading of the last chain.

## Not in the file

- Per-spot intensities or spectra.
- Per-spot raster index, x, y.
- Region names and region spot membership.
- Pixel size.
- Feature names, 1/K0 values, CCS values.
- Segment labels (`Labels.Labels` is NULL).
- Optical image bytes and registration transforms.

## Open points

| Point | Status |
| --- | --- |
| Meaning of `RasterXOffset` / `RasterYOffset` | inferred from the names. Not checked against spot coordinates |
| Match of leaf regions to rasters | unknown |
| Purpose of the spanning tree | inferred |
| ColumnType 5 encoding | no data in this file |
| `secure` row | opaque |
| Free page masks | purpose unknown |
| Other .slx versions | only this file (schema 12) is checked |

To export imzML you need the `.sbd` file or the `.d` folder in addition to this file.
