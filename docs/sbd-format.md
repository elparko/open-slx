# SCiLS Lab .sbd file format (HDF5)

Source file: Zenodo 21251439, `Romics_paper_SCils_file_case_2.sbd` (5,195,696,605 bytes).
HDF5 superblock version 3. Produced by SCiLS Lab 14.00.17781. Root dataset `Version` = 36.
The root attributes `nSpots` = 283631, `AxisType` = MRMS, `SpectraType` = Peaks, `uuid` are text arrays of length 1.

## Root tree

| Path | Shape / content |
| --- | --- |
| `SpectraGroup/` | all spectra (see below) |
| `SamplePositions/GlobalMassAxis/SamplePositions` | float64 (661027,). m/z axis, 245.69 to 1100.01 |
| `Registrations/0/` | spot coordinates (see below) |
| `MeasurementRasterIndices` | int32 (1, 283631). Raster number per spot. 41 rasters, sorted, one run each |
| `Regions/` | region tree. Each region has `SpotList` (uint64), `PolygonPoints`, `PolygonOffsets`, `MeanSpectra/` |
| `MeanSpectra/<id>/` | float32 (661027,): `Intensities`, `Skyline` (maximum), `Variances`. Axis-sized |
| `MzRanges/`, `OpticalImages/`, `AttributeRegions/` | imported feature lists, images, annotations |
| `Bookmarks`, `Labels`, `ClassificationModels`, `ComponentAnalyses`, `RecalibrationFiles`, `StorageSnapshots` | empty groups |
| `secured` | opaque string (2680 bytes) |

No dataset uses chunk compression filters, except where noted. No shuffle filter at the HDF5 level.

## SpectraGroup

| Dataset | Shape | dtype | Chunks | Meaning |
| --- | --- | --- | --- | --- |
| `compressedMz` | (613623218,) | int8 | 4096 | one zstd frame plus a 17-byte trailer per spot. Attribute `compressionType` = `zstd` |
| `mzOffsets` | (283631,) | uint64 | none | byte start of each spot in `compressedMz`. The end of the last spot is the array length |
| `intensities` | (217447635,) | float32 | 4096 | peak intensities of all spots, in spot order |
| `intensityOffsets` | (283631,) | uint64 | none | start of each spot in `intensities`. Example: 0, 426, 857, 1267 |
| `imageIntensities`, `imageIntensityOffsets`, `compressedImageSpots`, `imageSpotOffsets` | 217447512 / 661027 / 76185842 / 661027 | float32, uint64, int8, uint64 | | the same peaks sorted by axis position (per-m/z lists of spots). Not needed for export |
| `Normalizations/<id>/Values` | (283631,) float32 | | `sum(intensities) / TargetValue` for `Total Ion Count`. Also `Root Mean Square` |

Peak count of spot i = `intensityOffsets[i+1] - intensityOffsets[i]`.
Example: spot 0 has 426 peaks, spot 1 has 431.
The spectra are sparse peak lists. They are not profiles on the shared axis.

## m/z blob of one spot

Bytes `compressedMz[mzOffsets[i] : mzOffsets[i+1]]` have two parts.

| Part | Size | Content |
| --- | --- | --- |
| zstd frame | variable | decompresses to `k * n` bytes (n = peak count) |
| trailer | 17 bytes | `uint8 k`, `float64 scale`, `float64 mz0` (little-endian) |

Example trailers:

| Spot | Trailer (hex) | k | scale | mz0 |
| --- | --- | --- | --- | --- |
| 0 | `04 0000000080841e41 dae7cfa1c2217040` | 4 | 500000.0 | 258.1100176 |
| 1 | `03 0000000080841e41 526855afc2217040` | 3 | 500000.0 | 258.1100305 |
| 281771 | `04 0000000080841e41 6db2625bcf346f40` | 4 | 500000.0 | 249.6503121 |
| 281861 | `03 0000000080841e41 2ee496ed6d296f40` | 3 | 500000.0 | 249.2946690 |

Decoding:

1. Decompress the frame. Use `zstandard` with `decompressobj`. The unused data after the frame is the trailer.
2. Reshape the bytes to `(k, n)`. Row p holds byte p of every value (byte planes, low byte first).
3. Combine the planes to unsigned integers `d[0..n-1]`. `d[0]` is always 0.
4. `mz[j] = mz0 + cumsum(d)[j] / scale`.

The first peak is at `mz0`. Each value after it is the gap to the previous peak in units of 1/500000 m/z.
Peaks are sorted by m/z. The number of planes k is the smallest byte count that holds the largest gap.

Check: the region `Regions/1ec86fac.../b94375a1...` has 13 spots (281771 to 282040).
Its `Skyline` mean spectrum lets us find the axis index of single peaks by exact intensity match.
For 4 spots (166 to 278 matched peaks each), `axis[j] - cumsum(d)/500000` is constant with median absolute deviation 0.0007 to 0.0013.
This equals the axis spacing at that m/z. The constant equals the trailer `mz0` (249.65 or 249.295).

## Decodings ruled out

| Idea | Why it fails |
| --- | --- |
| `d` are indices on the axis | values reach 22.7 million. The axis has 661027 points |
| `d` are deltas of axis indices | the sum exceeds 661027 within a few peaks |
| float32 m/z, byte-shuffled | high plane holds 0, not 0x43 or 0x44 |
| XOR or float delta | gives denormal numbers |
| absolute m/z in the blob | the first value is 0 in every spot. The trailer holds the start |
| `mz0` constant for the file | it differs per spot (249.29, 249.65, 258.11, 273.04) |

## Spot coordinates (`Registrations/0`)

Attributes: `Initial` = true, `name` = Initial, `Dimensions` = 2.

| Dataset | Shape | dtype | Example (spots 0 to 4) |
| --- | --- | --- | --- |
| `Coordinates` | (3, 283631) | float32 | x: 1233.0082, 1133.0082, 1183.0082. y: 7355.8735, 7305.8735. z: 0. World coordinates in micrometers |
| `OriginalCoordinates` | (2, 283631) | uint64 | x: 414, 410, 412. y: 572, 574. Pixel numbers of the instrument raster |
| `OriginalSpectrumIndices` | (283631,) | uint64 | 0, 1, 2, 3, 4. Not monotone over the whole file |
| `SpotSizes` | (3, 283631) | float32 | 50, 50, 0 (raster 0). 25 for rasters of 25 micrometers |
| `Rasters/<r>/PixelToWorld` | (4, 4) | float64 | raster 0: x scale 50, y scale -50, translation (83.0082, 7355.87) |
| `Rasters/<r>/HasOriginalCoordinates` | (1,) | int32 | 1 |

There are 41 rasters (0 to 40). Rasters have 1420 to 16009 spots. Pixel size is 50 or 25 micrometers.

Pixel index of a spot in its raster:

- `px = (Coordinates[0] - M[0,3]) / M[0,0]`
- `py = (Coordinates[1] - M[1,3]) / M[1,1]`

where M = `PixelToWorld` of the raster. All 41 rasters give whole numbers (error 0.0).
Pixels are unique in each raster. Raster 0 has px 0 to 59 and py 0 to 35.
`OriginalCoordinates` step is 2 for 50 micrometer rasters and 1 for 25 micrometer rasters.

## Link to .slx spot ids

Spot k of the .sbd is spot k of the .slx (`Labels.SpotList` runs 0 to nSpots-1).
Empty spectra are skipped in both files. `OriginalSpectrumIndices[k]` is the spectrum number in the instrument data.
The .slx JSON gives the m/z axis type and size. The .sbd stores its own axis.

## Check of the decoding

First 200 spots (raster 0, 100677 peaks, x 10 to 50, y 1 to 10): the m/z range is 249.6509 to 1098.3752.
All intensities are above 0. For each spot, the sum of intensities equals `Total Ion Count Values * TargetValue`
(target 349941746.77; spot 0: 141148770 against 0.40334922 * 349941746.77) within 1e-4 relative.
This tests the intensity offsets and the peak counts. The m/z decoding is tested with the skyline match above.

## Converter (`open-slx`)

- `x = px + 1`, `y = py + 1` from the formula above. The pixel size is `abs(M[0,0])`, `abs(M[1,1])` in micrometers.
- One imzML file per raster, because the pixel size and the pixel grid differ between rasters.
  One raster: `OUT.imzML`. More rasters: `OUT_r00.imzML`, `OUT_r01.imzML`, and so on.
- `processed` mode (default) writes the peaks as they are: float64 m/z and float32 intensity per spot.
- `continuous` mode writes the global axis once. Each spot has 661027 float32 values (2.6 MB). Peaks go to the nearest axis point.
  Use it only with `--max-spots`.

## Open points

| Point | Status |
| --- | --- |
| .sbd of the private project | not available. Layout is checked on the Zenodo file only |
| Intensity normalization | `intensities` hold raw values. The normalization factors are separate |
| Region names for rasters | region tree `Regions/<file>/<file>/<region>` holds names and `SpotList` |
