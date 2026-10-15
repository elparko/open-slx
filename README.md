<p align="center"><img src="assets/logo.svg" width="96" alt="open-slx logo"></p>

# open-slx

open-slx converts SCiLS Lab projects to imzML.
It is tested on macOS (Apple silicon). It is written in Python and uses h5py, numpy, and zstandard. It is not tested on Linux or Windows.
It needs no SCiLS Lab installation, no SCiLS Python API, and no licence key.

imzML is the open file format for mass spectrometry imaging (MSI). MSI maps the mass spectra of a tissue section, spot by spot.
Many MSI instruments use matrix-assisted laser desorption/ionization (MALDI).

## Input files

A SCiLS Lab project has two files. open-slx needs both.

| File | Format | Content |
| --- | --- | --- |
| `PROJECT.slx` | SQLite | Feature lists (m/z low and high), region tree, labels, a reference to the raw data source, project settings. It holds no intensities and no spot coordinates. |
| `PROJECT.sbd` | HDF5 (Hierarchical Data Format 5) | Spectra, spot coordinates in micrometres, raster index per spot, global m/z axis, regions with spot lists, mean spectra, normalizations. |

The `.sbd` file normally has the same base name as the `.slx` file. It is in the same folder.
If no such file exists, open-slx uses the file name in the `.slx` storage settings (key `storage`, field `file`, shown in `OUT_metadata.json`). It looks for that name in the folder of the `.slx` file.

Common mistake: a user copies only the `.slx` file.
SCiLS Lab writes both files when it imports data.
Copy both files. The `.sbd` file is usually several gigabytes.

The original Bruker `.d` folder is not necessary.

## If you have only the .slx file

- SCiLS Lab writes the `.sbd` file at import. It is next to the `.slx` file and has the same base name. Look in the same folder on the computer that ran SCiLS Lab.
- The `.slx` table `DataSources` names the path of the original raw data. The command `open-slx PROJECT.slx` writes this path to `OUT_metadata.json`, under `data_source` (key `source_path`). It does this even when the `.sbd` file is missing.
- If the `.sbd` file is lost, use one of these options:
  - Import the Bruker raw data (`.d` folder) into SCiLS Lab again. This makes a new `.sbd` file.
  - Export imzML directly from SCiLS Lab, if you have a licence.
- open-slx does not read Bruker `.d` folders.

## Install

Use `uv`.

```
uv tool install git+https://github.com/elparko/open-slx
```

## Usage

```
open-slx [PROJECT.slx] [--sbd PATH_OR_URL] [-o OUT] [--mode processed|continuous] [--max-spots N]
```

| Option | Meaning |
| --- | --- |
| `--sbd PATH_OR_URL` | The data file. The default is the `.sbd` file next to the `.slx` file. If that file does not exist, open-slx uses the file name in the `storage.file` entry of the `.slx` project, in the same folder as the `.slx` file. An http or https URL is allowed. open-slx then reads the file in parts with byte ranges. |
| `-o OUT` | The base name of the output files, with an optional folder. open-slx removes a trailing `.imzML` or `.ibd` that you type. A dot in the name is safe: `-o 2023.05.12_kidney` gives `2023.05.12_kidney.imzML`. |
| `--mode processed` | Default. Each spot has its own m/z values and intensities. The spectra are centroid peak lists. |
| `--mode continuous` | Every spot uses the shared m/z axis. The files are much larger. Use it with `--max-spots`. |
| `--max-spots N` | Convert only the first N spots. N must be 1 or more. |

### Default output folder and name

If you do not give `-o`, open-slx writes to the current folder.
The name is the base name of the `.slx` file.
With `--sbd` and no `.slx` file, the name is the base name of the `.sbd` file or URL path.

### Output names for more than one raster

A raster is one scanned area of tissue. Each raster has its own pixel size and pixel grid.
open-slx writes one imzML file for each raster.

- One raster: `OUT.imzML` and `OUT.ibd`.
- More than one raster: `OUT_r00.imzML` and `OUT_r00.ibd`, then `OUT_r01.imzML` and `OUT_r01.ibd`, and so on.

The number after `_r` is the raster index in the `.sbd` file.

### Mode without a `.slx` file

You can run open-slx with `--sbd` only:

```
open-slx --sbd PROJECT.sbd -o PROJECT
```

In this mode open-slx writes the imzML and `.ibd` files and `OUT_metadata.json`.
It does not write `OUT_features.csv`. The feature list is only in the `.slx` file.
This `OUT_metadata.json` has three keys:

| Key | Content |
| --- | --- |
| `sbd` | The path or URL of the `.sbd` file. |
| `n_spots` | The number of spots in the `.sbd` file. |
| `attributes` | The root attributes of the `.sbd` file as text, for example `AxisType`, `SpectraType`, `nSpots`, `uuid`. |

The check of the spot count between the two files (exit code 3) does not run in this mode.

### Exit codes

| Code | Meaning |
| --- | --- |
| 0 | The conversion finished. |
| 1 | An unexpected error occurred. open-slx deletes the partial `.imzML` and `.ibd` files. |
| 2 | An input file is missing, or the command line is wrong (argparse also uses 2). |
| 3 | The `.sbd` file and the `.slx` file have different spot counts. |
| 4 | The version is not supported (see "Limits"). |

If the `.sbd` file is missing, open-slx does this, in this order:

- It writes `OUT_features.csv` and `OUT_metadata.json` from the `.slx` file.
- It prints a message that names the missing file.
- It exits with code 2.

### Example session

open-slx prints one line for each imzML file. The line shows the raster index, the spot count, the pixel size in micrometres, and the path.
This session is real output for a small test file with two rasters.

```
$ open-slx --sbd PROJECT.sbd -o PROJECT
raster 0: 3 spots, pixel size 50 x 50 um, PROJECT_r00.imzML
raster 1: 4 spots, pixel size 25 x 25 um, PROJECT_r01.imzML
```

Errors go to standard error. They start with `error:`.

## Output files

| File | Content |
| --- | --- |
| `OUT.imzML` | XML metadata (imzML 1.1): pixel size, maximum x and y, UUID, SHA-1 checksum of the `.ibd` file. |
| `OUT.ibd` | Binary data: m/z as float64, intensity as float32. |
| `OUT_features.csv` | The feature list from the `.slx` file. |
| `OUT_metadata.json` | Project metadata from the `.slx` file. |

### Coordinate to pixel rule

The `.sbd` file stores each spot as world coordinates in micrometres.
For each raster, the file has a 4 by 4 matrix `M` (`PixelToWorld`).
`M[0,3]` and `M[1,3]` are the translation (the world position of pixel 0). `M[0,0]` and `M[1,1]` are the pixel pitch (the pixel size, y is negative).
For a spot with world coordinates X and Y, open-slx writes:

- `x = round((X - M[0,3]) / M[0,0]) + 1`
- `y = round((Y - M[1,3]) / M[1,1]) + 1`

The `+ 1` makes the pixel numbers start at 1, as the imzML format requires.
The pixel size in the imzML file is `abs(M[0,0])` and `abs(M[1,1])`.

pyimzML opens the output (tested). The files follow imzML 1.1, so other imzML readers such as M2aia, Cardinal and METASPACE are expected to open them. These readers are not tested yet.

## How it works

The `.slx` file is an SQLite database. open-slx reads the feature lists, the region tree, and the project settings from it.

The `.sbd` file is an HDF5 file. For each spot, it stores the intensities as float32 values.
It stores the m/z values as gaps between neighbouring peaks. It compresses these gaps with zstd.
open-slx decodes the gaps and adds them to get the m/z values.

The full format notes are in these files:

- [docs/slx-format.md](docs/slx-format.md)
- [docs/sbd-format.md](docs/sbd-format.md)

## Testing and verification

- Test data: the Zenodo 21251439 `.sbd` file. SCiLS Lab 14.00.17781 wrote it. It has 283631 spots and centroid peak lists.
- Automated network test (`test_zenodo_first_200_spots` in `tests/test_cli.py`): we converted the first 200 spots. For every spot, the sum of intensities equals the stored SCiLS Total Ion Count normalization times its `TargetValue`, within 1e-4 relative. The test asserts this. The m/z range of these spots is 249.6509 to 1098.3752.
- Manual end-to-end run on macOS (Apple silicon): we ran the installed command with `--sbd <Zenodo URL> --max-spots 300`. pyimzML opens the output. Result: 300 spots, x 6 to 51, y 1 to 12, pixels of 50 x 50 um, m/z 249.2955 to 1098.3752, 426 peaks in spectrum 0.
- m/z check: the decoded peaks match a region Skyline (maximum) spectrum on the global axis within one axis step. The median absolute deviation is 0.0007 to 0.0013 in m/z, for 4 spots in one region.
- Unit tests: "Run the tests" gives the run counts.

## Limits

- open-slx reads SCiLS Lab 14 files (version 14.00.17781): `.sbd` `Version` 36 and `.slx` schema version 12.
- open-slx reads centroid data only: the `.sbd` attribute `SpectraType` must be `Peaks`.
- For a `.sbd` `Version` other than 36, a `.sbd` `SpectraType` other than `Peaks`, or an `.slx` schema version other than 12, open-slx prints an error that names the value and exits with code 4. It writes no imzML file.
- open-slx does not read or change the `secured` field in the `.sbd` file.
- open-slx only reads the input files. It never writes to them.
- A full conversion of a large dataset is not tested. Only 300 spots of the test data were converted.
- open-slx checks only the `.sbd` `Version` (36), the `.sbd` `SpectraType` (`Peaks`), and the `.slx` schema version (12). It does not check `AxisType`.
- We did not test resampled timsTOF projects. They can store spectra another way. If open-slx fails on such a project, open an issue with `Version`, `SpectraType`, and `AxisType`.
- Reading an `.sbd` file over http or https is slow for full conversions. Copy the file first.

## Reporting problems

Report problems at <https://github.com/elparko/open-slx/issues>.

## Development

1. Run `uv sync` to install the dependencies.
2. Run `uv run pytest` to run the tests.

### Run the tests

- `uv run pytest` runs the tests that use synthetic files.
- `SLX_TEST_FILE=/path/to/project.slx uv run pytest` adds the tests that read a private `.slx` file. This variable is optional. Without it, open-slx skips these tests.
- `SLX_NETWORK=1 uv run pytest` adds the tests that download from Zenodo. They read parts of a 5.2 GB file over the network. The time depends on the speed of Zenodo (1 to 12 minutes in our runs).

## Credits

- [slx2imzml](https://github.com/CeMOS-Mannheim/slx2imzml) by CeMOS (Mannheim University of Applied Sciences), licence GPL-2.0.
  It exports `.slx` files to imzML through the proprietary SCiLS Lab Python API.
  open-slx used it as the reference for the content of an export. open-slx does not contain code from it.
- Test data: Zenodo record 21251439, "RomicsProcessor use case 2: Histology guided MALDI-MSI", licence CC BY 4.0,
  DOI 10.5281/zenodo.21251439. The project is [RomicsProcessor](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor) from PNNL (Pacific Northwest National Laboratory).
- The imzML standard: <https://ms-imaging.org/imzml/>

## Disclaimer

open-slx is not affiliated with Bruker. Bruker does not endorse it.
SCiLS is a trademark of Bruker.
The format was decoded from public data for interoperability.

## Citation

If you use open-slx in published work, please cite it. This is a request. It is not a condition of the licence.
The file [CITATION.cff](CITATION.cff) has the citation data.

## Licence

MIT. Author: Parker Smith.
