# Xpid

[![PyPI version](https://img.shields.io/pypi/v/xpid)](https://pypi.org/project/xpid/)
[![Python 3.9+](https://img.shields.io/pypi/pyversions/xpid)](https://pypi.org/project/xpid/)
[![License](https://img.shields.io/github/license/SeanWang5868/xpid2)](https://github.com/SeanWang5868/xpid2/blob/main/LICENSE)

**Xpid** detects XH–π interactions in protein structures from PDB/mmCIF files using [Gemmi](https://gemmi.readthedocs.io/).

## Installation

```bash
pip install xpid
```

## Quick Start

Scan one structure:

```bash
xpid 1abc.cif --file-type csv
```

Scan a directory:

```bash
xpid ./structures --file-type json
```

Use a PDB code list with a local mirror:

```bash
xpid --pdb-list codes.txt --pdb-mirror /path/to/pdb/mirror
```

Prefer PDB-REDO structures when available:

```bash
xpid --pdb-list codes.txt --redo-mirror /path/to/pdb-redo/mirror --pdb-mirror /path/to/pdb/mirror
```

## Python API

```python
from xpid import detect

results = detect("structure.cif")
for hit in results:
    print(hit["pdb"], hit["pi_res"], hit["X_res"], hit["dist_X_Pi"])
```

Output column names currently remain ASCII for compatibility with existing scripts.

## Demo Notebook

`xpid_demo.ipynb` is a complete guided tutorial. It:

1. lists the required Python packages,
2. downloads PDB entry `5FJJ` beside the notebook,
3. verifies which CCP4 monomer library Xpid is using,
4. runs detection,
5. summarizes XH–π interactions with tables and plots,
6. demonstrates filtering and command-line usage.

## Geometric Criteria

Definitions: Cπ = ring centroid, **n** = ring normal, X = donor heavy atom, Xp = projection of X onto the π plane, H = hydrogen.

### Hudson System

| Parameter | Threshold |
| :--- | :--- |
| d(X–Cπ) | ≤ 4.5 Å |
| angle(X–H–**n**) | ≤ 40° |
| d(Xp–Cπ) | ≤ 1.6 Å for 5-membered rings; ≤ 2.0 Å for 6-membered rings |

### Plevin System

| Parameter | Threshold |
| :--- | :--- |
| d(X–Cπ) | ≤ donor-element cutoff |
| angle(X–H–Cπ) | ≥ 120° |
| angle(X–Cπ–**n**) | < 25° |

## Command-Line Options

### Input

| Argument | Description |
| :--- | :--- |
| `inputs` | PDB/CIF file(s) or directory path(s). |
| `--pdb-list` | Text file containing PDB codes, separated by commas or whitespace. |
| `--pdb-mirror` | Path to a local PDB mirror. |
| `--redo-mirror` | Path to a local PDB-REDO mirror, prioritized over `--pdb-mirror`. |

### Output

| Argument | Description |
| :--- | :--- |
| `--out-dir` | Override the output directory. |
| `--output-name` | Merged output filename stem. Default: `xpid_results`. |
| `--separate` | Write separate output files per structure. |
| `--file-type` | `json`, `csv`, or `parquet`. Default: `json`. |
| `-v`, `--verbose` | Include detailed geometric columns. |
| `--log` | Save a run log. |

### Processing

| Argument | Description |
| :--- | :--- |
| `--jobs N` | Number of worker processes. Default: 1. |
| `--h-mode N` | Gemmi hydrogen mode: 0=NoChange, 1=Shift, 2=Remove, 3=ReAdd, 4=ReAddButWater, 5=ReAddKnown. |
| `--model ID` | Model index to analyze, or `all`. |
| `--cone` | Enable implicit cone rescue for rotatable groups. This is the default. |
| `--no-cone` | Disable implicit cone rescue and use explicit hydrogens only. |
| `--sym-contacts` | Detect contacts across crystallographic symmetry mates. |
| `--include-water` | Include water molecules as potential donors. |
| `--max-b N` | Exclude contacts when any π-ring atom or X atom has B-factor above `N`. `0` disables this filter. |

### Filters

| Argument | Description |
| :--- | :--- |
| `--pi-res` | Limit π-acceptor residues, for example `TRP,TYR`. |
| `--donor-res` | Limit donor residues, for example `LYS,ARG`. |
| `--donor-atom` | Limit donor element symbols or exact atom names, for example `N,O,C` or `OG,NZ`. |
| `--min-occ N` | Minimum combined occupancy to report. Default: 0.0. |

## Output Data

Simple mode includes structure ID, resolution, donor/π-acceptor residue IDs, X atom, H atom, d(X–Cπ), Hudson/Plevin flags, remarks, and symmetry operation index.

Verbose mode adds secondary-structure annotations, π-center and X coordinates, angles, projection distance, sequence separation, and B-factors.

## Notes

- Same-residue donor/π-acceptor contacts are excluded.
- The automatic monomer-library download is stored in the user cache, not inside the installed package directory.
- Output column names currently remain ASCII (`pi_res`, `dist_X_Pi`) to avoid breaking existing scripts.

## Contact

Sean Wang, York Structural Biology Laboratory (YSBL), University of York

sean.wang@york.ac.uk
