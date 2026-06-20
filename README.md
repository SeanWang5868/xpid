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

By default, Xpid reports the Hudson/Plevin result set. A row is reported when
Hudson or Plevin is positive. The finite P-slab system is still available for
method-development comparisons, but it is opt-in rather than part of the
default output.

For directionality-background analysis, Xpid can also export X-H candidate rows
with `--xh-candidates`. In this mode, Xpid first applies only the Hudson/Plevin
X-position criteria around each π ring, then records every bonded X-H vector and
assigns the Hudson/Plevin labels post hoc. This mode is intended for downstream
annotation and null/background comparisons; it is not a new positive XH–π
definition.

Definitions: X = donor heavy atom, H = hydrogen, Xp = orthogonal projection of X
onto the aromatic plane, Hp = entry point where the X->H ray first intersects
the finite P slab, and P = finite aromatic π slab.

| System | Parameters | Threshold |
| :--- | :--- | :--- |
| Hudson | `dist_X_centroid`, `proj_dist`, `theta` | d(X, centroid) ≤ element cutoff; Xp inside P radius; `theta` ≤ 40° |
| Plevin | `dist_X_centroid`, `angle_XPCN`, `angle_XH_Pi` | d(X, centroid) < element cutoff; `angle_XPCN` < 25°; `angle_XH_Pi` ≥ 120° |
| P-slab optional | `dist_X_Pi`, `proj_dist`, `h_proj_dist`, `H_ray_t` | Enable with `--include-p-slab`; d(X, P plane) ≤ element cutoff; Xp inside P radius; X->H ray enters the finite P slab inside P |

Element cutoffs are ≤ 4.3 Å for N/O, ≤ 4.5 Å for C, and ≤ 4.8 Å for S.
P radius is 1.6 Å for 5-membered rings and 2.0 Å for 6-membered rings. When
the optional P-slab system is enabled, the P-slab half-thickness is 0.5 Å above
and below the aromatic plane.

The X->H ray is directional: H must lie between X and the finite P slab. This
prevents an infinite X-H line from being counted when the hydrogen points away
from P, while avoiding hard-boundary failures caused by treating the π region as
a zero-thickness disk.

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
| `--cone` | Enable implicit cone rescue for rotatable groups. |
| `--no-cone` | Disable implicit cone rescue and use explicit hydrogens only. This is the default. |
| `--include-p-slab`, `--p-slab` | Include the optional P-slab system, P-slab output columns, and P-slab summary counts. |
| `--xh-candidates` | Export all explicit X-H bonds passing Hudson/Plevin X-position filters, including direction-failed candidates. Cone virtual H is ignored in this mode. |
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

Simple mode includes structure ID, resolution, donor/π-acceptor residue IDs, X
atom, H atom, H source, the `is_hudson` and `is_plevin` labels, the main
Hudson/Plevin geometry columns, TRP 5-ring flag, T-shaped π-π flag, and
symmetry operation index. With `--include-p-slab`, simple mode also includes
`is_p_slab`, `P_radius`, `P_slab_half_thickness`, `h_proj_dist`, `H_ray_t`,
and related ray-entry geometry.

With `--xh-candidates`, simple mode additionally includes diagnostic columns
such as `is_xh_candidate`, `is_hudson_spatial`, `is_plevin_spatial`,
`hudson_direction_ok`, `plevin_direction_ok`, `xh_centroid_cos`,
`xh_lateral_inward_score`, `h_proj_dist`, `H_ray_t`, `h_plane_proj_dist`, and
`H_plane_t`. Unless `--include-p-slab` is also used, candidate output does not
include an `is_p_slab` label.

Verbose mode adds secondary-structure annotations, P center and X coordinates,
sequence separation, and B-factors.

Every default reported row satisfies Hudson or Plevin. The label columns are
integer 1/0 flags and can be used for downstream filtering. The command-line
summary prints Hudson-positive, Plevin-positive, Hudson/Plevin union, and
Hudson+Plevin overlap counts. With `--include-p-slab`, it additionally prints
P-slab-positive and all-three overlap counts.

In `--xh-candidates` mode, the total row count is the number of exported X-H
background candidates. `is_hudson` and `is_plevin` remain strict positive labels;
rows with both labels set to 0 are spatially nearby X-H vectors that failed the
directional filters.

## Notes

- Same-residue donor/π-acceptor contacts are excluded.
- The automatic monomer-library download is stored in the user cache, not inside the installed package directory.
- Output column names currently remain ASCII (`pi_res`, `dist_X_Pi`) to avoid breaking existing scripts. `dist_X_Pi` now means `d(X, P plane)`, not distance from X to the ring centroid.

## Contact

Sean Wang, York Structural Biology Laboratory (YSBL), University of York

sean.wang@york.ac.uk
