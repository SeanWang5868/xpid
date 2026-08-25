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

Scan one protein structure:

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

Analysis failures raise `XPIDError`; use `detect(..., on_error="empty")` only
when legacy empty-on-error behavior is explicitly desired.

Output column names currently remain ASCII for compatibility with existing scripts.

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
| `--include-coordinates` | Include absolute π-center, X, H coordinates, canonical π normal, and `X_side_of_pi`. |
| `--sasa` | Include solvent accessible surface area columns (`pi_sasa_avg`, `X_sasa`, `H_sasa`). |
| `--cooperativity` | Annotate hits with cooperativity metrics (donor counts per ring face). Enabled by default; use `--no-cooperativity` to disable. |
| `--provenance` | Write a `_metadata.json` companion file recording all run parameters for reproducibility. |
| `--log` | Save a run log. |

### Processing

| Argument | Description |
| :--- | :--- |
| `--jobs N` | Number of worker processes. Default: 1. |
| `--h-mode N` | Gemmi hydrogen mode: 0=NoChange, 1=Shift, 2=Remove, 3=ReAdd, 4=ReAddButWater, 5=ReAddKnown. |
| `--model ID` | Model index to analyze, or `all`. |
| `--no-cone` | Disable binary cone detection and use explicit hydrogen positions only. Rotatable Ser, Thr, Tyr, Cys and methyl groups use cone detection by default. Cation–π groups such as Lys/Arg are excluded by design. |
| `--include-p-slab`, `--p-slab` | Include the optional P-slab system, P-slab output columns, and P-slab summary counts. |
| `--xh-candidates` | Export all explicit X-H bonds passing Hudson/Plevin X-position filters, including direction-failed candidates. Cone virtual H is ignored in this mode. |
| `--sym-contacts` | Detect contacts across crystallographic symmetry mates. |
| `--include-water` | Include water molecules as potential donors. |
| `--max-b N` | Exclude contacts when any π-ring atom or X atom has B-factor above `N`. `0` disables this filter. |

### Filters

| Argument | Description |
| :--- | :--- |
| `--pi-res` | Limit π-acceptor residues, for example `TRP,TYR`. |
| `--donor-res` | Limit donor residues, for example `SER,THR`. Cation–π Lys/Arg groups are excluded by design. |
| `--donor-atom` | Limit donor element symbols or exact atom names, for example `N,O,C` or `OG,NZ`. |
| `--residue-pair SEL1 SEL2` | Restrict detection to XH–π interactions between two residue selections, for example `--residue-pair //A/12 //A/18`. Either selected residue can be the π acceptor; all donor atoms on the opposite residue are considered. |
| `--min-occ N` | Minimum combined occupancy to report. Default: 0.0. |

## Output Data

Simple mode includes structure ID, resolution, donor/π-acceptor residue IDs, X
atom, H atom, H source, the `is_hudson` and `is_plevin` labels, the main
Hudson/Plevin geometry columns, TRP 5-ring flag, T-shaped π-π flag, symmetry
operation index, and the full `symmetry_code` including lattice translation.
With `--include-p-slab`, simple mode also includes
`is_p_slab`, `P_radius`, `P_slab_half_thickness`, `h_proj_dist`, `H_ray_t`,
and related ray-entry geometry.

With `--xh-candidates`, simple mode additionally includes diagnostic columns
such as `is_xh_candidate`, `is_hudson_spatial`, `is_plevin_spatial`,
`hudson_direction_ok`, `plevin_direction_ok`, `xh_centroid_cos`,
`xh_lateral_inward_score`, `h_proj_dist`, `H_ray_t`, `h_plane_proj_dist`, and
`H_plane_t`. Unless `--include-p-slab` is also used, candidate output does not
include an `is_p_slab` label.

With `--include-coordinates`, simple mode additionally includes
`pi_center_x/y/z`, `pi_normal_x/y/z`, `X_xyz_x/y/z`, `H_xyz_x/y/z`, and
`X_side_of_pi`. The π normal is canonicalized so that its largest absolute
component is positive; `X_side_of_pi` is `1` when X lies on the positive side of
that normal, `-1` on the negative side, and `0` when X is effectively in the
π plane.


With `--cooperativity` (default), simple mode additionally includes
`coop_same_face_donors`, `coop_opp_face_donors`, `coop_total_donors`, and
`coop_bi_face`.  Verbose mode additionally includes conventional H-bond
competition columns (`hbond_competing`, `hbond_acceptor_*`, `hbond_HA_dist`,
`hbond_DHA_angle`, `hbond_vs_xhpi_score`) and `hbond_relation`, which
distinguishes same-hydrogen, same-conformer-other-hydrogen, and
alternative-conformer H-bond context.

Verbose mode adds secondary-structure annotations, P center and X coordinates,
sequence separation, and B-factors.

Every run also writes `<output-name>_diagnostics.json`.  It records the input
source and canonical structure ID, hydrogen-preparation status, missing monomer
components, incomplete aromatic rings, Cone groups with a missing parent atom,
ambiguous or altloc-incompatible Cone parents, and any per-structure failure.
Cone parent selection is conformer-aware: a labelled donor uses the matching
parent altloc or one shared blank parent, while unresolved or chemically
incompatible parent assignments are skipped and audited rather than selected
by atom order.  If Gemmi identifies one specific residue whose
atom names or topology are incompatible with its monomer definition, Xpid
marks hydrogen preparation as `partial`, retains that residue's experimental
heavy atoms in the steric environment, but disables its donor capability.  A
chemically complete aromatic ring in the same residue can still serve as the
pi system.  Hydrogens on adjacent polymer residues are also suppressed when a
temporary chain break could create false terminal donors.  The skipped and
protected residue identities are recorded in the diagnostics JSON; the run
audit prints skipped identities and aggregate counts.  Failures that cannot be
localized safely still fail the whole structure; they are never reported as a
zero-interaction result.

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

- Same-residue donor/π-acceptor contacts within the asymmetric unit are
  excluded. A crystallographic symmetry copy is treated as a distinct molecule.
- The automatic monomer-library download is stored in the user cache, not inside the installed package directory.
- Output column names currently remain ASCII (`pi_res`, `dist_X_Pi`) to avoid breaking existing scripts. `dist_X_Pi` now means `d(X, P plane)`, not distance from X to the ring centroid.
- Rotatable donor groups (Ser, Thr, Tyr, Cys and standard amino-acid methyl
  groups) use binary
  cone detection by default. The detector keeps the observed heavy atoms fixed
  and asks whether a chemically valid hydrogen conformer satisfies Hudson or
  Plevin. OH/SH groups are sampled as one hydrogen through 360 degrees. CH3
  groups are sampled as complete three-hydrogen conformers through their
  120-degree rotational period. Lys/Arg cationic donors remain excluded from
  production XH-pi detection. Group-specific CCP4 nuclear-position bond
  lengths and angles are used. Cys uses S-H = 1.338 Angstrom and
  C-beta-S-H = 97.543 degrees from the CCP4 CYS dictionary.
- Cone conformers are rejected for severe non-bonded clashes. A chemically
  valid H-bond contact is not treated as a clash. Potential conventional
  H-bonds are descriptive context and do not remove other sterically valid
  conformers from the default binary detector. This avoids treating the mere
  existence of a hypothetical H-bond rotamer as proof that the donor is
  committed to that direction.
- Every reported cone row uses one self-consistent conformer: its H
  coordinates, Hudson/Plevin labels, and geometric columns all refer to the
  same virtual hydrogen. Use `--no-cone` to rely on explicit hydrogens only
  (appropriate for neutron diffraction or structures whose H positions are
  otherwise known to be experimental).

## Contact

Sean Wang, York Structural Biology Laboratory (YSBL), University of York

sean.wang@york.ac.uk
