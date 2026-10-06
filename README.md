# Xpid

[![PyPI version](https://img.shields.io/pypi/v/xpid)](https://pypi.org/project/xpid/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22875176.svg)](https://doi.org/10.5281/zenodo.22875176)
[![Python 3.9+](https://img.shields.io/pypi/pyversions/xpid)](https://pypi.org/project/xpid/)
[![License](https://img.shields.io/github/license/SeanWang5868/xpid)](https://github.com/SeanWang5868/xpid/blob/main/LICENSE)

**Xpid** detects XH–π interactions in macromolecular structures. It reads
PDB and mmCIF coordinates and uses [Gemmi](https://gemmi.readthedocs.io/) with
the CCP4 monomer library for chemical perception and hydrogen placement.

Xpid:

- identifies aromatic rings from CCP4 monomer lib;
- evaluates rigid X–H vectors and chemically supported rotatable X–H groups;
- applies steric filtering and the Hudson/Plevin geometric criteria;
- can include crystallographic symmetry contacts, SASA and cooperativity;

## Installation

```bash
pip install xpid
```

Parquet output requires the optional dependencies:

```bash
pip install "xpid[parquet]"
```

## Quick start

Analyse one structure:

```bash
xpid model.cif --out-dir xpid_results --file-type csv
```

Analyse a directory in parallel:

```bash
xpid structures/ --out-dir xpid_results --file-type csv --jobs 8
```

Analyse a list of PDB entries from a local mirror:

```bash
xpid \
  --pdb-list pdb_ids.txt \
  --pdb-mirror /path/to/pdb/mmCIF \
  --out-dir xpid_results \
  --file-type csv \
  --jobs 8 \
  --provenance
```

Add `--redo-mirror /path/to/pdb-redo` to prefer PDB-REDO models, or
`--sym-contacts` to include crystallographic symmetry mates.

Run `xpid --help` for the complete command-line reference.

## Detection model

By default, Xpid reports a contact when either the Hudson or Plevin definition
is satisfied. The X-to-ring-centroid distance limits are 4.3 Å for N/O,
4.5 Å for C and 4.8 Å for S donors.

| Definition | Directional requirement |
| :--- | :--- |
| Hudson | X projects within the ring region and the X–H direction is within 40° of the ring-facing direction. |
| Plevin | `XPCN < 25°` and `XH–π ≥ 120°`. |

Important defaults:

- Hydrogen mode 4 (`ReAddButWater`) is used unless `--h-mode` is specified.
- Supported rotatable hydroxyl, thiol and methyl groups use Cone detection;
  `--no-cone` uses explicit hydrogen positions only.
- Cone conformers with severe non-bonded clashes are rejected.
- Conventional hydrogen bonds are reported as context and do not veto an
  otherwise valid XH–π contact.
- Lys/Arg cationic donors are excluded because Xpid does not report cation–π
  interactions.
- Water donors and crystallographic symmetry contacts are excluded unless
  explicitly requested. Water remains in the environment for steric filtering,
  conventional hydrogen-bond context and SASA.
- Xpid is a geometric detector; it does not calculate interaction energies.

When a component dictionary is available, only complete five- or six-membered
rings formed by `aromatic=y` bonds are treated as aromatic. If no dictionary is
available, built-in ring definitions are used only for PHE, TYR, TRP and HIS.

## Output

The default merged output is `xpid_results.json`; CSV and Parquet are selected
with `--file-type`. Each reported row identifies the aromatic ring, donor X and
hydrogen, and includes the main geometry together with `is_hudson` and
`is_plevin` labels.

`H_source` records `input_model` (input H/D coordinates), `added` (hydrogens
generated during preparation), `adjusted` (input H/D moved during preparation)
or `cone_virtual`. Input coordinates alone do not establish that a hydrogen
was experimentally observed. Virtual H inherit the active parent/X conformer.
Cone requires a terminal, covalently compatible group with available X–H
valence; declared cross-links and coordination are retained during this check.
Atoms below 0.10 occupancy are excluded from donor-H, Cone-parent and H-bond
environment checks. H-bond context also respects `--min-occ`, `--max-b` and the
detected alternate state, and identifies the acceptor's crystallographic image.
Cooperativity is computed separately for each ring and ring alternate state,
independently of coordinate export; an unknown face remains unassigned.

Useful optional output includes:

- `--verbose` for additional geometry and hydrogen-bond context;
- `--include-coordinates` for ring-centre, X and H coordinates;
- `--sasa` for solvent-accessible surface area;
- `--no-cooperativity` to disable the default cooperativity annotation;
- `--provenance` for a companion file containing the run configuration.

Every run writes `xpid_results_diagnostics.json`, including input resolution,
hydrogen-preparation status, missing monomer components and skipped structures
or residues. Localised chemistry problems are reported rather than silently
converted into zero-interaction results.

## Python API

```python
from xpid import detect

for hit in detect("model.cif"):
    print(hit["pi_res"], hit["X_res"], hit["dist_X_Pi"])
```

## Citation

If you use Xpid, please cite the software version used in your analysis:

> Wang, S. *et al.* (2026). *Xpid: A Detector for XH–π Interactions in
> Macromolecular Structures* (Version 2.1.11) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.22875176

This concept DOI represents all versions; use the version-specific DOI when
one is available for the release being cited.

## License and contact

Xpid is distributed under the [MIT License](LICENSE).

Shuai Wang, York Structural Biology Laboratory, University of York

[sean.wang@york.ac.uk](mailto:sean.wang@york.ac.uk)
