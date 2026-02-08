# Xpid

`Xpid` is a Gemmi-based tool to detect XH-π interactions in PDB/mmCIF files.

## Installation

Requires Python 3.9+.

```bash
pip install xpid
```

## Configuration

The detection of XH–π interactions depends on the positions of hydrogen atoms. To add hydrogens to the structure before detection, the path to a monomer library must be specified (e.g., the [CCP4 Monomer Library](https://github.com/MonomerLibrary/monomers)).


```bash
xpid --set-mon-lib /Users/abc123/monomers
```

## Quick Start

Scans a directory or PDB/mmCIF file and save results into a JSON file.

```bash
xpid ./data
```

> **Output**: `./data/xpid_output/xpid_results.json`

## Geometric Criteria (An XH–π interaction is recorded when at least one of the following set criteria is met.)

**Definitions:** * **C<sub>π</sub>**: Ring Centroid
* **n**: π-Ring Normal Vector
* **X**: Donor Atom
* **X<sub>proj</sub>**: The projection of Donor Atom onto the π plane
* **H**: Hydrogen

### [Hudson System](https://doi.org/10.1021/jacs.5b08424)

* **Distance (X &#x2d; C<sub>π</sub>)**: &le; 4.5 &#197;
* **Angle (X &#x2d; H, n)**: &le; 40&deg;
* **Distance (X<sub>proj</sub> &#x2d; C<sub>π</sub>)**: 
    * &le; 1.6 &#197; (for His, Trp-A)
    * &le; 2.0 &#197; (for Phe, Trp-B, Tyr)

### [Plevin System](https://doi.org/10.1038/nchem.650)

* **Distance (X &#x2d; C<sub>π</sub>)**: < 4.3 &#197;
* **Angle (X &#x2d; H &#x2d; C<sub>π</sub>)**: > 120&deg;
* **Angle (X &#x2d; C<sub>π</sub>, n)**: < 25&deg;


## Command Options

| Argument | Description |
| :--- | :--- |
| `inputs` | Input file (`.cif`, `.pdb`) or directory path. |
| `--out-dir` | Specify custom output directory. |
| `--separate` | Save results as separate files per PDB (Default: Merge). |
| `--file-type` | Output format: `json` (default) or `csv`. |
| `-v`, `--verbose` | Output detailed metrics (angles, coords, B-factors). |
| `--log` | Enable log file saving. |
| `--h-mode N` | Hydrogen handling mode (0=NoChange, 4=ReAddButWater). |
| `--jobs N` | Number of CPU cores to use (Default: 1). |
| `--model ID` | Model index to analyze (Default: `0`; use `all` for NMR). |
| `--mon-lib` | Custom Monomer Library path ("/Users/abc123/monomers"). |
| `--set-mon-lib` | Set default Monomer Library ("/Users/abc123/monomers"). |
| `--show-mon-lib-config` | Show current Monomer Library status. |
| `--pi-res` | Limit acceptor residues (e.g., `TRP,TYR`). |
| `--donor-res` | Limit donor residues (e.g., `HIS,ARG`). |
| `--donor-atom` | Limit donor element types (e.g., `N,O`). |

## Output Data

### Simple Mode (Default)
* **Structure Info**: PDB ID, Resolution
* **Residue Info**: Chain, Name, ID for X-donor and &pi; Residues
* **Metric**: Distance (X &#x2d; C<sub>&pi;</sub>)

### Detailed Mode (-v)
Includes all *Simple Mode* fields plus:

* **Secondary Structure**: Type (H/G/I/E/C) and Region IDs
* **Coordinates**: Flattened x, y, z for &pi;-center and X-atom
* **Geometric Parameters**:
    * Angle (X &#x2d; H, **n**)
    * Angle (X &#x2d; H &#x2d; C<sub>&pi;</sub>)
    * Angle (X &#x2d; C<sub>&pi;</sub>, **n**)
    * Distance (X<sub>proj</sub> &#x2d; C<sub>&pi;</sub>)
* **B-factors**: Average B-factor for ring atoms and X-atom

## Dependencies

  * `gemmi`
  * `numpy`

-----

## Contact

**Sean Wang** (sean.wang@york.ac.uk)

York Structural Biology Laboratory (YSBL), University of York
