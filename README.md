# Xpid

`Xpid` is a Gemmi-based tool designed to detect XH-π interactions in PDB/mmCIF files.

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

## Geometric Criteria

Definitions: $C_\pi$ (Ring Centroid), $\vec{n}$ (Ring Normal), $X$ (Donor Heavy Atom), $Xp$ (The projection of X onto the π plane), $H$ (Hydrogen).

### [Hudson System](https://doi.org/10.1021/jacs.5b08424)

$d_{X \text{--} C_\pi}$: $\le 4.5$ Å, $\angle X\text{--}H \text{--} \vec{n}$): $\le 40^\circ$. $d_{Xp \text{--} C_\pi}$: $\le 1.6$ Å $\text{for His, Trp-A}$, $\le 2.0$ Å $\text{for Phe, Trp-B, Tyr}$.

### [Plevin System](https://doi.org/10.1038/nchem.650)

$d_{X \text{--} C_\pi}$: $< 4.3$ Å, $\angle X\text{--}H \text{--} C_\pi$: $> 120^\circ$, $\angle X \text{--} C_\pi \text{--} \vec{n}$): $< 25^\circ$.


## Command Options

| Argument | Description |
| :--- | :--- |
| `inputs` | Input file (`.cif`, `.pdb`) or directory path. |
| `--out-dir` | Specify custom output directory. |
| `--separate` | Save results as separate files per PDB (Default: Merge). |
| `--file-type` | Output format: `json` (default) or `csv`. |
| `-v`, `--verbose` | Output detailed metrics (angles, coords, B-factors). |
| `--log` | Enable log file saving. |
| `--jobs N` | Number of CPU cores to use (Default: 1). |
| `--h-mode N` | Hydrogen handling mode (0=NoChange, 4=ReAddButWater). |
| `--model ID` | Model index to analyze (Default: `0`; use `all` for NMR). |
| `--pi-res` | Limit acceptor residues (e.g., `TRP,TYR`). |
| `--donor-res` | Limit donor residues (e.g., `HIS,ARG`). |
| `--donor-atom` | Limit donor element types (e.g., `N,O`). |

## Output Data

**Simple Mode (Default)**

  * PDB ID, Resolution
  * Chain, Name, ID for X-donor and $\pi$ Residues.
  * Distance ($d_{X \text{--} C_\pi}$)

**Detailed Mode (`-v`)**

  * **Includes all Simple fields plus:**
  * **Secondary structure**: Type (H/G/I/E/C) and Region IDs.
  * **Coordinates**: Flattened x, y, z for $\pi$-center and X-atom.
  * **Geometric parameters**: $\angle X\text{--}H \text{--} \vec{n}$, $\angle X\text{--}H \text{--} C_\pi$, $\angle X \text{--} C_\pi \text{--} \vec{n}$, $d_{Xp \text{--} C_\pi}$
  * **B-factors**: Average B-factor for ring atoms and X-atom.

## Dependencies

  * `gemmi`
  * `numpy`

-----

## Contact

**Sean Wang** (sean.wang@york.ac.uk)

York Structural Biology Laboratory (YSBL), University of York
