# Xpid

`xpid` is a Gemmi-based tool designed to detect XH-π interactions in PDB/mmCIF files.

## Installation

Requires Python 3.9+.

```bash
git clone https://github.com/SeanWang5868/xpid
cd xpid
pip install .
```

## Configuration

The detection of XH-π interactions depends on the position of H atoms. In order to add H to the structure before detecting, the path to the monomer library (e.g. CCP4 monomer library) needs to be specified.

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

Definitions: $C_\pi$ (Ring Centroid), $\vec{n}$ (Ring Normal), $X$ (Donor Heavy Atom), $Xp$ (The projection of X onto the \pi plane), $H$ (Hydrogen).

### [Hudson System](https://doi.org/10.1021/jacs.5b08424)

$d_{\text{X --} \, C_\pi}$: $\le 4.5\,\text{Å}$, 
$\angle\, \text{X -- H --} \vec{n}$: $\le 40^\circ$.  
$d_{\text{Xp --} \, C_\pi}$: $\le 1.6\,\text{Å}$ for \text{His} and \text{Trp\_A}, $\le 2.0\,\text{Å}$ for \text{Phe}, \text{Trp\_B}, \text{Tyr}

### [Plevin System](https://doi.org/10.1038/nchem.650)

$d_{\text{X --} \, C_\pi}$: $< 4.3\,\text{Å}$,  
$\angle\, \text{X -- H --} C_\pi$: $> 120^\circ$,  
$\angle\, \text{X --} C_\pi \text{ --} \vec{n}$: $< 25^\circ$.


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
