"""
sasa.py
Solvent Accessible Surface Area calculation using the Shrake-Rupley algorithm.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import gemmi
import numpy as np

# ---------------------------------------------------------------------------
# van der Waals radii (Å) — Bondi / Chothia values
# ---------------------------------------------------------------------------

VDW_RADII: Dict[str, float] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "H": 1.20,
    "D": 1.20,
}

DEFAULT_PROBE_RADIUS = 1.4
DEFAULT_N_POINTS = 256

# ---------------------------------------------------------------------------
# Golden-spiral unit-sphere point set (precomputed once per n_points)
# ---------------------------------------------------------------------------

_SPHERE_POINT_CACHE: Dict[int, np.ndarray] = {}


def _unit_sphere_points(n: int) -> np.ndarray:
    """Generate *n* approximately-even points on the unit sphere (golden spiral)."""
    if n in _SPHERE_POINT_CACHE:
        return _SPHERE_POINT_CACHE[n]

    indices = np.arange(n)
    y = 1.0 - (indices / (n - 1)) * 2.0
    radius_at_y = np.sqrt(1.0 - y * y)
    theta = math.pi * (1.0 + math.sqrt(5.0)) * indices

    x = np.cos(theta) * radius_at_y
    z = np.sin(theta) * radius_at_y

    points = np.column_stack((x, y, z))
    _SPHERE_POINT_CACHE[n] = points
    return points


def _atom_radius(element: str) -> float:
    return VDW_RADII.get(element.upper(), 1.70)


def compute_sasa(
    structure: gemmi.Structure,
    probe_radius: float = DEFAULT_PROBE_RADIUS,
    n_points: int = DEFAULT_N_POINTS,
    model_index: int = 0,
) -> Dict[Tuple[int, str, str, int], float]:
    """Compute per-atom SASA for one model of a structure.

    Returns a dict mapping ``(model_index, chain_name, residue_seqid, atom_index)``
    to the SASA value in Å².  Atom index is the position within the residue.
    """
    if len(structure) == 0 or model_index >= len(structure):
        return {}

    model = structure[model_index]
    sphere_points = _unit_sphere_points(n_points)

    atoms: List[Tuple[float, float, float, float, Tuple]] = []

    for chain in model:
        for residue in chain:
            for i, atom in enumerate(residue):
                elem = atom.element.name.upper() if atom.element.name else "C"
                r_vdw = _atom_radius(elem)
                R = r_vdw + probe_radius
                key = (model_index, chain.name, str(residue.seqid).strip(), i)
                atoms.append((atom.pos.x, atom.pos.y, atom.pos.z, R, key))

    if not atoms:
        return {}

    n_atoms = len(atoms)
    coords = np.array([(a[0], a[1], a[2]) for a in atoms])
    radii = np.array([a[3] for a in atoms])
    keys = [a[4] for a in atoms]

    area_per_point = 4.0 * math.pi / n_points
    sasa_values: Dict[Tuple, float] = {}

    for i in range(n_atoms):
        R_i = radii[i]
        test_points = coords[i] + sphere_points * R_i
        exposed = np.ones(n_points, dtype=bool)

        centre_dists = np.linalg.norm(coords - coords[i], axis=1)
        R_sum = R_i + radii
        neighbour_mask = (centre_dists < R_sum) & (np.arange(n_atoms) != i)
        neighbour_indices = np.where(neighbour_mask)[0]

        for j in neighbour_indices:
            R_j = radii[j]
            dists = np.linalg.norm(test_points - coords[j], axis=1)
            exposed &= (dists >= R_j)

        n_exposed = int(np.sum(exposed))
        sasa_values[keys[i]] = round(n_exposed * area_per_point * (R_i ** 2), 3)

    return sasa_values


def atom_sasa(
    sasa_map: Dict[Tuple, float],
    model_index: int,
    chain_name: str,
    residue_seqid: str,
    atom_index: int,
) -> Optional[float]:
    """Look up a single atom's SASA value from a precomputed map."""
    return sasa_map.get((model_index, chain_name, residue_seqid, atom_index))


def residue_atom_indices(residue: gemmi.Residue) -> List[int]:
    """Return the 0-based indices of heavy atoms in a residue."""
    return [
        i for i, atom in enumerate(residue)
        if atom.element.name.upper() not in ("H", "D", "")
    ]


def average_ring_sasa(
    sasa_map: Dict[Tuple, float],
    model_index: int,
    chain_name: str,
    residue: gemmi.Residue,
    ring_atom_indices: List[int],
) -> float:
    """Compute average SASA over a set of ring atoms."""
    values = []
    for idx in ring_atom_indices:
        v = atom_sasa(sasa_map, model_index, chain_name, str(residue.seqid).strip(), idx)
        if v is not None:
            values.append(v)
    if not values:
        return 0.0
    return round(sum(values) / len(values), 3)
