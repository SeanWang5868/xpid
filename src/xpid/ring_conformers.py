"""Aromatic-ring coordinate context and altloc resolution."""
from dataclasses import dataclass
import gemmi
from typing import List, Set, Any, Dict
import numpy as np

def _altloc(atom: gemmi.Atom) -> str:
    return "" if atom.altloc in ("", "\0") else atom.altloc


def _sorted_altlocs(altlocs: Set[str]) -> List[str]:
    return sorted(altlocs, key=lambda alt: (alt != "", alt))


def _atom_variants_for_names(residue, atom_names: Set[str]) -> List[tuple]:
    """Return complete atom-name sets for each compatible altloc state."""
    atoms_by_name = {name: [] for name in atom_names}
    for atom in residue:
        if atom.name in atoms_by_name:
            atoms_by_name[atom.name].append(atom)

    if any(not atoms for atoms in atoms_by_name.values()):
        return []

    altlocs = {""}
    for atoms in atoms_by_name.values():
        altlocs.update(_altloc(atom) for atom in atoms if _altloc(atom))

    variants = []
    seen = set()
    for alt in _sorted_altlocs(altlocs):
        selected = []
        for name in sorted(atom_names):
            atoms = atoms_by_name[name]
            exact = [atom for atom in atoms if _altloc(atom) == alt]
            blank = [atom for atom in atoms if _altloc(atom) == ""]
            if alt:
                atom = exact[0] if exact else (blank[0] if blank else None)
            else:
                atom = blank[0] if blank else None
            if atom is None:
                selected = []
                break
            selected.append(atom)

        if selected:
            signature = tuple(id(atom) for atom in selected)
            if signature not in seen:
                seen.add(signature)
                variants.append((alt, selected))
    return variants



@dataclass(frozen=True)
class RingContext:
    """Immutable context for one aromatic ring being analyzed."""
    pdb_name: str
    resolution: float
    model: gemmi.Model
    model_id: str
    chain: Any
    residue: Any
    ns: gemmi.NeighborSearch
    ss_index: Dict
    pi_center_arr: np.ndarray
    pi_normal: np.ndarray
    pi_b_mean: float
    pi_alt: str
    mode: str
    ring_size: int
    min_occ: float
    avg_pi_occ: float
    p_radius: float
    p_slab_half_thickness: float
