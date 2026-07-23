"""
rings.py
Aromatic ring context, altloc resolution, and residue-level helpers.
"""
import gemmi
from typing import List, Set, NamedTuple, Any, Dict
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



class _RingContext(NamedTuple):
    """Immutable context for one aromatic ring being analyzed."""
    pdb_name: str
    resolution: float
    model: gemmi.Model
    model_id: str
    chain: Any  # gemmi.Chain
    residue: Any  # gemmi.Residue
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

def select_best_altconf(structure: gemmi.Structure):
    """Select highest-occupancy altconf per residue; if tied, prefer alphabetically first (usually 'A')."""
    for model in structure:
        for chain in model:
            for residue in chain:
                altlocs = set()
                for atom in residue:
                    if atom.altloc != '\0':
                        altlocs.add(atom.altloc)
                if not altlocs:
                    continue
                if len(altlocs) == 1:
                    for atom in residue:
                        if atom.altloc != '\0':
                            atom.altloc = '\0'
                    continue
                occ_sum = {alt: 0.0 for alt in altlocs}
                occ_cnt = {alt: 0 for alt in altlocs}
                for atom in residue:
                    if atom.altloc in altlocs:
                        occ_sum[atom.altloc] += atom.occ
                        occ_cnt[atom.altloc] += 1
                avg_occ = {alt: occ_sum[alt] / occ_cnt[alt] if occ_cnt[alt] > 0 else 0.0
                           for alt in altlocs}
                best = min(altlocs, key=lambda x: (-avg_occ[x], x))
                to_remove = []
                for i in range(len(residue)):
                    atom = residue[i]
                    if atom.altloc != '\0' and atom.altloc != best:
                        to_remove.append(i)
                    elif atom.altloc == best:
                        atom.altloc = '\0'
                for i in reversed(to_remove):
                    del residue[i]

