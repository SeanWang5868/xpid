"""
config.py
Configuration constants, atom definitions, and geometric thresholds.
"""
import gemmi
import logging
from typing import Dict, Set, Optional, List

from collections import defaultdict
from . import monlib

logger = logging.getLogger("xpid.config")

# --- Defaults ---
DEFAULT_MON_LIB_PATH = monlib.get_monomer_library_path()
DEFAULT_H_CHANGE = 4 

FALLBACK_RINGS: Dict[str, List[Set[str]]] = {
    'TRP': [
        {'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},  # 6-ring
        {'CD1', 'CD2', 'NE1', 'CG', 'CE2'}           # 5-ring
    ],
    'TYR': [{'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'}],
    'PTR': [{'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'}],
    'PHE': [{'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'}],
    'HIS': [{'CE1', 'ND1', 'NE2', 'CG', 'CD2'}],
    'BER': [
        {'C1', 'N1', 'C3', 'C6', 'C8', 'C12'},
        {'C8', 'C12', 'C13', 'C15', 'C16', 'C18'},
        {'C2', 'C4', 'C5', 'C9', 'C11', 'C14'},
    ],
    '4PO': [
        {'N2', 'C6', 'C7', 'C8', 'C9', 'C10'}
    ]
}

# --- Cache ---
AROMATIC_RINGS_CACHE: Dict[tuple, List[Set[str]]] = {}
BONDED_HYDROGENS_CACHE: Dict[tuple, Set[str]] = {}

TARGET_ELEMENTS_X = {'C', 'N', 'O', 'S'}

# --- Aromatic Ring Detection (CIF plane restraints -> DFS cycle search -> fallback) ---
def get_aromatic_rings(res_name: str) -> List[Set[str]]:
    """
    Detect all possible 5/6-membered aromatic rings by direct CIF parsing.
    Priority:
      1. Plane restraints (most accurate when defined)
      2. Aromatic-bond DFS cycle search
      3. Manual fallback only if no library/no rings found
    Fully cached by residue name and monomer-library source.
    """
    res_upper = res_name.upper()
    cache_key = (res_upper, DEFAULT_MON_LIB_PATH)
    if cache_key in AROMATIC_RINGS_CACHE:
        return AROMATIC_RINGS_CACHE[cache_key]

    rings: List[Set[str]] = []
    seen: Set[frozenset] = set()
    cif_path = monlib.find_monomer_cif(res_upper)

    if cif_path:
        try:
            doc = gemmi.cif.read_file(str(cif_path))  # Handles .gz automatically

            block_name = f"comp_{res_upper}"
            block = doc.find_block(block_name)
            if block is None and len(doc.blocks) == 1:
                block = doc[0]

            if block:
                plane_atoms: Dict[str, List[str]] = defaultdict(list)
                for plane_id_item, atom_item in block.find([
                    '_chem_comp_plane_atom.plane_id',
                    '_chem_comp_plane_atom.atom_id',
                ]):
                    atom = atom_item.strip()
                    if atom and plane_id_item:
                        plane_atoms[plane_id_item].append(atom)

                for atoms_in_plane in plane_atoms.values():
                    if len(atoms_in_plane) in (5, 6):
                        ring_set = frozenset(atoms_in_plane)
                        if ring_set not in seen:
                            rings.append(set(atoms_in_plane))
                            seen.add(ring_set)

                aromatic_edges = []
                for row in block.find([
                    '_chem_comp_bond.atom_id_1',
                    '_chem_comp_bond.atom_id_2',
                    '_chem_comp_bond.aromatic',
                ]):
                    a1, a2, arom = row
                    a1_stripped = a1.strip()
                    a2_stripped = a2.strip()
                    if a1_stripped and a2_stripped and str(arom).lower() in {'y', 'yes'}:
                        aromatic_edges.append((a1_stripped, a2_stripped))

                if aromatic_edges and not rings:  # Only if no planes found (priority)
                    graph = defaultdict(set)
                    for a, b in aromatic_edges:
                        graph[a].add(b)
                        graph[b].add(a)

                    found_rings = set()
                    max_dfs_visits = 500
                    dfs_visits = 0

                    def dfs(path: List[str], start: str):
                        nonlocal dfs_visits
                        dfs_visits += 1
                        if dfs_visits > max_dfs_visits:
                            return
                        if len(path) > 8:
                            return
                        cur = path[-1]
                        for nb in graph[cur]:
                            if nb == start and len(path) in {5, 6}:
                                found_rings.add(tuple(sorted(path)))
                            elif nb not in path:
                                dfs(path + [nb], start)

                    for atom in list(graph.keys()):
                        dfs([atom], atom)

                    for ring_tuple in found_rings:
                        ring_set = frozenset(ring_tuple)
                        if ring_set not in seen:
                            rings.append(set(ring_tuple))
                            seen.add(ring_set)

        except Exception as e:
            logger.warning(f"Failed to parse CIF for {res_name} at {cif_path}: {e}")

    if not rings and res_upper in FALLBACK_RINGS:
        rings = FALLBACK_RINGS[res_upper]

    AROMATIC_RINGS_CACHE[cache_key] = rings
    return rings


def get_bonded_hydrogens(res_name: str, atom_name: str) -> Set[str]:
    """Return dictionary-defined H/D atoms bonded to atom_name, if available."""
    res_upper = res_name.upper()
    atom_upper = atom_name.upper()
    cache_key = (res_upper, atom_upper, DEFAULT_MON_LIB_PATH)
    if cache_key in BONDED_HYDROGENS_CACHE:
        return BONDED_HYDROGENS_CACHE[cache_key]

    hydrogens: Set[str] = set()
    cif_path = monlib.find_monomer_cif(res_upper)

    if cif_path:
        try:
            doc = gemmi.cif.read_file(str(cif_path))
            block_name = f"comp_{res_upper}"
            block = doc.find_block(block_name)
            if block is None and len(doc.blocks) == 1:
                block = doc[0]

            if block:
                atom_elements: Dict[str, str] = {}
                for atom_id, elem in block.find([
                    '_chem_comp_atom.atom_id',
                    '_chem_comp_atom.type_symbol',
                ]):
                    atom_elements[atom_id.strip().upper()] = elem.strip().upper()

                for a1, a2 in block.find([
                    '_chem_comp_bond.atom_id_1',
                    '_chem_comp_bond.atom_id_2',
                ]):
                    b1 = a1.strip().upper()
                    b2 = a2.strip().upper()
                    if b1 == atom_upper and atom_elements.get(b2) in {'H', 'D'}:
                        hydrogens.add(b2)
                    elif b2 == atom_upper and atom_elements.get(b1) in {'H', 'D'}:
                        hydrogens.add(b1)
        except Exception as e:
            logger.warning(f"Failed to parse bonded hydrogens for {res_name}/{atom_name}: {e}")

    BONDED_HYDROGENS_CACHE[cache_key] = hydrogens
    return hydrogens

ROTATABLE_MAPPING: Dict[str, Dict[str, str]] = {
    'ALA': {'CB': 'CA'},
    'VAL': {'CG1': 'CB', 'CG2': 'CB'},
    'LEU': {'CD1': 'CG', 'CD2': 'CG'},
    'ILE': {'CD1': 'CG1', 'CG2': 'CB'}, 
    'MET': {'CE': 'SD'},
    'MSE': {'CE': 'SE'},
    'THR': {'CG2': 'CB', 'OG1': 'CB'},
    'SER': {'OG': 'CB'},
    'TYR': {'OH': 'CZ'},
    'CYS': {'SG': 'CB'},
    'LYS': {'NZ': 'CE'}
}

# Flexible donors: low rotational barrier (~1-2 kcal/mol), continuous 360° scan
FLEXIBLE_DONORS = {'OG', 'OG1', 'OH', 'SG'}

# Rigid rotors: high three-fold torsional barrier (~3 kcal/mol), discrete staggered states
RIGID_DONORS = {'CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE', 'NZ'} 

# Grandparent atom mapping: defines the reference plane for computing
# staggered rotamer positions of methyl / ammonium groups.
GRANDPARENT_MAPPING: Dict[str, Dict[str, str]] = {
    'ALA': {'CB': 'N'},
    'VAL': {'CG1': 'CA', 'CG2': 'CA'},
    'LEU': {'CD1': 'CB', 'CD2': 'CB'},
    'ILE': {'CD1': 'CB', 'CG2': 'CA'},
    'MET': {'CE': 'CG'},
    'MSE': {'CE': 'CG'},
    'THR': {'CG2': 'CA', 'OG1': 'CA'},
    'SER': {'OG': 'CA'},
    'TYR': {'OH': 'CE1'},
    'CYS': {'SG': 'CA'},
    'LYS': {'NZ': 'CD'}
}


BOND_LENGTHS = {
    'C': 1.09, 
    'N': 1.01, 
    'O': 0.96, 
    'S': 1.33
}
TETRAHEDRAL_ANGLE = 109.5 

DIST_SEARCH_LIMIT = 6.0
THRESHOLDS = {
    'N': 4.3,
    'O': 4.3,
    'C': 4.5,
    'S': 4.8, 
    'default': 4.5
}
MIN_COVALENT_XH = 0.5
DIST_CUTOFF_H = 1.3

# Cation-π donor atoms (positively charged groups)
CATION_DONORS = {
    ('LYS', 'NZ'),
    ('ARG', 'NH1'),
    ('ARG', 'NH2'),
    ('ARG', 'NE'),
}

PI_PI_DIST_MAX = 5.5
PI_PI_ANGLE_TSHAPED_MIN = 60.0
