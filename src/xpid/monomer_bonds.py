"""Covalent hydrogen information from CCP4 monomer dictionaries."""
from __future__ import annotations

import logging
from typing import Dict, Optional, Set

import gemmi

from . import monlib

logger = logging.getLogger("xpid.monomer_bonds")
_BONDED_HYDROGEN_CACHE: Dict[tuple, Optional[Set[str]]] = {}
_BONDED_HEAVY_CACHE: Dict[tuple, Optional[Set[str]]] = {}


def clear_cache() -> None:
    _BONDED_HYDROGEN_CACHE.clear()
    _BONDED_HEAVY_CACHE.clear()


def get_bonded_heavy_atom_names(res_name: str, atom_name: str) -> Optional[Set[str]]:
    """Dictionary heavy neighbours, including modified component chemistry."""
    path = monlib.find_monomer_cif(res_name.upper())
    key = (res_name.upper(), atom_name.upper(), str(path) if path else None)
    if key in _BONDED_HEAVY_CACHE:
        return _BONDED_HEAVY_CACHE[key]
    neighbours = None
    if path:
        try:
            doc = gemmi.cif.read_file(str(path))
            block = doc.find_block(f"comp_{key[0]}")
            if block is None and len(doc) == 1:
                block = doc[0]
            if block is not None:
                elements = {
                    name.upper(): element.upper() for name, element in block.find([
                        "_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"])}
                if key[1] in elements:
                    neighbours = set()
                    for left, right in block.find([
                            "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2"]):
                        left, right = left.upper(), right.upper()
                        other = right if left == key[1] else left if right == key[1] else None
                        if other and elements.get(other) not in {"H", "D", None}:
                            neighbours.add(other)
        except Exception as exc:
            logger.warning("Failed to parse heavy bonds for %s/%s: %s",
                           res_name, atom_name, exc)
    _BONDED_HEAVY_CACHE[key] = neighbours
    return neighbours


def canonical_hydrogen_name(atom_name: str) -> str:
    """Normalize neutron D atom names to equivalent dictionary H names."""
    name = atom_name.strip().upper()
    return "H" + name[1:] if name.startswith("D") else name


def get_bonded_hydrogen_names(res_name: str,
                              atom_name: str) -> Optional[Set[str]]:
    """Return dictionary-defined hydrogen names.

    ``None`` means that no usable dictionary was available. An empty set
    means the dictionary explicitly defines no hydrogen bonded to this atom.
    """
    res_upper = res_name.upper()
    atom_upper = atom_name.upper()
    cif_path = monlib.find_monomer_cif(res_upper)
    cache_key = (
        res_upper, atom_upper, str(cif_path) if cif_path else None)
    if cache_key in _BONDED_HYDROGEN_CACHE:
        return _BONDED_HYDROGEN_CACHE[cache_key]
    if not cif_path:
        _BONDED_HYDROGEN_CACHE[cache_key] = None
        return None

    hydrogens: Set[str] = set()
    try:
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.find_block(f"comp_{res_upper}")
        if block is None and len(doc) == 1:
            block = doc[0]
        if block is None:
            _BONDED_HYDROGEN_CACHE[cache_key] = None
            return None

        atom_elements: Dict[str, str] = {}
        for atom_id, element in block.find([
            "_chem_comp_atom.atom_id",
            "_chem_comp_atom.type_symbol",
        ]):
            atom_elements[atom_id.strip().upper()] = element.strip().upper()

        for atom_1, atom_2 in block.find([
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
        ]):
            atom_1 = atom_1.strip().upper()
            atom_2 = atom_2.strip().upper()
            if atom_1 == atom_upper and atom_elements.get(atom_2) in {"H", "D"}:
                hydrogens.add(atom_2)
            elif atom_2 == atom_upper and atom_elements.get(atom_1) in {"H", "D"}:
                hydrogens.add(atom_1)
    except Exception as exc:
        logger.warning(
            "Failed to parse bonded hydrogens for %s/%s: %s",
            res_name, atom_name, exc)
        _BONDED_HYDROGEN_CACHE[cache_key] = None
        return None

    _BONDED_HYDROGEN_CACHE[cache_key] = hydrogens
    return hydrogens
