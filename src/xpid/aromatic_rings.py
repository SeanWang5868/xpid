"""Aromatic-ring definitions derived from CCP4 monomer dictionaries."""
from __future__ import annotations

from collections import defaultdict
import logging
from typing import Dict, List, Set

import gemmi

from . import monlib

logger = logging.getLogger("xpid.aromatic_rings")

STANDARD_AROMATIC_RINGS: Dict[str, List[Set[str]]] = {
    "TRP": [
        {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
        {"CD1", "CD2", "NE1", "CG", "CE2"},
    ],
    "TYR": [{"CD1", "CD2", "CE1", "CE2", "CZ", "CG"}],
    "PHE": [{"CD1", "CD2", "CE1", "CE2", "CZ", "CG"}],
    "HIS": [{"CE1", "ND1", "NE2", "CG", "CD2"}],
}

_CACHE: Dict[tuple, List[Set[str]]] = {}


def clear_cache() -> None:
    """Clear parsed component-ring definitions (primarily for tests)."""
    _CACHE.clear()


def _canonical_cycle(cycle: List[str]) -> tuple[str, ...]:
    variants = []
    for ordered in (cycle, list(reversed(cycle))):
        for offset in range(len(ordered)):
            variants.append(tuple(ordered[offset:] + ordered[:offset]))
    return min(variants)


def _is_chordless_cycle(cycle: tuple[str, ...],
                        graph: Dict[str, Set[str]]) -> bool:
    cycle_size = len(cycle)
    for index, atom in enumerate(cycle):
        adjacent = {
            cycle[(index - 1) % cycle_size],
            cycle[(index + 1) % cycle_size],
        }
        if any(other in graph[atom]
               for other in cycle
               if other != atom and other not in adjacent):
            return False
    return True


def _minimal_aromatic_cycles(
        aromatic_edges: List[tuple[str, str]]) -> List[Set[str]]:
    """Extract complete chordless five- and six-membered aromatic cycles."""
    graph: Dict[str, Set[str]] = defaultdict(set)
    for atom_1, atom_2 in aromatic_edges:
        if atom_1 == atom_2:
            continue
        graph[atom_1].add(atom_2)
        graph[atom_2].add(atom_1)

    cycles: Set[tuple[str, ...]] = set()

    def visit(start: str, path: List[str]) -> None:
        if len(path) > 6:
            return
        current = path[-1]
        for neighbor in sorted(graph[current]):
            if neighbor == start:
                if len(path) in {5, 6}:
                    canonical = _canonical_cycle(path)
                    if _is_chordless_cycle(canonical, graph):
                        cycles.add(canonical)
                continue
            if neighbor in path or len(path) == 6:
                continue
            visit(start, path + [neighbor])

    for atom in sorted(graph):
        visit(atom, [atom])

    ordered_cycles = sorted(
        cycles, key=lambda cycle: (len(cycle), tuple(sorted(cycle))))
    return [set(cycle) for cycle in ordered_cycles]


def _order_standard_rings(res_name: str,
                          rings: List[Set[str]]) -> List[Set[str]]:
    fallback = STANDARD_AROMATIC_RINGS.get(res_name)
    if not fallback:
        return rings
    rank = {
        frozenset(atom_names): index
        for index, atom_names in enumerate(fallback)
    }
    return sorted(
        rings,
        key=lambda atom_names: (
            rank.get(frozenset(atom_names), len(rank)),
            len(atom_names),
            tuple(sorted(atom_names)),
        ),
    )


def _standard_rings_are_complete(res_name: str,
                                 rings: List[Set[str]]) -> bool:
    fallback = STANDARD_AROMATIC_RINGS.get(res_name)
    if fallback is None:
        return True
    observed = {frozenset(atom_names) for atom_names in rings}
    expected = {frozenset(atom_names) for atom_names in fallback}
    return observed == expected


def get_aromatic_rings(res_name: str) -> List[Set[str]]:
    """Return dictionary-backed five- and six-membered aromatic rings.

    Only complete chordless cycles made entirely from
    ``_chem_comp_bond.aromatic = y`` bonds are accepted. PHE, TYR, TRP and
    HIS use built-in definitions when their dictionary is missing or unusable.
    """
    res_upper = res_name.upper()
    cif_path = monlib.find_monomer_cif(res_upper)
    cache_key = (res_upper, str(cif_path) if cif_path else None)
    if cache_key in _CACHE:
        return _CACHE[cache_key]

    rings: List[Set[str]] = []
    if cif_path:
        try:
            doc = gemmi.cif.read_file(str(cif_path))
            block = doc.find_block(f"comp_{res_upper}")
            if block is None and len(doc.blocks) == 1:
                block = doc[0]

            if block is not None:
                aromatic_edges = []
                for atom_1, atom_2, aromatic in block.find([
                    "_chem_comp_bond.atom_id_1",
                    "_chem_comp_bond.atom_id_2",
                    "_chem_comp_bond.aromatic",
                ]):
                    atom_1 = atom_1.strip()
                    atom_2 = atom_2.strip()
                    if (atom_1 and atom_2 and
                            str(aromatic).lower() in {"y", "yes"}):
                        aromatic_edges.append((atom_1, atom_2))
                rings = _minimal_aromatic_cycles(aromatic_edges)
        except Exception as exc:
            logger.warning(
                "Failed to parse aromatic bonds for %s at %s: %s",
                res_name, cif_path, exc)

    if (res_upper in STANDARD_AROMATIC_RINGS and
            not _standard_rings_are_complete(res_upper, rings)):
        rings = [
            set(atom_names)
            for atom_names in STANDARD_AROMATIC_RINGS[res_upper]
        ]
    else:
        rings = _order_standard_rings(res_upper, rings)

    _CACHE[cache_key] = rings
    return rings
