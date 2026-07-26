"""Chemically conservative conventional hydrogen-bond acceptor typing.

This module does not identify aromatic π acceptors.  Aromatic π-system
recognition lives in :mod:`xpid.aromatic_rings`.
"""
from __future__ import annotations

from typing import Optional

import gemmi


BACKBONE_ACCEPTORS = {"O", "OXT"}

SIDECHAIN_ACCEPTORS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "ASN": {"OD1"},
    "GLN": {"OE1"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "TYR": {"OH"},
    "MET": {"SD"},
    "MSE": {"SE"},
}

# Neutral hydroxyl oxygens can both donate and accept hydrogen bonds.  Their
# explicit H therefore does not make them non-acceptors.
HYDROXYL_ACCEPTORS = {
    ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"),
}

NON_ACCEPTOR_NITROGEN = {
    ("LYS", "NZ"),
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("ASN", "ND2"), ("GLN", "NE2"),
}


def _has_bonded_hydrogen(residue: gemmi.Residue, atom: gemmi.Atom) -> bool:
    """Use explicit local H positions as conservative protonation evidence."""
    for candidate in residue:
        if candidate.element.name.upper() not in {"H", "D"}:
            continue
        if atom.pos.dist(candidate.pos) <= 1.45:
            return True
    return False


def is_hbond_acceptor(residue: gemmi.Residue, atom: gemmi.Atom) -> bool:
    """Return whether *atom* is a defensible conventional H-bond acceptor.

    This deliberately favours false negatives over allowing chemically
    impossible atoms to lock a rotatable donor away from a pi interaction.
    """
    name = atom.name.upper()
    res_name = residue.name.upper()
    element = atom.element.name.upper()

    if name in BACKBONE_ACCEPTORS and element == "O":
        return True
    if (res_name, name) in HYDROXYL_ACCEPTORS:
        return True
    if name in SIDECHAIN_ACCEPTORS.get(res_name, set()):
        return not _has_bonded_hydrogen(residue, atom)
    if (res_name, name) in NON_ACCEPTOR_NITROGEN:
        return False

    if res_name == "HIS" and name in {"ND1", "NE2"}:
        return not _has_bonded_hydrogen(residue, atom)

    if residue.is_water() and element == "O":
        return True

    # Unknown O/S atoms are accepted only when not explicitly protonated.
    # Unknown nitrogens require valence-aware monomer typing and are therefore
    # not allowed to impose a deterministic cone constraint.
    if element in {"O", "S"}:
        return not _has_bonded_hydrogen(residue, atom)
    return False


def altlocs_compatible(left: Optional[str], right: Optional[str]) -> bool:
    left = "" if left in (None, "", "\0") else left
    right = "" if right in (None, "", "\0") else right
    return not left or not right or left == right
