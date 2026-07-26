"""Chemical definitions for rotatable X-H groups.

The cone detector keeps the observed heavy-atom coordinates fixed and rotates
only terminal hydrogens around the parent-X bond.  Definitions are explicit so
that chemically different groups do not inherit a single tetrahedral model.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional, Tuple

import gemmi


class RotatableGroupKind(str, Enum):
    ROTATABLE_SINGLE_H = "rotatable_single_h"
    ROTATABLE_NH3 = "rotatable_nh3"
    ROTATABLE_CH3 = "rotatable_ch3"


@dataclass(frozen=True)
class RotatableGroupDefinition:
    parent_atom_name: str
    kind: RotatableGroupKind
    element: str
    hydrogen_count: int
    xh_bond_length: float
    parent_x_h_angle: float
    rotation_period: int


def _single(parent: str, element: str, length: float,
            angle: float) -> RotatableGroupDefinition:
    return RotatableGroupDefinition(
        parent, RotatableGroupKind.ROTATABLE_SINGLE_H, element, 1,
        length, angle, 360,
    )


def _three(parent: str, kind: RotatableGroupKind, element: str, length: float,
           angle: float) -> RotatableGroupDefinition:
    return RotatableGroupDefinition(
        parent, kind, element, 3, length, angle, 120,
    )


# Bond lengths and parent-X-H angles follow the standard CCP4 monomer
# restraints.  Small residue-specific differences are retained where they
# materially define the cone (notably Cys S-H).
ROTATABLE_GROUPS: Dict[Tuple[str, str], RotatableGroupDefinition] = {
    ("SER", "OG"): _single("CB", "O", 0.972, 108.539),
    ("THR", "OG1"): _single("CB", "O", 0.972, 109.544),
    ("TYR", "OH"): _single("CZ", "O", 0.966, 109.970),
    ("CYS", "SG"): _single("CB", "S", 1.212, 108.4),
    ("LYS", "NZ"): _three(
        "CE", RotatableGroupKind.ROTATABLE_NH3, "N", 1.018, 109.659),

    ("ALA", "CB"): _three(
        "CA", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.742),
    ("VAL", "CG1"): _three(
        "CB", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("VAL", "CG2"): _three(
        "CB", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("LEU", "CD1"): _three(
        "CG", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("LEU", "CD2"): _three(
        "CG", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("ILE", "CD1"): _three(
        "CG1", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("ILE", "CG2"): _three(
        "CB", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("MET", "CE"): _three(
        "SD", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("MSE", "CE"): _three(
        "SE", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.5),
    ("THR", "CG2"): _three(
        "CB", RotatableGroupKind.ROTATABLE_CH3, "C", 1.092, 109.532),
}


def get_rotatable_group(
        residue_name: str,
        atom_name: str) -> Optional[RotatableGroupDefinition]:
    return ROTATABLE_GROUPS.get((residue_name.upper(), atom_name.upper()))


def is_rotatable(residue_name: str, atom_name: str) -> bool:
    return get_rotatable_group(residue_name, atom_name) is not None


def find_parent_atom(
        residue: gemmi.Residue,
        definition: RotatableGroupDefinition) -> Optional[gemmi.Atom]:
    return next(
        (atom for atom in residue if atom.name == definition.parent_atom_name),
        None,
    )
