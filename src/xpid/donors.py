"""Chemical definitions for rotatable X-H donor groups.

The cone detector keeps the observed heavy-atom coordinates fixed and rotates
only terminal hydrogens around the parent-X bond.  Definitions are explicit so
that chemically different groups do not inherit a single tetrahedral model.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional, Tuple

import gemmi


class DonorKind(str, Enum):
    ROTATABLE_SINGLE_H = "rotatable_single_h"
    ROTATABLE_NH3 = "rotatable_nh3"
    ROTATABLE_CH3 = "rotatable_ch3"


@dataclass(frozen=True)
class DonorDefinition:
    parent_atom_name: str
    kind: DonorKind
    element: str
    hydrogen_count: int
    xh_bond_length: float
    parent_x_h_angle: float
    rotation_period: int
    allow_hbond_constraint: bool


def _single(parent: str, element: str, length: float, angle: float) -> DonorDefinition:
    return DonorDefinition(
        parent, DonorKind.ROTATABLE_SINGLE_H, element, 1,
        length, angle, 360, True,
    )


def _three(parent: str, kind: DonorKind, element: str, length: float,
           angle: float, hbond: bool) -> DonorDefinition:
    return DonorDefinition(
        parent, kind, element, 3, length, angle, 120, hbond,
    )


# Bond lengths and parent-X-H angles follow the standard CCP4 monomer
# restraints.  Small residue-specific differences are retained where they
# materially define the cone (notably Cys S-H).
ROTATABLE_DONORS: Dict[Tuple[str, str], DonorDefinition] = {
    ("SER", "OG"): _single("CB", "O", 0.972, 108.539),
    ("THR", "OG1"): _single("CB", "O", 0.972, 109.544),
    ("TYR", "OH"): _single("CZ", "O", 0.966, 109.970),
    ("CYS", "SG"): _single("CB", "S", 1.212, 108.4),
    ("LYS", "NZ"): _three(
        "CE", DonorKind.ROTATABLE_NH3, "N", 1.018, 109.659, True),

    ("ALA", "CB"): _three(
        "CA", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.742, False),
    ("VAL", "CG1"): _three(
        "CB", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("VAL", "CG2"): _three(
        "CB", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("LEU", "CD1"): _three(
        "CG", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("LEU", "CD2"): _three(
        "CG", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("ILE", "CD1"): _three(
        "CG1", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("ILE", "CG2"): _three(
        "CB", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("MET", "CE"): _three(
        "SD", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("MSE", "CE"): _three(
        "SE", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.5, False),
    ("THR", "CG2"): _three(
        "CB", DonorKind.ROTATABLE_CH3, "C", 1.092, 109.532, False),
}


def get_definition(residue_name: str, atom_name: str) -> Optional[DonorDefinition]:
    return ROTATABLE_DONORS.get((residue_name.upper(), atom_name.upper()))


def is_rotatable(residue_name: str, atom_name: str) -> bool:
    return get_definition(residue_name, atom_name) is not None


def parent_atom(residue: gemmi.Residue,
                definition: DonorDefinition) -> Optional[gemmi.Atom]:
    return next(
        (atom for atom in residue if atom.name == definition.parent_atom_name),
        None,
    )
