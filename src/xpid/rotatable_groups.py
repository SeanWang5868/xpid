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


def normalized_altloc(atom: gemmi.Atom) -> str:
    """Return Gemmi's blank altloc marker as an ordinary empty string."""
    return "" if atom.altloc in ("", "\0") else atom.altloc


@dataclass(frozen=True)
class DonorConformer:
    """One chemically self-consistent rotatable donor conformer."""

    model_index: int
    chain_name: str
    residue_name: str
    residue_seqid: str
    x_atom: gemmi.Atom
    x_altloc: str
    x_occupancy: float
    parent_atom: gemmi.Atom
    parent_altloc: str
    parent_occupancy: float
    parent_selection: str

    @property
    def active_altloc(self) -> str:
        """Conformer label used for other atoms in the donor residue."""
        return self.x_altloc or self.parent_altloc

    @property
    def occupancy(self) -> float:
        return min(self.x_occupancy, self.parent_occupancy)


@dataclass(frozen=True)
class DonorConformerResolution:
    conformers: Tuple[DonorConformer, ...] = ()
    issue: Optional[str] = None


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


# X-H bond lengths represent nuclear positions and follow the CCP4
# ``value_dist_nucleus`` convention.  Parent-X-H angles follow the standard
# CCP4 geometry, including the non-tetrahedral Cys C-beta-S-H angle.
ROTATABLE_GROUPS: Dict[Tuple[str, str], RotatableGroupDefinition] = {
    ("SER", "OG"): _single("CB", "O", 0.972, 108.539),
    ("THR", "OG1"): _single("CB", "O", 0.972, 109.544),
    ("TYR", "OH"): _single("CZ", "O", 0.966, 109.970),
    ("CYS", "SG"): _single("CB", "S", 1.338, 97.543),
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
        definition: RotatableGroupDefinition,
        x_atom: Optional[gemmi.Atom] = None) -> Optional[gemmi.Atom]:
    """Return a unique compatible parent, never an order-dependent first hit.

    New code should use :func:`resolve_donor_conformers`, which can represent
    a shared blank-altloc X atom paired with multiple alternate parent states.
    This wrapper remains for compatibility with callers expecting one atom.
    """
    parents = [
        atom for atom in residue
        if atom.name == definition.parent_atom_name
    ]
    if x_atom is None:
        return parents[0] if len(parents) == 1 else None
    resolution = resolve_donor_conformers(
        residue, x_atom, definition, model_index=0, chain_name="")
    if len(resolution.conformers) != 1:
        return None
    return resolution.conformers[0].parent_atom


def resolve_donor_conformers(
        residue: gemmi.Residue,
        x_atom: gemmi.Atom,
        definition: RotatableGroupDefinition,
        model_index: int,
        chain_name: str,
        ) -> DonorConformerResolution:
    """Resolve parent atoms without mixing deposited alternate conformers.

    A labelled X atom requires a parent with the same label, falling back only
    to one shared blank parent.  A blank X atom first uses a blank parent; if
    none exists, each unique labelled parent defines a separate self-consistent
    conformer.  Duplicate atoms within one altloc are chemically ambiguous and
    are rejected rather than resolved by file order.
    """
    parents = [
        atom for atom in residue
        if atom.name == definition.parent_atom_name
    ]
    if not parents:
        return DonorConformerResolution(issue="missing_parent")

    by_altloc: Dict[str, list[gemmi.Atom]] = {}
    for parent in parents:
        by_altloc.setdefault(normalized_altloc(parent), []).append(parent)

    duplicate_altlocs = sorted(
        altloc or "<blank>"
        for altloc, atoms in by_altloc.items()
        if len(atoms) != 1
    )
    if duplicate_altlocs:
        return DonorConformerResolution(
            issue="duplicate_parent_altloc:" + ",".join(duplicate_altlocs))

    x_altloc = normalized_altloc(x_atom)
    selected: list[tuple[gemmi.Atom, str]] = []
    if x_altloc:
        if x_altloc in by_altloc:
            selected = [(by_altloc[x_altloc][0], "matching_altloc")]
        elif "" in by_altloc:
            selected = [(by_altloc[""][0], "shared_blank_parent")]
        else:
            available = ",".join(sorted(by_altloc)) or "none"
            return DonorConformerResolution(
                issue=f"incompatible_parent_altloc:{available}")
    elif "" in by_altloc:
        selected = [(by_altloc[""][0], "matching_blank_altloc")]
    else:
        selected = [
            (by_altloc[altloc][0], "alternate_parent_for_shared_x")
            for altloc in sorted(by_altloc)
        ]

    conformers = tuple(
        DonorConformer(
            model_index=model_index,
            chain_name=chain_name,
            residue_name=residue.name,
            residue_seqid=str(residue.seqid).strip(),
            x_atom=x_atom,
            x_altloc=x_altloc,
            x_occupancy=float(x_atom.occ),
            parent_atom=parent,
            parent_altloc=normalized_altloc(parent),
            parent_occupancy=float(parent.occ),
            parent_selection=selection,
        )
        for parent, selection in selected
    )
    return DonorConformerResolution(conformers=conformers)
