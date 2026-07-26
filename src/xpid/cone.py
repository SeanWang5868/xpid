"""Binary conformer-existence detector for rotatable X-H groups."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import gemmi
import numpy as np

from . import hbond_acceptors
from . import rotatable_groups
from . import xhpi_criteria


VDW_RADII = {
    "H": 1.20, "D": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
    "S": 1.80, "P": 1.80, "SE": 1.90,
}
DEFAULT_VDW_RADIUS = 1.70
CLASH_SCALE = 0.75
ABSOLUTE_MIN_HA = 1.30
STRONG_HBOND_HA_MAX = 2.50
STRONG_HBOND_DHA_MIN = 140.0
HBOND_CONTACT_HA_MAX = 2.80
HBOND_CONTACT_DHA_MIN = 120.0
ACCEPTOR_MIN_OCCUPANCY = 0.20


@dataclass(frozen=True)
class HydrogenConformer:
    phi: float
    hydrogen_positions: Tuple[np.ndarray, ...]


@dataclass(frozen=True)
class EnvironmentAtom:
    position: np.ndarray
    atom: gemmi.Atom
    residue: gemmi.Residue
    chain_name: str
    image_idx: int


@dataclass(frozen=True)
class PositiveEvidence:
    conformer: HydrogenConformer
    hydrogen_index: int
    metrics: Dict[str, Any]
    hbond_relation: str = "none"


def _orthonormal_frame(parent_pos: np.ndarray,
                       x_pos: np.ndarray) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    axis = x_pos - parent_pos
    norm = np.linalg.norm(axis)
    if norm == 0:
        return None
    u = axis / norm
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(u, reference))) > 0.99:
        reference = np.array([0.0, 1.0, 0.0])
    v = np.cross(u, reference)
    v /= np.linalg.norm(v)
    w = np.cross(u, v)
    return u, v, w


def generate_conformers(parent_pos: np.ndarray, x_pos: np.ndarray,
                        definition: rotatable_groups.RotatableGroupDefinition,
                        step_degrees: int = 1) -> List[HydrogenConformer]:
    """Generate complete one- or three-hydrogen rotational conformers."""
    frame = _orthonormal_frame(parent_pos, x_pos)
    if frame is None or step_degrees <= 0:
        return []
    u, v, w = frame

    beta = np.radians(180.0 - definition.parent_x_h_angle)
    axial = definition.xh_bond_length * np.cos(beta)
    radial = definition.xh_bond_length * np.sin(beta)

    conformers: List[HydrogenConformer] = []
    for phi in range(0, definition.rotation_period, step_degrees):
        offsets = (0.0,) if definition.hydrogen_count == 1 else (0.0, 120.0, 240.0)
        positions = []
        for offset in offsets:
            angle = np.radians(phi + offset)
            position = (
                x_pos + axial * u +
                radial * np.cos(angle) * v +
                radial * np.sin(angle) * w
            )
            positions.append(position)
        conformers.append(HydrogenConformer(float(phi), tuple(positions)))
    return conformers


def _dha_angle(x_pos: np.ndarray, h_pos: np.ndarray,
               acceptor_pos: np.ndarray) -> Optional[float]:
    h_to_x = x_pos - h_pos
    h_to_a = acceptor_pos - h_pos
    denominator = np.linalg.norm(h_to_x) * np.linalg.norm(h_to_a)
    if denominator == 0:
        return None
    cosine = np.clip(np.dot(h_to_x, h_to_a) / denominator, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def _hbond_geometry(x_pos: np.ndarray, h_pos: np.ndarray,
                    env: EnvironmentAtom) -> Tuple[bool, bool]:
    """Return (chemically valid H-bond contact, strong direction constraint)."""
    if env.atom.occ < ACCEPTOR_MIN_OCCUPANCY:
        return False, False
    if not hbond_acceptors.is_hbond_acceptor(env.residue, env.atom):
        return False, False
    distance = float(np.linalg.norm(env.position - h_pos))
    if not (ABSOLUTE_MIN_HA <= distance <= HBOND_CONTACT_HA_MAX):
        return False, False
    angle = _dha_angle(x_pos, h_pos, env.position)
    if angle is None or angle < HBOND_CONTACT_DHA_MIN:
        return False, False
    strong = (
        distance <= STRONG_HBOND_HA_MAX and
        angle >= STRONG_HBOND_DHA_MIN
    )
    return True, strong


def _hydrogen_state(x_pos: np.ndarray, h_pos: np.ndarray,
                    environment: Sequence[EnvironmentAtom]) -> Tuple[bool, bool]:
    """Return (sterically_valid, has_strong_hbond)."""
    has_strong_hbond = False
    for env in environment:
        distance = float(np.linalg.norm(env.position - h_pos))
        hbond_contact, strong_hbond = _hbond_geometry(
            x_pos, h_pos, env)
        if hbond_contact:
            has_strong_hbond = has_strong_hbond or strong_hbond
            continue

        if distance < ABSOLUTE_MIN_HA:
            return False, has_strong_hbond
        env_radius = VDW_RADII.get(
            env.atom.element.name.upper(), DEFAULT_VDW_RADIUS)
        if distance < CLASH_SCALE * (VDW_RADII["H"] + env_radius):
            return False, has_strong_hbond
    return True, has_strong_hbond


def _strong_hbond_directions(
    x_pos: np.ndarray,
    h_pos: np.ndarray,
    environment: Sequence[EnvironmentAtom],
) -> List[np.ndarray]:
    """Return unit H→acceptor directions for strong conventional H-bonds."""
    directions = []
    for env in environment:
        _, strong = _hbond_geometry(x_pos, h_pos, env)
        if not strong:
            continue
        direction = env.position - h_pos
        norm = np.linalg.norm(direction)
        if norm:
            directions.append(direction / norm)
    return directions


def _candidate_hbond_relation(
    conformer: HydrogenConformer,
    hydrogen_index: int,
    x_pos: np.ndarray,
    environment: Sequence[EnvironmentAtom],
    alternative_hbond_exists: bool,
) -> str:
    selected_strong = False
    other_strong = False
    for index, h_pos in enumerate(conformer.hydrogen_positions):
        strong = bool(_strong_hbond_directions(x_pos, h_pos, environment))
        if index == hydrogen_index:
            selected_strong = strong
        else:
            other_strong = other_strong or strong
    if selected_strong and other_strong:
        return "multiple"
    if selected_strong:
        return "same_hydrogen"
    if other_strong:
        return "same_conformer_other_hydrogen"
    if alternative_hbond_exists:
        return "alternative_conformer"
    return "none"


def classify_conformers(
    conformers: Sequence[HydrogenConformer],
    x_pos: np.ndarray,
    environment: Sequence[EnvironmentAtom],
) -> Tuple[List[HydrogenConformer], List[HydrogenConformer]]:
    """Return sterically allowed and H-bond-capable conformers."""
    sterically_allowed: List[HydrogenConformer] = []
    hbond_capable: List[HydrogenConformer] = []

    for conformer in conformers:
        conformer_valid = True
        conformer_hbond = False
        for h_pos in conformer.hydrogen_positions:
            h_valid, h_hbond = _hydrogen_state(x_pos, h_pos, environment)
            if not h_valid:
                conformer_valid = False
                break
            conformer_hbond = conformer_hbond or h_hbond
        if conformer_valid:
            sterically_allowed.append(conformer)
            if conformer_hbond:
                hbond_capable.append(conformer)
    return sterically_allowed, hbond_capable


def _evidence_rank(evidence: PositiveEvidence) -> tuple:
    metrics = evidence.metrics
    both = int(metrics["is_hudson"] and metrics["is_plevin"])
    theta = metrics.get("theta")
    xh_pi = metrics.get("angle_XH_Pi")
    return (
        both,
        int(metrics["is_hudson"]) + int(metrics["is_plevin"]),
        -(theta if theta is not None else 999.0),
        xh_pi if xh_pi is not None else -999.0,
        -evidence.conformer.phi,
        -evidence.hydrogen_index,
    )


def evaluate_binary(
    rctx,
    definition: rotatable_groups.RotatableGroupDefinition,
    parent_pos: np.ndarray,
    x_pos: np.ndarray,
    x_element: str,
    environment: Sequence[EnvironmentAtom],
    dist_x_plane: Optional[float],
    dist_x_centroid: float,
    proj_dist: Optional[float],
    include_p_slab: bool = False,
) -> Optional[PositiveEvidence]:
    """Return one self-consistent positive H conformer, or ``None``.

    The production default is a steric-only allowed set.  H-bond geometry is
    retained as descriptive context and never changes the binary result.
    """
    conformers = generate_conformers(parent_pos, x_pos, definition)
    sterically_allowed, hbond_capable = classify_conformers(
        conformers, x_pos, environment)
    if not sterically_allowed:
        return None

    alternative_hbond_exists = bool(hbond_capable)

    positive: List[PositiveEvidence] = []
    for conformer in sterically_allowed:
        for h_index, h_pos in enumerate(conformer.hydrogen_positions):
            metrics = xhpi_criteria.evaluate_xhpi_geometry(
                rctx, x_element, x_pos, h_pos,
                dist_x_plane, dist_x_centroid, proj_dist,
            )
            if metrics["is_hudson"] or metrics["is_plevin"] or (
                    include_p_slab and metrics["is_p_slab"]):
                relation = _candidate_hbond_relation(
                    conformer, h_index, x_pos, environment,
                    alternative_hbond_exists)
                positive.append(PositiveEvidence(
                    conformer, h_index, metrics, relation))

    return max(positive, key=_evidence_rank) if positive else None
