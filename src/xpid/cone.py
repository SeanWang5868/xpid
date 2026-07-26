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
    vdw_radius: Optional[float] = None
    is_hbond_acceptor: Optional[bool] = None
    occupancy: Optional[float] = None


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


def _environment_properties(
        env: EnvironmentAtom) -> Tuple[float, bool, float]:
    """Return cached vdW radius, acceptor typing and occupancy."""
    radius = (
        env.vdw_radius if env.vdw_radius is not None else
        VDW_RADII.get(
            env.atom.element.name.upper(), DEFAULT_VDW_RADIUS)
    )
    is_acceptor = (
        env.is_hbond_acceptor
        if env.is_hbond_acceptor is not None else
        hbond_acceptors.is_hbond_acceptor(env.residue, env.atom)
    )
    occupancy = (
        env.occupancy if env.occupancy is not None else env.atom.occ)
    return radius, is_acceptor, occupancy


def _hbond_geometry(x_pos: np.ndarray, h_pos: np.ndarray,
                    env: EnvironmentAtom) -> Tuple[bool, bool]:
    """Return (chemically valid H-bond contact, strong direction constraint)."""
    _, is_acceptor, occupancy = _environment_properties(env)
    if occupancy < ACCEPTOR_MIN_OCCUPANCY:
        return False, False
    if not is_acceptor:
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
        env_radius, _, _ = _environment_properties(env)
        if distance < CLASH_SCALE * (VDW_RADII["H"] + env_radius):
            return False, has_strong_hbond
    return True, has_strong_hbond


def _candidate_hbond_relation(
    hydrogen_index: int,
    strong_flags: np.ndarray,
    alternative_hbond_exists: bool,
) -> str:
    selected_strong = bool(strong_flags[hydrogen_index])
    other_strong = bool(np.any(np.delete(strong_flags, hydrogen_index)))
    if selected_strong and other_strong:
        return "multiple"
    if selected_strong:
        return "same_hydrogen"
    if other_strong:
        return "same_conformer_other_hydrogen"
    if alternative_hbond_exists:
        return "alternative_conformer"
    return "none"


def _classify_conformer_arrays(
    conformers: Sequence[HydrogenConformer],
    x_pos: np.ndarray,
    environment: Sequence[EnvironmentAtom],
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized steric validity and per-hydrogen strong-H-bond flags."""
    if not conformers:
        return np.zeros(0, dtype=bool), np.zeros((0, 0), dtype=bool)

    h_positions = np.asarray([
        conformer.hydrogen_positions for conformer in conformers
    ], dtype=float)
    n_conformers, n_hydrogens, _ = h_positions.shape
    if not environment:
        return (
            np.ones(n_conformers, dtype=bool),
            np.zeros((n_conformers, n_hydrogens), dtype=bool),
        )

    env_positions = np.asarray(
        [env.position for env in environment], dtype=float)
    properties = [_environment_properties(env) for env in environment]
    env_radii = np.asarray([item[0] for item in properties])
    acceptor_mask = np.asarray([
        item[1] and item[2] >= ACCEPTOR_MIN_OCCUPANCY
        for item in properties
    ], dtype=bool)

    h_to_env = (
        env_positions[None, None, :, :] - h_positions[:, :, None, :])
    distances = np.linalg.norm(h_to_env, axis=-1)

    h_to_x = x_pos[None, None, :] - h_positions
    xh_norms = np.linalg.norm(h_to_x, axis=-1)
    denominators = xh_norms[:, :, None] * distances
    dot_products = np.einsum("chd,ched->che", h_to_x, h_to_env)
    cosines = np.divide(
        dot_products, denominators,
        out=np.ones_like(dot_products),
        where=denominators != 0,
    )
    angles = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))

    hbond_contacts = (
        acceptor_mask[None, None, :] &
        (distances >= ABSOLUTE_MIN_HA) &
        (distances <= HBOND_CONTACT_HA_MAX) &
        (angles >= HBOND_CONTACT_DHA_MIN)
    )
    strong_hbonds = (
        hbond_contacts &
        (distances <= STRONG_HBOND_HA_MAX) &
        (angles >= STRONG_HBOND_DHA_MIN)
    )

    clash_thresholds = (
        CLASH_SCALE * (VDW_RADII["H"] + env_radii))
    clashes = (
        ~hbond_contacts &
        (
            (distances < ABSOLUTE_MIN_HA) |
            (distances < clash_thresholds[None, None, :])
        )
    )
    hydrogen_valid = ~np.any(clashes, axis=2)
    conformer_valid = np.all(hydrogen_valid, axis=1)
    hydrogen_strong = np.any(strong_hbonds, axis=2)
    return conformer_valid, hydrogen_strong


def classify_conformers(
    conformers: Sequence[HydrogenConformer],
    x_pos: np.ndarray,
    environment: Sequence[EnvironmentAtom],
) -> Tuple[List[HydrogenConformer], List[HydrogenConformer]]:
    """Return sterically allowed and H-bond-capable conformers."""
    valid_mask, strong_flags = _classify_conformer_arrays(
        conformers, x_pos, environment)
    sterically_allowed = [
        conformer for index, conformer in enumerate(conformers)
        if valid_mask[index]
    ]
    hbond_capable = [
        conformer for index, conformer in enumerate(conformers)
        if valid_mask[index] and np.any(strong_flags[index])
    ]
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
    valid_mask, strong_flags = _classify_conformer_arrays(
        conformers, x_pos, environment)
    valid_indices = np.flatnonzero(valid_mask)
    if valid_indices.size == 0:
        return None

    alternative_hbond_exists = bool(np.any(strong_flags[valid_mask]))

    spatial = xhpi_criteria.prepare_spatial_criteria(
        rctx, x_element, x_pos, dist_x_centroid, proj_dist)

    # P-slab is optional and relatively uncommon. Preserve its exact scalar
    # path, while the default Hudson/Plevin production path is vectorized.
    if include_p_slab:
        positive: List[PositiveEvidence] = []
        for conformer_index in valid_indices:
            conformer = conformers[int(conformer_index)]
            for h_index, h_pos in enumerate(conformer.hydrogen_positions):
                metrics = xhpi_criteria.evaluate_xhpi_geometry(
                    rctx, x_element, x_pos, h_pos,
                    dist_x_plane, dist_x_centroid, proj_dist,
                    spatial=spatial,
                )
                if (metrics["is_hudson"] or metrics["is_plevin"] or
                        metrics["is_p_slab"]):
                    relation = _candidate_hbond_relation(
                        h_index, strong_flags[conformer_index],
                        alternative_hbond_exists)
                    positive.append(PositiveEvidence(
                        conformer, h_index, metrics, relation))
        return max(positive, key=_evidence_rank) if positive else None

    candidate_map = [
        (int(conformer_index), h_index)
        for conformer_index in valid_indices
        for h_index in range(
            len(conformers[int(conformer_index)].hydrogen_positions))
    ]
    h_positions = np.asarray([
        conformers[conformer_index].hydrogen_positions[h_index]
        for conformer_index, h_index in candidate_map
    ])
    theta, xh_pi, is_hudson, is_plevin = (
        xhpi_criteria.evaluate_hudson_plevin_batch(
            rctx, x_pos, h_positions, spatial))
    positive_indices = np.flatnonzero(is_hudson | is_plevin)
    if positive_indices.size == 0:
        return None

    def rank(candidate_index: int) -> tuple:
        conformer_index, h_index = candidate_map[candidate_index]
        hudson = int(is_hudson[candidate_index])
        plevin = int(is_plevin[candidate_index])
        theta_value = (
            float(theta[candidate_index])
            if np.isfinite(theta[candidate_index]) else None)
        xh_pi_value = (
            float(xh_pi[candidate_index])
            if np.isfinite(xh_pi[candidate_index]) else None)
        return (
            int(hudson and plevin),
            hudson + plevin,
            -(theta_value if theta_value is not None else 999.0),
            xh_pi_value if xh_pi_value is not None else -999.0,
            -conformers[conformer_index].phi,
            -h_index,
        )

    selected_index = max(
        (int(index) for index in positive_indices), key=rank)
    conformer_index, h_index = candidate_map[selected_index]
    conformer = conformers[conformer_index]
    metrics = xhpi_criteria.evaluate_xhpi_geometry(
        rctx, x_element, x_pos, conformer.hydrogen_positions[h_index],
        dist_x_plane, dist_x_centroid, proj_dist,
        spatial=spatial, include_direction_metrics=False,
    )
    relation = _candidate_hbond_relation(
        h_index, strong_flags[conformer_index],
        alternative_hbond_exists)
    return PositiveEvidence(conformer, h_index, metrics, relation)
