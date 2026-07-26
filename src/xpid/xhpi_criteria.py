"""Geometric evaluation of Hudson, Plevin, and P-slab XH–π criteria."""
from dataclasses import dataclass
import numpy as np
from typing import Dict, Any, Optional, Tuple
from . import config
from . import geometry


@dataclass(frozen=True)
class SpatialCriteria:
    p_dmax: float
    xpcn_angle: Optional[float]
    hudson_dist_ok: int
    hudson_proj_ok: int
    is_hudson_spatial: int
    plevin_dist_ok: int
    plevin_xpcn_ok: int
    is_plevin_spatial: int


def prepare_spatial_criteria(
        rctx, x_elem: str, x_pos: np.ndarray,
        dist_x_centroid: float,
        proj_dist: Optional[float]) -> SpatialCriteria:
    """Precompute X/π geometry that is invariant across H conformers."""
    p_dmax = config.P_PLANE_DMAX.get(
        x_elem, config.P_PLANE_DMAX["default"])
    xpcn_angle = geometry.calculate_xpcn_angle(
        x_pos, rctx.pi_center_arr, rctx.pi_normal)
    hudson_dist_ok = int(dist_x_centroid <= p_dmax)
    hudson_proj_ok = int(
        proj_dist is not None and proj_dist <= rctx.p_radius)
    plevin_dist_ok = int(dist_x_centroid < p_dmax)
    plevin_xpcn_ok = int(
        xpcn_angle is not None and
        xpcn_angle < config.PLEVIN_XPCN_MAX)
    return SpatialCriteria(
        p_dmax=p_dmax,
        xpcn_angle=xpcn_angle,
        hudson_dist_ok=hudson_dist_ok,
        hudson_proj_ok=hudson_proj_ok,
        is_hudson_spatial=int(hudson_dist_ok and hudson_proj_ok),
        plevin_dist_ok=plevin_dist_ok,
        plevin_xpcn_ok=plevin_xpcn_ok,
        is_plevin_spatial=int(plevin_dist_ok and plevin_xpcn_ok),
    )


def evaluate_hudson_plevin_batch(
        rctx,
        x_pos: np.ndarray,
        h_positions: np.ndarray,
        spatial: SpatialCriteria,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized binary Hudson/Plevin direction evaluation."""
    count = len(h_positions)
    theta = np.full(count, np.nan)
    xh_pi_angle = np.full(count, np.nan)
    if count == 0:
        empty = np.zeros(0, dtype=bool)
        return theta, xh_pi_angle, empty, empty

    v_xh = h_positions - x_pos
    norm_xh = np.linalg.norm(v_xh, axis=1)

    v_x_pi = rctx.pi_center_arr - x_pos
    normal = rctx.pi_normal
    norm_normal = np.linalg.norm(normal)
    toward_ring = np.einsum("ij,j->i", v_xh, v_x_pi) > 0
    hudson_valid = toward_ring & (norm_xh != 0) & (norm_normal != 0)
    if np.any(hudson_valid):
        cosines = (
            np.einsum("ij,j->i", v_xh[hudson_valid], normal) /
            (norm_xh[hudson_valid] * norm_normal)
        )
        angles = np.degrees(
            np.arccos(np.clip(cosines, -1.0, 1.0)))
        angles = np.where(angles >= 90.0, 180.0 - angles, angles)
        theta[hudson_valid] = angles

    v_hx = -v_xh
    v_hc = rctx.pi_center_arr - h_positions
    norm_hc = np.linalg.norm(v_hc, axis=1)
    plevin_valid = (norm_xh != 0) & (norm_hc != 0)
    if np.any(plevin_valid):
        cosines = (
            np.einsum(
                "ij,ij->i", v_hx[plevin_valid], v_hc[plevin_valid]) /
            (norm_xh[plevin_valid] * norm_hc[plevin_valid])
        )
        xh_pi_angle[plevin_valid] = np.degrees(
            np.arccos(np.clip(cosines, -1.0, 1.0)))

    is_hudson = (
        bool(spatial.is_hudson_spatial) &
        np.isfinite(theta) &
        (theta <= config.HUDSON_THETA_MAX)
    )
    is_plevin = (
        bool(spatial.is_plevin_spatial) &
        np.isfinite(xh_pi_angle) &
        (xh_pi_angle >= config.PLEVIN_XH_PI_MIN)
    )
    return theta, xh_pi_angle, is_hudson, is_plevin


def _cosine_between(vec_a: np.ndarray, vec_b: np.ndarray) -> Optional[float]:
    norm_a = np.linalg.norm(vec_a)
    norm_b = np.linalg.norm(vec_b)
    if norm_a == 0 or norm_b == 0:
        return None
    return float(np.dot(vec_a, vec_b) / (norm_a * norm_b))


def _ray_entry_distance(x_pos: np.ndarray, h_pos: np.ndarray, ray_t: Optional[float]) -> Optional[float]:
    if ray_t is None:
        return None
    return float(ray_t * np.linalg.norm(h_pos - x_pos))


def _canonical_unit_normal(normal: np.ndarray) -> Optional[np.ndarray]:
    norm = np.linalg.norm(normal)
    if norm == 0:
        return None
    unit = normal / norm
    pivot = int(np.argmax(np.abs(unit)))
    if unit[pivot] < 0:
        unit = -unit
    return unit


def _side_of_plane(point: np.ndarray, plane_point: np.ndarray, canonical_normal: np.ndarray) -> int:
    signed = float(np.dot(point - plane_point, canonical_normal))
    if abs(signed) < 1e-6:
        return 0
    return 1 if signed > 0 else -1


def _calculate_direction_metrics(rctx: "ring_conformers.RingContext", x_pos: np.ndarray, h_pos: np.ndarray,
                                 proj_dist: Optional[float]) -> Dict[str, Any]:
    v_xh = h_pos - x_pos
    v_x_centroid = rctx.pi_center_arr - x_pos
    xh_centroid_cos = _cosine_between(v_xh, v_x_centroid)

    xh_lateral_inward_score = None
    x_projection = geometry.project_point_to_plane(x_pos, rctx.pi_center_arr, rctx.pi_normal)
    norm_n = np.linalg.norm(rctx.pi_normal)
    if x_projection is not None and norm_n != 0:
        unit_normal = rctx.pi_normal / norm_n
        lateral_xh = v_xh - np.dot(v_xh, unit_normal) * unit_normal
        inward = rctx.pi_center_arr - x_projection
        xh_lateral_inward_score = _cosine_between(lateral_xh, inward)

    h_proj_dist = None
    h_ray_t = None
    h_ray_entry_dist = None
    slab_entry = geometry.calculate_xh_ray_p_slab_entry(
        x_pos, h_pos, rctx.pi_center_arr, rctx.pi_normal, rctx.p_slab_half_thickness)
    if slab_entry is not None:
        h_hit_pos, h_ray_t = slab_entry
        h_proj_dist = geometry.calculate_projection_dist(rctx.pi_normal, rctx.pi_center_arr, h_hit_pos)
        h_ray_entry_dist = _ray_entry_distance(x_pos, h_pos, h_ray_t)

    h_plane_proj_dist = None
    h_plane_t = None
    h_plane_entry_dist = None
    plane_entry = geometry.calculate_xh_ray_plane_entry(
        x_pos, h_pos, rctx.pi_center_arr, rctx.pi_normal)
    if plane_entry is not None:
        h_plane_pos, h_plane_t = plane_entry
        h_plane_proj_dist = geometry.calculate_projection_dist(
            rctx.pi_normal, rctx.pi_center_arr, h_plane_pos)
        h_plane_entry_dist = _ray_entry_distance(x_pos, h_pos, h_plane_t)

    delta_h_proj_dist = None
    if h_plane_proj_dist is not None and proj_dist is not None:
        delta_h_proj_dist = h_plane_proj_dist - proj_dist

    return {
        'xh_centroid_cos': xh_centroid_cos,
        'xh_lateral_inward_score': xh_lateral_inward_score,
        'h_proj_dist': h_proj_dist,
        'H_ray_t': h_ray_t,
        'H_ray_entry_dist': h_ray_entry_dist,
        'h_plane_proj_dist': h_plane_proj_dist,
        'H_plane_t': h_plane_t,
        'H_plane_entry_dist': h_plane_entry_dist,
        'delta_h_proj_dist': delta_h_proj_dist,
    }


def evaluate_xhpi_geometry(
        rctx: "ring_conformers.RingContext", x_elem: str, x_pos: np.ndarray,
        h_pos: np.ndarray, dist_x_plane: Optional[float],
        dist_x_centroid: float,
        proj_dist: Optional[float],
        spatial: Optional[SpatialCriteria] = None,
        include_direction_metrics: bool = True) -> Dict[str, Any]:
    """Calculate Hudson, Plevin, and P-slab metrics for one X-H candidate."""
    theta = geometry.calculate_hudson_theta(rctx.pi_center_arr, x_pos, h_pos, rctx.pi_normal)
    xh_pi_angle = geometry.calculate_xh_picenter_angle(rctx.pi_center_arr, x_pos, h_pos)

    if spatial is None:
        spatial = prepare_spatial_criteria(
            rctx, x_elem, x_pos, dist_x_centroid, proj_dist)
    xpcn_angle = spatial.xpcn_angle
    direction_metrics = (
        _calculate_direction_metrics(rctx, x_pos, h_pos, proj_dist)
        if include_direction_metrics else {}
    )

    hudson_dist_ok = spatial.hudson_dist_ok
    hudson_proj_ok = spatial.hudson_proj_ok
    hudson_direction_ok = int(theta is not None and theta <= config.HUDSON_THETA_MAX)
    is_hudson_spatial = spatial.is_hudson_spatial

    plevin_dist_ok = spatial.plevin_dist_ok
    plevin_xpcn_ok = spatial.plevin_xpcn_ok
    plevin_direction_ok = int(
        xh_pi_angle is not None and xh_pi_angle >= config.PLEVIN_XH_PI_MIN)
    is_plevin_spatial = spatial.is_plevin_spatial

    is_hudson = int(is_hudson_spatial and hudson_direction_ok)
    is_plevin = int(is_plevin_spatial and plevin_direction_ok)
    is_p_slab = int(
        include_direction_metrics and
        dist_x_plane is not None and dist_x_plane <= spatial.p_dmax and
        proj_dist is not None and proj_dist <= rctx.p_radius and
        direction_metrics.get('h_proj_dist') is not None and
        direction_metrics['h_proj_dist'] <= rctx.p_radius
    )

    metrics = {
        'dist_X_centroid': dist_x_centroid,
        'theta': theta,
        'angle_XPCN': xpcn_angle,
        'angle_XH_Pi': xh_pi_angle,
        'hudson_dist_ok': hudson_dist_ok,
        'hudson_proj_ok': hudson_proj_ok,
        'hudson_direction_ok': hudson_direction_ok,
        'is_hudson_spatial': is_hudson_spatial,
        'plevin_dist_ok': plevin_dist_ok,
        'plevin_xpcn_ok': plevin_xpcn_ok,
        'plevin_direction_ok': plevin_direction_ok,
        'is_plevin_spatial': is_plevin_spatial,
        'is_hudson': is_hudson,
        'is_plevin': is_plevin,
        'is_p_slab': is_p_slab,
    }
    metrics.update(direction_metrics)
    return metrics
