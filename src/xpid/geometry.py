"""
geometry.py
Geometric calculations for the P-model and cone alignment.
"""
import numpy as np
import gemmi
from typing import Tuple, Optional, List
from . import config

EPSILON = 1e-8

def get_pi_info(atoms: List[gemmi.Atom]) -> Tuple[gemmi.Position, np.ndarray, np.ndarray, float]:
    positions = np.array([[atom.pos.x, atom.pos.y, atom.pos.z] for atom in atoms])
    center_array = np.mean(positions, axis=0)
    pi_center = gemmi.Position(*center_array)
    b_mean = sum(atom.b_iso for atom in atoms) / len(atoms)
    
    centered_pos = positions - center_array
    _, _, vh = np.linalg.svd(centered_pos)
    normal_vector = vh[2, :] 
    
    return pi_center, center_array, normal_vector, b_mean

def calculate_planarity_deviation(atoms: List[gemmi.Atom]) -> float:
    """Calculate maximum deviation from the best-fit plane (in Å)."""
    if len(atoms) < 3:
        return 999.0
    
    _, center_arr, normal, _ = get_pi_info(atoms)
    norm_normal = np.linalg.norm(normal)
    if norm_normal == 0:
        return 999.0
    
    normal = normal / norm_normal  # Normalize
    
    max_dev = 0.0
    for atom in atoms:
        pos_arr = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
        vec = pos_arr - center_arr
        dev = np.abs(np.dot(normal, vec))
        if dev > max_dev:
            max_dev = dev
    
    return max_dev

def calculate_distance(pos1_array: np.ndarray, pos2_array: np.ndarray) -> float:
    return np.linalg.norm(pos1_array - pos2_array)


def calculate_plane_distance(point: np.ndarray, plane_point: np.ndarray, normal: np.ndarray) -> Optional[float]:
    """Perpendicular distance from a point to the P plane."""
    norm_n = np.linalg.norm(normal)
    if norm_n == 0:
        return None
    return abs(float(np.dot(point - plane_point, normal) / norm_n))


def project_point_to_plane(point: np.ndarray, plane_point: np.ndarray, normal: np.ndarray) -> Optional[np.ndarray]:
    """Orthogonally project a point onto the P plane."""
    denominator = float(np.dot(normal, normal))
    if denominator == 0:
        return None
    t = float(np.dot(point - plane_point, normal) / denominator)
    return point - t * normal


def calculate_p_offset(point_on_plane: np.ndarray, p_center: np.ndarray) -> float:
    """Distance from a point on the P plane to the center of the finite P region."""
    return calculate_distance(point_on_plane, p_center)


def point_is_in_p_region(point_on_plane: np.ndarray, p_center: np.ndarray, p_radius: float) -> bool:
    """Return True when a projected point lies inside the finite P region."""
    return calculate_p_offset(point_on_plane, p_center) <= p_radius


def calculate_xh_ray_p_slab_entry(
    x_pos: np.ndarray,
    h_pos: np.ndarray,
    plane_point: np.ndarray,
    normal: np.ndarray,
    half_thickness: float,
    min_t: float = 1.0,
) -> Optional[Tuple[np.ndarray, float]]:
    """Intersect the directional X->H ray with the near surface of the P slab.

    The finite P slab is centered on the aromatic plane and extends
    +/- half_thickness along the plane normal. For a donor X outside the slab,
    the relevant entry surface is the face on the same side of the aromatic
    plane as X. Requiring t > 1 keeps the same X->H directionality: H must lie
    between X and the P slab.
    """
    norm_n = np.linalg.norm(normal)
    if norm_n == 0:
        return None

    unit_normal = normal / norm_n
    z_x = float(np.dot(x_pos - plane_point, unit_normal))
    if abs(z_x) <= half_thickness:
        return None

    direction = h_pos - x_pos
    z_dir = float(np.dot(direction, unit_normal))
    if abs(z_dir) < EPSILON:
        return None

    entry_z = np.copysign(half_thickness, z_x)
    t = (entry_z - z_x) / z_dir
    if t <= min_t:
        return None

    return x_pos + t * direction, float(t)

def calculate_projection_dist(normal: np.ndarray, pi_center: np.ndarray, x_pos: np.ndarray) -> Optional[float]:
    projection_point = project_point_to_plane(x_pos, pi_center, normal)
    if projection_point is None:
        return None
    return calculate_p_offset(projection_point, pi_center)

def check_hbond_locked(x_pos: np.ndarray, 
                       orig_h_positions: list, 
                       acceptor_coords: np.ndarray, 
                       dist_cutoff: float = 3.5, 
                       angle_cutoff_deg: float = 120.0) -> bool:
    """Check if a donor's polar hydrogen is locked by a strong hydrogen bond.

    Returns True if any D-H...A angle >= 120° and H...A distance <= 3.5 Å.
    """
    if len(acceptor_coords) == 0 or len(orig_h_positions) == 0:
        return False
        
    for h_pos in orig_h_positions:
        vec_xh = h_pos - x_pos
        norm_xh = np.linalg.norm(vec_xh)
        if norm_xh == 0: continue
        vec_xh_norm = vec_xh / norm_xh
        
        # Compute distances from H to all potential acceptors
        dists = np.linalg.norm(acceptor_coords - h_pos, axis=1)
        valid_acceptors = acceptor_coords[dists <= dist_cutoff]
        
        for acc in valid_acceptors:
            vec_ha = acc - h_pos
            norm_ha = np.linalg.norm(vec_ha)
            if norm_ha == 0: continue
            
            # D-H...A angle via dot product
            cos_theta = np.dot(-vec_xh_norm, vec_ha / norm_ha)
            angle = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
            
            if angle >= angle_cutoff_deg:
                return True  # Locked by a strong hydrogen bond
                
    return False

def generate_rotated_hydrogens(parent_pos: np.ndarray, 
                               x_pos: np.ndarray, 
                               element: str, 
                               env_coords: np.ndarray = None, 
                               clash_cutoff: float = 2.0,
                               num_samples: int = 72) -> list:
    """Vectorized cone hydrogen generator with steric clash filtering."""
    axis = x_pos - parent_pos
    norm_axis = np.linalg.norm(axis)
    if norm_axis == 0: 
        return []
    u = axis / norm_axis 
    
    arbitrary_vec = np.array([1.0, 0.0, 0.0])
    if np.abs(np.dot(u, arbitrary_vec)) > 0.99:
        arbitrary_vec = np.array([0.0, 1.0, 0.0])
    
    v = np.cross(u, arbitrary_vec)
    v = v / np.linalg.norm(v)
    w = np.cross(u, v)
    
    bond_length = config.BOND_LENGTHS.get(element, 1.09)
    theta_rad = np.radians(config.TETRAHEDRAL_ANGLE)
    
    h_proj_u = bond_length * np.cos(np.pi - theta_rad)
    h_radius = bond_length * np.sin(np.pi - theta_rad)
    
    # Vectorized: generate all H positions at once
    angles = np.linspace(0, 2 * np.pi, num_samples, endpoint=False)
    cos_phi = np.cos(angles)  # (N,)
    sin_phi = np.sin(angles)  # (N,)
    
    # h_positions shape: (N, 3)
    h_positions = (x_pos + h_proj_u * u + 
                   h_radius * cos_phi[:, None] * v + 
                   h_radius * sin_phi[:, None] * w)
    
    # Vectorized clash check
    if env_coords is not None and len(env_coords) > 0:
        # diffs shape: (N, M, 3) where M = number of env atoms
        diffs = h_positions[:, None, :] - env_coords[None, :, :]
        dists = np.linalg.norm(diffs, axis=2)  # (N, M)
        min_dists = dists.min(axis=1)  # (N,)
        mask = min_dists >= clash_cutoff
        h_positions = h_positions[mask]
    
    return [h_positions[i] for i in range(len(h_positions))]


def calculate_pi_pi_geometry(center1: np.ndarray, normal1: np.ndarray,
                             center2: np.ndarray, normal2: np.ndarray
                             ) -> Tuple[float, float, float]:
    """Return centroid distance, inter-normal angle, and lateral offset."""
    vec = center2 - center1
    dist = np.linalg.norm(vec)

    n1_norm = np.linalg.norm(normal1)
    n2_norm = np.linalg.norm(normal2)
    if n1_norm == 0 or n2_norm == 0:
        return dist, 90.0, dist

    n1 = normal1 / n1_norm
    n2 = normal2 / n2_norm

    cos_angle = np.clip(np.abs(np.dot(n1, n2)), 0.0, 1.0)
    angle = np.degrees(np.arccos(cos_angle))

    proj_along = abs(np.dot(vec, n1)) if dist > 0 else 0.0
    offset = np.sqrt(max(dist**2 - proj_along**2, 0.0))

    return dist, angle, offset
