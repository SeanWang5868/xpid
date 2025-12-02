"""
geometry.py
Geometric calculations for distance, angles, and vector projections.
"""
import numpy as np
import gemmi
from typing import Tuple, Optional, List

def get_pi_info(atoms: List[gemmi.Atom]) -> Tuple[gemmi.Position, np.ndarray, np.ndarray, float]:
    """
    Calculates Pi-system center, normal vector (via SVD), and mean B-factor.
    """
    positions = np.array([atom.pos.tolist() for atom in atoms])
    
    # 1. Geometric Center
    center_array = np.mean(positions, axis=0)
    pi_center = gemmi.Position(*center_array)
    
    # 2. Mean B-factor
    b_mean = sum(atom.b_iso for atom in atoms) / len(atoms)
    
    # 3. Normal Vector (SVD Fitting)
    centered_pos = positions - center_array
    _, _, vh = np.linalg.svd(centered_pos)
    normal_vector = vh[2, :] # The vector corresponding to the smallest singular value
    
    return pi_center, center_array, normal_vector, b_mean

def calculate_distance(pos1_array: np.ndarray, pos2_array: np.ndarray) -> float:
    return np.linalg.norm(pos1_array - pos2_array)

def calculate_xpcn_angle(x_pos: np.ndarray, pi_center: np.ndarray, pi_normal: np.ndarray) -> Optional[float]:
    """Plevin: Angle between X-PiCenter vector and Normal vector."""
    v_x_pi = pi_center - x_pos
    norm_v = np.linalg.norm(v_x_pi)
    norm_n = np.linalg.norm(pi_normal)
    
    if norm_v == 0 or norm_n == 0: return None

    dot_product = np.dot(v_x_pi, pi_normal)
    cos_theta = np.clip(dot_product / (norm_v * norm_n), -1.0, 1.0)
    angle = np.degrees(np.arccos(cos_theta))
    
    if angle > 90:
        angle = 180 - angle
    return angle

def calculate_xh_picenter_angle(pi_center: np.ndarray, x_pos: np.ndarray, h_pos: np.ndarray) -> Optional[float]:
    """Plevin: Angle between X-H vector and X-PiCenter vector."""
    v_hx = h_pos - x_pos
    v_hc = h_pos - pi_center 
    
    norm_hx = np.linalg.norm(v_hx)
    norm_hc = np.linalg.norm(v_hc)
    
    if norm_hx == 0 or norm_hc == 0: return None
    
    cos_theta = np.clip(np.dot(v_hx, v_hc) / (norm_hx * norm_hc), -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

def calculate_hudson_theta(pi_center: np.ndarray, x_pos: np.ndarray, h_pos: np.ndarray, normal: np.ndarray) -> Optional[float]:
    """Hudson: Angle between X-H vector and Normal vector (checking projection direction)."""
    v_x_pi = pi_center - x_pos
    v_xh = h_pos - x_pos
    
    norm_xpi = np.linalg.norm(v_x_pi)
    if norm_xpi == 0: return None

    # Projection check: H must point towards the ring
    proj_len = np.dot(v_xh, v_x_pi) / norm_xpi
    
    if proj_len > 0:
        norm_n = np.linalg.norm(normal)
        norm_xh = np.linalg.norm(v_xh)
        if norm_n == 0 or norm_xh == 0: return None
        
        cos_angle = np.clip(np.dot(normal, v_xh) / (norm_n * norm_xh), -1.0, 1.0)
        angle = np.degrees(np.arccos(cos_angle))
        
        if angle >= 90:
            angle = 180 - angle
        return angle
    return None

def calculate_projection_dist(normal: np.ndarray, pi_center: np.ndarray, x_pos: np.ndarray) -> Optional[float]:
    """Hudson: Distance from the Pi center to the projection of X onto the plane."""
    numerator = np.dot(normal, pi_center - x_pos)
    denominator = np.dot(normal, normal)
    if denominator == 0: return None
    
    t = numerator / denominator
    projection_point = x_pos + t * normal
    return np.linalg.norm(projection_point - pi_center)