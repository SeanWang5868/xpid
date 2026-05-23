import numpy as np
import pytest
from xpid import geometry

def test_distance():
    a = np.array([0., 0., 0.])
    b = np.array([3., 4., 0.])
    assert geometry.calculate_distance(a, b) == 5.0


def test_p_plane_distance_and_projection():
    p_center = np.array([0., 0., 0.])
    normal = np.array([0., 0., 2.])
    point = np.array([1.5, -0.5, 3.0])

    assert np.isclose(geometry.calculate_plane_distance(point, p_center, normal), 3.0)
    projected = geometry.project_point_to_plane(point, p_center, normal)

    assert np.allclose(projected, np.array([1.5, -0.5, 0.0]))
    assert np.isclose(geometry.calculate_p_offset(projected, p_center), np.sqrt(2.5))


def test_legacy_hudson_plevin_angles_are_directional():
    p_center = np.array([0., 0., 0.])
    normal = np.array([0., 0., 1.])
    x_pos = np.array([0., 0., 3.])
    h_pos = np.array([0., 0., 2.])

    assert np.isclose(geometry.calculate_xpcn_angle(x_pos, p_center, normal), 0.0)
    assert np.isclose(geometry.calculate_xh_picenter_angle(p_center, x_pos, h_pos), 180.0)
    assert np.isclose(geometry.calculate_hudson_theta(p_center, x_pos, h_pos, normal), 0.0)

    h_away = np.array([0., 0., 4.])
    assert geometry.calculate_hudson_theta(p_center, x_pos, h_away, normal) is None


def test_xh_ray_p_slab_entry_uses_near_surface():
    p_center = np.array([0., 0., 0.])
    normal = np.array([0., 0., 1.])
    x_pos = np.array([0.83, 0., 3.97])
    h_pos = np.array([1.13, 0., 2.9775])

    entry = geometry.calculate_xh_ray_p_slab_entry(
        x_pos, h_pos, p_center, normal, half_thickness=0.5)

    assert entry is not None
    hit, t = entry
    assert np.isclose(hit[2], 0.5)
    assert t > 1.0
    assert np.isclose(
        geometry.calculate_projection_dist(normal, p_center, hit),
        1.879, atol=1e-3)

    assert geometry.calculate_xh_ray_p_slab_entry(
        x_pos, np.array([0.83, 0., 4.5]), p_center, normal, 0.5) is None
