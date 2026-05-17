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


def test_xh_ray_p_intersection_requires_h_pointing_to_plane():
    p_center = np.array([0., 0., 0.])
    normal = np.array([0., 0., 1.])
    x_pos = np.array([0., 0., 3.0])
    h_pos = np.array([0., 0., 2.0])

    intersection = geometry.calculate_xh_ray_p_intersection(x_pos, h_pos, p_center, normal)

    assert intersection is not None
    hit, t = intersection
    assert np.allclose(hit, p_center)
    assert np.isclose(t, 3.0)

    assert geometry.calculate_xh_ray_p_intersection(
        x_pos, np.array([0., 0., 4.0]), p_center, normal) is None
    assert geometry.calculate_xh_ray_p_intersection(
        x_pos, np.array([1., 0., 3.0]), p_center, normal) is None
