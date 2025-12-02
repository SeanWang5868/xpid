import numpy as np
import pytest
from xpid import geometry

def test_distance():
    a = np.array([0., 0., 0.])
    b = np.array([3., 4., 0.])
    assert geometry.calculate_distance(a, b) == 5.0

def test_xpcn_angle():
    # X on Z axis, Center at Origin, Normal is Z
    pi_center = np.array([0., 0., 0.])
    pi_normal = np.array([0., 0., 1.])
    x_pos = np.array([0., 0., 5.])
    
    # Vector X->Pi is [0,0,-5]. Dot with Normal [0,0,1] is -5.
    # Angle is 180. Logic converts >90 to 180-angle -> 0.
    angle = geometry.calculate_xpcn_angle(x_pos, pi_center, pi_normal)
    assert np.isclose(angle, 0.0)

def test_hudson_theta():
    # Setup where H points towards ring
    pi_center = np.array([0., 0., 0.])
    normal = np.array([0., 0., 1.])
    x_pos = np.array([5., 0., 0.])
    h_pos = np.array([4., 0., 0.]) # X-H points to center
    
    # X-H vector: [-1, 0, 0]. Normal: [0, 0, 1]. Dot is 0. Angle 90.
    angle = geometry.calculate_hudson_theta(pi_center, x_pos, h_pos, normal)
    assert np.isclose(angle, 90.0) # Hudson logic allows <=90 normalization