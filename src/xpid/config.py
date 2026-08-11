"""Numerical thresholds and small static chemical exclusions."""
from typing import Dict

from .aromatic_rings import (
    STANDARD_AROMATIC_RINGS as FALLBACK_RINGS,
    _CACHE as AROMATIC_RINGS_CACHE,
    get_aromatic_rings,
)
from .monomer_bonds import (
    _BONDED_HYDROGEN_CACHE as BONDED_HYDROGENS_CACHE,
    get_bonded_hydrogen_names,
)

# --- Defaults ---
DEFAULT_H_CHANGE = 4 

TARGET_ELEMENTS_X = {'C', 'N', 'O', 'S'}


def get_bonded_hydrogens(res_name: str, atom_name: str) -> set[str]:
    """Backward-compatible set-only dictionary lookup."""
    return get_bonded_hydrogen_names(res_name, atom_name) or set()

DIST_SEARCH_LIMIT = 6.0

# Unified P-model thresholds. P is the finite aromatic pi-plane region used
# for X projection and X-H ray intersection tests.
P_PLANE_DMAX = {
    'N': 4.3,
    'O': 4.3,
    'C': 4.5,
    'S': 4.8, 
    'default': 4.5
}
HUDSON_THETA_MAX = 40.0
PLEVIN_XPCN_MAX = 25.0
PLEVIN_XH_PI_MIN = 120.0
P_RADIUS_BY_RING_SIZE = {
    5: 1.6,
    6: 2.0,
}
# Half-thickness of the finite P slab used for the X-H directional ray test.
# X-position is still filtered against the zero-thickness P-plane projection.
P_SLAB_HALF_THICKNESS = 0.5

MIN_COVALENT_XH = 0.5
# Element-specific upper bounds for an explicit covalent X-H/D bond.  A
# single 1.3 A cutoff is not chemically valid for thiols and caused genuine
# neutron S-D positions to disappear in ``--no-cone`` mode.
COVALENT_XH_MAX = {
    'C': 1.25,
    'N': 1.25,
    'O': 1.20,
    'S': 1.55,
}
DIST_CUTOFF_H = max(COVALENT_XH_MAX.values())

# Cation-π donor atoms (positively charged groups)
CATION_DONORS = {
    ('LYS', 'NZ'),
    ('ARG', 'NH1'),
    ('ARG', 'NH2'),
    ('ARG', 'NE'),
}

PI_PI_DIST_MAX = 5.5
PI_PI_ANGLE_TSHAPED_MIN = 60.0
