"""
schema.py
Centralised registry of all XH-π hit output columns.

Every column that can appear in a result row is defined here exactly once.
``core.py`` uses this to validate hits at build time; ``output.py`` uses it
to derive column lists and type coercions automatically.

Adding a new column means adding one entry here — no need to touch
multiple files or remember which list to update.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Field descriptor
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Field:
    """Descriptor for one column in an XH-π hit row."""

    name: str
    dtype: str = "float"  # "float", "int", "str"
    group: str = "base"  # "base", "p_slab", "coords", "sasa", "coop", "hbond", "candidate"
    default: object = None
    description: str = ""


# ---------------------------------------------------------------------------
# All columns — define once, use everywhere
# ---------------------------------------------------------------------------

FIELDS: Tuple[Field, ...] = (
    # -- Identity -----------------------------------------------------------
    Field("pdb", "str", "base", "", "PDB code or filename stem"),
    Field("model", "str", "base", "", "Model identifier"),
    Field("resolution", "float", "base", 0.0, "Nominal resolution (Å)"),
    # -- π acceptor ---------------------------------------------------------
    Field("pi_chain", "str", "base", "", "Chain of π-ring residue"),
    Field("pi_res", "str", "base", "", "Residue name of π acceptor"),
    Field("pi_id", "str", "base", "", "Sequence identifier of π residue"),
    Field("pi_ring_id", "str", "base", "", "Stable ring index within the π residue"),
    Field("pi_ring_size", "int", "base", 0, "Number of atoms in the π ring"),
    Field("pi_altloc", "str", "base", "", "Alternate-conformation label of the π ring"),
    Field("pi_ss_type", "str", "base", "", "Secondary-structure type (H/E/C)"),
    Field("pi_ss_id", "str", "base", "", "SS element identifier"),
    Field("pi_avg_b", "float", "base", 0.0, "Mean B-factor of π-ring atoms"),
    Field("is_trp_5ring_acceptor", "int", "base", 0, "1 if TRP 5-ring is the π acceptor"),
    # -- Donor --------------------------------------------------------------
    Field("X_chain", "str", "base", "", "Chain of donor residue"),
    Field("X_res", "str", "base", "", "Residue name of donor"),
    Field("X_id", "str", "base", "", "Sequence identifier of donor residue"),
    Field("X_atom", "str", "base", "", "Donor heavy-atom name"),
    Field("X_altloc", "str", "base", "", "Alternate-conformation label of donor X"),
    Field("X_b", "float", "base", 0.0, "B-factor of donor heavy atom"),
    Field("X_ss_type", "str", "base", "", "Donor SS type"),
    Field("X_ss_id", "str", "base", "", "Donor SS identifier"),
    # -- Hydrogen -----------------------------------------------------------
    Field("H_atom", "str", "base", "", "Hydrogen atom name (or 'virt' for cone)"),
    Field("H_source", "str", "base", "", "experimental / added / cone_virtual"),
    Field("H_altloc", "str", "base", "", "Alternate-conformation label of explicit H"),
    Field("combined_occupancy", "float", "base", 1.0,
          "Minimum occupancy across the ring, donor X and explicit H"),
    # -- Geometry (always present) ------------------------------------------
    Field("dist_X_Pi", "float", "base", None, "Perpendicular distance X → π plane (Å)"),
    Field("dist_X_centroid", "float", "base", None, "Distance X → ring centroid (Å)"),
    Field("proj_dist", "float", "base", None, "Lateral offset of X projection from ring centre (Å)"),
    Field("theta", "float", "base", None, "Hudson angle between X–H and ring normal (°)"),
    Field("angle_XPCN", "float", "base", None, "Plevin X→centroid→normal angle (°)"),
    Field("angle_XH_Pi", "float", "base", None, "Plevin X–H→centroid angle (°)"),
    # -- System labels ------------------------------------------------------
    Field("is_hudson", "int", "base", 0, "1 if Hudson-positive"),
    Field("is_plevin", "int", "base", 0, "1 if Plevin-positive"),
    # -- Additional flags ---------------------------------------------------
    Field("is_pi_pi_tshaped", "int", "base", 0, "1 if donor is also part of a T-shaped π-π stack"),
    Field("sym_op", "int", "base", 0, "Symmetry operation index (0 = asymmetric unit)"),
    Field("symmetry_code", "str", "base", "1_555",
          "Full crystallographic symmetry code including lattice translation"),
    Field("seq_sep", "int", "base", 0, "Sequence separation (π residue - donor residue)"),

    # -- P-slab (conditional) -----------------------------------------------
    Field("is_p_slab", "int", "p_slab", 0, "1 if P-slab-positive"),
    Field("P_radius", "float", "p_slab", None, "P-region radius (Å)"),
    Field("P_slab_half_thickness", "float", "p_slab", None, "Half-thickness of P slab (Å)"),
    Field("h_proj_dist", "float", "p_slab", None, "H-ray entry point lateral offset (Å)"),
    Field("H_ray_t", "float", "p_slab", None, "H-ray parameter t at slab entry"),
    Field("H_ray_entry_dist", "float", "p_slab", None, "Distance X → slab entry along H ray (Å)"),
    Field("h_plane_proj_dist", "float", "p_slab", None, "H-plane entry lateral offset (Å)"),
    Field("H_plane_t", "float", "p_slab", None, "H-ray parameter t at plane entry"),
    Field("H_plane_entry_dist", "float", "p_slab", None, "Distance X → plane entry (Å)"),
    Field("delta_h_proj_dist", "float", "p_slab", None, "h_plane_proj_dist − proj_dist (Å)"),

    # -- Coordinates (conditional) ------------------------------------------
    Field("pi_center_x", "float", "coords", None, "π-ring centroid x (Å)"),
    Field("pi_center_y", "float", "coords", None, "π-ring centroid y (Å)"),
    Field("pi_center_z", "float", "coords", None, "π-ring centroid z (Å)"),
    Field("pi_normal_x", "float", "coords", None, "Canonical π normal x"),
    Field("pi_normal_y", "float", "coords", None, "Canonical π normal y"),
    Field("pi_normal_z", "float", "coords", None, "Canonical π normal z"),
    Field("X_xyz_x", "float", "coords", None, "Donor X atom x (Å)"),
    Field("X_xyz_y", "float", "coords", None, "Donor X atom y (Å)"),
    Field("X_xyz_z", "float", "coords", None, "Donor X atom z (Å)"),
    Field("H_xyz_x", "float", "coords", None, "Hydrogen x (Å) — always computed internally"),
    Field("H_xyz_y", "float", "coords", None, "Hydrogen y (Å) — always computed internally"),
    Field("H_xyz_z", "float", "coords", None, "Hydrogen z (Å) — always computed internally"),
    Field("X_side_of_pi", "int", "coords", 0, "+1 / -1 / 0 for donor side of π plane"),

    # -- SASA (conditional) -------------------------------------------------
    Field("pi_sasa_avg", "float", "sasa", None, "Mean SASA of π-ring residue heavy atoms (Å²)"),
    Field("X_sasa", "float", "sasa", None, "SASA of donor heavy atom (Å²)"),
    Field("H_sasa", "float", "sasa", None, "SASA of hydrogen atom (Å²)"),

    # -- Cooperativity (always computed if coordinates present) --------------
    Field("coop_same_face_donors", "int", "coop", 1, "Distinct donors on the same π face"),
    Field("coop_opp_face_donors", "int", "coop", 0, "Distinct donors on the opposite π face"),
    Field("coop_total_donors", "int", "coop", 1, "Total distinct donors for this π ring"),
    Field("coop_bi_face", "int", "coop", 0, "1 if donors on both faces"),

    # -- H-bond competition (always computed, verbose output) ---------------
    Field("hbond_competing", "int", "hbond", 0, "1 if a competing H-bond acceptor was found"),
    Field("hbond_relation", "str", "hbond", "none",
          "none / same_hydrogen / same_conformer_other_hydrogen / alternative_conformer / multiple"),
    Field("hbond_also_hbonded", "int", "hbond", 0,
          "1 if an explicit donor H also satisfies conventional H-bond geometry"),
    Field("hbond_acceptor_atom", "str", "hbond", None, "Atom name of strongest competing acceptor"),
    Field("hbond_acceptor_res", "str", "hbond", None, "Residue of strongest competing acceptor"),
    Field("hbond_acceptor_chain", "str", "hbond", None, "Chain of strongest competing acceptor"),
    Field("hbond_HA_dist", "float", "hbond", None, "H···A distance (Å)"),
    Field("hbond_DHA_angle", "float", "hbond", None, "D–H···A angle (°)"),
    Field("hbond_vs_xhpi_score", "float", "hbond", None, "Legacy name for pure H-bond geometry score (0–1)"),

    # -- XH-candidate diagnostics (conditional) -----------------------------
    Field("is_xh_candidate", "int", "candidate", 0, "1 if exported as X-H candidate"),
    Field("is_hudson_spatial", "int", "candidate", 0, "1 if Hudson spatial criteria pass"),
    Field("is_plevin_spatial", "int", "candidate", 0, "1 if Plevin spatial criteria pass"),
    Field("hudson_dist_ok", "int", "candidate", 0, "Hudson distance check passed"),
    Field("hudson_proj_ok", "int", "candidate", 0, "Hudson projection check passed"),
    Field("hudson_direction_ok", "int", "candidate", 0, "Hudson direction check passed"),
    Field("plevin_dist_ok", "int", "candidate", 0, "Plevin distance check passed"),
    Field("plevin_xpcn_ok", "int", "candidate", 0, "Plevin XPCN angle check passed"),
    Field("plevin_direction_ok", "int", "candidate", 0, "Plevin direction check passed"),
    Field("xh_centroid_cos", "float", "candidate", None, "Cosine of X-H → centroid angle"),
    Field("xh_lateral_inward_score", "float", "candidate", None, "Cosine of lateral X-H vs inward direction"),
)


# ---------------------------------------------------------------------------
# Derived lookup helpers (built once at import time)
# ---------------------------------------------------------------------------

_FIELD_MAP: Dict[str, Field] = OrderedDict((f.name, f) for f in FIELDS)
_FLOAT_NAMES: List[str] = [f.name for f in FIELDS if f.dtype == "float"]
_INT_NAMES: List[str] = [f.name for f in FIELDS if f.dtype == "int"]

# Columns that always appear in simple (non-verbose) output
SIMPLE_NAMES: List[str] = [f.name for f in FIELDS if f.group == "base"]

# Grouped column name lists for conditional output
P_SLAB_NAMES: List[str] = [f.name for f in FIELDS if f.group == "p_slab"]
COORDS_NAMES: List[str] = [f.name for f in FIELDS if f.group == "coords"]
SASA_NAMES: List[str] = [f.name for f in FIELDS if f.group == "sasa"]
COOP_NAMES: List[str] = [f.name for f in FIELDS if f.group == "coop"]
HBOND_NAMES: List[str] = [f.name for f in FIELDS if f.group == "hbond"]
CANDIDATE_NAMES: List[str] = [f.name for f in FIELDS if f.group == "candidate"]


def field(name: str) -> Optional[Field]:
    """Look up a column descriptor by name."""
    return _FIELD_MAP.get(name)


def validate_hit(hit: Dict) -> List[str]:
    """Check that every key in *hit* is a known column.

    Returns a list of unexpected keys (empty means valid).
    """
    return [k for k in hit if k not in _FIELD_MAP and not k.startswith("_")]
