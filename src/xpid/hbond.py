"""
hbond.py
Conventional hydrogen-bond competition analysis for XH–π hits.

For each detected XH–π interaction, searches for nearby H-bond acceptors
(O, N, S) that could form a competing conventional hydrogen bond with the
same donor group.  Reports the strongest competitor with geometric metrics
and a competition score.
"""
from __future__ import annotations

import math
from typing import Dict, List, Any, Optional, Tuple

import gemmi
import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HBOND_SEARCH_RADIUS = 3.0       # Å — search for acceptors around H
HBOND_HA_MAX = 2.8              # Å — max H···A distance for competing H-bond
HBOND_DHA_MIN = 120.0           # ° — min D–H···A angle for competing H-bond
HBOND_ACCEPTOR_ELEMENTS = {"O", "N", "S"}

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def _competition_score(
    d_ha: float,
    angle_dha: float,
    dist_x_pi: float,
    angle_xh_pi: Optional[float],
) -> float:
    """Compute a competition score comparing H-bond vs XH-π geometry fit.

    Positive → conventional H-bond is geometrically preferred.
    Negative → XH-π is geometrically preferred.
    Near zero → ambiguous.
    """
    hbond_term = (d_ha / HBOND_HA_MAX) * (angle_dha / 180.0)
    xh_pi_term = (dist_x_pi / 4.5)
    if angle_xh_pi is not None:
        xh_pi_term *= (angle_xh_pi / 180.0)
    else:
        xh_pi_term *= 1.0

    return hbond_term - xh_pi_term


# ---------------------------------------------------------------------------
# Main annotation function
# ---------------------------------------------------------------------------


def annotate_hbond_competition(
    hits: List[Dict[str, Any]],
    structure: gemmi.Structure,
    model_index: int = 0,
    ns: Optional[gemmi.NeighborSearch] = None,
) -> List[Dict[str, Any]]:
    """Annotate each hit with conventional H-bond competition metrics.

    Requires ``include_coordinates=True`` so that H atom positions are available
    (``H_xyz_x/y/z`` columns). For each hit, searches for nearby O/N/S atoms
    that could compete for the donor H and reports the strongest competitor.

    Appends the following columns:

    ``hbond_competing``
        0/1 — whether a competing H-bond acceptor was found within cutoffs.
    ``hbond_acceptor_atom``
        Atom name of the strongest competing acceptor.
    ``hbond_acceptor_res``
        Residue name of the strongest competing acceptor.
    ``hbond_acceptor_chain``
        Chain of the strongest competing acceptor.
    ``hbond_HA_dist``
        H···A distance (Å) to the strongest competitor.
    ``hbond_DHA_angle``
        D–H···A angle (°) to the strongest competitor.
    ``hbond_vs_xhpi_score``
        Competition score: positive → H-bond preferred,
        negative → XH-π preferred.

    Returns the modified list (in-place).
    """
    if len(structure) == 0 or model_index >= len(structure):
        return hits

    if ns is None:
        model = structure[model_index]
        ns = gemmi.NeighborSearch(model, structure.cell, 5.0)
        ns.populate(include_h=False)
    else:
        model = structure[model_index]

    for hit in hits:
        # Extract H atom position
        hx = hit.get("H_xyz_x")
        hy = hit.get("H_xyz_y")
        hz = hit.get("H_xyz_z")

        if hx is None or hy is None or hz is None:
            _set_no_competition(hit)
            continue

        h_pos = gemmi.Position(hx, hy, hz)

        # Extract X atom position for angle calculation
        xx = hit.get("X_xyz_x")
        xy = hit.get("X_xyz_y")
        xz = hit.get("X_xyz_z")

        if xx is None:
            _set_no_competition(hit)
            continue

        x_pos_arr = np.array([xx, xy, xz])
        h_pos_arr = np.array([hx, hy, hz])

        # Get pi-plane distance for scoring
        dist_x_pi = float(hit.get("dist_X_Pi", 4.5))
        angle_xh_pi = hit.get("angle_XH_Pi")

        # Exclude atoms from the donor residue itself
        donor_chain = hit.get("X_chain", "")
        donor_res = hit.get("X_res", "")
        donor_id = str(hit.get("X_id", ""))

        # Search for nearby atoms
        best_acceptor = None
        best_score = None

        marks = ns.find_atoms(h_pos, radius=HBOND_SEARCH_RADIUS)

        for mark in marks:
            cra = mark.to_cra(model)
            acceptor = cra.atom

            # Skip non-acceptor elements
            elem = acceptor.element.name.upper()
            if elem not in HBOND_ACCEPTOR_ELEMENTS:
                continue

            # Skip atoms from the same residue as the donor
            if (cra.chain.name == donor_chain and
                cra.residue.name == donor_res and
                str(cra.residue.seqid).strip() == donor_id):
                continue

            # Compute H-bond geometry
            a_pos_arr = np.array([mark.pos.x, mark.pos.y, mark.pos.z])
            d_ha = float(np.linalg.norm(a_pos_arr - h_pos_arr))

            if d_ha > HBOND_HA_MAX:
                continue

            vec_hd = x_pos_arr - h_pos_arr
            vec_ha = a_pos_arr - h_pos_arr
            norm_hd = np.linalg.norm(vec_hd)
            norm_ha = np.linalg.norm(vec_ha)

            if norm_hd == 0 or norm_ha == 0:
                continue

            cos_angle = np.dot(vec_hd, vec_ha) / (norm_hd * norm_ha)
            angle_dha = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))

            if angle_dha < HBOND_DHA_MIN:
                continue

            score = _competition_score(d_ha, angle_dha, dist_x_pi, angle_xh_pi)

            if best_score is None or score > best_score:
                best_score = score
                best_acceptor = (cra, d_ha, angle_dha, score)

        if best_acceptor is not None:
            cra, d_ha, angle_dha, score = best_acceptor
            hit["hbond_competing"] = 1
            hit["hbond_acceptor_atom"] = cra.atom.name
            hit["hbond_acceptor_res"] = cra.residue.name
            hit["hbond_acceptor_chain"] = cra.chain.name
            hit["hbond_HA_dist"] = round(d_ha, 3)
            hit["hbond_DHA_angle"] = round(angle_dha, 2)
            hit["hbond_vs_xhpi_score"] = round(score, 4)
        else:
            _set_no_competition(hit)

    return hits




def annotate_rigid_hbond(hits: List[Dict[str, Any]], structure: gemmi.Structure,
                         model_index: int = 0,
                         ns: Optional[gemmi.NeighborSearch] = None) -> List[Dict[str, Any]]:
    """Annotate rigid-donor hits with a binary H-bond flag only.

    For rigid donors (backbone N-H, Trp NE1, His, Arg, Asn/Gln sidechains),
    the hydrogen cannot rotate, so there is no "competition" — it can
    simultaneously participate in XH-π and a conventional H-bond.
    This function only records whether a H-bond acceptor is present.
    """
    if len(structure) == 0 or model_index >= len(structure):
        return hits

    if ns is None:
        model = structure[model_index]
        ns = gemmi.NeighborSearch(model, structure.cell, 5.0)
        ns.populate(include_h=False)
    else:
        model = structure[model_index]

    for hit in hits:
        # Exclude cone_virtual hits (they go through the competition path)
        if hit.get("H_source") == "cone_virtual":
            hit["hbond_also_hbonded"] = 0
            continue

        hx = hit.get("H_xyz_x")
        hy = hit.get("H_xyz_y")
        hz = hit.get("H_xyz_z")
        if hx is None:
            hit["hbond_also_hbonded"] = 0
            continue

        h_pos = gemmi.Position(hx, hy, hz)
        xx = hit.get("X_xyz_x")
        if xx is None:
            hit["hbond_also_hbonded"] = 0
            continue

        x_pos_arr = np.array([xx, hit.get("X_xyz_y", 0), hit.get("X_xyz_z", 0)])
        h_pos_arr = np.array([hx, hy, hz])

        donor_chain = hit.get("X_chain", "")
        donor_res = hit.get("X_res", "")
        donor_id = str(hit.get("X_id", ""))

        marks = ns.find_atoms(h_pos, radius=HBOND_SEARCH_RADIUS)
        found = False

        for mark in marks:
            cra = mark.to_cra(model)
            elem = cra.atom.element.name.upper()
            if elem not in HBOND_ACCEPTOR_ELEMENTS:
                continue
            if (cra.chain.name == donor_chain and
                cra.residue.name == donor_res and
                str(cra.residue.seqid).strip() == donor_id):
                continue

            a_pos_arr = np.array([mark.pos.x, mark.pos.y, mark.pos.z])
            d_ha = float(np.linalg.norm(a_pos_arr - h_pos_arr))
            if d_ha > HBOND_HA_MAX:
                continue

            vec_hd = x_pos_arr - h_pos_arr
            vec_ha = a_pos_arr - h_pos_arr
            norm_hd = np.linalg.norm(vec_hd)
            norm_ha = np.linalg.norm(vec_ha)
            if norm_hd == 0 or norm_ha == 0:
                continue

            cos_angle = np.dot(vec_hd, vec_ha) / (norm_hd * norm_ha)
            angle_dha = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
            if angle_dha >= HBOND_DHA_MIN:
                found = True
                break

        hit["hbond_also_hbonded"] = 1 if found else 0

    return hits

def _set_no_competition(hit: Dict[str, Any]) -> None:
    """Fill competition columns with fallback values."""
    hit["hbond_competing"] = 0
    hit["hbond_acceptor_atom"] = None
    hit["hbond_acceptor_res"] = None
    hit["hbond_acceptor_chain"] = None
    hit["hbond_HA_dist"] = None
    hit["hbond_DHA_angle"] = None
    hit["hbond_vs_xhpi_score"] = None
