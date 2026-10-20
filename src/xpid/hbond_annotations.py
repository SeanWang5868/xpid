"""Conventional hydrogen-bond geometry consistent with a detected conformer.

Annotations are descriptive and never create or veto an XH–π detection.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np

from . import cone, config, hbond_acceptors

HBOND_SEARCH_RADIUS = 3.0
HBOND_HA_MAX = cone.HBOND_CONTACT_HA_MAX
HBOND_DHA_MIN = cone.HBOND_CONTACT_DHA_MIN


def _nearest_mark_position(
        ns: gemmi.NeighborSearch, reference: gemmi.Position,
        mark: gemmi.NeighborSearch.Mark,
        base_position: gemmi.Position) -> Tuple[gemmi.Position, str]:
    position = ns.grid_cell.find_nearest_pbc_position(
        reference, base_position, mark.image_idx)
    nearest = ns.grid_cell.find_nearest_pbc_image(
        reference, base_position, mark.image_idx)
    return position, nearest.symmetry_code()


def _competition_score(d_ha: float, angle_dha: float,
                       dist_x_pi: Optional[float] = None,
                       angle_xh_pi: Optional[float] = None) -> float:
    """Legacy API: pure geometry score, not an energy or XH–π competition."""
    return (1.0 - d_ha / HBOND_HA_MAX) * angle_dha / 180.0


def _hbond_contacts(hit: Dict[str, Any], model: gemmi.Model,
                    ns: gemmi.NeighborSearch, min_occ: float, max_b: float,
                    protonation=None, model_index=0):
    """Yield acceptors in the same occupied alternate state as X, H and π."""
    h_xyz = [hit.get(f"H_xyz_{axis}") for axis in "xyz"]
    x_xyz = [hit.get(f"X_xyz_{axis}") for axis in "xyz"]
    if any(value is None for value in h_xyz + x_xyz):
        return
    h_arr, x_arr = np.asarray(h_xyz), np.asarray(x_xyz)
    h_pos = gemmi.Position(*h_xyz)
    hit_occ = hit.get("combined_occupancy", 1.0)
    if hit_occ < config.MIN_ATOM_OCCUPANCY or hit_occ < min_occ:
        return
    labels = [hit.get(field, "") for field in ("pi_altloc", "X_altloc", "H_altloc")]
    active_alt = next((label for label in labels if label not in ("", "\0", None)), "")
    if any(not hbond_acceptors.altlocs_compatible(active_alt, label) for label in labels):
        return
    vec_hd = x_arr - h_arr
    norm_hd = np.linalg.norm(vec_hd)
    if norm_hd == 0:
        return
    for mark in ns.find_atoms(h_pos, alt=active_alt or "\0", radius=HBOND_SEARCH_RADIUS):
        cra = mark.to_cra(model)
        atom = cra.atom
        if atom.occ < config.MIN_ATOM_OCCUPANCY or min(hit_occ, atom.occ) < min_occ:
            continue
        if max_b > 0 and atom.b_iso > max_b:
            continue
        if not hbond_acceptors.altlocs_compatible(active_alt, atom.altloc):
            continue
        acceptor = (protonation.is_acceptor(model_index, cra.chain.name, cra.residue, atom, active_alt)
                    if protonation is not None else
                    hbond_acceptors.is_hbond_acceptor(cra.residue, atom, altloc=active_alt))
        if not acceptor:
            continue
        # Resolve around the H being annotated, not around the ring centroid.
        position, code = _nearest_mark_position(ns, h_pos, mark, atom.pos)
        if (code == hit.get("symmetry_code", "1_555")
                and cra.chain.name == hit.get("X_chain", "")
                and cra.residue.name == hit.get("X_res", "")
                and str(cra.residue.seqid).strip() == str(hit.get("X_id", ""))):
            continue
        vec_ha = np.array([position.x, position.y, position.z]) - h_arr
        distance = float(np.linalg.norm(vec_ha))
        if not cone.ABSOLUTE_MIN_HA <= distance <= HBOND_HA_MAX:
            continue
        cosine = np.dot(vec_hd, vec_ha) / (norm_hd * distance)
        angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
        if angle >= HBOND_DHA_MIN:
            yield cra, code, distance, angle, _competition_score(distance, angle)


def _context_model(structure, model_index, ns):
    if model_index < 0 or model_index >= len(structure):
        return None, None
    model = structure[model_index]
    if ns is None:
        ns = gemmi.NeighborSearch(model, structure.cell, 5.0)
        ns.populate(include_h=False)
    return model, ns


def annotate_cone_hbond_context(
        hits: List[Dict[str, Any]], structure: gemmi.Structure,
        model_index: int = 0, ns: Optional[gemmi.NeighborSearch] = None,
        min_occ: float = 0.0, max_b: float = 0.0, protonation=None) -> List[Dict[str, Any]]:
    """Report the strongest nearby H-bond geometry for every selected H."""
    model, ns = _context_model(structure, model_index, ns)
    if model is None:
        return hits
    for hit in hits:
        best = max(_hbond_contacts(hit, model, ns, min_occ, max_b, protonation, model_index),
                   key=lambda item: item[4], default=None)
        if best is None:
            _set_no_competition(hit)
            continue
        cra, code, distance, angle, score = best
        hit.update({
            "hbond_competing": 1,
            "hbond_acceptor_atom": cra.atom.name,
            "hbond_acceptor_res": cra.residue.name,
            "hbond_acceptor_chain": cra.chain.name,
            "hbond_acceptor_id": str(cra.residue.seqid).strip(),
            "hbond_acceptor_altloc": "" if cra.atom.altloc == "\0" else cra.atom.altloc,
            "hbond_acceptor_symmetry_code": code,
            "hbond_HA_dist": round(distance, 3),
            "hbond_DHA_angle": round(angle, 2),
            "hbond_vs_xhpi_score": round(score, 4),
        })
    return hits


def annotate_explicit_hbond_context(
        hits: List[Dict[str, Any]], structure: gemmi.Structure,
        model_index: int = 0, ns: Optional[gemmi.NeighborSearch] = None,
        min_occ: float = 0.0, max_b: float = 0.0, protonation=None) -> List[Dict[str, Any]]:
    """Report whether an explicit H also has conventional H-bond geometry."""
    model, ns = _context_model(structure, model_index, ns)
    if model is None:
        return hits
    for hit in hits:
        hit["hbond_also_hbonded"] = int(
            hit.get("H_source") != "cone_virtual"
            and next(_hbond_contacts(hit, model, ns, min_occ, max_b, protonation, model_index), None) is not None)
    return hits


def _set_no_competition(hit: Dict[str, Any]) -> None:
    hit["hbond_competing"] = 0
    for field in ("hbond_acceptor_atom", "hbond_acceptor_res", "hbond_acceptor_chain",
                  "hbond_acceptor_id", "hbond_acceptor_altloc", "hbond_acceptor_symmetry_code",
                  "hbond_HA_dist", "hbond_DHA_angle", "hbond_vs_xhpi_score"):
        hit[field] = None
