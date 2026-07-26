"""
hits.py
Hit recording, deduplication, and helper functions.
"""
import gemmi
import logging
import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Set

from . import config
from . import geometry
from . import ss
from . import sasa
from . import rings
from . import systems
from . import schema


def _round_float(value, ndigits: int):
    return round(float(value), ndigits)



def _score_float(hit: Dict[str, Any], key: str, default: float = 999.0) -> float:
    value = hit.get(key)
    return default if value is None else float(value)



def _dedup_key(hit: Dict[str, Any]) -> tuple:
    return (
        hit.get("pdb"), hit.get("model"), hit.get("_pi_ring_key"),
        hit.get("pi_chain"), hit.get("pi_res"), hit.get("pi_id"),
        hit.get("X_chain"), hit.get("X_res"), hit.get("X_id"),
        hit.get("X_atom"), hit.get("H_atom"), hit.get("sym_op"),
    )



def _dedup_score(hit: Dict[str, Any]) -> tuple:
    return (
        hit.get("_combined_occ", 0.0),
        int(hit.get("is_hudson", 0)) + int(hit.get("is_plevin", 0)) + int(hit.get("is_p_slab", 0)),
        int(hit.get("is_p_slab", 0)),
        int(hit.get("is_hudson", 0)) + int(hit.get("is_plevin", 0)),
        -_score_float(hit, "dist_X_Pi"),
        -_score_float(hit, "h_proj_dist"),
        -_score_float(hit, "proj_dist"),
    )



def _candidate_dedup_score(hit: Dict[str, Any]) -> tuple:
    return (
        hit.get("_combined_occ", 0.0),
        -_score_float(hit, "dist_X_Pi"),
        -_score_float(hit, "dist_X_centroid"),
    )



def _candidate_dedup_key(hit: Dict[str, Any]) -> tuple:
    return _dedup_key(hit) + (
        hit.get("_pi_alt", ""),
        hit.get("_X_alt", ""),
        hit.get("_H_alt", ""),
    )



def _h_source(h_atom: Optional[gemmi.Atom], is_cone: bool) -> str:
    if is_cone:
        return "cone_virtual"
    if h_atom is not None and h_atom.flag == "G":
        return "added"
    return "experimental"



def _deduplicate_hits(hits: List[Dict[str, Any]], prefer_directional: bool = True) -> List[Dict[str, Any]]:
    selected: Dict[tuple, Dict[str, Any]] = {}
    order: List[tuple] = []
    score_fn = _dedup_score if prefer_directional else _candidate_dedup_score
    key_fn = _dedup_key if prefer_directional else _candidate_dedup_key

    for hit in hits:
        key = key_fn(hit)
        if key not in selected:
            selected[key] = hit
            order.append(key)
            continue
        if score_fn(hit) > score_fn(selected[key]):
            selected[key] = hit

    deduped = []
    for key in order:
        hit = selected[key]
        hit.pop("_combined_occ", None)
        hit.pop("_pi_ring_key", None)
        hit.pop("_pi_alt", None)
        hit.pop("_X_alt", None)
        hit.pop("_H_alt", None)
        deduped.append(hit)
    return deduped



def _record_hit(hits: List[Dict[str, Any]], rctx: "rings._RingContext",
                x_cra, x_atom, h_name: str, dist: float,
                dist_x_centroid: float, x_pos: np.ndarray, proj: Optional[float],
                metrics: Dict[str, Any],
                is_cone: bool = False,
                combined_occ: float = 1.0, sym_op: int = 0,
                h_atom: Optional[gemmi.Atom] = None,
                include_p_slab: bool = False,
                include_candidate_metrics: bool = False,
                include_coordinates: bool = False,
                sasa_map: Optional[Dict] = None,
                h_pos: Optional[np.ndarray] = None,
                hbond_relation: str = "none"):
    
    if combined_occ < rctx.min_occ:
        return
    
    pi_ss_type, pi_ss_uid = ss.get_info(rctx.chain.name, rctx.residue.seqid.num, rctx.ss_index)
    x_ss_type, x_ss_uid = ss.get_info(x_cra.chain.name, x_cra.residue.seqid.num, rctx.ss_index)

    seq_sep = 0
    if rctx.chain.name == x_cra.chain.name:
        seq_sep = rctx.residue.seqid.num - x_cra.residue.seqid.num

    h_source = _h_source(h_atom, is_cone)
    is_trp_5ring = int(rctx.residue.name == 'TRP' and rctx.ring_size == 5)
    is_pi_pi_tshaped = 0

    donor_rings = config.get_aromatic_rings(x_cra.residue.name)
    for d_ring_atoms in donor_rings:
        d_pi_atoms = [a for a in x_cra.residue if a.name in d_ring_atoms]
        if len(d_pi_atoms) != len(d_ring_atoms):
            continue
        _, d_center, d_normal, _ = geometry.get_pi_info(d_pi_atoms)
        pp_dist, pp_angle, _ = geometry.calculate_pi_pi_geometry(
            rctx.pi_center_arr, rctx.pi_normal, d_center, d_normal)
        if 3.0 <= pp_dist <= config.PI_PI_DIST_MAX and pp_angle >= config.PI_PI_ANGLE_TSHAPED_MIN:
            is_pi_pi_tshaped = 1
            break

    hit = {
        'pdb': rctx.pdb_name,
        'model': rctx.model_id,
        'resolution': rctx.resolution,
        'pi_chain': rctx.chain.name,
        'pi_res': rctx.residue.name,
        'pi_id': str(rctx.residue.seqid),
        'pi_ring_id': rctx.mode,
        'pi_ring_size': rctx.ring_size,
        'X_chain': x_cra.chain.name,
        'X_res': x_cra.residue.name,
        'X_id': str(x_cra.residue.seqid),
        'X_atom': x_atom.name,
        'H_atom': h_name,
        'dist_X_Pi': _round_float(dist, 3),
        'dist_X_centroid': _round_float(dist_x_centroid, 3),
        'H_source': h_source,
        'is_hudson': metrics['is_hudson'],
        'is_plevin': metrics['is_plevin'],
        'is_trp_5ring_acceptor': is_trp_5ring,
        'is_pi_pi_tshaped': is_pi_pi_tshaped,
        'pi_ss_type': pi_ss_type,
        'pi_ss_id': pi_ss_uid,
        'X_ss_type': x_ss_type,
        'X_ss_id': x_ss_uid,
        'pi_avg_b': _round_float(rctx.pi_b_mean, 2),
        'pi_center_x': _round_float(rctx.pi_center_arr[0], 3),
        'pi_center_y': _round_float(rctx.pi_center_arr[1], 3),
        'pi_center_z': _round_float(rctx.pi_center_arr[2], 3),
        'X_b': _round_float(x_atom.b_iso, 2),
        'X_xyz_x': _round_float(x_pos[0], 3),
        'X_xyz_y': _round_float(x_pos[1], 3),
        'X_xyz_z': _round_float(x_pos[2], 3),
        'seq_sep': seq_sep,
        'proj_dist': _round_float(proj, 3) if proj is not None else None,
        'theta': _round_float(metrics['theta'], 2) if metrics['theta'] is not None else None,
        'angle_XPCN': _round_float(metrics['angle_XPCN'], 2) if metrics['angle_XPCN'] is not None else None,
        'angle_XH_Pi': _round_float(metrics['angle_XH_Pi'], 2) if metrics['angle_XH_Pi'] is not None else None,
        'sym_op': sym_op,
        '_combined_occ': float(combined_occ),
        '_pi_ring_key': rctx.mode,
        '_pi_alt': rctx.pi_alt,
        '_X_alt': rings._altloc(x_atom),
        '_H_alt': rings._altloc(h_atom) if h_atom is not None else "",
        'hbond_relation': hbond_relation,
    }

    if include_candidate_metrics:
        hit.update({
            'is_xh_candidate': 1,
            'is_hudson_spatial': metrics['is_hudson_spatial'],
            'is_plevin_spatial': metrics['is_plevin_spatial'],
            'hudson_dist_ok': metrics['hudson_dist_ok'],
            'hudson_proj_ok': metrics['hudson_proj_ok'],
            'hudson_direction_ok': metrics['hudson_direction_ok'],
            'plevin_dist_ok': metrics['plevin_dist_ok'],
            'plevin_xpcn_ok': metrics['plevin_xpcn_ok'],
            'plevin_direction_ok': metrics['plevin_direction_ok'],
            'xh_centroid_cos': _round_float(metrics['xh_centroid_cos'], 4) if metrics['xh_centroid_cos'] is not None else None,
            'xh_lateral_inward_score': _round_float(metrics['xh_lateral_inward_score'], 4) if metrics['xh_lateral_inward_score'] is not None else None,
        })

    # Always store H coordinates (needed for hbond competition)
    hit.update({
        'H_xyz_x': _round_float(h_pos[0], 3) if h_pos is not None else None,
        'H_xyz_y': _round_float(h_pos[1], 3) if h_pos is not None else None,
        'H_xyz_z': _round_float(h_pos[2], 3) if h_pos is not None else None,
    })

    if include_coordinates:
        canonical_normal = systems._canonical_unit_normal(rctx.pi_normal)
        hit.update({
            'pi_normal_x': _round_float(canonical_normal[0], 6) if canonical_normal is not None else None,
            'pi_normal_y': _round_float(canonical_normal[1], 6) if canonical_normal is not None else None,
            'pi_normal_z': _round_float(canonical_normal[2], 6) if canonical_normal is not None else None,
            'X_side_of_pi': systems._side_of_plane(x_pos, rctx.pi_center_arr, canonical_normal) if canonical_normal is not None else 0,
        })

    if sasa_map:
        pi_ring_indices = sasa.residue_atom_indices(rctx.residue)
        hit['pi_sasa_avg'] = _round_float(
            sasa.average_ring_sasa(sasa_map, 0, rctx.chain.name, rctx.residue, pi_ring_indices), 2)
        x_atom_idx = None
        for idx, atom in enumerate(x_cra.residue):
            if atom is x_atom:
                x_atom_idx = idx
                break
        hit['X_sasa'] = _round_float(
            sasa.atom_sasa(sasa_map, 0, x_cra.chain.name, str(x_cra.residue.seqid).strip(), x_atom_idx), 2
        ) if x_atom_idx is not None else None
        if h_atom is not None:
            h_atom_idx = None
            for idx, atom in enumerate(x_cra.residue):
                if atom is h_atom:
                    h_atom_idx = idx
                    break
            hit['H_sasa'] = _round_float(
                sasa.atom_sasa(sasa_map, 0, x_cra.chain.name, str(x_cra.residue.seqid).strip(), h_atom_idx), 2
            ) if h_atom_idx is not None else None
        else:
            hit['H_sasa'] = None


    if include_p_slab or include_candidate_metrics:
        p_geometry = {
            'P_radius': _round_float(rctx.p_radius, 3),
            'P_slab_half_thickness': _round_float(rctx.p_slab_half_thickness, 3),
            'h_proj_dist': _round_float(metrics['h_proj_dist'], 3) if metrics['h_proj_dist'] is not None else None,
            'H_ray_t': _round_float(metrics['H_ray_t'], 3) if metrics['H_ray_t'] is not None else None,
            'H_ray_entry_dist': _round_float(metrics['H_ray_entry_dist'], 3) if metrics['H_ray_entry_dist'] is not None else None,
            'h_plane_proj_dist': _round_float(metrics['h_plane_proj_dist'], 3) if metrics['h_plane_proj_dist'] is not None else None,
            'H_plane_t': _round_float(metrics['H_plane_t'], 3) if metrics['H_plane_t'] is not None else None,
            'H_plane_entry_dist': _round_float(metrics['H_plane_entry_dist'], 3) if metrics['H_plane_entry_dist'] is not None else None,
            'delta_h_proj_dist': _round_float(metrics['delta_h_proj_dist'], 3) if metrics['delta_h_proj_dist'] is not None else None,
        }
        if include_p_slab:
            p_geometry['is_p_slab'] = metrics['is_p_slab']
        hit.update(p_geometry)

    hits.append(hit)


def _merge_cone_system_metrics(legacy_best: Optional[Dict[str, Any]],
                               p_slab_best: Optional[Dict[str, Any]],
                               needs_legacy: bool,
                               needs_p_slab: bool) -> Optional[Dict[str, Any]]:
    """Merge independently selected cone candidates into one reported virtual hit."""
    if legacy_best is None and p_slab_best is None:
        return None

    base = dict(legacy_best if legacy_best is not None else p_slab_best)

    base['is_hudson'] = int(needs_legacy and legacy_best is not None and legacy_best['is_hudson'])
    base['is_plevin'] = int(needs_legacy and legacy_best is not None and legacy_best['is_plevin'])
    base['is_p_slab'] = int(needs_p_slab and p_slab_best is not None and p_slab_best['is_p_slab'])

    if p_slab_best is not None:
        base['h_proj_dist'] = p_slab_best['h_proj_dist']
        base['H_ray_t'] = p_slab_best['H_ray_t']
        base['H_ray_entry_dist'] = p_slab_best['H_ray_entry_dist']
        base['h_plane_proj_dist'] = p_slab_best['h_plane_proj_dist']
        base['H_plane_t'] = p_slab_best['H_plane_t']
        base['H_plane_entry_dist'] = p_slab_best['H_plane_entry_dist']
        base['delta_h_proj_dist'] = p_slab_best['delta_h_proj_dist']

    if legacy_best is not None:
        base['theta'] = legacy_best['theta']
        base['angle_XPCN'] = legacy_best['angle_XPCN']
        base['angle_XH_Pi'] = legacy_best['angle_XH_Pi']

    return base


def _pos_to_arr(pos: gemmi.Position) -> np.ndarray:
    """Convert gemmi.Position to numpy array without intermediate list."""
    return np.array([pos.x, pos.y, pos.z])
