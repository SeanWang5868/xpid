"""
core.py
Core logic for detecting XH-π interactions with Hudson/Plevin labels by default.
The P-slab system is retained as an opt-in method-development mode.
Track 1: Explicit H (Default)
Track 2: Implicit/Cone rescue (disabled by default; enable with use_cone=True)
"""
import gemmi
import logging
import numpy as np
from typing import List, Dict, Any, Optional, Union, Set, NamedTuple
from . import config
from . import geometry
from . import sasa
from . import cooperativity
from . import hbond
from . import rings
from . import systems
from . import hits as _hits
from . import ss

logger = logging.getLogger("xpid.core")

ResiduePairKeys = Optional[tuple[Set[tuple], Set[tuple]]]

def _pos_to_arr(pos: gemmi.Position) -> np.ndarray:
    """Convert gemmi.Position to numpy array without intermediate list."""
    return np.array([pos.x, pos.y, pos.z])

def _residue_match_key(chain_name: str, residue: gemmi.Residue) -> tuple:
    return (chain_name, str(residue.seqid).strip())


def _resolve_residue_selector_keys(model: gemmi.Model, selector_text: str) -> Set[tuple]:
    selection = gemmi.Selection(selector_text)
    keys: Set[tuple] = set()
    for chain in selection.chains(model):
        for residue in selection.residues(chain):
            keys.add(_residue_match_key(chain.name, residue))
    return keys


def _resolve_residue_pair_keys(model: gemmi.Model,
                               residue_pair: Optional[tuple[str, str]]) -> ResiduePairKeys:
    if not residue_pair:
        return None
    left, right = residue_pair
    return (
        _resolve_residue_selector_keys(model, left),
        _resolve_residue_selector_keys(model, right),
    )


def _residue_in_pair_side(pair_side: Set[tuple], chain_name: str, residue: gemmi.Residue) -> bool:
    return _residue_match_key(chain_name, residue) in pair_side


def _pair_allows_pi_residue(pair_keys: ResiduePairKeys,
                            chain_name: str, residue: gemmi.Residue) -> bool:
    if pair_keys is None:
        return True
    left, right = pair_keys
    return (
        _residue_in_pair_side(left, chain_name, residue) or
        _residue_in_pair_side(right, chain_name, residue)
    )


def _pair_allows_donor_residue(pair_keys: ResiduePairKeys,
                               pi_chain_name: str, pi_residue: gemmi.Residue,
                               donor_chain_name: str, donor_residue: gemmi.Residue) -> bool:
    if pair_keys is None:
        return True
    left, right = pair_keys
    pi_key = _residue_match_key(pi_chain_name, pi_residue)
    donor_key = _residue_match_key(donor_chain_name, donor_residue)
    if pi_key in left:
        return donor_key in right
    if pi_key in right:
        return donor_key in left
    return False

def _is_cation_pi_donor(res_name: str, atom_name: str) -> bool:
    return (res_name, atom_name) in config.CATION_DONORS
BLOCKING_METALS = {
    'ZN', 'FE', 'CU', 'MN', 'MG', 'CO', 'NI', 'CA', 'CD', 'HG',
    'NA', 'K', 'PT', 'AU', 'AG', 'FE2', 'FE3'
}

def detect_interactions_in_structure(structure: gemmi.Structure, 
                                     pdb_name: str,
                                     filter_pi: Optional[List[str]] = None,
                                     filter_donor: Optional[List[str]] = None,
                                     filter_donor_atom: Optional[List[str]] = None,
                                     model_mode: Union[str, int] = 0,
                                     cone_mode: str = "auto",  # "auto" (cone for rotatable) | "none" (explicit H only)
                                     include_p_slab: bool = False,
                                     report_xh_candidates: bool = False,
                                     include_coordinates: bool = False,
                                     residue_pair: Optional[tuple[str, str]] = None,
                                     min_occ: float = 0.0,
                                     external_ss_index: Optional[Dict] = None,
                                     sym_contacts: bool = False,
                                     include_water: bool = False,
                                     max_b: float = 0.0, compute_sasa: bool = False, annotate_cooperativity: bool = True) -> List[Dict[str, Any]]:
    results = []
    if not structure or len(structure) == 0: return []

    if not include_water:
        structure.remove_waters()
    structure.remove_empty_chains()

    models_with_ids = [] 
    if model_mode == 'all':
        for i, m in enumerate(structure):
            m_name = getattr(m, 'name', str(i+1))
            models_with_ids.append((m, m_name))
    else:
        try:
            idx = int(model_mode)
            if 0 <= idx < len(structure):
                m = structure[idx]
                m_name = getattr(m, 'name', str(idx+1))
                models_with_ids = [(m, m_name)]
            else: 
                return []
        except ValueError:
            m = structure[0]
            m_name = getattr(m, 'name', "1")
            models_with_ids = [(m, m_name)]

    resolution = structure.resolution if structure.resolution else 0.0

    ss_index = external_ss_index if external_ss_index else ss.build_index(structure)

    sasa_map = {}
    if compute_sasa:
        sasa_map = sasa.compute_sasa(structure)


    if sym_contacts:
        structure.setup_cell_images()

    for model, model_id in models_with_ids:
        ns = gemmi.NeighborSearch(model, structure.cell, config.DIST_SEARCH_LIMIT)
        ns.populate(include_h=True)
        pair_keys = _resolve_residue_pair_keys(model, residue_pair)

        for chain in model:
            for residue in chain:
                res_name = residue.name
                if filter_pi and res_name not in filter_pi: continue
                if not _pair_allows_pi_residue(pair_keys, chain.name, residue):
                    continue

                aromatic_rings = config.get_aromatic_rings(res_name)
                if not aromatic_rings: continue

                for ring_idx, target_atoms in enumerate(aromatic_rings):
                    ring_size = len(target_atoms)
                    mode = f'ring{ring_idx+1}'

                    for pi_alt, pi_atoms in rings._atom_variants_for_names(residue, target_atoms):
                        results.extend(_detect_residue(
                            pdb_name, resolution, model, model_id, chain, residue, ns, ss_index,
                            pi_atoms, pi_alt, mode, filter_donor, filter_donor_atom, cone_mode,
                            include_p_slab, report_xh_candidates, include_coordinates,
                            pair_keys, ring_size, sasa_map, min_occ,
                            sym_contacts=sym_contacts, max_b=max_b
                        ))
    hits = _hits._deduplicate_hits(results, prefer_directional=not report_xh_candidates)
    if annotate_cooperativity:
        hits = cooperativity.annotate_cooperativity(hits)
    # Annotate H-bond competition
    # Build a single NeighborSearch shared by both annotation passes
    _hbond_ns = gemmi.NeighborSearch(structure[0], structure.cell, 5.0)
    _hbond_ns.populate(include_h=False)
    hbond.annotate_rigid_hbond(hits, structure, ns=_hbond_ns)
    hbond.annotate_hbond_competition(hits, structure, ns=_hbond_ns)
    return hits

def _is_donor_blocked(x_atom: gemmi.Atom, model: gemmi.Model, ns: gemmi.NeighborSearch,
                      x_pos: Optional[gemmi.Position] = None) -> bool:
    radius = 2.6
    search_pos = x_pos if x_pos is not None else x_atom.pos
    neighbors = ns.find_atoms(search_pos, radius=radius)
    
    x_elem = x_atom.element.name.upper()
    
    for mark in neighbors:
        dist = mark.pos.dist(search_pos)
        if dist < 0.01: continue
        
        cra = mark.to_cra(model)
        neighbor_atom = cra.atom
        neighbor_el = neighbor_atom.element.name.upper()
        
        if x_elem == 'S' and neighbor_el == 'S' and neighbor_atom.name == 'SG':
            if 1.8 <= dist <= 2.2:
                return True
        
        if neighbor_el in BLOCKING_METALS:
            if dist <= 2.6:
                return True
    
    return False

def _detect_residue(pdb_name, resolution, model, model_id, chain, residue, ns, ss_index,
                    pi_atoms: List[gemmi.Atom], pi_alt: str, mode: str, filter_donor: Optional[List[str]],
                    filter_donor_atom: Optional[List[str]], cone_mode: str,
                    include_p_slab: bool, report_xh_candidates: bool,
                    include_coordinates: bool, pair_keys: ResiduePairKeys,
                    ring_size: int, sasa_map: Dict,
                    min_occ: float, sym_contacts: bool = False, max_b: float = 0.0):
    hits = []

    max_planar_dev = geometry.calculate_planarity_deviation(pi_atoms)
    if max_planar_dev > 0.5:
        return []

    pi_occs = [atom.occ for atom in pi_atoms]
    avg_pi_occ = sum(pi_occs) / len(pi_occs)
    if avg_pi_occ < 0.10:
        return []
    if max_b > 0 and any(atom.b_iso > max_b for atom in pi_atoms):
        return []

    pi_center, pi_center_arr, pi_normal, pi_b_mean = geometry.get_pi_info(pi_atoms)
    p_radius = config.P_RADIUS_BY_RING_SIZE.get(ring_size, config.P_RADIUS_BY_RING_SIZE[6])
    
    rctx = rings._RingContext(
        pdb_name=pdb_name, resolution=resolution, model=model, model_id=model_id,
        chain=chain, residue=residue, ns=ns, ss_index=ss_index,
        pi_center_arr=pi_center_arr, pi_normal=pi_normal, pi_b_mean=pi_b_mean,
        pi_alt=pi_alt, mode=mode, ring_size=ring_size, min_occ=min_occ,
        avg_pi_occ=avg_pi_occ, p_radius=p_radius,
        p_slab_half_thickness=config.P_SLAB_HALF_THICKNESS,
    )
    
    x_candidates = ns.find_atoms(pi_center, alt=pi_alt or "\0", radius=config.DIST_SEARCH_LIMIT)
    
    for x_mark in x_candidates:
        x_cra = x_mark.to_cra(model)
        x_atom = x_cra.atom
        x_res = x_cra.residue
        x_res_name = x_res.name
        x_elem = x_atom.element.name.upper()
        if x_elem not in config.TARGET_ELEMENTS_X:
            continue
        if not _pair_allows_donor_residue(pair_keys, chain.name, residue, x_cra.chain.name, x_res):
            continue
        
        is_sym_mate = (x_mark.image_idx != 0)
        if is_sym_mate and not sym_contacts:
            continue
        
        if filter_donor and x_res_name not in filter_donor: continue
        if filter_donor_atom and x_elem not in filter_donor_atom and x_atom.name.upper() not in filter_donor_atom:
            continue

        if x_cra.chain.name == chain.name and x_res.seqid == residue.seqid and x_res.name == residue.name:
            continue

        if (x_res_name in ('ASP', 'GLU') and x_atom.name in ('OD1', 'OD2', 'OE1', 'OE2')) or \
           (x_atom.name == 'OXT'):
            continue

        if _is_cation_pi_donor(x_res_name, x_atom.name):
            continue
        
        # Classify donor: rotatable vs rigid
        is_rotatable_donor = config.is_rotatable(x_res_name, x_atom.name)

        if is_sym_mate:
            if _is_donor_blocked(x_atom, model, ns, x_pos=x_mark.pos):
                continue
        else:
            if _is_donor_blocked(x_atom, model, ns):
                continue

        if x_atom.occ < 0.10: continue
        if max_b > 0 and x_atom.b_iso > max_b: continue

        if pi_alt and x_atom.altloc and pi_alt != x_atom.altloc:
            continue
        
        combined_occ = min(avg_pi_occ, x_atom.occ)

        if is_sym_mate:
            x_pos_arr = _pos_to_arr(x_mark.pos)
        else:
            x_pos_arr = _pos_to_arr(x_atom.pos)
        dist_x_pi = geometry.calculate_plane_distance(x_pos_arr, pi_center_arr, pi_normal)
        dist_x_centroid = geometry.calculate_distance(x_pos_arr, pi_center_arr)

        max_plane_dist = config.P_PLANE_DMAX.get(x_elem, config.P_PLANE_DMAX['default'])

        if dist_x_centroid > max_plane_dist and (
            dist_x_pi is None or dist_x_pi > max_plane_dist
        ):
            continue

        proj_dist = geometry.calculate_projection_dist(pi_normal, pi_center_arr, x_pos_arr)
        
        sym_op = x_mark.image_idx if is_sym_mate else 0

        # --- Detection path ---
        if is_rotatable_donor and cone_mode == "auto" and not report_xh_candidates:
            # Pre-filter: skip cone if even the closest possible H position
            # (X shifted towards ring by bond length) exceeds the element cutoff.
            _bond = config.BOND_LENGTHS.get(x_elem, 1.09)
            if dist_x_centroid - _bond > max_plane_dist:
                continue
            # Rotatable donor in auto mode: skip explicit-H track entirely.
            # Crystal-structure hydrogens on rotatable groups are riding
            # hydrogens with no experimental support at typical resolutions.
            # Cone rescue scans the full rotational space.
            _run_cone_track(
                rctx, x_cra, x_atom, x_mark, x_res, x_res_name, x_elem,
                x_pos_arr, is_sym_mate, dist_x_pi, dist_x_centroid, proj_dist,
                combined_occ, sym_op, include_p_slab, include_coordinates,
                sasa_map, hits, report_xh_candidates
            )
        else:
            # Rigid donor, or rotatable in --no-cone mode:
            # use explicit hydrogens from the structure file.
            _run_explicit_track(
                rctx, x_cra, x_atom, x_mark, x_pos_arr, is_sym_mate,
                x_elem, dist_x_pi, dist_x_centroid, proj_dist, combined_occ,
                sym_op, include_p_slab, report_xh_candidates, include_coordinates, sasa_map, hits
            )

    return hits





def _run_explicit_track(rctx: "rings._RingContext", x_cra, x_atom, x_mark,
                        x_pos_arr, is_sym_mate,
                        x_elem, dist_x_pi, dist_x_centroid, proj_dist,
                        combined_occ, sym_op, include_p_slab: bool,
                        report_xh_candidates: bool, include_coordinates: bool,
                        sasa_map: Dict,
                        hits) -> tuple:
    """Track 1: Explicit hydrogen geometry. Returns (found_systems, orig_h_positions)."""
    found_systems: Set[str] = set()
    orig_h_positions = []
    
    h_search_pos = x_mark.pos if is_sym_mate else x_atom.pos
    h_candidates = rctx.ns.find_atoms(h_search_pos, alt=x_atom.altloc, radius=config.DIST_CUTOFF_H)

    for h_mark in h_candidates:
        if is_sym_mate and h_mark.image_idx != x_mark.image_idx:
            continue

        h_cra = h_mark.to_cra(rctx.model)
        h_atom = h_cra.atom
        
        if h_atom.element.name.upper() not in {'H', 'D'}: 
            continue
                
        h_pos_arr = _pos_to_arr(h_mark.pos) if is_sym_mate else _pos_to_arr(h_atom.pos)
        orig_h_positions.append(h_pos_arr)

        xh_dist = float(np.linalg.norm(h_pos_arr - x_pos_arr))
        if xh_dist <= config.MIN_COVALENT_XH or xh_dist > config.DIST_CUTOFF_H:
            continue

        bonded_hydrogens = config.get_bonded_hydrogens(x_cra.residue.name, x_atom.name)
        if bonded_hydrogens and h_atom.name.upper() not in bonded_hydrogens:
            continue
        
        if h_atom.altloc and x_atom.altloc and h_atom.altloc != x_atom.altloc:
            continue
        
        h_combined_occ = min(combined_occ, h_atom.occ)
        metrics = systems._evaluate_systems(
            rctx, x_elem, x_pos_arr, h_pos_arr, dist_x_pi, dist_x_centroid, proj_dist)
        is_p_slab_hit = include_p_slab and metrics['is_p_slab']
        is_candidate_hit = (
            report_xh_candidates and
            (metrics['is_hudson_spatial'] or metrics['is_plevin_spatial'])
        )
        if not (metrics['is_hudson'] or metrics['is_plevin'] or is_p_slab_hit or is_candidate_hit):
            continue

        if metrics['is_hudson']:
            found_systems.add("hudson")
        if metrics['is_plevin']:
            found_systems.add("plevin")
        if is_p_slab_hit:
            found_systems.add("p_slab")

        _hits._record_hit(hits, rctx, x_cra, x_atom, h_atom.name, dist_x_pi,
                    dist_x_centroid, x_pos_arr, proj_dist, metrics,
                    is_cone=False, combined_occ=h_combined_occ, sym_op=sym_op,
                    h_atom=h_atom, include_p_slab=include_p_slab,
                    include_candidate_metrics=report_xh_candidates,
                    include_coordinates=include_coordinates,
                    sasa_map=sasa_map,
                    h_pos=h_pos_arr)
    
    return found_systems, orig_h_positions


def _run_cone_track(rctx: "rings._RingContext", x_cra, x_atom, x_mark, x_res, x_res_name, x_elem,
                    x_pos_arr, is_sym_mate, dist_x_pi, dist_x_centroid, proj_dist,
                    combined_occ, sym_op, include_p_slab: bool,
                    include_coordinates: bool, sasa_map: Dict, hits,
                    report_xh_candidates: bool = False):
    """Cone rescue: scan full rotational space for rotatable donors.

    Generates 72 virtual hydrogen positions around the X–parent bond axis,
    filters by steric clash and H-bond competition, then evaluates
    Hudson/Plevin/P-slab geometry on surviving conformers.
    """
    parent_name = config.ROTATABLE_MAPPING[x_res_name].get(x_atom.name)
    if not parent_name:
        return
    parent_atom = next((a for a in x_res if a.name == parent_name), None)
    if not parent_atom:
        return

    # Resolve parent position (sym mates need transformed coords)
    if is_sym_mate:
        parent_mark_candidates = rctx.ns.find_atoms(x_mark.pos, radius=2.0)
        parent_pos_arr = None
        for pm in parent_mark_candidates:
            if pm.image_idx == x_mark.image_idx:
                pm_cra = pm.to_cra(rctx.model)
                if (pm_cra.atom.name == parent_name and
                    pm_cra.residue.seqid == x_res.seqid and
                    pm_cra.chain.name == x_cra.chain.name):
                    parent_pos_arr = _pos_to_arr(pm.pos)
                    break
        if parent_pos_arr is None:
            return
    else:
        parent_pos_arr = _pos_to_arr(parent_atom.pos)

    # Gather environment atoms for steric clash check
    cone_search_pos = gemmi.Position(x_pos_arr[0], x_pos_arr[1], x_pos_arr[2])
    neighbors = rctx.ns.find_atoms(cone_search_pos, radius=4.0)

    env_coords_list = []
    for n_mark in neighbors:
        if n_mark.pos.dist(cone_search_pos) < 0.01:
            continue
        n_cra = n_mark.to_cra(rctx.model)
        if n_cra.residue.seqid == x_res.seqid and n_cra.chain.name == x_cra.chain.name:
            continue
        n_elem = n_cra.atom.element.name.upper()
        if n_elem in ('H', 'D', ''):
            continue
        n_pos_arr = _pos_to_arr(n_mark.pos)
        dist = np.linalg.norm(n_pos_arr - x_pos_arr)
        if dist <= 4.0:
            env_coords_list.append(n_pos_arr)

    env_coords = np.array(env_coords_list) if env_coords_list else np.empty((0, 3))

    # Generate 72 virtual H positions around parent→X axis (full 360°)
    h_candidates_cone = geometry.generate_rotated_hydrogens(
        parent_pos_arr, x_pos_arr, x_elem,
        env_coords=env_coords, clash_cutoff=2.0, num_samples=72
    )

    if not h_candidates_cone:
        return

    # --- H-bond competition gate ---
    # Filter out conformers where the H would be captured by a nearby acceptor.
    # Reuses the existing NeighborSearch (rctx.ns) instead of building a new one.
    from . import hbond as _hbond
    hbond_candidates = []

    for h_pos_np in h_candidates_cone:
        h_pos_gemmi = gemmi.Position(h_pos_np[0], h_pos_np[1], h_pos_np[2])
        blocked = False
        for a_mark in rctx.ns.find_atoms(h_pos_gemmi, radius=3.0):
            a_cra = a_mark.to_cra(rctx.model)
            a_elem = a_cra.atom.element.name.upper()
            if a_elem not in ('O', 'N', 'S'):
                continue
            if (a_cra.chain.name == x_cra.chain.name and
                a_cra.residue.seqid == x_res.seqid):
                continue
            a_pos_arr = _pos_to_arr(a_mark.pos)
            d_ha = float(np.linalg.norm(a_pos_arr - h_pos_np))
            if d_ha > 2.8:
                continue
            vec_hd = x_pos_arr - h_pos_np
            vec_ha = a_pos_arr - h_pos_np
            norm_hd = np.linalg.norm(vec_hd)
            norm_ha = np.linalg.norm(vec_ha)
            if norm_hd == 0 or norm_ha == 0:
                continue
            cos_angle = np.dot(vec_hd, vec_ha) / (norm_hd * norm_ha)
            angle_dha = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
            if angle_dha < 120.0:
                continue
            score = _hbond._competition_score(d_ha, angle_dha, dist_x_pi if dist_x_pi is not None else 4.5, None)
            if score > 0.2:
                blocked = True
                break
        if not blocked:
            hbond_candidates.append(h_pos_np)

    if not hbond_candidates:
        return

    # --- Geometric evaluation on surviving conformers ---
    legacy_best = None
    legacy_best_score = None
    legacy_best_h_pos = None
    p_slab_best = None
    p_slab_best_score = None
    p_slab_best_h_pos = None

    h_pos_array = np.array(hbond_candidates)  # (N, 3)
    batch_metrics = systems._evaluate_systems_batch(
        rctx, x_elem, x_pos_arr, h_pos_array, dist_x_pi, dist_x_centroid, proj_dist)

    for i, metrics in enumerate(batch_metrics):
        h_pos_np = hbond_candidates[i]

        if metrics['is_hudson'] or metrics['is_plevin']:
            legacy_score = metrics['angle_XH_Pi'] if metrics['angle_XH_Pi'] is not None else -1.0
            if legacy_best_score is None or legacy_score > legacy_best_score:
                legacy_best_score = legacy_score
                legacy_best = metrics
                legacy_best_h_pos = h_pos_np

        if include_p_slab and metrics['is_p_slab']:
            h_proj_score = metrics['h_proj_dist'] if metrics['h_proj_dist'] is not None else 999.0
            h_ray_score = metrics['H_ray_t'] if metrics['H_ray_t'] is not None else 999.0
            p_slab_score = (-h_proj_score, -h_ray_score)
            if p_slab_best_score is None or p_slab_score > p_slab_best_score:
                p_slab_best_score = p_slab_score
                p_slab_best = metrics
                p_slab_best_h_pos = h_pos_np

    combined = _hits._merge_cone_system_metrics(legacy_best, p_slab_best, True, include_p_slab)
    if combined is not None:
        selected_h_pos = legacy_best_h_pos if legacy_best_h_pos is not None else p_slab_best_h_pos
        _hits._record_hit(hits, rctx, x_cra, x_atom, "virt", dist_x_pi,
                    dist_x_centroid, x_pos_arr, proj_dist, combined,
                    is_cone=True, combined_occ=combined_occ, sym_op=sym_op,
                    h_atom=None, include_p_slab=include_p_slab,
                    include_coordinates=include_coordinates,
                    sasa_map=sasa_map,
                    h_pos=selected_h_pos)


