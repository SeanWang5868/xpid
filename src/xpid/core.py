"""
core.py
Core logic for detecting XH-π interactions with Hudson, Plevin, and P-slab labels.
Track 1: Explicit H (Default)
Track 2: Implicit/Cone rescue (enabled by default via use_cone=True)
"""
import gemmi
import logging
import numpy as np
from typing import List, Dict, Any, Optional, Union, Set, NamedTuple
from . import config
from . import geometry
from . import ss

logger = logging.getLogger("xpid.core")

def _pos_to_arr(pos: gemmi.Position) -> np.ndarray:
    """Convert gemmi.Position to numpy array without intermediate list."""
    return np.array([pos.x, pos.y, pos.z])


def _round_float(value, ndigits: int):
    return round(float(value), ndigits)


def _score_float(hit: Dict[str, Any], key: str, default: float = 999.0) -> float:
    value = hit.get(key)
    return default if value is None else float(value)


def _altloc(atom: gemmi.Atom) -> str:
    return "" if atom.altloc in ("", "\0") else atom.altloc


def _sorted_altlocs(altlocs: Set[str]) -> List[str]:
    return sorted(altlocs, key=lambda alt: (alt != "", alt))


def _atom_variants_for_names(residue, atom_names: Set[str]) -> List[tuple]:
    """Return complete atom-name sets for each compatible altloc state."""
    atoms_by_name = {name: [] for name in atom_names}
    for atom in residue:
        if atom.name in atoms_by_name:
            atoms_by_name[atom.name].append(atom)

    if any(not atoms for atoms in atoms_by_name.values()):
        return []

    altlocs = {""}
    for atoms in atoms_by_name.values():
        altlocs.update(_altloc(atom) for atom in atoms if _altloc(atom))

    variants = []
    seen = set()
    for alt in _sorted_altlocs(altlocs):
        selected = []
        for name in sorted(atom_names):
            atoms = atoms_by_name[name]
            exact = [atom for atom in atoms if _altloc(atom) == alt]
            blank = [atom for atom in atoms if _altloc(atom) == ""]
            if alt:
                atom = exact[0] if exact else (blank[0] if blank else None)
            else:
                atom = blank[0] if blank else None
            if atom is None:
                selected = []
                break
            selected.append(atom)

        if selected:
            signature = tuple(id(atom) for atom in selected)
            if signature not in seen:
                seen.add(signature)
                variants.append((alt, selected))
    return variants


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


def _h_source(h_atom: Optional[gemmi.Atom], is_cone: bool) -> str:
    if is_cone:
        return "cone_virtual"
    if h_atom is not None and h_atom.flag == "G":
        return "added"
    return "experimental"


def _is_cation_pi_donor(res_name: str, atom_name: str) -> bool:
    return (res_name, atom_name) in config.CATION_DONORS


def _evaluate_systems(rctx: "_RingContext", x_elem: str, x_pos: np.ndarray, h_pos: np.ndarray,
                      dist_x_plane: Optional[float], dist_x_centroid: float,
                      proj_dist: Optional[float]) -> Dict[str, Any]:
    """Calculate Hudson, Plevin, and P-slab metrics for one X-H candidate."""
    theta = geometry.calculate_hudson_theta(rctx.pi_center_arr, x_pos, h_pos, rctx.pi_normal)
    xpcn_angle = geometry.calculate_xpcn_angle(x_pos, rctx.pi_center_arr, rctx.pi_normal)
    xh_pi_angle = geometry.calculate_xh_picenter_angle(rctx.pi_center_arr, x_pos, h_pos)

    p_dmax = config.P_PLANE_DMAX.get(x_elem, config.P_PLANE_DMAX['default'])
    h_proj_dist = None
    h_ray_t = None
    slab_entry = geometry.calculate_xh_ray_p_slab_entry(
        x_pos, h_pos, rctx.pi_center_arr, rctx.pi_normal, rctx.p_slab_half_thickness)
    if slab_entry is not None:
        h_hit_pos, h_ray_t = slab_entry
        h_proj_dist = geometry.calculate_projection_dist(rctx.pi_normal, rctx.pi_center_arr, h_hit_pos)

    is_hudson = int(
        dist_x_centroid <= p_dmax and
        proj_dist is not None and proj_dist <= rctx.p_radius and
        theta is not None and theta <= config.HUDSON_THETA_MAX
    )
    is_plevin = int(
        dist_x_centroid < p_dmax and
        xpcn_angle is not None and xpcn_angle < config.PLEVIN_XPCN_MAX and
        xh_pi_angle is not None and xh_pi_angle >= config.PLEVIN_XH_PI_MIN
    )
    is_p_slab = int(
        dist_x_plane is not None and dist_x_plane <= p_dmax and
        proj_dist is not None and proj_dist <= rctx.p_radius and
        h_proj_dist is not None and h_proj_dist <= rctx.p_radius
    )

    return {
        'dist_X_centroid': dist_x_centroid,
        'theta': theta,
        'angle_XPCN': xpcn_angle,
        'angle_XH_Pi': xh_pi_angle,
        'h_proj_dist': h_proj_dist,
        'H_ray_t': h_ray_t,
        'is_hudson': is_hudson,
        'is_plevin': is_plevin,
        'is_p_slab': is_p_slab,
    }


def _deduplicate_hits(hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    selected: Dict[tuple, Dict[str, Any]] = {}
    order: List[tuple] = []

    for hit in hits:
        key = _dedup_key(hit)
        if key not in selected:
            selected[key] = hit
            order.append(key)
            continue
        if _dedup_score(hit) > _dedup_score(selected[key]):
            selected[key] = hit

    deduped = []
    for key in order:
        hit = selected[key]
        hit.pop("_combined_occ", None)
        hit.pop("_pi_ring_key", None)
        deduped.append(hit)
    return deduped


class _RingContext(NamedTuple):
    """Immutable context for one aromatic ring being analyzed."""
    pdb_name: str
    resolution: float
    model: gemmi.Model
    model_id: str
    chain: Any  # gemmi.Chain
    residue: Any  # gemmi.Residue
    ns: gemmi.NeighborSearch
    ss_index: Dict
    pi_center_arr: np.ndarray
    pi_normal: np.ndarray
    pi_b_mean: float
    pi_alt: str
    mode: str
    ring_size: int
    min_occ: float
    avg_pi_occ: float
    p_radius: float
    p_slab_half_thickness: float

BLOCKING_METALS = {
    'ZN', 'FE', 'CU', 'MN', 'MG', 'CO', 'NI', 'CA', 'CD', 'HG',
    'NA', 'K', 'PT', 'AU', 'AG', 'FE2', 'FE3'
}

def select_best_altconf(structure: gemmi.Structure):
    """Select highest-occupancy altconf per residue; if tied, prefer alphabetically first (usually 'A')."""
    for model in structure:
        for chain in model:
            for residue in chain:
                altlocs = set()
                for atom in residue:
                    if atom.altloc != '\0':
                        altlocs.add(atom.altloc)
                if not altlocs:
                    continue
                if len(altlocs) == 1:
                    for atom in residue:
                        if atom.altloc != '\0':
                            atom.altloc = '\0'
                    continue
                occ_sum = {alt: 0.0 for alt in altlocs}
                occ_cnt = {alt: 0 for alt in altlocs}
                for atom in residue:
                    if atom.altloc in altlocs:
                        occ_sum[atom.altloc] += atom.occ
                        occ_cnt[atom.altloc] += 1
                avg_occ = {alt: occ_sum[alt] / occ_cnt[alt] if occ_cnt[alt] > 0 else 0.0
                           for alt in altlocs}
                best = min(altlocs, key=lambda x: (-avg_occ[x], x))
                to_remove = []
                for i in range(len(residue)):
                    atom = residue[i]
                    if atom.altloc != '\0' and atom.altloc != best:
                        to_remove.append(i)
                    elif atom.altloc == best:
                        atom.altloc = '\0'
                for i in reversed(to_remove):
                    del residue[i]

def detect_interactions_in_structure(structure: gemmi.Structure, 
                                     pdb_name: str,
                                     filter_pi: Optional[List[str]] = None,
                                     filter_donor: Optional[List[str]] = None,
                                     filter_donor_atom: Optional[List[str]] = None,
                                     model_mode: Union[str, int] = 0,
                                     use_cone: bool = True,
                                     min_occ: float = 0.0,
                                     external_ss_index: Optional[Dict] = None,
                                     sym_contacts: bool = False,
                                     include_water: bool = False,
                                     max_b: float = 0.0) -> List[Dict[str, Any]]:
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

    if sym_contacts:
        structure.setup_cell_images()

    for model, model_id in models_with_ids:
        ns = gemmi.NeighborSearch(model, structure.cell, config.DIST_SEARCH_LIMIT)
        ns.populate(include_h=True)

        for chain in model:
            for residue in chain:
                res_name = residue.name
                if filter_pi and res_name not in filter_pi: continue

                rings = config.get_aromatic_rings(res_name)
                if not rings: continue

                for ring_idx, target_atoms in enumerate(rings):
                    ring_size = len(target_atoms)
                    mode = f'ring{ring_idx+1}'

                    for pi_alt, pi_atoms in _atom_variants_for_names(residue, target_atoms):
                        results.extend(_detect_residue(
                            pdb_name, resolution, model, model_id, chain, residue, ns, ss_index,
                            pi_atoms, pi_alt, mode, filter_donor, filter_donor_atom, use_cone,
                            ring_size, min_occ, sym_contacts=sym_contacts, max_b=max_b
                        ))
    return _deduplicate_hits(results)

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
                    filter_donor_atom: Optional[List[str]], use_cone: bool, ring_size: int,
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
    
    rctx = _RingContext(
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
            
        allow_cone_scan = True
        if (x_res_name == 'ARG' and x_atom.name in ('NH1', 'NH2', 'NE')) or \
           (x_res_name in ('ASN', 'GLN') and x_atom.name in ('ND2', 'NE2')) or \
           (x_res_name == 'HIS' and x_atom.name in ('ND1', 'NE2')) or \
           (x_res_name == 'TRP' and x_atom.name == 'NE1'):
            allow_cone_scan = False

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

        # --- Track 1: Explicit H ---
        found_systems, orig_h_positions = _run_explicit_track(
            rctx, x_cra, x_atom, x_mark, x_pos_arr, is_sym_mate,
            x_elem, dist_x_pi, dist_x_centroid, proj_dist, combined_occ, sym_op, hits
        )

        # --- Track 2: Cone rescue ---
        # Cone rescue is evaluated per system family. This preserves the old
        # single-system behavior when Hudson/Plevin or P-slab need rescue even
        # if another system already had an explicit-H hit.
        cone_needed_systems = set()
        if not (found_systems & {"hudson", "plevin"}):
            cone_needed_systems.update({"hudson", "plevin"})
        if "p_slab" not in found_systems:
            cone_needed_systems.add("p_slab")

        if cone_needed_systems and use_cone and allow_cone_scan:
            _run_cone_track(
                rctx, x_cra, x_atom, x_mark, x_res, x_res_name, x_elem,
                x_pos_arr, is_sym_mate, dist_x_pi, dist_x_centroid, proj_dist,
                combined_occ, orig_h_positions, sym_op, cone_needed_systems, hits
            )

    return hits


def _run_explicit_track(rctx: _RingContext, x_cra, x_atom, x_mark,
                        x_pos_arr, is_sym_mate,
                        x_elem, dist_x_pi, dist_x_centroid, proj_dist,
                        combined_occ, sym_op, hits) -> tuple:
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
        metrics = _evaluate_systems(
            rctx, x_elem, x_pos_arr, h_pos_arr, dist_x_pi, dist_x_centroid, proj_dist)
        if not (metrics['is_hudson'] or metrics['is_plevin'] or metrics['is_p_slab']):
            continue

        if metrics['is_hudson']:
            found_systems.add("hudson")
        if metrics['is_plevin']:
            found_systems.add("plevin")
        if metrics['is_p_slab']:
            found_systems.add("p_slab")

        _record_hit(hits, rctx, x_cra, x_atom, h_atom.name, dist_x_pi,
                    dist_x_centroid, x_pos_arr, proj_dist, metrics,
                    is_cone=False, combined_occ=h_combined_occ, sym_op=sym_op, h_atom=h_atom)
    
    return found_systems, orig_h_positions


def _run_cone_track(rctx: _RingContext, x_cra, x_atom, x_mark, x_res, x_res_name, x_elem,
                    x_pos_arr, is_sym_mate, dist_x_pi, dist_x_centroid, proj_dist,
                    combined_occ, orig_h_positions, sym_op, needed_systems: Set[str], hits):
    """Track 2: Cone rescue for rotatable groups."""
    if x_res_name not in config.ROTATABLE_MAPPING:
        return
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
    
    # Extract local heavy atoms and polar acceptors for steric/hbond checks
    cone_search_pos = gemmi.Position(x_pos_arr[0], x_pos_arr[1], x_pos_arr[2])
    neighbors = rctx.ns.find_atoms(cone_search_pos, radius=4.0)
    
    env_coords_list = []
    acceptor_coords_list = []
    
    for n_mark in neighbors:
        if n_mark.pos.dist(cone_search_pos) < 0.01: continue
        n_cra = n_mark.to_cra(rctx.model)
        
        if n_cra.residue.seqid == x_res.seqid and n_cra.chain.name == x_cra.chain.name:
            continue
            
        n_elem = n_cra.atom.element.name.upper()
        if n_elem in ('H', 'D', ''): continue
        
        n_pos_arr = _pos_to_arr(n_mark.pos)
        dist = np.linalg.norm(n_pos_arr - x_pos_arr)
        
        if dist <= 4.0:
            env_coords_list.append(n_pos_arr)
        if dist <= 3.5 and n_elem in ('O', 'N', 'S'):
            acceptor_coords_list.append(n_pos_arr)
            
    env_coords = np.array(env_coords_list) if env_coords_list else np.empty((0, 3))
    acceptor_coords = np.array(acceptor_coords_list) if acceptor_coords_list else np.empty((0, 3))
    
    is_locked = geometry.check_hbond_locked(x_pos_arr, orig_h_positions, acceptor_coords)
    
    h_candidates_cone = []
    
    if not is_locked:
        flexible_donors = {('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH'), ('CYS', 'SG')}
        
        if (x_res_name, x_atom.name) in flexible_donors:
            h_candidates_cone = geometry.generate_rotated_hydrogens(
                parent_pos_arr, x_pos_arr, x_elem, 
                env_coords=env_coords, clash_cutoff=2.0, num_samples=72
            )
        else:
            axis = x_pos_arr - parent_pos_arr
            axis_norm = np.linalg.norm(axis)
            if axis_norm > 1e-5:
                axis = axis / axis_norm
                wobble_angles_deg = [a for a in range(-20, 21, 5) if a != 0]
                
                for h_pos_orig in orig_h_positions:
                    vec_xh = h_pos_orig - x_pos_arr
                    for angle_deg in wobble_angles_deg:
                        theta_rad = np.radians(angle_deg)
                        cos_t = np.cos(theta_rad)
                        sin_t = np.sin(theta_rad)
                        
                        cross_prod = np.cross(axis, vec_xh)
                        dot_prod = np.dot(axis, vec_xh)
                        
                        vec_xh_rotated = (vec_xh * cos_t + 
                                          cross_prod * sin_t + 
                                          axis * dot_prod * (1 - cos_t))
                        
                        h_pos_wobbled = x_pos_arr + vec_xh_rotated
                        
                        if len(env_coords) > 0:
                            min_d = np.min(np.linalg.norm(env_coords - h_pos_wobbled, axis=1))
                            if min_d < 2.0: continue
                                
                        h_candidates_cone.append(h_pos_wobbled)

    # Score and select best cone candidates independently for the legacy
    # Hudson/Plevin family and for the P-slab model. Their historical
    # single-system behavior used different scoring rules.
    legacy_best = None
    legacy_best_score = None
    p_slab_best = None
    p_slab_best_score = None
    needs_legacy = bool(needed_systems & {"hudson", "plevin"})
    needs_p_slab = "p_slab" in needed_systems
    
    for h_pos_np in h_candidates_cone:
        metrics = _evaluate_systems(
            rctx, x_elem, x_pos_arr, h_pos_np, dist_x_pi, dist_x_centroid, proj_dist)

        if needs_legacy and (metrics['is_hudson'] or metrics['is_plevin']):
            legacy_score = metrics['angle_XH_Pi'] if metrics['angle_XH_Pi'] is not None else -1.0
            if legacy_best_score is None or legacy_score > legacy_best_score:
                legacy_best_score = legacy_score
                legacy_best = metrics

        if needs_p_slab and metrics['is_p_slab']:
            h_proj_score = metrics['h_proj_dist'] if metrics['h_proj_dist'] is not None else 999.0
            h_ray_score = metrics['H_ray_t'] if metrics['H_ray_t'] is not None else 999.0
            p_slab_score = (-h_proj_score, -h_ray_score)
            if p_slab_best_score is None or p_slab_score > p_slab_best_score:
                p_slab_best_score = p_slab_score
                p_slab_best = metrics

    combined = _merge_cone_system_metrics(legacy_best, p_slab_best, needs_legacy, needs_p_slab)
    if combined is not None:
        _record_hit(hits, rctx, x_cra, x_atom, "virt", dist_x_pi,
                    dist_x_centroid, x_pos_arr, proj_dist, combined,
                    is_cone=True, combined_occ=combined_occ, sym_op=sym_op, h_atom=None)


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

    if legacy_best is not None:
        base['theta'] = legacy_best['theta']
        base['angle_XPCN'] = legacy_best['angle_XPCN']
        base['angle_XH_Pi'] = legacy_best['angle_XH_Pi']

    return base


def _record_hit(hits: List[Dict[str, Any]], rctx: _RingContext,
                x_cra, x_atom, h_name: str, dist: float,
                dist_x_centroid: float, x_pos: np.ndarray, proj: Optional[float],
                metrics: Dict[str, Any],
                is_cone: bool = False,
                combined_occ: float = 1.0, sym_op: int = 0,
                h_atom: Optional[gemmi.Atom] = None):
    
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

    hits.append({
        'pdb': rctx.pdb_name,
        'model': rctx.model_id,
        'resolution': rctx.resolution,
        'pi_chain': rctx.chain.name,
        'pi_res': rctx.residue.name,
        'pi_id': str(rctx.residue.seqid),
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
        'is_p_slab': metrics['is_p_slab'],
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
        'P_radius': _round_float(rctx.p_radius, 3),
        'P_slab_half_thickness': _round_float(rctx.p_slab_half_thickness, 3),
        'proj_dist': _round_float(proj, 3) if proj is not None else None,
        'theta': _round_float(metrics['theta'], 2) if metrics['theta'] is not None else None,
        'angle_XPCN': _round_float(metrics['angle_XPCN'], 2) if metrics['angle_XPCN'] is not None else None,
        'angle_XH_Pi': _round_float(metrics['angle_XH_Pi'], 2) if metrics['angle_XH_Pi'] is not None else None,
        'h_proj_dist': _round_float(metrics['h_proj_dist'], 3) if metrics['h_proj_dist'] is not None else None,
        'H_ray_t': _round_float(metrics['H_ray_t'], 3) if metrics['H_ray_t'] is not None else None,
        'sym_op': sym_op,
        '_combined_occ': float(combined_occ),
        '_pi_ring_key': rctx.mode,
    })
