"""
core.py
Core logic for detecting XH-π interactions with Hudson/Plevin labels by default.
The P-slab system is retained as an opt-in method-development mode.
Track 1: Explicit H (Default)
Track 2: Binary conformer-existence test for chemically rotatable groups
"""
import gemmi
import logging
import numpy as np
from typing import List, Dict, Any, Optional, Union, Set, NamedTuple
from . import config
from . import aromatic_rings
from . import monomer_bonds
from . import geometry
from . import sasa
from . import cooperativity
from . import cone
from . import rotatable_groups
from . import hbond_acceptors
from . import hbond_annotations
from . import ring_conformers
from . import xhpi_criteria
from . import hits as _hits
from . import ss
from . import schema

logger = logging.getLogger("xpid.core")

ResiduePairKeys = Optional[tuple[Set[tuple], Set[tuple]]]

def _pos_to_arr(pos: gemmi.Position) -> np.ndarray:
    """Convert gemmi.Position to numpy array without intermediate list."""
    return np.array([pos.x, pos.y, pos.z])


def _nearest_mark_position(
        reference: gemmi.Position,
        mark: gemmi.NeighborSearch.Mark,
        cell: gemmi.UnitCell,
        base_position: gemmi.Position,
        ) -> tuple[gemmi.Position, bool, str]:
    """Resolve a NeighborSearch mark into the nearest crystallographic copy.

    ``Mark.pos`` can be a grid-wrapped coordinate and ``image_idx`` does not
    encode the lattice translation (P1 translations still use image index 0),
    so the atom's original ASU position is required explicitly.
    """
    position = cell.find_nearest_pbc_position(
        reference, base_position, mark.image_idx)
    nearest = cell.find_nearest_image(
        reference, base_position, gemmi.Asu.Any)
    is_symmetry_mate = (
        mark.image_idx != 0 or position.dist(base_position) > 1e-5)
    symmetry_code = nearest.symmetry_code() if is_symmetry_mate else "1_555"
    return position, is_symmetry_mate, symmetry_code

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
                                     max_b: float = 0.0, compute_sasa: bool = False,
                                     annotate_cooperativity: bool = True,
                                     diagnostics: Optional[Dict[str, Any]] = None
                                     ) -> List[Dict[str, Any]]:
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

                residue_rings = aromatic_rings.get_aromatic_rings(res_name)
                if not residue_rings: continue

                for ring_idx, target_atoms in enumerate(residue_rings):
                    ring_size = len(target_atoms)
                    mode = f'ring{ring_idx+1}'

                    ring_variants = ring_conformers._atom_variants_for_names(
                        residue, target_atoms)
                    if not ring_variants and diagnostics is not None:
                        diagnostics.setdefault(
                            "incomplete_aromatic_rings", set()).add(
                                f"{model_id}:{chain.name}:"
                                f"{str(residue.seqid).strip()}:"
                                f"{res_name}:{mode}")
                    for pi_alt, pi_atoms in ring_variants:
                        results.extend(_detect_residue(
                            pdb_name, resolution, model, model_id, chain, residue, ns, ss_index,
                            pi_atoms, pi_alt, mode, filter_donor, filter_donor_atom, cone_mode,
                            include_p_slab, report_xh_candidates, include_coordinates,
                            pair_keys, ring_size, sasa_map, min_occ,
                            sym_contacts=sym_contacts, max_b=max_b,
                            diagnostics=diagnostics,
                        ))
    hits = _hits._deduplicate_hits(results, prefer_directional=not report_xh_candidates)
    if annotate_cooperativity:
        hits = cooperativity.annotate_cooperativity(hits)
    # Annotate each model against its own coordinates.  Reusing model 0 for
    # ``--model all`` can assign acceptors from the wrong conformer/model.
    for model_index, model in enumerate(structure):
        model_name = getattr(model, "name", str(model_index + 1))
        model_hits = [hit for hit in hits if str(hit.get("model")) == str(model_name)]
        if not model_hits:
            continue
        hbond_ns = gemmi.NeighborSearch(model, structure.cell, 5.0)
        hbond_ns.populate(include_h=False)
        hbond_annotations.annotate_explicit_hbond_context(
            model_hits, structure, model_index=model_index, ns=hbond_ns)
        hbond_annotations.annotate_cone_hbond_context(
            model_hits, structure, model_index=model_index, ns=hbond_ns)
    for hit in hits:
        unexpected_fields = schema.validate_hit(hit)
        if unexpected_fields:
            raise RuntimeError(
                "Internal hit schema mismatch: "
                + ", ".join(sorted(unexpected_fields)))
    return hits

def _is_donor_blocked(x_atom: gemmi.Atom, model: gemmi.Model, ns: gemmi.NeighborSearch,
                      x_pos: Optional[gemmi.Position] = None) -> bool:
    radius = 2.6
    search_pos = x_pos if x_pos is not None else x_atom.pos
    neighbors = ns.find_atoms(search_pos, radius=radius)

    x_elem = x_atom.element.name.upper()

    for mark in neighbors:
        cra = mark.to_cra(model)
        neighbor_atom = cra.atom
        neighbor_pos = ns.grid_cell.find_nearest_pbc_position(
            search_pos, neighbor_atom.pos, mark.image_idx)
        dist = neighbor_pos.dist(search_pos)
        if dist < 0.01: continue
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
                    min_occ: float, sym_contacts: bool = False,
                    max_b: float = 0.0,
                    diagnostics: Optional[Dict[str, Any]] = None):
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

    rctx = ring_conformers.RingContext(
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

        x_effective_pos, is_sym_mate, symmetry_code = (
            _nearest_mark_position(
                pi_center, x_mark, ns.grid_cell, x_atom.pos))
        if is_sym_mate and not sym_contacts:
            continue

        if filter_donor and x_res_name not in filter_donor: continue
        if filter_donor_atom and x_elem not in filter_donor_atom and x_atom.name.upper() not in filter_donor_atom:
            continue

        # A residue in the asymmetric unit must not interact with its own ring,
        # but the same residue in a crystallographic symmetry mate is a
        # physically distinct molecular copy and is a valid contact.
        if (not is_sym_mate and
                x_cra.chain.name == chain.name and
                x_res.seqid == residue.seqid and
                x_res.name == residue.name):
            continue

        if (x_res_name in ('ASP', 'GLU') and x_atom.name in ('OD1', 'OD2', 'OE1', 'OE2')) or \
           (x_atom.name == 'OXT'):
            continue

        if _is_cation_pi_donor(x_res_name, x_atom.name):
            continue

        donor_definition = rotatable_groups.get_rotatable_group(
            x_res_name, x_atom.name)
        is_rotatable_donor = donor_definition is not None

        if is_sym_mate:
            if _is_donor_blocked(
                    x_atom, model, ns, x_pos=x_effective_pos):
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
            x_pos_arr = _pos_to_arr(x_effective_pos)
        else:
            x_pos_arr = _pos_to_arr(x_atom.pos)
        dist_x_pi = geometry.calculate_plane_distance(x_pos_arr, pi_center_arr, pi_normal)
        dist_x_centroid = geometry.calculate_distance(x_pos_arr, pi_center_arr)
        proj_dist = geometry.calculate_projection_dist(
            pi_normal, pi_center_arr, x_pos_arr)

        max_plane_dist = config.P_PLANE_DMAX.get(x_elem, config.P_PLANE_DMAX['default'])

        xpcn_angle = geometry.calculate_xpcn_angle(
            x_pos_arr, pi_center_arr, pi_normal)
        hudson_spatial = (
            dist_x_centroid <= max_plane_dist and
            proj_dist is not None and proj_dist <= p_radius
        )
        plevin_spatial = (
            dist_x_centroid < max_plane_dist and
            xpcn_angle is not None and xpcn_angle < config.PLEVIN_XPCN_MAX
        )
        p_slab_spatial = (
            include_p_slab and dist_x_pi is not None and
            dist_x_pi <= max_plane_dist
        )
        if not (hudson_spatial or plevin_spatial or p_slab_spatial):
            continue

        p_slab_spatial = (
            include_p_slab and dist_x_pi is not None and
            dist_x_pi <= max_plane_dist and
            proj_dist is not None and proj_dist <= p_radius
        )
        if not (hudson_spatial or plevin_spatial or p_slab_spatial):
            continue

        sym_op = x_mark.image_idx if is_sym_mate else 0

        # --- Detection path ---
        if is_rotatable_donor and cone_mode == "auto":
            # Background-candidate export requires an observed/fixed H
            # direction.  A cone direction selected by the detector would
            # contaminate that background by construction.
            if report_xh_candidates:
                continue

            # Rotatable donor in auto mode: skip explicit-H track entirely.
            # Crystal-structure hydrogens on rotatable groups are riding
            # hydrogens with no experimental support at typical resolutions.
            # Cone rescue scans the full rotational space.
            _run_cone_track(
                rctx, x_cra, x_atom, x_mark, x_res, x_res_name, x_elem,
                x_pos_arr, is_sym_mate, dist_x_pi, dist_x_centroid, proj_dist,
                combined_occ, sym_op, include_p_slab, include_coordinates,
                sasa_map, hits, report_xh_candidates,
                donor_definition, diagnostics, symmetry_code,
            )
        else:
            # Rigid donor, or rotatable in --no-cone mode:
            # use explicit hydrogens from the structure file.
            _run_explicit_track(
                rctx, x_cra, x_atom, x_mark, x_pos_arr, is_sym_mate,
                x_elem, dist_x_pi, dist_x_centroid, proj_dist, combined_occ,
                sym_op, include_p_slab, report_xh_candidates,
                include_coordinates, sasa_map, hits, symmetry_code,
            )

    return hits





def _run_explicit_track(rctx: "ring_conformers.RingContext", x_cra, x_atom, x_mark,
                        x_pos_arr, is_sym_mate,
                        x_elem, dist_x_pi, dist_x_centroid, proj_dist,
                        combined_occ, sym_op, include_p_slab: bool,
                        report_xh_candidates: bool, include_coordinates: bool,
                        sasa_map: Dict,
                        hits, symmetry_code: str = "1_555") -> tuple:
    """Track 1: Explicit hydrogen geometry. Returns (found_systems, orig_h_positions)."""
    found_systems: Set[str] = set()
    orig_h_positions: List[np.ndarray] = []
    spatial = xhpi_criteria.prepare_spatial_criteria(
        rctx, x_elem, x_pos_arr, dist_x_centroid, proj_dist)

    h_search_pos = gemmi.Position(*x_pos_arr)
    xh_max = config.COVALENT_XH_MAX.get(x_elem)
    if xh_max is None:
        return found_systems, orig_h_positions
    h_candidates = rctx.ns.find_atoms(
        h_search_pos, alt=x_atom.altloc, radius=xh_max)

    for h_mark in h_candidates:
        h_cra = h_mark.to_cra(rctx.model)
        h_atom = h_cra.atom

        if h_atom.element.name.upper() not in {'H', 'D'}:
            continue

        if (h_cra.chain.name != x_cra.chain.name or
                h_cra.residue.seqid != x_cra.residue.seqid or
                h_cra.residue.name != x_cra.residue.name):
            continue

        h_effective_pos = rctx.ns.grid_cell.find_nearest_pbc_position(
            h_search_pos, h_atom.pos, h_mark.image_idx)
        h_pos_arr = _pos_to_arr(h_effective_pos)
        orig_h_positions.append(h_pos_arr)

        xh_dist = float(np.linalg.norm(h_pos_arr - x_pos_arr))
        if xh_dist <= config.MIN_COVALENT_XH or xh_dist > xh_max:
            continue

        bonded_hydrogens = monomer_bonds.get_bonded_hydrogen_names(
            x_cra.residue.name, x_atom.name)
        if bonded_hydrogens is not None:
            canonical_h_name = monomer_bonds.canonical_hydrogen_name(
                h_atom.name)
            canonical_bonded = {
                monomer_bonds.canonical_hydrogen_name(name)
                for name in bonded_hydrogens
            }
            if canonical_h_name not in canonical_bonded:
                continue
        elif not _is_unambiguous_covalent_hydrogen(
                x_cra.residue, x_atom, h_atom):
            continue

        if h_atom.altloc and x_atom.altloc and h_atom.altloc != x_atom.altloc:
            continue

        h_combined_occ = min(combined_occ, h_atom.occ)
        metrics = xhpi_criteria.evaluate_xhpi_geometry(
            rctx, x_elem, x_pos_arr, h_pos_arr, dist_x_pi,
            dist_x_centroid, proj_dist, spatial=spatial,
            include_direction_metrics=(
                include_p_slab or report_xh_candidates))
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
                    symmetry_code=symmetry_code,
                    h_atom=h_atom, include_p_slab=include_p_slab,
                    include_candidate_metrics=report_xh_candidates,
                    include_coordinates=include_coordinates,
                    sasa_map=sasa_map,
                    h_pos=h_pos_arr)

    return found_systems, orig_h_positions


def _is_unambiguous_covalent_hydrogen(
        residue: gemmi.Residue,
        x_atom: gemmi.Atom,
        h_atom: gemmi.Atom) -> bool:
    """Conservative fallback when no monomer dictionary is available."""
    max_xh = config.COVALENT_XH_MAX
    x_element = x_atom.element.name.upper()
    cutoff = max_xh.get(x_element)
    if cutoff is None:
        return False

    def compatible_altloc(atom: gemmi.Atom) -> bool:
        x_alt = ring_conformers._altloc(x_atom)
        atom_alt = ring_conformers._altloc(atom)
        return not x_alt or not atom_alt or x_alt == atom_alt

    if not compatible_altloc(h_atom):
        return False

    h_pos = _pos_to_arr(h_atom.pos)
    candidates = []
    for atom in residue:
        element = atom.element.name.upper()
        atom_cutoff = max_xh.get(element)
        if atom_cutoff is None or not compatible_altloc(atom):
            continue
        distance = float(np.linalg.norm(_pos_to_arr(atom.pos) - h_pos))
        if config.MIN_COVALENT_XH < distance <= atom_cutoff:
            candidates.append((distance / atom_cutoff, atom))

    if not candidates:
        return False
    candidates.sort(key=lambda item: item[0])
    if candidates[0][1] is not x_atom:
        return False
    if len(candidates) == 1:
        return True

    # Reject ambiguous ownership instead of assigning a nearby hydrogen to the
    # wrong heteroatom/carbon in an unknown ligand.
    return candidates[1][0] - candidates[0][0] >= 0.15


def _run_cone_track(rctx: "ring_conformers.RingContext", x_cra, x_atom, x_mark, x_res, x_res_name, x_elem,
                    x_pos_arr, is_sym_mate, dist_x_pi, dist_x_centroid, proj_dist,
                    combined_occ, sym_op, include_p_slab: bool,
                    include_coordinates: bool, sasa_map: Dict, hits,
                    report_xh_candidates: bool = False,
                    donor_definition: Optional[
                        rotatable_groups.RotatableGroupDefinition] = None,
                    diagnostics: Optional[Dict[str, Any]] = None,
                    symmetry_code: str = "1_555"):
    """Evaluate a chemically complete binary cone for a rotatable donor."""
    if donor_definition is None:
        donor_definition = rotatable_groups.get_rotatable_group(
            x_res_name, x_atom.name)
    if donor_definition is None:
        return

    parent_name = donor_definition.parent_atom_name
    parent_atom = rotatable_groups.find_parent_atom(
        x_res, donor_definition)
    if not parent_atom:
        if diagnostics is not None:
            diagnostics.setdefault("cone_missing_parent_groups", set()).add(
                f"{rctx.model_id}:{x_cra.chain.name}:"
                f"{str(x_res.seqid).strip()}:{x_res_name}:"
                f"{x_atom.name}:{parent_name}")
        return

    # Resolve parent position (sym mates need transformed coords)
    if is_sym_mate:
        parent_effective_pos = (
            rctx.ns.grid_cell.find_nearest_pbc_position(
                gemmi.Position(*x_pos_arr), parent_atom.pos,
                x_mark.image_idx))
        parent_pos_arr = _pos_to_arr(parent_effective_pos)
    else:
        parent_pos_arr = _pos_to_arr(parent_atom.pos)

    # Gather a model/image/altloc-compatible heavy-atom environment.  Only the
    # donor X and its directly bonded parent are excluded; other same-residue
    # atoms remain capable of causing a real non-bonded clash.
    cone_search_pos = gemmi.Position(x_pos_arr[0], x_pos_arr[1], x_pos_arr[2])
    neighbors = rctx.ns.find_atoms(cone_search_pos, radius=4.0)
    environment = []
    for n_mark in neighbors:
        n_cra = n_mark.to_cra(rctx.model)
        n_elem = n_cra.atom.element.name.upper()
        if n_elem in ('H', 'D', ''):
            continue
        same_residue = (
            n_cra.chain.name == x_cra.chain.name and
            n_cra.residue.seqid == x_res.seqid
        )
        n_effective_pos = (
            rctx.ns.grid_cell.find_nearest_pbc_position(
                cone_search_pos, n_cra.atom.pos, n_mark.image_idx))
        if same_residue:
            if (n_cra.atom.name == x_atom.name and
                    n_effective_pos.dist(cone_search_pos) < 0.01):
                continue
            if (n_cra.atom.name == parent_name and
                    np.linalg.norm(
                        _pos_to_arr(n_effective_pos) - parent_pos_arr) < 0.01):
                continue
        if not hbond_acceptors.altlocs_compatible(
                x_atom.altloc, n_cra.atom.altloc):
            continue
        environment.append(cone.EnvironmentAtom(
            _pos_to_arr(n_effective_pos), n_cra.atom, n_cra.residue,
            n_cra.chain.name, n_mark.image_idx,
            cone.VDW_RADII.get(n_elem, cone.DEFAULT_VDW_RADIUS),
            hbond_acceptors.is_hbond_acceptor(
                n_cra.residue, n_cra.atom),
            n_cra.atom.occ,
        ))

    evidence = cone.evaluate_binary(
        rctx, donor_definition, parent_pos_arr, x_pos_arr, x_elem,
        environment, dist_x_pi, dist_x_centroid, proj_dist,
        include_p_slab=include_p_slab,
    )
    if evidence is not None:
        selected_h_pos = evidence.conformer.hydrogen_positions[
            evidence.hydrogen_index]
        _hits._record_hit(hits, rctx, x_cra, x_atom, "virt", dist_x_pi,
                    dist_x_centroid, x_pos_arr, proj_dist, evidence.metrics,
                    is_cone=True, combined_occ=combined_occ, sym_op=sym_op,
                    symmetry_code=symmetry_code,
                    h_atom=None, include_p_slab=include_p_slab,
                    include_candidate_metrics=report_xh_candidates,
                    include_coordinates=include_coordinates,
                    sasa_map=sasa_map,
                    h_pos=selected_h_pos,
                    hbond_relation=evidence.hbond_relation)
