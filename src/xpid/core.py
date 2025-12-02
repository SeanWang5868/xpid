"""
core.py
Core logic for detecting XH-pi interactions based on geometric and topological criteria.
Implements the Hudson and Plevin criteria.
"""
import gemmi
import numpy as np
from typing import List, Dict, Any, Optional, Union
from . import config
from . import geometry
from . import residue_ss

def detect_interactions_in_structure(structure: gemmi.Structure, 
                                     pdb_name: str,
                                     filter_pi: Optional[List[str]] = None,
                                     filter_donor: Optional[List[str]] = None,
                                     filter_donor_atom: Optional[List[str]] = None,
                                     model_mode: Union[str, int] = 0) -> List[Dict[str, Any]]:
    """
    Scans the structure for interactions.
    Supports multi-model analysis (e.g., NMR).
    """
    results = []
    
    if not structure or len(structure) == 0:
        return []

    # Determine which models to process
    models_with_ids = [] # List of tuples (model, model_id_str)
    
    if model_mode == 'all':
        for i, m in enumerate(structure):
            # [Bug Fix] Safety check for name attribute
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
                # Fallback to first model if index out of bounds, or handle empty
                pass 
        except ValueError:
            m = structure[0]
            m_name = getattr(m, 'name', "1")
            models_with_ids = [(m, m_name)]

    resolution = structure.resolution if structure.resolution else 0.0

    # Iterate selected models
    for model, model_id in models_with_ids:
        # 1. Neighbor Search Grid
        ns = gemmi.NeighborSearch(model, structure.cell, config.DIST_SEARCH_LIMIT)
        ns.populate(include_h=True)

        # 2. Build Secondary Structure Index
        ss_index = residue_ss.build_index(structure)
        
        for chain in model:
            for residue in chain:
                res_name = residue.name
                
                # Apply Pi Residue Filter
                if filter_pi and res_name not in filter_pi:
                    continue

                # Scan Main Pi Systems
                if res_name in config.RING_ATOMS:
                    results.extend(_detect_residue(
                        pdb_name, resolution, model, model_id, chain, residue, ns, ss_index,
                        config.RING_ATOMS[res_name], 'main', filter_donor, filter_donor_atom
                    ))
                
                # Scan Trp A Systems
                if res_name in config.TRP_A_ATOMS:
                    results.extend(_detect_residue(
                        pdb_name, resolution, model, model_id, chain, residue, ns, ss_index,
                        config.TRP_A_ATOMS[res_name], 'trpA', filter_donor, filter_donor_atom
                    ))
    return results

def _detect_residue(pdb_name, resolution, model, model_id, chain, residue, ns, ss_index, 
                    target_atoms, mode, filter_donor, filter_donor_atom):
    hits = []
    
    # --- Pi Side ---
    pi_atoms = [atom for atom in residue if atom.name in target_atoms]
    if len(pi_atoms) != len(target_atoms): return []
    
    pi_center, pi_center_arr, pi_normal, pi_b_mean = geometry.get_pi_info(pi_atoms)
    alt_pi = pi_atoms[0].altloc
    
    # --- Search X Donors ---
    x_candidates = ns.find_atoms(pi_center, alt=alt_pi, radius=config.DIST_SEARCH_LIMIT)
    
    for x_mark in x_candidates:
        x_cra = x_mark.to_cra(model)
        x_atom = x_cra.atom
        x_res_name = x_cra.residue.name
        
        # Element Check
        if x_atom.element not in config.TARGET_ELEMENTS_X: continue
        
        # Filters
        if filter_donor and x_res_name not in filter_donor: continue
        if filter_donor_atom and x_atom.element.name not in filter_donor_atom: continue

        x_pos_arr = np.array(x_atom.pos.tolist())
        dist_x_pi = geometry.calculate_distance(x_pos_arr, pi_center_arr)
        if dist_x_pi > config.DIST_HUDSON_MAX: continue
        
        # --- Search H ---
        h_candidates = ns.find_atoms(x_atom.pos, alt=x_atom.altloc, radius=config.DIST_CUTOFF_H)
        
        for h_mark in h_candidates:
            h_cra = h_mark.to_cra(model)
            h_atom = h_cra.atom
            
            if h_atom.element not in config.TARGET_ELEMENTS_H: continue
            
            # --- Geometric Checks ---
            h_pos_arr = np.array(h_atom.pos.tolist())
            xpcn_angle = geometry.calculate_xpcn_angle(x_pos_arr, pi_center_arr, pi_normal)
            xh_pi_angle = geometry.calculate_xh_picenter_angle(pi_center_arr, x_pos_arr, h_pos_arr)
            theta = geometry.calculate_hudson_theta(pi_center_arr, x_pos_arr, h_pos_arr, pi_normal)
            
            if xh_pi_angle is None or theta is None or xpcn_angle is None: continue
            
            # --- Criteria ---
            plevin = 0
            if (dist_x_pi < config.DIST_PLEVIN_MAX and 
                xh_pi_angle > 120.0 and 
                xpcn_angle < 25.0):
                plevin = 1
            
            proj_threshold = None
            if mode == 'trpA': proj_threshold = 1.6
            elif residue.name == 'HIS': proj_threshold = 1.6
            elif residue.name in ['TRP', 'TYR', 'PHE']: proj_threshold = 2.0
            
            hudson = 0
            proj_dist = None
            if proj_threshold is not None:
                proj_dist = geometry.calculate_projection_dist(pi_normal, pi_center_arr, x_pos_arr)
                if (theta <= 40.0 and 
                    dist_x_pi <= config.DIST_HUDSON_MAX and 
                    proj_dist is not None and 
                    proj_dist <= proj_threshold):
                    hudson = 1

            if plevin == 0 and hudson == 0: continue
            
            # --- Result Construction ---
            pi_ss_type, pi_ss_uid = residue_ss.get_info(chain.name, residue.seqid.num, ss_index)
            x_ss_type, x_ss_uid = residue_ss.get_info(x_cra.chain.name, x_cra.residue.seqid.num, ss_index)

            if chain.name == x_cra.chain.name:
                try:
                    seq_sep = abs(residue.seqid.num - x_cra.residue.seqid.num)
                except:
                    seq_sep = -1
            else:
                seq_sep = -1

            remark = ""
            if residue.name == "TRP":
                remark = "6-ring" if mode == "main" else "5-ring"

            pi_cx, pi_cy, pi_cz = np.round(pi_center_arr, 3)
            x_cx, x_cy, x_cz = np.round(x_pos_arr, 3)

            hits.append({
                # Metadata
                'pdb': pdb_name,
                'model': model_id, # [Modified] Passed explicitly as string
                'resolution': resolution,
                
                # Pi
                'pi_chain': chain.name,
                'pi_res': residue.name,
                'pi_id': residue.seqid.num,
                
                # X
                'X_chain': x_cra.chain.name,
                'X_res': x_res_name,
                'X_id': x_cra.residue.seqid.num,
                'X_atom': x_atom.name,
                'H_atom': h_atom.name,
                'dist_X_Pi': round(dist_x_pi, 3),
                
                # Validation
                'is_plevin': plevin,
                'is_hudson': hudson,
                'remark': remark,
                
                # Detailed
                'pi_ss_type': pi_ss_type,
                'pi_ss_id': pi_ss_uid,
                'X_ss_type': x_ss_type,
                'X_ss_id': x_ss_uid,
                'pi_avg_b': round(pi_b_mean, 2),
                'pi_center_x': pi_cx,
                'pi_center_y': pi_cy,
                'pi_center_z': pi_cz,
                'X_b': round(x_atom.b_iso, 2),
                'X_xyz_x': x_cx,
                'X_xyz_y': x_cy,
                'X_xyz_z': x_cz,
                'seq_sep': seq_sep,
                'theta': round(theta, 2),
                'angle_XH_Pi': round(xh_pi_angle, 2),
                'angle_XPCN': round(xpcn_angle, 2),
                'proj_dist': round(proj_dist, 3) if proj_dist is not None else None,
            })
            
    return hits