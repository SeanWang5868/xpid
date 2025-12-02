"""
prep.py
Handles structure preparation, including Monomer Library loading and Hydrogen addition.
"""
import gemmi
import os
import logging
from typing import Optional, Set
from . import config

logger = logging.getLogger("xpid.prep")

# Process-level cache for Monomer Library
_CACHED_MONLIB: Optional[gemmi.MonLib] = None
_CACHED_LIB_PATH: Optional[str] = None

def _get_shared_monlib(mon_lib_path: Optional[str]) -> gemmi.MonLib:
    """Retrieves or initializes the cached Monomer Library."""
    global _CACHED_MONLIB, _CACHED_LIB_PATH
    
    if _CACHED_MONLIB is None or _CACHED_LIB_PATH != mon_lib_path:
        monlib = gemmi.MonLib()
        _CACHED_MONLIB = monlib
        _CACHED_LIB_PATH = mon_lib_path
        
    return _CACHED_MONLIB

def add_hydrogens_memory(structure: gemmi.Structure, 
                         mon_lib_path: Optional[str] = None, 
                         h_change_val: int = config.DEFAULT_H_CHANGE) -> Optional[gemmi.Structure]:
    """
    Adds hydrogens to the structure in-memory using Gemmi topology.
    Skips processing if h_change_val is 0 (e.g., for Neutron structures).
    """
    try:
        # Mode 0: Return original structure (safe for Neutron data)
        if h_change_val == 0:
            return structure

        # 1. Identify used residues
        all_codes = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name.strip():
                        all_codes.add(residue.name)
        
        codes_to_load = list(all_codes)

        # 2. Get MonLib (Cached)
        monlib = _get_shared_monlib(mon_lib_path)
        
        # 3. Incrementally load missing monomers
        missing_codes = [code for code in codes_to_load if code not in monlib.monomers]
        
        if missing_codes:
            if mon_lib_path and os.path.exists(mon_lib_path):
                try:
                    monlib.read_monomer_lib(mon_lib_path, missing_codes)
                except Exception as e:
                    logger.warning(f"Failed to load monomers {missing_codes}: {e}")
            else:
                # Fallback to Gemmi's internal library
                pass

        # 4. Prepare Topology
        gemmi.prepare_topology(
            structure, 
            monlib, 
            model_index=0, 
            h_change=h_change_val, 
            reorder=False, 
            ignore_unknown_links=True 
        )
        
        return structure

    except Exception as e:
        logger.error(f"Topology preparation failed: {e}")
        return None