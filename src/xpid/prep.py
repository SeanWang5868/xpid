"""
prep.py
Handles structure preparation and Hydrogen addition.
"""
import gemmi
import re
import logging
from pathlib import Path
from typing import Optional

from . import monlib

logger = logging.getLogger("xpid.prep")
_CACHED_MONLIB: Optional[gemmi.MonLib] = None
_CACHED_LIB_PATH: Optional[str] = None

def _get_shared_monlib(residue_names: set) -> gemmi.MonLib:
    global _CACHED_MONLIB, _CACHED_LIB_PATH
    source_path = monlib.get_monomer_library_path()
    cache_key = source_path

    if _CACHED_MONLIB is None or _CACHED_LIB_PATH != cache_key:
        _CACHED_MONLIB = gemmi.MonLib()
        _CACHED_LIB_PATH = cache_key

    gemmi_monlib = _CACHED_MONLIB

    if source_path and monlib.is_monomer_library(Path(source_path)):
        missing = [c for c in residue_names if c not in gemmi_monlib.monomers]
        if missing:
            try:
                gemmi_monlib.read_monomer_lib(source_path, missing)
            except Exception as exc:
                logger.warning(f"Could not read local monomer library {source_path}: {exc}")
                for code in missing:
                    cif_path = monlib.find_monomer_cif(code)
                    if not cif_path:
                        continue
                    try:
                        gemmi_monlib.read_monomer_cif(str(cif_path))
                    except Exception as cif_exc:
                        logger.warning(f"Could not read monomer {code} from {cif_path}: {cif_exc}")

    return gemmi_monlib

def add_hydrogens_memory(structure: gemmi.Structure,
                         h_change_val: int = 4) -> gemmi.Structure:
    try:
        if h_change_val == 0: return structure

        all_codes = set()
        for model in structure:
            for chain in model:
                for residue in chain: all_codes.add(residue.name)
        
        monlib = _get_shared_monlib(all_codes)

        max_attempts = 10  # Maximum auto-removal retries for problematic residues
        for attempt in range(max_attempts):
            try:
                # Build topology and add hydrogens
                gemmi.prepare_topology(structure, monlib, model_index=0, h_change=h_change_val, reorder=False, ignore_unknown_links=True)
                break  # Success: exit retry loop
                
            except Exception as topo_err:
                err_msg = str(topo_err)
                
                # Strategy 1: Clear explicit links on link-related errors
                if "link" in err_msg.lower():
                    logger.warning("  -> Bad explicit link detected. Clearing connections and retrying...")
                    structure.connections.clear()
                    continue
                    
                # Strategy 2: Surgically remove distorted residues
                # Parse error format, e.g.: "bonded to V/NAG 2/O4 failed"
                match = re.search(r"bonded to ([^/]+)/([^ ]+) ([^/]+)/([^ ]+) failed", err_msg)
                if match:
                    bad_chain = match.group(1)
                    bad_seqid = match.group(3)
                    
                    removed = False
                    for model in structure:
                        for chain in model:
                            if chain.name == bad_chain:
                                for i in range(len(chain)):
                                    if str(chain[i].seqid).strip() == bad_seqid.strip():
                                        del chain[i]  # Remove the problematic residue
                                        removed = True
                                        logger.warning(f"  -> Removed twisted residue {bad_chain}/{bad_seqid} to save the rest of the structure.")
                                        break
                            if removed: break
                        if removed: break
                        
                    if removed:
                        continue  # Retry after successful removal
                
                # Strategy 3: Give up if error is unrecognized
                logger.warning(f"  -> Topology incomplete after {attempt} retries: {err_msg}.")
                break  # Continue with partial hydrogen placement
        # ---------------------------------------------------------

        structure.setup_cell_images()
        return structure
        
    except Exception as e:
        logger.error(f"Critical error in prep: {e}")
        # Never return None here — return the original structure so the pipeline continues
        return structure
