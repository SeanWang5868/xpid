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

def _residue_key(model_idx: int, chain: gemmi.Chain, residue: gemmi.Residue) -> tuple:
    return (model_idx, chain.name, str(residue.seqid).strip(), residue.name)


def _has_hydrogen(residue: gemmi.Residue) -> bool:
    return any(atom.element.name.upper() in {"H", "D"} for atom in residue)


def _residues_needing_hydrogen(structure: gemmi.Structure,
                               h_change_val: int) -> set:
    residue_keys = set()
    for model_idx, model in enumerate(structure):
        for chain in model:
            for residue in chain:
                if _has_hydrogen(residue):
                    continue
                if h_change_val == 4 and residue.is_water():
                    continue
                residue_keys.add(_residue_key(model_idx, chain, residue))
    return residue_keys


def _residue_lookup(structure: gemmi.Structure) -> dict:
    lookup = {}
    for model_idx, model in enumerate(structure):
        for chain in model:
            for residue in chain:
                lookup[_residue_key(model_idx, chain, residue)] = residue
    return lookup


def _merge_new_hydrogens(target: gemmi.Structure,
                         source: gemmi.Structure,
                         residue_keys: set) -> None:
    target_residues = _residue_lookup(target)
    source_residues = _residue_lookup(source)

    for key in residue_keys:
        target_residue = target_residues.get(key)
        source_residue = source_residues.get(key)
        if target_residue is None or source_residue is None:
            continue
        if _has_hydrogen(target_residue):
            continue
        for atom in source_residue:
            if atom.element.name.upper() in {"H", "D"}:
                h_atom = atom.clone()
                h_atom.calc_flag = gemmi.CalcFlag.Calculated
                h_atom.flag = "G"
                target_residue.add_atom(h_atom)


def _prepare_topology_with_retries(structure: gemmi.Structure,
                                   monlib: gemmi.MonLib,
                                   h_change_val: int) -> None:
    max_attempts = 10
    for attempt in range(max_attempts):
        try:
            gemmi.prepare_topology(
                structure, monlib, model_index=0, h_change=h_change_val,
                reorder=False, ignore_unknown_links=True)
            return
        except Exception as topo_err:
            err_msg = str(topo_err)

            if "link" in err_msg.lower():
                logger.warning("  -> Bad explicit link detected. Clearing connections and retrying...")
                structure.connections.clear()
                continue

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
                                    del chain[i]
                                    removed = True
                                    logger.warning(
                                        f"  -> Removed twisted residue {bad_chain}/{bad_seqid} "
                                        "from hydrogen-placement copy.")
                                    break
                        if removed:
                            break
                    if removed:
                        break

                if removed:
                    continue

            logger.warning(f"  -> Topology incomplete after {attempt} retries: {err_msg}.")
            return


def add_hydrogens_memory(structure: gemmi.Structure,
                         h_change_val: int = 4) -> gemmi.Structure:
    try:
        if h_change_val == 0: return structure

        residue_keys = _residues_needing_hydrogen(structure, h_change_val)
        if not residue_keys:
            structure.setup_cell_images()
            return structure

        all_codes = set()
        for model in structure:
            for chain in model:
                for residue in chain: all_codes.add(residue.name)
        
        monlib = _get_shared_monlib(all_codes)

        working = structure.clone()
        _prepare_topology_with_retries(working, monlib, h_change_val)
        _merge_new_hydrogens(structure, working, residue_keys)

        structure.setup_cell_images()
        return structure
        
    except Exception as e:
        logger.error(f"Critical error in prep: {e}")
        # Never return None here — return the original structure so the pipeline continues
        return structure
