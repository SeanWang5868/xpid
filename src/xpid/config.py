"""
config.py
Configuration constants and atom definitions.
"""
import gemmi
import os
import json
from pathlib import Path
from typing import Dict, Set, Optional

# --- Configuration Management ---
CONFIG_FILE = Path.home() / ".xpid_config.json"

def load_saved_mon_lib() -> Optional[str]:
    """Loads the saved Monomer Library path from local config."""
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r') as f:
                data = json.load(f)
                return data.get("monomer_library_path", None)
        except Exception:
            return None
    return None

def validate_monomer_library(path: str) -> Optional[str]:
    """
    Validates if the path is a valid CCP4 Monomer Library.
    Checks for 'list/mon_lib_list.cif'.
    Tries to auto-correct path if 'monomers' subdir exists.
    
    Returns:
        Absolute valid path (str) if valid.
        None if invalid.
    """
    p = Path(path).resolve()
    
    # Check direct path
    if (p / "list" / "mon_lib_list.cif").exists():
        return str(p)
    
    # Check strict 'monomers' subdirectory (Common user mistake)
    if (p / "monomers" / "list" / "mon_lib_list.cif").exists():
        return str(p / "monomers")
    
    return None

def save_mon_lib_path(path: str) -> bool: # 注意：返回值改为 bool 以指示成功/失败
    """Saves Monomer Library path to local config ONLY if valid."""
    valid_path = validate_monomer_library(path)
    
    if not valid_path:
        return False
    
    data = {}
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r') as f:
                data = json.load(f)
        except Exception:
            data = {}
            
    data["monomer_library_path"] = valid_path
    
    with open(CONFIG_FILE, 'w') as f:
        json.dump(data, f, indent=4)
    
    return True

def clear_mon_lib_path() -> bool:
    """Removes the Monomer Library path from local config."""
    if not CONFIG_FILE.exists():
        return True
    
    try:
        with open(CONFIG_FILE, 'r') as f:
            data = json.load(f)
        
        if "monomer_library_path" in data:
            del data["monomer_library_path"]
            
            with open(CONFIG_FILE, 'w') as f:
                json.dump(data, f, indent=4)
        return True
    except Exception:
        return False
    
# --- Defaults ---
_env_path = os.environ.get("GEMMI_MON_LIB_PATH", None)
_saved_path = load_saved_mon_lib()
DEFAULT_MON_LIB_PATH = _env_path if _env_path else _saved_path

# H-Change Mode: 4 = ReAddButWater (Default)
DEFAULT_H_CHANGE = 4 

# --- Atom Definitions ---
RING_ATOMS: Dict[str, Set[str]] = {
    'TRP': {'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'TYR': {'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'},
    'PHE': {'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'},
    'HIS': {'CE1', 'ND1', 'NE2', 'CG', 'CD2'}
}

TRP_A_ATOMS: Dict[str, Set[str]] = {
    'TRP': {'CD1', 'CD2', 'NE1', 'CG', 'CE2'}
}

TARGET_ELEMENTS_X = {gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')}

# Support both Hydrogen (H) and Deuterium (D) for neutron structures
TARGET_ELEMENTS_H = {gemmi.Element('H'), gemmi.Element('D')}

# --- Geometric Thresholds ---
DIST_SEARCH_LIMIT = 6.0
DIST_PLEVIN_MAX = 4.3
DIST_HUDSON_MAX = 4.5
DIST_CUTOFF_H = 1.3