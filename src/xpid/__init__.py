import gemmi
from pathlib import Path
from typing import List, Dict, Any, Optional, Union

from . import prep
from . import core
from . import config

def detect(
    file_path: Union[str, Path],
    mon_lib_path: Optional[str] = None,
    h_mode: int = 4,
    filter_pi: Optional[List[str]] = None,
    filter_donor: Optional[List[str]] = None,
    filter_donor_atom: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """
    API entry point to run analysis on a single file from Python code.
    """
    path_obj = Path(file_path)
    # Handle filename extraction
    if path_obj.name.count('.') > 1:
         pdb_name = path_obj.name.split('.')[0]
    else:
         pdb_name = path_obj.stem

    try:
        structure = gemmi.read_structure(str(path_obj))
        structure = prep.add_hydrogens_memory(structure, mon_lib_path, h_change_val=h_mode)
        
        if not structure:
            return []

        results = core.detect_interactions_in_structure(
            structure, 
            pdb_name, 
            filter_pi=filter_pi, 
            filter_donor=filter_donor,
            filter_donor_atom=filter_donor_atom
        )
        return results

    except Exception as e:
        print(f"Analysis error: {e}")
        return []