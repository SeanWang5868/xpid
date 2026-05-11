import gemmi
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Union

from . import prep
from . import core

logger = logging.getLogger("xpid")

def detect(
    file_path: Union[str, Path],
    h_mode: int = 4,
    filter_pi: Optional[List[str]] = None,
    filter_donor: Optional[List[str]] = None,
    filter_donor_atom: Optional[List[str]] = None,
    model_mode: Union[str, int] = 0,
    use_cone: bool = True,
    min_occ: float = 0.0,
    sym_contacts: bool = False,
    include_water: bool = False,
    max_b: float = 0.0
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
        structure = prep.add_hydrogens_memory(structure, h_change_val=h_mode)
        
        if not structure:
            return []

        results = core.detect_interactions_in_structure(
            structure, 
            pdb_name, 
            filter_pi=filter_pi, 
            filter_donor=filter_donor,
            filter_donor_atom=filter_donor_atom,
            model_mode=model_mode,
            use_cone=use_cone,
            min_occ=min_occ,
            sym_contacts=sym_contacts,
            include_water=include_water,
            max_b=max_b
        )
        return results

    except Exception as e:
        logger.error(f"Analysis error: {e}")
        return []
