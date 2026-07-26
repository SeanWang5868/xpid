import logging
from pathlib import Path
from typing import List, Dict, Any, Literal, Optional, Union

from . import hydrogen_prep
from . import detector
from . import structure_io

logger = logging.getLogger("xpid")


class XPIDError(RuntimeError):
    """Raised when an XPID analysis cannot be completed."""


def detect(
    file_path: Union[str, Path],
    h_mode: int = 4,
    filter_pi: Optional[List[str]] = None,
    filter_donor: Optional[List[str]] = None,
    filter_donor_atom: Optional[List[str]] = None,
    model_mode: Union[str, int] = 0,
    cone_mode: str = "auto",
    include_p_slab: bool = False,
    report_xh_candidates: bool = False,
    include_coordinates: bool = False,
    compute_sasa: bool = False,
    annotate_cooperativity: bool = True,
    residue_pair: Optional[tuple[str, str]] = None,
    min_occ: float = 0.0,
    sym_contacts: bool = False,
    include_water: bool = False,
    max_b: float = 0.0,
    on_error: Literal["raise", "empty"] = "raise",
) -> List[Dict[str, Any]]:
    """
    API entry point to run analysis on a single file from Python code.
    """
    if on_error not in {"raise", "empty"}:
        raise ValueError("on_error must be 'raise' or 'empty'")
    path_obj = Path(file_path)
    # Handle filename extraction
    if path_obj.name.count('.') > 1:
         pdb_name = path_obj.name.split('.')[0]
    else:
         pdb_name = path_obj.stem

    try:
        structure = structure_io.read_structure(path_obj)
        structure = hydrogen_prep.add_hydrogens_memory(
            structure, h_change_val=h_mode)

        if not structure:
            return []

        results = detector.detect_interactions_in_structure(
            structure,
            pdb_name,
            filter_pi=filter_pi,
            filter_donor=filter_donor,
            filter_donor_atom=filter_donor_atom,
            model_mode=model_mode,
            cone_mode=cone_mode,
            include_p_slab=include_p_slab,
            report_xh_candidates=report_xh_candidates,
            include_coordinates=include_coordinates,
            compute_sasa=compute_sasa,
            annotate_cooperativity=annotate_cooperativity,
            residue_pair=residue_pair,
            min_occ=min_occ,
            sym_contacts=sym_contacts,
            include_water=include_water,
            max_b=max_b
        )
        return results

    except Exception as exc:
        if on_error == "empty":
            logger.error("Analysis error: %s", exc)
            return []
        raise XPIDError(f"Analysis failed for {path_obj}: {exc}") from exc
