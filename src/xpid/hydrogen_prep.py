"""Gemmi-backed hydrogen preparation with auditable failure semantics."""
from dataclasses import dataclass
import gemmi
import logging
from pathlib import Path
from typing import Optional, Tuple

from . import monlib

logger = logging.getLogger("xpid.prep")
_CACHED_MONLIB: Optional[gemmi.MonLib] = None
_CACHED_LIB_PATH: Optional[str] = None


@dataclass(frozen=True)
class TopologyPreparationReport:
    model_index: int
    status: str
    attempts: int
    connections_cleared: bool = False
    message: Optional[str] = None


@dataclass(frozen=True)
class HydrogenPreparationReport:
    status: str
    h_mode: int
    model_reports: Tuple[TopologyPreparationReport, ...] = ()
    missing_monomer_components: Tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "h_mode": self.h_mode,
            "missing_monomer_components": list(
                self.missing_monomer_components),
            "models": [
                {
                    "model_index": item.model_index,
                    "status": item.status,
                    "attempts": item.attempts,
                    "connections_cleared": item.connections_cleared,
                    "message": item.message,
                }
                for item in self.model_reports
            ],
        }


@dataclass(frozen=True)
class HydrogenPreparationResult:
    structure: gemmi.Structure
    report: HydrogenPreparationReport


class HydrogenPreparationError(RuntimeError):
    """Raised when Gemmi cannot prepare a chemically consistent topology."""

    def __init__(self, message: str, report: HydrogenPreparationReport):
        super().__init__(message)
        self.report = report

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


def _residue_lookup(structure: gemmi.Structure) -> dict:
    lookup = {}
    for model_idx, model in enumerate(structure):
        for chain in model:
            for residue in chain:
                lookup[_residue_key(model_idx, chain, residue)] = residue
    return lookup


def _synchronize_hydrogens(target: gemmi.Structure,
                           source: gemmi.Structure) -> None:
    """Replace target H/D atoms with Gemmi's prepared result.

    Heavy atoms remain untouched so topology-recovery removals made in the
    working copy cannot delete experimental coordinates from the caller's
    structure.
    """
    target_residues = _residue_lookup(target)
    source_residues = _residue_lookup(source)

    for key, target_residue in target_residues.items():
        source_residue = source_residues.get(key)
        if source_residue is None:
            continue

        hydrogen_indices = [
            index for index, atom in enumerate(target_residue)
            if atom.element.name.upper() in {"H", "D"}
        ]
        for index in reversed(hydrogen_indices):
            del target_residue[index]

        for atom in source_residue:
            if atom.element.name.upper() in {"H", "D"}:
                h_atom = atom.clone()
                target_residue.add_atom(h_atom)


def _prepare_topology_with_retries(structure: gemmi.Structure,
                                   monlib: gemmi.MonLib,
                                   h_change: gemmi.HydrogenChange,
                                   model_index: int
                                   ) -> TopologyPreparationReport:
    """Prepare one model, allowing one explicitly reported link recovery.

    A bad explicit link can be retried without the connection table, but the
    recovery is marked ``partial`` because covalent links can affect hydrogen
    placement.  Residues are never deleted to force topology preparation.
    """
    connections_cleared = False
    for attempt in range(1, 3):
        try:
            gemmi.prepare_topology(
                structure, monlib, model_index=model_index, h_change=h_change,
                reorder=False, ignore_unknown_links=True)
            return TopologyPreparationReport(
                model_index=model_index,
                status="partial" if connections_cleared else "success",
                attempts=attempt,
                connections_cleared=connections_cleared,
            )
        except Exception as topo_err:
            err_msg = str(topo_err)

            if "link" in err_msg.lower() and not connections_cleared:
                logger.warning(
                    "Bad explicit link detected in model %d; retrying "
                    "without explicit connections.", model_index)
                structure.connections.clear()
                connections_cleared = True
                continue
            return TopologyPreparationReport(
                model_index=model_index,
                status="failed",
                attempts=attempt,
                connections_cleared=connections_cleared,
                message=err_msg,
            )

    return TopologyPreparationReport(
        model_index=model_index,
        status="failed",
        attempts=2,
        connections_cleared=connections_cleared,
        message="Topology preparation exhausted all attempts.",
    )


def prepare_hydrogens_memory(
        structure: gemmi.Structure,
        h_change_val: int = 4) -> HydrogenPreparationResult:
    """Prepare hydrogens and return both the structure and an audit report.

    A failed topology raises :class:`HydrogenPreparationError`; it is never
    converted into an apparently successful zero-hit structure.
    """
    try:
        hydrogen_changes = {
            0: gemmi.HydrogenChange.NoChange,
            1: gemmi.HydrogenChange.Shift,
            2: gemmi.HydrogenChange.Remove,
            3: gemmi.HydrogenChange.ReAdd,
            4: gemmi.HydrogenChange.ReAddButWater,
            5: gemmi.HydrogenChange.ReAddKnown,
        }
        if h_change_val not in hydrogen_changes:
            raise ValueError(f"Unsupported Gemmi hydrogen mode: {h_change_val}")
        if h_change_val == 0:
            structure.setup_cell_images()
            return HydrogenPreparationResult(
                structure,
                HydrogenPreparationReport("not_requested", h_change_val),
            )

        all_codes = set()
        for model in structure:
            for chain in model:
                for residue in chain: all_codes.add(residue.name)

        monlib = _get_shared_monlib(all_codes)
        missing_monomers = tuple(sorted(
            code for code in all_codes if code not in monlib.monomers))

        working = structure.clone()
        h_change = hydrogen_changes[h_change_val]
        model_reports = []
        for model_index in range(len(working)):
            model_report = _prepare_topology_with_retries(
                working, monlib, h_change, model_index)
            model_reports.append(model_report)

        failed = [item for item in model_reports if item.status == "failed"]
        if failed:
            report = HydrogenPreparationReport(
                "failed", h_change_val, tuple(model_reports),
                missing_monomers)
            details = "; ".join(
                f"model {item.model_index}: {item.message}"
                for item in failed)
            raise HydrogenPreparationError(
                f"Hydrogen topology preparation failed ({details})", report)

        _synchronize_hydrogens(structure, working)

        structure.setup_cell_images()
        status = (
            "partial"
            if missing_monomers or any(
                item.status == "partial" for item in model_reports)
            else "success"
        )
        return HydrogenPreparationResult(
            structure,
            HydrogenPreparationReport(
                status, h_change_val, tuple(model_reports),
                missing_monomers),
        )

    except HydrogenPreparationError:
        raise
    except Exception as exc:
        report = HydrogenPreparationReport("failed", h_change_val)
        raise HydrogenPreparationError(
            f"Critical error in hydrogen preparation: {exc}", report
        ) from exc


def add_hydrogens_memory(structure: gemmi.Structure,
                         h_change_val: int = 4) -> gemmi.Structure:
    """Backward-compatible structure-only wrapper with strict failures."""
    return prepare_hydrogens_memory(
        structure, h_change_val=h_change_val).structure
