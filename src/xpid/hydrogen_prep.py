"""Gemmi-backed hydrogen preparation with auditable failure semantics."""
from dataclasses import dataclass
import gemmi
import logging
from pathlib import Path
import re
from typing import Dict, Optional, Set, Tuple

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
class ResiduePreparationExclusion:
    """A residue whose generated donor hydrogens cannot be trusted."""

    model_index: int
    chain: str
    seqid: str
    residue_name: str
    reason: str

    @property
    def key(self) -> tuple[int, str, str, str]:
        return (self.model_index, self.chain, self.seqid, self.residue_name)

    def to_dict(self) -> dict:
        return {
            "model_index": self.model_index,
            "chain": self.chain,
            "seqid": self.seqid,
            "residue_name": self.residue_name,
            "reason": self.reason,
        }


@dataclass(frozen=True)
class HydrogenPreparationReport:
    status: str
    h_mode: int
    model_reports: Tuple[TopologyPreparationReport, ...] = ()
    missing_monomer_components: Tuple[str, ...] = ()
    skipped_residues: Tuple[ResiduePreparationExclusion, ...] = ()
    protected_neighbor_residues: Tuple[ResiduePreparationExclusion, ...] = ()

    def excluded_donor_keys(self) -> Set[tuple[int, str, str, str]]:
        return {
            item.key
            for item in self.skipped_residues + self.protected_neighbor_residues
        }

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "h_mode": self.h_mode,
            "missing_monomer_components": list(
                self.missing_monomer_components),
            "skipped_residues": [
                item.to_dict() for item in self.skipped_residues
            ],
            "protected_neighbor_residues": [
                item.to_dict() for item in self.protected_neighbor_residues
            ],
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


def _synchronize_hydrogens(
        target: gemmi.Structure,
        source: gemmi.Structure,
        excluded_residue_keys: Optional[Set[tuple[int, str, str, str]]] = None,
        ) -> None:
    """Replace target H/D atoms with Gemmi's prepared result.

    Heavy atoms remain untouched so topology-recovery removals made in the
    working copy cannot delete experimental coordinates from the caller's
    structure.
    """
    excluded_residue_keys = excluded_residue_keys or set()
    target_residues = _residue_lookup(target)
    source_residues = _residue_lookup(source)

    for key, target_residue in target_residues.items():
        if key in excluded_residue_keys:
            hydrogen_indices = [
                index for index, atom in enumerate(target_residue)
                if atom.element.name.upper() in {"H", "D"}
            ]
            for index in reversed(hydrogen_indices):
                del target_residue[index]
            continue

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
    placement.  This low-level retry does not delete residues.
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


_RESIDUE_IN_ERROR = re.compile(
    r"(?P<chain>[^/\s]*)/(?P<resname>[^/\s]+)\s+"
    r"(?P<seqid>[^/\s]+)/(?P<atom>[^\s:]+)"
)


def _find_failed_residue(structure: gemmi.Structure, model_index: int,
                         message: Optional[str]):
    """Locate a unique residue explicitly identified by a Gemmi error."""
    if not message or not (0 <= model_index < len(structure)):
        return None
    model = structure[model_index]
    for match in _RESIDUE_IN_ERROR.finditer(message):
        chain_name = match.group("chain")
        residue_name = match.group("resname")
        seqid = match.group("seqid")
        for chain_index, chain in enumerate(model):
            if chain.name != chain_name:
                continue
            matches = [
                (residue_index, residue)
                for residue_index, residue in enumerate(chain)
                if residue.name == residue_name
                and str(residue.seqid).strip() == seqid
            ]
            if len(matches) == 1:
                residue_index, residue = matches[0]
                return chain_index, residue_index, residue, match.group("atom")
    return None


def _is_polymer_residue(residue: gemmi.Residue,
                        monomer_library: gemmi.MonLib) -> bool:
    if residue.entity_type == gemmi.EntityType.Polymer:
        return True
    if residue.name not in monomer_library.monomers:
        return False
    component = monomer_library.monomers[residue.name]
    return component.group in {
        gemmi.ChemComp.Group.Peptide,
        gemmi.ChemComp.Group.PPeptide,
        gemmi.ChemComp.Group.MPeptide,
        gemmi.ChemComp.Group.Dna,
        gemmi.ChemComp.Group.Rna,
        gemmi.ChemComp.Group.DnaRna,
    }


def _has_polymer_context(chain: gemmi.Chain, residue_index: int,
                         monomer_library: gemmi.MonLib) -> bool:
    """Treat an unknown residue between polymer residues conservatively."""
    residue = chain[residue_index]
    if _is_polymer_residue(residue, monomer_library):
        return True
    if residue.entity_type != gemmi.EntityType.Unknown:
        return False
    return any(
        _is_polymer_residue(chain[index], monomer_library)
        for index in (residue_index - 1, residue_index + 1)
        if 0 <= index < len(chain)
    )


def _remove_connections_for_residue(
        structure: gemmi.Structure,
        chain_name: str,
        residue: gemmi.Residue,
        ) -> int:
    """Remove only explicit links to a residue removed from the working copy."""
    seqid = str(residue.seqid).strip()

    def matches(partner: gemmi.AtomAddress) -> bool:
        return (
            partner.chain_name == chain_name
            and partner.res_id.name == residue.name
            and str(partner.res_id.seqid).strip() == seqid
        )

    removed = 0
    for index in reversed(range(len(structure.connections))):
        connection = structure.connections[index]
        if matches(connection.partner1) or matches(connection.partner2):
            structure.connections.pop(index)
            removed += 1
    return removed


def _prepare_model_with_quarantine(
        structure: gemmi.Structure,
        monomer_library: gemmi.MonLib,
        h_change: gemmi.HydrogenChange,
        model_index: int,
        ) -> tuple[
            TopologyPreparationReport,
            Tuple[ResiduePreparationExclusion, ...],
            Tuple[ResiduePreparationExclusion, ...],
        ]:
    """Prepare one model, quarantining only localizable failing residues."""
    skipped: Dict[tuple, ResiduePreparationExclusion] = {}
    protected: Dict[tuple, ResiduePreparationExclusion] = {}
    total_attempts = 0
    connections_cleared = False

    while True:
        report = _prepare_topology_with_retries(
            structure, monomer_library, h_change, model_index)
        total_attempts += report.attempts
        connections_cleared = (
            connections_cleared or report.connections_cleared)

        if report.status != "failed":
            partial = bool(skipped) or report.status == "partial"
            message = report.message
            if skipped:
                message = (
                    f"Skipped {len(skipped)} residue(s) after localized "
                    "topology failures; adjacent polymer donor hydrogens "
                    "were protected."
                )
            return (
                TopologyPreparationReport(
                    model_index=model_index,
                    status="partial" if partial else "success",
                    attempts=total_attempts,
                    connections_cleared=connections_cleared,
                    message=message,
                ),
                tuple(skipped.values()),
                tuple(protected.values()),
            )

        located = _find_failed_residue(
            structure, model_index, report.message)
        if located is None:
            return (
                TopologyPreparationReport(
                    model_index=model_index,
                    status="failed",
                    attempts=total_attempts,
                    connections_cleared=connections_cleared,
                    message=report.message,
                ),
                tuple(skipped.values()),
                tuple(protected.values()),
            )

        chain_index, residue_index, residue, atom_name = located
        chain = structure[model_index][chain_index]
        residue_key = _residue_key(model_index, chain, residue)
        if residue_key in skipped:
            return (
                TopologyPreparationReport(
                    model_index=model_index,
                    status="failed",
                    attempts=total_attempts,
                    connections_cleared=connections_cleared,
                    message=report.message,
                ),
                tuple(skipped.values()),
                tuple(protected.values()),
            )

        reason = report.message or "Localized topology preparation failure"
        skipped[residue_key] = ResiduePreparationExclusion(
            model_index, chain.name, str(residue.seqid).strip(),
            residue.name, reason)
        protected.pop(residue_key, None)
        logger.warning(
            "Skipping H preparation for model %d %s/%s %s (atom %s): %s",
            model_index, chain.name, residue.name,
            str(residue.seqid).strip(), atom_name, reason)

        if _has_polymer_context(chain, residue_index, monomer_library):
            for neighbor_index in (residue_index - 1, residue_index + 1):
                if not (0 <= neighbor_index < len(chain)):
                    continue
                neighbor = chain[neighbor_index]
                neighbor_key = _residue_key(model_index, chain, neighbor)
                if neighbor_key in skipped:
                    continue
                protected[neighbor_key] = ResiduePreparationExclusion(
                    model_index, chain.name, str(neighbor.seqid).strip(),
                    neighbor.name,
                    f"Adjacent to skipped polymer residue "
                    f"{chain.name}/{residue.name} "
                    f"{str(residue.seqid).strip()}; generated boundary "
                    "hydrogens are not trusted.",
                )

        removed_connections = _remove_connections_for_residue(
            structure, chain.name, residue)
        if removed_connections:
            logger.warning(
                "Removed %d explicit connection(s) to quarantined residue "
                "%s/%s %s from the H-preparation working copy.",
                removed_connections, chain.name, residue.name,
                str(residue.seqid).strip())
        del chain[residue_index]


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
        skipped_residues = []
        protected_neighbor_residues = []
        for model_index in range(len(working)):
            model_report, skipped, protected = _prepare_model_with_quarantine(
                working, monlib, h_change, model_index)
            model_reports.append(model_report)
            skipped_residues.extend(skipped)
            protected_neighbor_residues.extend(protected)

        failed = [item for item in model_reports if item.status == "failed"]
        if failed:
            report = HydrogenPreparationReport(
                "failed", h_change_val, tuple(model_reports),
                missing_monomers, tuple(skipped_residues),
                tuple(protected_neighbor_residues))
            details = "; ".join(
                f"model {item.model_index}: {item.message}"
                for item in failed)
            raise HydrogenPreparationError(
                f"Hydrogen topology preparation failed ({details})", report)

        excluded_residue_keys = {
            item.key
            for item in skipped_residues + protected_neighbor_residues
        }
        _synchronize_hydrogens(
            structure, working,
            excluded_residue_keys=excluded_residue_keys)

        structure.setup_cell_images()
        status = (
            "partial"
            if (missing_monomers or skipped_residues or
                protected_neighbor_residues or any(
                item.status == "partial" for item in model_reports)
            )
            else "success"
        )
        return HydrogenPreparationResult(
            structure,
            HydrogenPreparationReport(
                status, h_change_val, tuple(model_reports),
                missing_monomers, tuple(skipped_residues),
                tuple(protected_neighbor_residues)),
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
