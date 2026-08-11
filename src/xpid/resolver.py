"""
resolver.py
Input file resolution: PDB list parsing, local mirror lookup, directory scanning.
"""
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger("xpid.resolver")


@dataclass(frozen=True)
class InputResolution:
    """Resolved structure files and auditable PDB-list source counts."""

    files: List[Path]
    pdb_list_entries: int = 0
    pdb_redo: int = 0
    standard_pdb: int = 0
    missing: int = 0
    missing_codes: tuple[str, ...] = ()
    identities: Dict[Path, Tuple[str, str]] = field(default_factory=dict)

    def provenance_counts(self) -> Dict[str, int]:
        return {
            "pdb_list_entries": self.pdb_list_entries,
            "pdb_redo_mirror": self.pdb_redo,
            "standard_pdb_mirror": self.standard_pdb,
            "missing": self.missing,
        }

    def identity_for(self, path: Path) -> Tuple[str, str]:
        """Return the canonical structure identifier and input source."""
        resolved = path.resolve()
        return self.identities.get(
            resolved, (_identifier_from_path(resolved), "direct"))


def _identifier_from_path(path: Path) -> str:
    """Derive a stable identifier for a direct, non-mirror input."""
    name = path.name
    lowered = name.lower()
    for suffix in (
            ".cif.gz", ".mmcif.gz", ".pdb.gz", ".ent.gz",
            ".cif", ".mmcif", ".pdb", ".ent"):
        if lowered.endswith(suffix):
            name = name[:-len(suffix)]
            break
    match = re.fullmatch(r"([0-9][A-Za-z0-9]{3})_final", name)
    if match:
        return match.group(1).lower()
    return name.lower()


def parse_pdb_list_file(list_path: Path) -> Set[str]:
    """Read a text file and extract 4-character PDB codes (comma or newline separated)."""
    codes: Set[str] = set()
    try:
        content = list_path.read_text(encoding='utf-8')
        tokens = re.split(r'[,\s\n\t]+', content)
        for t in tokens:
            clean = t.strip()
            if len(clean) == 4 and clean.isalnum():
                codes.add(clean.lower())
    except Exception as e:
        logger.error(f"Failed to read PDB list file: {e}")
        sys.exit(1)
    return codes


def resolve_mirror_path(mirror_root: Path, pdb_code: str) -> Optional[Path]:
    """Resolve a PDB code to a file using the standard divided directory layout.

    Example: ``1fn0`` -> ``{mirror_root}/fn/1fn0.cif.gz``
    """
    middle = pdb_code[1:3].lower()

    candidates = [
        mirror_root / middle / f"{pdb_code}.cif.gz",
        mirror_root / middle / f"{pdb_code}.cif",
        mirror_root / f"{pdb_code}.cif.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def resolve_redo_path(mirror_root: Path, pdb_code: str) -> Optional[Path]:
    """Resolve a PDB code to a PDB-REDO mirror file.

    Tries official rsync layout, simplified layout, and standard PDB layout.
    """
    middle = pdb_code[1:3].lower()

    candidates = [
        # Official PDB-REDO rsync layout: /ab/1abc/1abc_final.cif
        mirror_root / middle / pdb_code / f"{pdb_code}_final.cif",
        mirror_root / middle / pdb_code / f"{pdb_code}_final.cif.gz",
        # Simplified layout: /ab/1abc_final.cif
        mirror_root / middle / f"{pdb_code}_final.cif",
        mirror_root / middle / f"{pdb_code}_final.cif.gz",
        # Standard PDB layout: /ab/1abc.cif
        mirror_root / middle / f"{pdb_code}.cif",
        mirror_root / middle / f"{pdb_code}.cif.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def resolve_inputs(inputs: List[str], pdb_list: str,
                   pdb_mirror: str, redo_mirror: str) -> InputResolution:
    """Resolve inputs and retain source counts for logging and provenance.

    PDB-REDO mirror is prioritized over the standard PDB mirror when both
    are provided.
    """
    final_files: Set[Path] = set()
    identities: Dict[Path, Tuple[str, str]] = {}

    # 1. Direct file or directory inputs
    if inputs:
        structure_suffixes = (
            ".cif", ".mmcif", ".pdb", ".ent",
            ".cif.gz", ".mmcif.gz", ".pdb.gz", ".ent.gz",
        )
        for inp in inputs:
            path = Path(inp)
            if path.is_file():
                resolved = path.resolve()
                final_files.add(resolved)
                identities[resolved] = (
                    _identifier_from_path(resolved), "direct")
            elif path.is_dir():
                for p in path.rglob("*"):
                    if p.is_file() and p.name.lower().endswith(structure_suffixes):
                        resolved = p.resolve()
                        final_files.add(resolved)
                        identities[resolved] = (
                            _identifier_from_path(resolved), "directory")

    # 2. PDB list + mirror resolution
    if pdb_list:
        if not pdb_mirror and not redo_mirror:
            logger.error("--pdb-mirror or --redo-mirror is REQUIRED when using --pdb-list.")
            sys.exit(1)

        pdb_root = Path(pdb_mirror).resolve() if pdb_mirror else None
        redo_root = Path(redo_mirror).resolve() if redo_mirror else None

        if pdb_root and not pdb_root.exists():
            logger.error(f"PDB mirror directory not found: {pdb_root}")
            sys.exit(1)
        if redo_root and not redo_root.exists():
            logger.error(f"PDB-REDO mirror directory not found: {redo_root}")
            sys.exit(1)

        codes = parse_pdb_list_file(Path(pdb_list))
        found_redo = 0
        found_pdb = 0
        missing: List[str] = []

        for code in codes:
            fpath = None
            source = None

            # Priority 1: PDB-REDO
            if redo_root:
                fpath = resolve_redo_path(redo_root, code)
                if fpath:
                    found_redo += 1
                    source = "pdb_redo"

            # Priority 2: Standard PDB
            if not fpath and pdb_root:
                fpath = resolve_mirror_path(pdb_root, code)
                if fpath:
                    found_pdb += 1
                    source = "pdb_mirror"

            if fpath:
                resolved = fpath.resolve()
                final_files.add(resolved)
                identities[resolved] = (code, source or "pdb_mirror")
            else:
                missing.append(code)

        return InputResolution(
            files=sorted(final_files),
            pdb_list_entries=len(codes),
            pdb_redo=found_redo,
            standard_pdb=found_pdb,
            missing=len(missing),
            missing_codes=tuple(sorted(missing)),
            identities=identities,
        )

    return InputResolution(
        files=sorted(final_files), identities=identities)


def gather_inputs(inputs: List[str], pdb_list: str,
                  pdb_mirror: str, redo_mirror: str) -> List[Path]:
    """Backward-compatible file-list-only input resolver."""
    return resolve_inputs(inputs, pdb_list, pdb_mirror, redo_mirror).files
