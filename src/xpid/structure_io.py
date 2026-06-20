"""
Robust structure reading helpers.
"""
from __future__ import annotations

import logging
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Union

import gemmi

logger = logging.getLogger("xpid.structure_io")


def structure_cache_dir() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "xpid" / "structures"
    return Path.home() / ".cache" / "xpid" / "structures"


def _pdb_code_from_path(path: Path) -> Optional[str]:
    name = path.name.lower()
    for suffix in (".cif.gz", ".mmcif.gz", ".cif", ".mmcif", ".pdb.gz", ".pdb"):
        if name.endswith(suffix):
            code = name[:-len(suffix)]
            break
    else:
        code = path.stem.split(".")[0].lower()

    if len(code) == 4 and code.isalnum():
        return code
    return None


def _looks_like_bad_gzip(path: Path, exc: Exception) -> bool:
    if path.suffix.lower() != ".gz":
        return False

    message = str(exc).lower()
    return any(fragment in message for fragment in (
        "cannot determine uncompressed size",
        "unexpected end of file",
        "compressed file ended before the end-of-stream",
        "not a gzipped file",
    ))


def download_rcsb_mmcif(pdb_code: str, timeout: float = 60.0) -> Path:
    code = pdb_code.lower()
    target = structure_cache_dir() / f"{code}.cif.gz"
    target.parent.mkdir(parents=True, exist_ok=True)

    url = f"https://files.rcsb.org/download/{code.upper()}.cif.gz"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            target.write_bytes(response.read())
    except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
        raise RuntimeError(f"Could not download {code.upper()} from RCSB: {exc}") from exc

    logger.warning("Downloaded replacement mmCIF for %s to %s", code.upper(), target)
    return target


def read_structure(path: Union[str, Path], allow_remote_recovery: bool = True) -> gemmi.Structure:
    structure_path = Path(path)
    try:
        return gemmi.read_structure(str(structure_path))
    except RuntimeError as exc:
        if not _looks_like_bad_gzip(structure_path, exc):
            raise

        pdb_code = _pdb_code_from_path(structure_path)
        if not allow_remote_recovery:
            raise RuntimeError(
                f"Could not read {structure_path}: the local gzip file appears corrupted. "
                "Automatic RCSB recovery is disabled for batch inputs; replace the local "
                "file or run xpid on this single file to allow one-file recovery."
            ) from exc

        if not pdb_code:
            raise RuntimeError(
                f"Could not read {structure_path}: the gzip file appears corrupted, "
                "and no 4-character PDB code could be inferred for automatic recovery."
            ) from exc

        try:
            replacement = download_rcsb_mmcif(pdb_code)
            return gemmi.read_structure(str(replacement))
        except Exception as fallback_exc:
            raise RuntimeError(
                f"Could not read {structure_path}: the local gzip file appears corrupted "
                f"and automatic RCSB recovery for {pdb_code.upper()} failed: {fallback_exc}"
            ) from exc
