"""
CCP4 monomer-library discovery and download helpers.
"""
from __future__ import annotations

import logging
import os
import shutil
import sys
import urllib.error
import urllib.request
import zipfile
from functools import lru_cache
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger("xpid.monlib")

MONOMER_REPO_ZIPS = (
    ("master", "https://github.com/MonomerLibrary/monomers/archive/refs/heads/master.zip"),
    ("main", "https://github.com/MonomerLibrary/monomers/archive/refs/heads/main.zip"),
)


def _candidate_paths() -> List[Path]:
    paths: List[Path] = []

    for env_name in ("CLIBD_MON", "GEMMI_MON_LIB_PATH"):
        val = os.environ.get(env_name)
        if val:
            paths.append(Path(val).expanduser())

    clibd = os.environ.get("CLIBD")
    if clibd:
        paths.append(Path(clibd).expanduser() / "monomers")

    ccp4 = os.environ.get("CCP4")
    if ccp4:
        ccp4_root = Path(ccp4).expanduser()
        paths.extend([
            ccp4_root / "lib" / "data" / "monomers",
            ccp4_root / "lib" / "data" / "monomer",
        ])

    if sys.platform == "darwin":
        paths.extend([
            Path("/Applications/ccp4-9/lib/data/monomers"),
            Path("/Applications/ccp4-8/lib/data/monomers"),
            Path("/Applications/ccp4-7/lib/data/monomers"),
        ])

    paths.append(managed_monomer_dir())

    unique: List[Path] = []
    seen = set()
    for path in paths:
        key = str(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def managed_monomer_dir() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "xpid" / "ccp4-monomers"

    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base / "xpid" / "ccp4-monomers"


def is_monomer_library(path: Optional[Path]) -> bool:
    if not path or not path.is_dir():
        return False
    return (
        (path / "ener_lib.cif").is_file()
        and any((path / name).is_dir() for name in ("a", "p", "list"))
        and any((path / "a").glob("ALA.cif"))
    )


@lru_cache(maxsize=1)
def find_installed_monomer_library() -> Optional[str]:
    for path in _candidate_paths():
        if is_monomer_library(path):
            return str(path.resolve())
    return None


@lru_cache(maxsize=1)
def get_monomer_library_path() -> str:
    local = find_installed_monomer_library()
    if local:
        return local

    downloaded = download_ccp4_monomer_library()
    if downloaded and is_monomer_library(downloaded):
        return str(downloaded.resolve())

    raise RuntimeError(
        "Could not find or download the CCP4 monomer library. "
        "Check your network connection or install CCP4."
    )


def download_ccp4_monomer_library(timeout: float = 120.0) -> Optional[Path]:
    """Download the official MonomerLibrary/monomers repository into the user cache."""
    target = managed_monomer_dir()
    if is_monomer_library(target):
        return target

    parent = target.parent
    parent.mkdir(parents=True, exist_ok=True)
    tmp_extract = parent / "monomers-download"

    try:
        last_error: Optional[Exception] = None
        extracted: Optional[Path] = None
        archive = parent / "monomers.zip"

        for branch, url in MONOMER_REPO_ZIPS:
            try:
                with urllib.request.urlopen(url, timeout=timeout) as response:
                    archive.write_bytes(response.read())

                if tmp_extract.exists():
                    shutil.rmtree(tmp_extract)
                tmp_extract.mkdir(parents=True)

                with zipfile.ZipFile(archive) as zf:
                    zf.extractall(tmp_extract)

                candidate = tmp_extract / f"monomers-{branch}"
                if is_monomer_library(candidate):
                    extracted = candidate
                    break
                last_error = RuntimeError(f"{url} does not look like a CCP4 monomer library")
            except (OSError, RuntimeError, zipfile.BadZipFile, urllib.error.URLError, urllib.error.HTTPError) as exc:
                last_error = exc

        if extracted is None:
            raise RuntimeError(str(last_error) if last_error else "No downloadable archive was available")

        if target.exists():
            shutil.rmtree(target)
        extracted.replace(target)
        logger.info(f"Downloaded CCP4 monomer library to {target}")
        return target

    except (OSError, RuntimeError, zipfile.BadZipFile, urllib.error.URLError, urllib.error.HTTPError) as exc:
        logger.error(f"Could not download CCP4 monomer library: {exc}")
        return None
    finally:
        try:
            (parent / "monomers.zip").unlink(missing_ok=True)
        except OSError:
            pass
        if tmp_extract.exists():
            shutil.rmtree(tmp_extract, ignore_errors=True)


@lru_cache(maxsize=None)
def find_monomer_cif(code: str) -> Optional[Path]:
    root = Path(get_monomer_library_path())
    code = code.upper()
    for candidate in (
        root / code[0].lower() / f"{code}.cif",
        root / code[0].lower() / f"{code}.cif.gz",
        root / f"{code}.cif",
        root / f"{code}.cif.gz",
    ):
        if candidate.is_file():
            return candidate
    return None


def clear_path_caches() -> None:
    """Clear monomer discovery caches after environment/path changes."""
    find_monomer_cif.cache_clear()
    get_monomer_library_path.cache_clear()
    find_installed_monomer_library.cache_clear()
