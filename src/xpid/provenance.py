"""
provenance.py
Write a companion _metadata.json file recording all run-time parameters
for reproducibility.
"""
from __future__ import annotations

import json
import sys
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, List, Mapping, Optional

import gemmi

from ._version import __version__


def build_metadata(
    args: object,
    output_path: Path,
    monomer_lib_path: str,
    file_count: int,
    input_resolution: Optional[Mapping[str, int]] = None,
) -> Dict[str, Any]:
    """Build a dictionary of provenance information from the argparse namespace."""
    metadata = {
        "tool": "xpid",
        "version": __version__,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "gemmi_version": gemmi.__version__ if hasattr(gemmi, "__version__") else "unknown",
        "monomer_library": monomer_lib_path,
        "output": str(output_path.resolve()),
        "file_count": file_count,
        "parameters": {
            "h_mode": getattr(args, "h_mode", 4),
            "file_type": getattr(args, "file_type", "json"),
            "verbose": getattr(args, "verbose", False),
            "include_coordinates": getattr(args, "include_coordinates", False),
            "sasa": getattr(args, "sasa", False),
            "cone_mode": (
                "none" if getattr(args, "no_cone", False) else "auto"),
            "include_p_slab": getattr(args, "include_p_slab", False),
            "xh_candidates": getattr(args, "report_xh_candidates", False),
            "include_water": getattr(args, "include_water", False),
            "sym_contacts": getattr(args, "sym_contacts", False),
            "max_b": getattr(args, "max_b", 0.0),
            "min_occ": getattr(args, "min_occ", 0.0),
            "model": getattr(args, "model", "0"),
            "jobs": getattr(args, "jobs", 1),
        },
    }
    if input_resolution is not None:
        metadata["input_resolution"] = dict(input_resolution)
    return metadata


def write_metadata(metadata: Dict[str, Any], output_dir: Path, stem: str) -> Optional[Path]:
    """Write *metadata* as ``{stem}_metadata.json`` in *output_dir*.

    Returns the path written, or ``None`` on failure.
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / f"{stem}_metadata.json"
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(metadata, fh, indent=2, default=str)
        return path
    except OSError:
        return None


def write_diagnostics(records: List[Dict[str, Any]], output_dir: Path,
                      stem: str) -> Optional[Path]:
    """Atomically write per-structure chemical-preparation diagnostics."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / f"{stem}_diagnostics.json"
        temporary = output_dir / f".{stem}_diagnostics.json.partial"
        payload = {
            "tool": "xpid",
            "version": __version__,
            "structures": records,
        }
        with open(temporary, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, default=str)
        temporary.replace(path)
        return path
    except OSError:
        return None
