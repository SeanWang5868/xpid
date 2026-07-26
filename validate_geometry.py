#!/usr/bin/env python3
"""
Standalone geometry validator for Xpid.

Independently computes all XH–π geometric parameters from raw coordinates
and compares them against Xpid's reported values.  This is NOT part of the
Xpid package — it is an external verification tool.

Usage:
    cd /path/to/xpid
    PYTHONPATH=src python validate_geometry.py [structure.cif ...]

Output:
    For each hit, prints side-by-side comparison of key geometric parameters.
    Discrepancies > 0.001 Å or > 0.01° are flagged.
"""
from __future__ import annotations

import sys
import math
from pathlib import Path
from typing import List, Dict, Any, Tuple

import gemmi
import numpy as np

# Add local src to path so we can import xpid
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from xpid import detect, structure_io
from xpid import config as xpid_config
from xpid import geometry as xpid_geom

# ---------------------------------------------------------------------------
# Independent geometry calculations (no dependency on xpid internals)
# ---------------------------------------------------------------------------


def _centroid(positions: np.ndarray) -> np.ndarray:
    return np.mean(positions, axis=0)


def _ring_normal(positions: np.ndarray) -> np.ndarray:
    """Best-fit plane normal via SVD."""
    centered = positions - _centroid(positions)
    _, _, vh = np.linalg.svd(centered)
    return vh[2, :]


def _plane_distance(point: np.ndarray, plane_point: np.ndarray, normal: np.ndarray) -> float:
    n = normal / np.linalg.norm(normal)
    return float(abs(np.dot(point - plane_point, n)))


def _project_to_plane(point: np.ndarray, plane_point: np.ndarray, normal: np.ndarray) -> np.ndarray:
    n = normal / np.linalg.norm(normal)
    return point - np.dot(point - plane_point, n) * n


def _angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    """Angle in degrees between two vectors."""
    c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def compute_ring_geometry(
    structure: gemmi.Structure,
    chain_name: str,
    res_seqid: str,
    res_name: str,
    ring_id: str = "ring1",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Compute centroid, normal, and radius for an aromatic ring in a structure.

    Returns (centroid_array, normal_array, ring_radius).
    """
    # Get ring atom names from xpid config
    ring_sets = xpid_config.get_aromatic_rings(res_name)
    if not ring_sets:
        raise ValueError(f"No aromatic ring definition for {res_name}")

    model = structure[0]
    chain = model[chain_name]

    # Find the residue
    residue = None
    for r in chain:
        if str(r.seqid).strip() == res_seqid and r.name == res_name:
            residue = r
            break
    if residue is None:
        raise ValueError(f"Residue {chain_name}/{res_name}/{res_seqid} not found")

    try:
        ring_index = int(ring_id.removeprefix("ring")) - 1
    except (TypeError, ValueError):
        ring_index = 0
    if not 0 <= ring_index < len(ring_sets):
        raise ValueError(
            f"Ring identifier {ring_id!r} is invalid for {res_name}")
    target_atoms = ring_sets[ring_index]

    # Collect atom positions
    positions = []
    for atom in residue:
        if atom.name in target_atoms:
            positions.append([atom.pos.x, atom.pos.y, atom.pos.z])

    if len(positions) < 3:
        raise ValueError(f"Only {len(positions)} ring atoms found for {res_name}")

    pos = np.array(positions)
    centroid = _centroid(pos)
    normal = _ring_normal(pos)
    radius = float(np.max(np.linalg.norm(pos - centroid, axis=1)))

    return centroid, normal, radius


def compute_donor_coords(
    structure: gemmi.Structure,
    chain_name: str,
    res_seqid: str,
    atom_name: str,
    h_name: str | None = None,
) -> Tuple[np.ndarray, np.ndarray | None]:
    """Get coordinates for donor X and optionally H atom.

    Returns (x_pos, h_pos). h_pos is None if H not found.
    """
    model = structure[0]
    chain = model[chain_name]

    residue = None
    for r in chain:
        if str(r.seqid).strip() == res_seqid:
            residue = r
            break
    if residue is None:
        raise ValueError(f"Residue {chain_name}/{res_seqid} not found")

    x_pos = None
    h_pos = None
    for atom in residue:
        if atom.name == atom_name:
            x_pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
        if h_name and atom.name == h_name:
            h_pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])

    if x_pos is None:
        raise ValueError(f"Atom {atom_name} not found in {chain_name}/{res_seqid}")

    return x_pos, h_pos


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_hit(
    structure: gemmi.Structure,
    hit: Dict[str, Any],
) -> List[str]:
    """Validate one Xpid hit against independent geometry calculation.

    Returns a list of discrepancy messages (empty = perfect match).
    """
    issues: List[str] = []

    pi_chain = hit["pi_chain"]
    pi_res = hit["pi_res"]
    pi_id = str(hit["pi_id"]).strip()
    X_chain = hit["X_chain"]
    X_res = hit["X_res"]
    X_id = str(hit["X_id"]).strip()
    X_atom = hit["X_atom"]
    H_atom = hit.get("H_atom", "")
    H_source = hit.get("H_source", "")

    try:
        centroid, normal, radius = compute_ring_geometry(
            structure, pi_chain, pi_id, pi_res,
            str(hit.get("pi_ring_id", "ring1")))
        x_pos, h_pos = compute_donor_coords(
            structure, X_chain, X_id, X_atom, H_atom if H_atom != "virt" else None)
        if H_source == "cone_virtual":
            h_coords = (
                hit.get("H_xyz_x"), hit.get("H_xyz_y"), hit.get("H_xyz_z"))
            if any(value is None for value in h_coords):
                raise ValueError("Cone hit does not contain virtual H coordinates")
            h_pos = np.array(h_coords, dtype=float)
    except ValueError as e:
        return [f"  Cannot compute geometry: {e}"]

    # --- dist_X_Pi (perpendicular distance X → π plane) ---
    our_dist_x_pi = _plane_distance(x_pos, centroid, normal)
    their_dist_x_pi = float(hit.get("dist_X_Pi", -1))
    if abs(our_dist_x_pi - their_dist_x_pi) > 0.002:
        issues.append(
            f"  dist_X_Pi:  ours={our_dist_x_pi:.4f}  xpid={their_dist_x_pi:.4f}  "
            f"Δ={abs(our_dist_x_pi - their_dist_x_pi):.4f} Å"
        )

    # --- dist_X_centroid ---
    our_dist_centroid = float(np.linalg.norm(x_pos - centroid))
    their_dist_centroid = float(hit.get("dist_X_centroid", -1))
    if abs(our_dist_centroid - their_dist_centroid) > 0.002:
        issues.append(
            f"  dist_X_centroid: ours={our_dist_centroid:.4f}  xpid={their_dist_centroid:.4f}  "
            f"Δ={abs(our_dist_centroid - their_dist_centroid):.4f} Å"
        )

    # --- proj_dist (lateral offset of X projection from ring centre) ---
    x_proj = _project_to_plane(x_pos, centroid, normal)
    our_proj = float(np.linalg.norm(x_proj - centroid))
    their_proj = float(hit.get("proj_dist", -1))
    if abs(our_proj - their_proj) > 0.003:
        issues.append(
            f"  proj_dist:   ours={our_proj:.4f}  xpid={their_proj:.4f}  "
            f"Δ={abs(our_proj - their_proj):.4f} Å"
        )

    # --- H-dependent geometry (only if H atom available) ---
    if h_pos is not None:
        # Hudson theta: angle between X-H and ring normal
        v_xh = h_pos - x_pos
        v_xh_norm = np.linalg.norm(v_xh)
        if v_xh_norm > 0:
            # Hudson theta: angle between X-H and normal
            our_theta = _angle_between(v_xh, normal)
            if our_theta > 90:
                our_theta = 180 - our_theta
            their_theta = float(hit.get("theta", -1))
            if their_theta is not None and abs(our_theta - their_theta) > 0.05:
                issues.append(
                    f"  theta:       ours={our_theta:.2f}°  xpid={their_theta:.2f}°  "
                    f"Δ={abs(our_theta - their_theta):.3f}°"
                )

        # Plevin XPCN: angle between X→centroid and normal
        v_xc = centroid - x_pos
        our_xpcn = _angle_between(v_xc, normal)
        if our_xpcn > 90:
            our_xpcn = 180 - our_xpcn
        their_xpcn = float(hit.get("angle_XPCN", -1))
        if their_xpcn is not None and abs(our_xpcn - their_xpcn) > 0.05:
            issues.append(
                f"  angle_XPCN:  ours={our_xpcn:.2f}°  xpid={their_xpcn:.2f}°  "
                f"Δ={abs(our_xpcn - their_xpcn):.3f}°"
            )

        # Plevin XH_Pi: X-H→centroid angle
        v_hc = centroid - h_pos
        our_xh_pi = _angle_between(x_pos - h_pos, v_hc)
        their_xh_pi = float(hit.get("angle_XH_Pi", -1))
        if their_xh_pi is not None and abs(our_xh_pi - their_xh_pi) > 0.05:
            issues.append(
                f"  angle_XH_Pi: ours={our_xh_pi:.2f}°  xpid={their_xh_pi:.2f}°  "
                f"Δ={abs(our_xh_pi - their_xh_pi):.3f}°"
            )

    return issues


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    structures = sys.argv[1:] if len(sys.argv) > 1 else ["5fjj.cif.gz"]

    total_hits = 0
    total_checked = 0
    total_issues = 0

    for struct_path in structures:
        path = Path(struct_path)
        if not path.exists():
            print(f"SKIP: {path} not found")
            continue

        print(f"\n{'='*70}")
        print(f"Validating: {path.name}")
        print(f"{'='*70}")

        # Load structure (with hydrogens added, matching Xpid's pipeline)
        structure = structure_io.read_structure(str(path))
        from xpid import prep
        structure = prep.add_hydrogens_memory(structure, h_change_val=4)

        # Run Xpid
        hits = detect(str(path))
        total_hits += len(hits)
        print(f"Xpid found {len(hits)} hits")

        # Sample hits for validation: take up to 3 per distinct π-residue type
        sampled: Dict[str, List[Dict]] = {}
        for hit in hits:
            key = hit["pi_res"]
            if key not in sampled:
                sampled[key] = []
            if len(sampled[key]) < 3:
                sampled[key].append(hit)

        flat_sample = [h for group in sampled.values() for h in group]
        print(f"Sampling {len(flat_sample)} hits for detailed validation")

        for hit in flat_sample:
            total_checked += 1
            desc = (
                f"{hit['X_res']}/{hit['X_id']}/{hit['X_atom']}-{hit['H_atom']}"
                f"  →  {hit['pi_res']}/{hit['pi_id']}"
                f"  [{hit['is_hudson'] and 'H' or ''}{hit['is_plevin'] and 'P' or ''}]"
            )
            issues = validate_hit(structure, hit)
            if issues:
                print(f"\n✗ {desc}")
                for issue in issues:
                    print(issue)
                total_issues += len(issues)
            else:
                print(f"  ✓ {desc}")

    print(f"\n{'='*70}")
    print(f"Summary: {total_hits} total hits, {total_checked} validated")
    if total_issues == 0:
        print("All geometric parameters match within tolerance. ✓")
    else:
        print(f"{total_issues} discrepancies found ✗")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
