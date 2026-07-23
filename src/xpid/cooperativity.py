"""
cooperativity.py
Post-processing analysis of cooperative XH–π interactions.

Groups hits by π-ring and annotates each with donor counts per face,
total donors, and a bi-face cooperativity flag.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Any, Optional, Tuple


def _ring_key(hit: Dict[str, Any]) -> Tuple:
    """Unique key for a π-ring across all hits."""
    return (
        hit.get("pdb", ""),
        hit.get("model", ""),
        hit.get("pi_chain", ""),
        hit.get("pi_res", ""),
        hit.get("pi_id", ""),
    )


def _donor_key(hit: Dict[str, Any]) -> Tuple:
    """Unique key for a donor within a ring's hits."""
    return (
        hit.get("X_chain", ""),
        hit.get("X_res", ""),
        hit.get("X_id", ""),
        hit.get("X_atom", ""),
        hit.get("sym_op", 0),
    )


def _face(hit: Dict[str, Any]) -> int:
    """Return which face the donor is on: +1, -1, or 0 (unknown)."""
    side = hit.get("X_side_of_pi")
    if side is not None:
        return int(side)
    # Fallback: use signed plane distance via dist_X_Pi
    # Without pi_normal direction we can't determine sign;
    # treat as unknown.
    return 0


def annotate_cooperativity(
    hits: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Annotate each hit with cooperativity metrics.

    Appends the following columns to each hit:

    ``coop_same_face_donors``
        Number of distinct donors on the same face (including self).
    ``coop_opp_face_donors``
        Number of distinct donors on the opposite face.
    ``coop_total_donors``
        Total distinct donors interacting with this π-ring.
    ``coop_bi_face``
        1 if donors occupy both faces, 0 otherwise.

    Returns the modified list (in-place).
    """
    # Group hits by ring
    ring_groups: Dict[Tuple, List[Dict[str, Any]]] = defaultdict(list)
    for hit in hits:
        ring_groups[_ring_key(hit)].append(hit)

    # Pre-compute per-ring donor counts by face
    ring_stats: Dict[Tuple, Dict[str, Any]] = {}
    for rkey, ring_hits in ring_groups:
        positive_donors: set = set()
        negative_donors: set = set()
        unknown_donors: set = set()

        for hit in ring_hits:
            dkey = _donor_key(hit)
            f = _face(hit)
            if f > 0:
                positive_donors.add(dkey)
            elif f < 0:
                negative_donors.add(dkey)
            else:
                # Unknown face: count in positive as a conservative default
                positive_donors.add(dkey)

        total = len(positive_donors | negative_donors | unknown_donors)
        pos_count = len(positive_donors)
        neg_count = len(negative_donors)

        ring_stats[rkey] = {
            "pos_count": pos_count,
            "neg_count": neg_count,
            "total": total,
            "bi_face": 1 if (pos_count > 0 and neg_count > 0) else 0,
        }

    # Annotate hits
    for hit in hits:
        rkey = _ring_key(hit)
        stats = ring_stats.get(rkey, {})
        f = _face(hit)

        if f > 0:
            same = stats.get("pos_count", 1)
            opp = stats.get("neg_count", 0)
        elif f < 0:
            same = stats.get("neg_count", 1)
            opp = stats.get("pos_count", 0)
        else:
            same = stats.get("pos_count", 1)
            opp = stats.get("neg_count", 0)

        hit["coop_same_face_donors"] = same
        hit["coop_opp_face_donors"] = opp
        hit["coop_total_donors"] = stats.get("total", 1)
        hit["coop_bi_face"] = stats.get("bi_face", 0)

    return hits
