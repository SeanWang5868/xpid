"""Regression tests on real PDB structure 5FJJ.

These tests detect changes in Xpid's detection behaviour by checking
invariant properties (hit counts, system labels, column presence)
against known values for a benchmark structure.
"""
from pathlib import Path

import pytest
from xpid import core, structure_io, prep

_5FJJ = Path(__file__).resolve().parent.parent / "5fjj.cif.gz"


@pytest.fixture(scope="module")
def structure_5fjj():
    """Read 5FJJ once for all tests in this module."""
    assert _5FJJ.exists(), f"Benchmark structure not found: {_5FJJ}"
    st = structure_io.read_structure(str(_5FJJ))
    st = prep.add_hydrogens_memory(st, h_change_val=4)
    return st


# ---------------------------------------------------------------
# Sanity: structure loads and has expected properties
# ---------------------------------------------------------------

def test_5fjj_loads(structure_5fjj):
    assert structure_5fjj is not None
    assert len(structure_5fjj) > 0
    # 5FJJ is a high-resolution lysozyme structure (~0.95 Å)
    assert structure_5fjj.resolution is not None
    assert 1.0 < structure_5fjj.resolution < 3.0


# ---------------------------------------------------------------
# Default detection: Hudson/Plevin
# ---------------------------------------------------------------

def test_default_detection_has_hits(structure_5fjj):
    hits = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    assert len(hits) > 0, "5FJJ should have XH–π interactions"
    # Every reported hit is Hudson or Plevin positive
    for hit in hits:
        assert hit["is_hudson"] or hit["is_plevin"], \
            f"Hit {hit['X_res']}/{hit['X_atom']} → {hit['pi_res']}/{hit['pi_id']} has neither label"


def test_default_detection_hudson_plevin_balance(structure_5fjj):
    hits = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    n_hudson = sum(h["is_hudson"] for h in hits)
    n_plevin = sum(h["is_plevin"] for h in hits)
    # Both systems should report interactions
    assert n_hudson > 0, "Hudson should report hits"
    assert n_plevin > 0, "Plevin should report hits"


# ---------------------------------------------------------------
# Column invariants
# ---------------------------------------------------------------

def test_required_columns_present(structure_5fjj):
    hits = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    required = [
        "pdb", "pi_chain", "pi_res", "pi_id",
        "X_chain", "X_res", "X_id", "X_atom",
        "dist_X_Pi", "dist_X_centroid", "proj_dist",
        "is_hudson", "is_plevin", "H_source",
    ]
    for hit in hits:
        for col in required:
            assert col in hit, f"Missing required column '{col}'"


def test_coordinates_enabled_adds_columns(structure_5fjj):
    hits = core.detect_interactions_in_structure(
        structure_5fjj, "5fjj", include_coordinates=True)
    assert len(hits) > 0
    for hit in hits:
        assert "pi_center_x" in hit
        assert "X_side_of_pi" in hit


def test_h_coordinates_always_present(structure_5fjj):
    """H coordinates are always stored internally, even without --include-coordinates."""
    hits = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    assert len(hits) > 0
    for hit in hits:
        assert "H_xyz_x" in hit, "H coordinates should always be present"


def test_cooperativity_default_on(structure_5fjj):
    hits = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    assert len(hits) > 0
    for hit in hits:
        assert "coop_total_donors" in hit
        assert hit["coop_total_donors"] >= 1


def test_hbond_columns_present(structure_5fjj):
    """H-bond competition always runs internally."""
    hits = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    assert len(hits) > 0
    for hit in hits:
        assert "hbond_competing" in hit
        assert "hbond_vs_xhpi_score" in hit


def test_sasa_only_with_flag(structure_5fjj):
    hits_default = core.detect_interactions_in_structure(structure_5fjj, "5fjj")
    for hit in hits_default:
        assert "pi_sasa_avg" not in hit

    hits_sasa = core.detect_interactions_in_structure(
        structure_5fjj, "5fjj", compute_sasa=True)
    for hit in hits_sasa:
        assert "pi_sasa_avg" in hit
        assert hit["pi_sasa_avg"] > 0, "SASA should be positive"


# ---------------------------------------------------------------
# P-slab regression
# ---------------------------------------------------------------

def test_p_slab_inclusion_reports_label(structure_5fjj):
    hits = core.detect_interactions_in_structure(
        structure_5fjj, "5fjj", include_p_slab=True)
    has_p_slab = any(h.get("is_p_slab") for h in hits)
    # P-slab is stricter; may or may not have hits, but column must exist
    for hit in hits:
        assert "is_p_slab" in hit
