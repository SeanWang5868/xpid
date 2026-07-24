"""Tests for sasa, cooperativity, and hbond modules."""

import math
import numpy as np
import gemmi
from xpid import sasa, cooperativity, hbond, core


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _atom(name, element, xyz, b_iso=10.0, occ=1.0, altloc="\0"):
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element(element)
    atom.pos = gemmi.Position(*xyz)
    atom.occ = occ
    atom.b_iso = b_iso
    atom.altloc = altloc
    return atom


def _seqid(num):
    return gemmi.SeqId(str(num))


def _make_ala_structure():
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    m = gemmi.Model("1")
    c = gemmi.Chain("A")
    r = gemmi.Residue()
    r.name = "ALA"
    r.seqid = _seqid(1)
    r.add_atom(_atom("CA", "C", (0, 0, 0)))
    r.add_atom(_atom("CB", "C", (1.5, 0, 0)))
    c.add_residue(r)
    m.add_chain(c)
    st.add_model(m)
    return st


# TODO: Add parent atoms (CB/CA) to synthetic residues so cone_mode="auto" works.
# Then remove explicit cone_mode="none" from all test calls.

def _make_phe_ser_structure():
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    phe = gemmi.Residue()
    phe.name = "PHE"
    phe.seqid = _seqid(1)
    ring = {
        "CG": (1.4, 0.0, 0.0), "CD1": (0.7, 1.212, 0.0),
        "CE1": (-0.7, 1.212, 0.0), "CZ": (-1.4, 0.0, 0.0),
        "CE2": (-0.7, -1.212, 0.0), "CD2": (0.7, -1.212, 0.0),
    }
    for name, xyz in ring.items():
        phe.add_atom(_atom(name, "C", xyz))
    chain.add_residue(phe)

    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 2.0)))
    chain.add_residue(ser)

    model.add_chain(chain)
    st.add_model(model)
    return st


# ===================================================================
# SASA tests
# ===================================================================


def test_sasa_isolated_atom_full_exposure():
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    m = gemmi.Model('1')
    c = gemmi.Chain('A')
    r = gemmi.Residue()
    r.name = 'ALA'
    r.seqid = _seqid(1)
    r.add_atom(_atom('CA', 'C', (0, 0, 0)))
    c.add_residue(r)
    m.add_chain(c)
    st.add_model(m)
    result = sasa.compute_sasa(st, n_points=256)
    assert len(result) == 1
    for k, v in result.items():
        r_vdw = 1.70
        R = r_vdw + 1.4
        max_sasa = 4 * math.pi * R * R
        assert v > 0.95 * max_sasa, f'Atom should be nearly fully exposed, got {v} vs {max_sasa}'




def test_sasa_two_contacting_atoms_partial_burial():
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    m = gemmi.Model("1")
    c = gemmi.Chain("A")
    r = gemmi.Residue()
    r.name = "ALA"
    r.seqid = _seqid(1)
    # Two carbons placed 2.0 Å apart (vdW radius 1.7, so significant overlap)
    r.add_atom(_atom("C1", "C", (0, 0, 0)))
    r.add_atom(_atom("C2", "C", (2.0, 0, 0)))
    c.add_residue(r)
    m.add_chain(c)
    st.add_model(m)

    result = sasa.compute_sasa(st, n_points=512)
    for k, v in result.items():
        r_vdw = 1.70
        R = r_vdw + 1.4
        max_sasa = 4 * math.pi * R * R
        # Should be significantly buried
        assert v < 0.85 * max_sasa, f"Atom should be partially buried, got {v} vs {max_sasa}"


def test_atom_sasa_lookup():
    st = _make_ala_structure()
    sasa_map = sasa.compute_sasa(st)
    v = sasa.atom_sasa(sasa_map, 0, "A", "1", 0)
    assert v is not None
    assert v > 0


def test_average_ring_sasa():
    st = _make_phe_ser_structure()
    sasa_map = sasa.compute_sasa(st)
    ring_indices = sasa.residue_atom_indices(st[0]["A"][0])
    avg = sasa.average_ring_sasa(sasa_map, 0, "A", st[0]["A"][0], ring_indices)
    assert avg > 0


def test_sasa_compute_returns_empty_for_empty_structure():
    st = gemmi.Structure()
    result = sasa.compute_sasa(st)
    assert result == {}


def test_sasa_cli_flag_integration():
    st = _make_phe_ser_structure()
    hits = core.detect_interactions_in_structure(st, "test", compute_sasa=True)
    for hit in hits:
        assert "pi_sasa_avg" in hit
        assert "X_sasa" in hit
        assert "H_sasa" in hit
        assert hit["pi_sasa_avg"] > 0


# ===================================================================
# Cooperativity tests
# ===================================================================


def test_cooperativity_single_donor():
    hits = [{
        "pdb": "test", "model": "1",
        "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
        "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG",
        "sym_op": 0, "X_side_of_pi": 1,
    }]
    cooperativity.annotate_cooperativity(hits)
    assert hits[0]["coop_same_face_donors"] == 1
    assert hits[0]["coop_opp_face_donors"] == 0
    assert hits[0]["coop_total_donors"] == 1
    assert hits[0]["coop_bi_face"] == 0


def test_cooperativity_same_face_two_donors():
    hits = [
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG", "sym_op": 0, "X_side_of_pi": 1},
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "LYS", "X_id": "3", "X_atom": "NZ", "sym_op": 0, "X_side_of_pi": 1},
    ]
    cooperativity.annotate_cooperativity(hits)
    for hit in hits:
        assert hit["coop_same_face_donors"] == 2
        assert hit["coop_opp_face_donors"] == 0
        assert hit["coop_total_donors"] == 2
        assert hit["coop_bi_face"] == 0


def test_cooperativity_bi_face():
    hits = [
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG", "sym_op": 0, "X_side_of_pi": 1},
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "B", "X_res": "LYS", "X_id": "3", "X_atom": "NZ", "sym_op": 0, "X_side_of_pi": -1},
    ]
    cooperativity.annotate_cooperativity(hits)
    assert hits[0]["coop_same_face_donors"] == 1
    assert hits[0]["coop_opp_face_donors"] == 1
    assert hits[0]["coop_total_donors"] == 2
    assert hits[0]["coop_bi_face"] == 1
    assert hits[1]["coop_same_face_donors"] == 1
    assert hits[1]["coop_opp_face_donors"] == 1
    assert hits[1]["coop_bi_face"] == 1


def test_cooperativity_separate_rings():
    hits = [
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG", "sym_op": 0, "X_side_of_pi": 1},
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "TYR", "pi_id": "5",
         "X_chain": "B", "X_res": "LYS", "X_id": "3", "X_atom": "NZ", "sym_op": 0, "X_side_of_pi": 1},
    ]
    cooperativity.annotate_cooperativity(hits)
    # Each ring has only one donor
    for hit in hits:
        assert hit["coop_total_donors"] == 1


def test_cooperativity_dedup_same_donor():
    """Same donor from different altloc should count as one."""
    hits = [
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG", "sym_op": 0, "X_side_of_pi": 1},
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG", "sym_op": 0, "X_side_of_pi": 1},
    ]
    cooperativity.annotate_cooperativity(hits)
    for hit in hits:
        assert hit["coop_total_donors"] == 1


def test_cooperativity_unknown_face_defaults_to_positive():
    hits = [
        {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
         "X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG", "sym_op": 0, "X_side_of_pi": 0},
    ]
    cooperativity.annotate_cooperativity(hits)
    assert hits[0]["coop_same_face_donors"] == 1
    assert hits[0]["coop_opp_face_donors"] == 0
    assert hits[0]["coop_bi_face"] == 0


def test_cooperativity_cli_flag_integration():
    st = _make_phe_ser_structure()
    hits = core.detect_interactions_in_structure(st, "test", annotate_cooperativity=True)
    for hit in hits:
        assert "coop_same_face_donors" in hit
        assert "coop_bi_face" in hit


# ===================================================================
# H-bond competition tests
# ===================================================================


def test_hbond_no_competing_acceptor():
    st = _make_phe_ser_structure()
    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", include_coordinates=True)
    # In this simple structure the only O near H is the donor's own OG
    for hit in hits:
        assert "hbond_competing" in hit
        # Should be 0 since no other acceptor nearby
        assert hit["hbond_competing"] == 0


def test_hbond_with_nearby_carbonyl():
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    phe = gemmi.Residue()
    phe.name = "PHE"
    phe.seqid = _seqid(1)
    for n, p in [("CG", (1.4,0,0)), ("CD1", (0.7,1.212,0)), ("CE1", (-0.7,1.212,0)),
                  ("CZ", (-1.4,0,0)), ("CE2", (-0.7,-1.212,0)), ("CD2", (0.7,-1.212,0))]:
        phe.add_atom(_atom(n, "C", p))
    chain.add_residue(phe)

    # Donor SER with OG and HG
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0, 0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0, 0, 2.0)))
    chain.add_residue(ser)

    # Nearby carbonyl oxygen (competing H-bond acceptor)
    ala = gemmi.Residue()
    ala.name = "ALA"
    ala.seqid = _seqid(3)
    ala.add_atom(_atom("O", "O", (0.0, 0, 1.0)))  # 2.0 Å from HG
    chain.add_residue(ala)

    model.add_chain(chain)
    st.add_model(model)

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", include_coordinates=True)

    assert len(hits) > 0
    hit = hits[0]
    assert hit["hbond_competing"] == 1
    assert hit["hbond_acceptor_res"] == "ALA"
    assert hit["hbond_acceptor_atom"] == "O"
    assert hit["hbond_HA_dist"] is not None
    assert hit["hbond_HA_dist"] <= 3.0  # within search radius
    assert hit["hbond_DHA_angle"] >= 120.0
    assert hit["hbond_vs_xhpi_score"] is not None


def test_hbond_ignores_same_residue_acceptor():
    """Atoms in the same residue as the donor should not count as competing."""
    st = _make_phe_ser_structure()
    # Add a CB atom to SER close to HG — should be ignored
    ser_res = st[0]["A"][1]
    ser_res.add_atom(_atom("CB", "C", (0.5, 0, 2.0)))

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", include_coordinates=True)

    assert len(hits) > 0
    assert hits[0]["hbond_competing"] == 0


def _old_test_hbond_competition_score_positive_for_close_linear_hbond():
    """When H is close to a linear acceptor, score should be positive (H-bond preferred)."""
    score = hbond._competition_score(d_ha=1.9, angle_dha=178.0, dist_x_pi=4.5, angle_xh_pi=120.0)
    assert score > 0, f"Expected positive score, got {score}"


def _old_test_hbond_competition_score_negative_for_xhpi():
    """When H is far from acceptor and XH-pi geometry is good, score should be negative."""
    score = hbond._competition_score(d_ha=2.8, angle_dha=121.0, dist_x_pi=3.0, angle_xh_pi=170.0)
    assert score < 0, f"Expected negative score, got {score}"


def test_hbond_cli_flag_integration():
    st = _make_phe_ser_structure()
    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", include_coordinates=True)
    for hit in hits:
        assert "hbond_competing" in hit
        assert "hbond_vs_xhpi_score" in hit
