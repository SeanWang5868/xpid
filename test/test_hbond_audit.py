"""Audit semantics and strict invariance of the production detector."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

from xpid import api, cli, cone, core, hbond_audit as audit, rotatable_groups
from test_scientific_fixes import atom, residue, structure, phe, ser


@pytest.mark.parametrize("pi,hb,same,other,disjoint", [
    ([[1], [0]], [[1], [0]], 1, 0, False),
    ([[1], [0]], [[0], [1]], 0, 0, True),
    ([[1, 0, 0]], [[0, 1, 0]], 0, 1, False),
    ([[1, 0, 0]], [[1, 1, 0]], 1, 1, False),
    ([[0], [0]], [[1], [0]], 0, 0, False),
    ([[1], [0]], [[0], [0]], 0, 0, False),
])
def test_relations_are_independent(pi, hb, same, other, disjoint):
    result = audit.relations(pi, hb)
    assert result["same_hydrogen_conformers"] == same
    assert result["same_conformer_other_hydrogen"] == other
    assert result["disjoint_only"] == disjoint


def test_contact_mask_ignores_clashing_conformers_and_wrong_direction():
    positions = np.array([[[1., 0, 0]], [[1., 0, 0]], [[-1., 0, 0]]])
    mask, distance, angle = audit.contact_mask([0, 0, 0], positions, [3, 0, 0], np.array([1, 0, 1], bool))
    assert mask.tolist() == [[True], [False], [False]]
    assert distance[0, 0] == 2 and angle[0, 0] == 180


@pytest.mark.parametrize("name,atom_name,element,expected", [
    ("ALA", "O", "O", "yes"), ("ALA", "N", "N", "no"),
    ("ALA", "OXT", "O", "unknown"), ("TRP", "NE1", "N", "no"),
    ("HIS", "ND1", "N", "unknown"), ("ASP", "OD1", "O", "unknown"),
    ("LIG", "O", "O", "unknown"), ("LIG", "S1", "S", "unknown"),
    ("HOH", "O", "O", "yes"), ("SER", "OG", "O", "yes"),
    ("ASN", "OD1", "O", "yes"), ("ASN", "ND2", "N", "no"),
])
def test_audit_acceptor_chemistry(name, atom_name, element, expected):
    r = residue(name, 1, [atom(atom_name, element, [0, 0, 0])])
    assert audit.classify_acceptor(r, r[0])[0] == expected


def test_his_zero_occupancy_is_not_a_deprotonation_assignment():
    r = residue("HIS", 1, [atom("ND1", "N", [0, 0, 0]), atom("HD1", "H", [1, 0, 0], occ=0)])
    assert audit.classify_acceptor(r, r[0])[0] == "unknown"
    r[1].occ = 1
    assert audit.classify_acceptor(r, r[0])[0] == "no"


@pytest.mark.parametrize("name,atom_name", [
    ("ALA", "N"), ("SER", "N"), ("PRO", "N"),
    ("ASN", "ND2"), ("GLN", "NE2"), ("TRP", "NE1"),
])
def test_standard_nonacceptor_nitrogen_stays_excluded_with_covalent_links(name, atom_name):
    r = residue(name, 1, [atom(atom_name, "N", [0, 0, 0])])
    assert audit.classify_acceptor(r, r[0], special=True)[0] == "no"


def test_neutral_lysine_microstate_is_not_overridden_by_amide_exclusion():
    r = residue("LYS", 1, [atom("NZ", "N", [0, 0, 0])])
    site = SimpleNamespace(status="resolved", source="propka", acceptors={"NZ": True}, reason="")
    assert audit.classify_acceptor(r, r[0], site=site)[0] == "yes"


def test_ph_source_and_uncertainty_are_preserved():
    r = residue("HIS", 1, [atom("ND1", "N", [0, 0, 0])])
    site = SimpleNamespace(status="resolved", source="propka", acceptors={"ND1": True}, reason="")
    assert audit.classify_acceptor(r, r[0], site=site) == ("yes", "assigned_microstate", "propka")
    site.status = "unresolved"
    assert audit.classify_acceptor(r, r[0], site=site)[0] == "unknown"
    site.status = "resolved"
    assert audit.classify_acceptor(r, r[0], site=site, special=True)[0] == "unknown"


def fixture(water_alt="\0", parent_alt="\0", water_occ=1, water_b=10):
    definition = rotatable_groups.get_rotatable_group("SER", "OG")
    x, parent = np.array([10, 10, 13]), np.array([11.42, 10, 13])
    confs = cone.generate_conformers(parent, x, definition)
    h = max((c.hydrogen_positions[0] for c in confs), key=lambda p: p[2])
    oxygen = h + 1.8 * (h - x) / np.linalg.norm(h - x)
    return structure(phe(), ser(parent_alt=parent_alt),
                     residue("HOH", 3, [atom("O", "O", oxygen, altloc=water_alt, occ=water_occ, b=water_b)]))


@pytest.mark.parametrize("slab,sym", [(False, False), (True, False), (False, True)])
def test_complete_audit_does_not_change_hits_or_model(slab, sym):
    s, diagnostics = fixture(), {}
    before = s.make_mmcif_document().as_string()
    baseline = core.detect_interactions_in_structure(s, "test", include_p_slab=slab, sym_contacts=sym)
    observed = core.detect_interactions_in_structure(s, "test", hbond_audit=True,
                                                    diagnostics=diagnostics, include_p_slab=slab, sym_contacts=sym)
    assert baseline == observed and baseline
    assert before == s.make_mmcif_document().as_string()
    report = diagnostics["hbond_audit"]
    assert report["pairs_seen"] == len(report["pairs"]) == 1
    pair = report["pairs"][0]
    assert pair["typed_contacts"]["disjoint_only"]
    assert pair["uncertain_contacts"]["hbond_conformers"] == 0
    assert pair["acceptors"][0]["residue_name"] == "HOH"
    assert pair["acceptors"][0]["symmetry_code"] == "1_555"
    mask = np.unpackbits(np.frombuffer(bytes.fromhex(pair["pi_mask"]), dtype=np.uint8))
    assert int(mask.sum()) >= pair["pi_positive_conformers"] > 0
    json.dumps(report, allow_nan=False)
    core.detect_interactions_in_structure(s, "test", diagnostics=diagnostics)
    assert "hbond_audit" not in diagnostics


@pytest.mark.parametrize("alt,occ,b,expected", [("B", 1, 10, False), ("A", 0, 10, False),
                                               ("A", 1, 80, False), ("A", 1, 10, True)])
def test_filters_and_cross_residue_altloc_uncertainty(alt, occ, b, expected):
    s, d = fixture(alt, "A", occ, b), {}
    core.detect_interactions_in_structure(s, "test", hbond_audit=True, diagnostics=d, max_b=30)
    pair = d["hbond_audit"]["pairs"][0]
    assert bool(pair["acceptors"]) == expected
    if expected:
        assert pair["acceptors"][0]["status"] == "unknown"
        assert pair["typed_contacts"]["hbond_conformers"] == 0


def test_audit_truncation_never_truncates_detection(monkeypatch):
    monkeypatch.setattr(audit, "MAX_PAIRS", 0)
    s, d = fixture(), {}
    baseline = core.detect_interactions_in_structure(s, "test")
    assert baseline == core.detect_interactions_in_structure(s, "test", diagnostics=d, hbond_audit=True)
    assert d["hbond_audit"]["pairs_omitted"] == 1
    assert d["hbond_audit"]["pairs"] == []


def test_negative_pair_and_unknown_contacts_are_retained():
    c = cone.HydrogenConformer(0, (np.array([1., 0, 0]),))
    raw = {"conformers": [c], "valid_mask": np.array([True]), "pi_mask": np.array([[False]])}
    r = audit.summarize(raw, [0, 0, 0], [{"position": [3, 0, 0], "status": "unknown"}], {}, "C")
    assert r["pi_positive_conformers"] == 0
    assert r["uncertain_contacts"]["hbond_conformers"] == 1
    assert r["typed_contacts"]["hbond_conformers"] == 0
    assert r["donor_class"] == "weak_CH_contact_only"


@pytest.mark.parametrize("kwargs", [{}, {"diagnostics": {}, "cone_mode": "none"},
                                   {"diagnostics": {}, "report_xh_candidates": True}])
def test_api_requires_a_destination_and_supported_mode(kwargs):
    with pytest.raises(ValueError, match="hbond_audit"):
        api.detect("unused.cif", hbond_audit=True, **kwargs)


def test_cli_has_opt_in_switch():
    assert not cli._build_parser().parse_args([]).hbond_audit
    assert cli._build_parser().parse_args(["--hbond-audit"]).hbond_audit


def test_audit_resolves_periodic_image_around_the_donor():
    s, d = fixture(), {}
    s[0][0][2][0].pos.x -= 30
    baseline = core.detect_interactions_in_structure(s, "test")
    assert baseline == core.detect_interactions_in_structure(s, "test", hbond_audit=True, diagnostics=d)
    assert d["hbond_audit"]["pairs"][0]["acceptors"][0]["symmetry_code"] == "1_655"


def test_coordinated_water_is_not_a_definite_competitor():
    s, d = fixture(), {}
    water = s[0][0][2][0].pos.tolist()
    s[0][0].add_residue(residue("ZN", 4, [atom("ZN", "Zn", np.array(water) + [2.4, 0, 0])]))
    core.detect_interactions_in_structure(s, "test", hbond_audit=True, diagnostics=d)
    water_records = [a for p in d["hbond_audit"]["pairs"] for a in p["acceptors"] if a["residue_name"] == "HOH"]
    assert water_records and all(a["status"] == "unknown" for a in water_records)


def test_peptide_link_does_not_introduce_uncertain_backbone_n_competition():
    import gemmi
    s, d = fixture(), {}
    nitrogen_residue = s[0][0][2]
    nitrogen_residue.name = "SER"
    nitrogen_residue[0].name = "N"
    nitrogen_residue[0].element = gemmi.Element("N")
    position = np.array(nitrogen_residue[0].pos.tolist())
    s[0][0].add_residue(residue("ALA", 4, [atom("C", "C", position + [0, 0, 1.33])]))
    baseline = core.detect_interactions_in_structure(s, "test")
    assert baseline == core.detect_interactions_in_structure(s, "test", diagnostics=d, hbond_audit=True)
    assert d["hbond_audit"]["pairs"]
    assert not any(a["atom"] == "N" for p in d["hbond_audit"]["pairs"] for a in p["acceptors"])


def test_sh_search_covers_hydrogen_length_plus_contact_radius():
    x, parent = np.array([10., 10, 13]), np.array([11.812, 10, 13])
    definition = rotatable_groups.get_rotatable_group("CYS", "SG")
    confs = cone.generate_conformers(parent, x, definition)
    h = max((c.hydrogen_positions[0] for c in confs), key=lambda p: p[2])
    oxygen = h + 2.78 * (h - x) / np.linalg.norm(h - x)
    assert np.linalg.norm(oxygen - x) > 4
    donor = residue("CYS", 2, [atom("SG", "S", x), atom("CB", "C", parent), atom("HG", "H", h)])
    s, d = structure(phe(), donor, residue("HOH", 3, [atom("O", "O", oxygen)])), {}
    baseline = core.detect_interactions_in_structure(s, "test")
    assert baseline == core.detect_interactions_in_structure(s, "test", hbond_audit=True, diagnostics=d)
    assert any(a["residue_name"] == "HOH" for p in d["hbond_audit"]["pairs"] for a in p["acceptors"])


def test_all_models_receive_separate_audit_identities():
    s, d = fixture(), {}
    second = s[0].clone()
    second.num = 2
    s.add_model(second)
    baseline = core.detect_interactions_in_structure(s, "test", model_mode="all")
    assert baseline == core.detect_interactions_in_structure(s, "test", model_mode="all", hbond_audit=True, diagnostics=d)
    assert {p["model_index"] for p in d["hbond_audit"]["pairs"]} == {0, 1}
