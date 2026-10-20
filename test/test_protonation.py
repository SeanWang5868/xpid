"""Chemical-state regressions: pH never silently creates an occupied microstate."""
import gzip
import json
from pathlib import Path

import gemmi
import numpy as np
import pytest

from xpid import api, cli, core, hydrogen_prep, protonation, protonation_geometry, protonation_pka
from xpid.output import ResultStreamer
from test_scientific_fixes import atom, residue, structure, phe, thiol


def his(label="", occupancy=1.0):
    # Ideal pentagon, with the two N-H directions pointing away from its centre.
    atoms = []
    for index, name in enumerate(["CG", "ND1", "CE1", "NE2", "CD2"]):
        theta = 2 * np.pi * index / 5
        direction = np.array([np.cos(theta), np.sin(theta), 0.0])
        xyz = np.array([10., 10., 10.]) + 1.15 * direction
        atoms.append(atom(name, "N" if name.startswith("N") else "C", xyz,
                          altloc=label or "\0", occ=occupancy))
        if name in {"ND1", "NE2"}:
            atoms.append(atom("HD1" if name == "ND1" else "HE2", "H", xyz + direction,
                              altloc=label or "\0", occ=0))
    r = residue("HIS", 12, atoms)
    r.seqid = gemmi.SeqId(12, "A")
    return r


def prepared_context(s, ph=8, **kwargs):
    return protonation.prepare(s, "unused.cif", ph, method="reference", **kwargs)


def heavy_snapshot(s):
    return [(ci, str(r.seqid), r.name, a.name, a.altloc, a.pos.tolist(), a.occ, a.b_iso, a.charge)
            for m in s for ci, c in enumerate(m) for r in c for a in r if not a.is_hydrogen()]


@pytest.mark.parametrize("ph", ["nan", "inf", "6-7", "7.0(1)", "unknown", -1, 15])
def test_invalid_ph_rejected(ph):
    with pytest.raises(ValueError):
        protonation_pka.numeric_ph(ph)


@pytest.mark.parametrize("mode", [0, 1, 2])
def test_ph_requires_rebuild(mode):
    with pytest.raises(ValueError, match="rebuilding"):
        protonation_pka.validate_options(7, "reference", mode)


@pytest.mark.parametrize("method,tag", [
    ("X-RAY DIFFRACTION", "_exptl_crystal_grow.pH"),
    ("ELECTRON MICROSCOPY", "_em_buffer.pH"),
    ("SOLUTION NMR", "_pdbx_nmr_exptl_sample_conditions.pH")])
def test_auto_ph_uses_experimental_method(tmp_path, method, tag):
    path = tmp_path / "input.cif.gz"
    with gzip.open(path, "wt") as f:
        f.write(f"data_x\n_exptl.method '{method}'\n{tag} 6.5\n")
    assert protonation_pka.resolve_ph(path, "auto") == (6.5, tag)


def test_auto_ph_does_not_pick_irrelevant_or_conflicting_metadata(tmp_path):
    path = tmp_path / "input.cif"
    path.write_text("data_x\n_exptl.method 'ELECTRON MICROSCOPY'\n_exptl_crystal_grow.pH 6\n")
    with pytest.raises(ValueError, match="unique"):
        protonation_pka.resolve_ph(path, "auto")
    path.write_text("data_x\n_exptl.method 'X-RAY DIFFRACTION'\nloop_\n_exptl_crystal_grow.pH\n6\n7\n")
    with pytest.raises(ValueError, match="unique"):
        protonation_pka.resolve_ph(path, "auto")
    assert protonation_pka.resolve_ph(path, 8) == (8, "user")


def test_pdb_ph(tmp_path):
    path = tmp_path / "x.pdb"
    path.write_text("REMARK 200  PH                             : 7.5\n")
    assert protonation_pka.resolve_ph(path, "auto")[0] == 7.5


def test_neutral_his_without_evidence_is_not_double_protonated():
    s = structure(his())
    before = heavy_snapshot(s)
    context = prepared_context(s)
    site = context.sites[0]
    assert site.charge == 0 and site.status == "unresolved"
    assert site.reason == "ambiguous_or_clashing_microstate"
    assert not any(a.is_hydrogen() for a in s[0][0][0])
    assert context.pi_allowed(0, "A", s[0][0][0], "")
    assert before == heavy_snapshot(s)


@pytest.mark.parametrize("state,h_name,acceptor", [("delta", "HD1", "NE2"), ("epsilon", "HE2", "ND1")])
def test_explicit_his_microstate_is_not_zero_occupancy_or_renamed(state, h_name, acceptor):
    s = structure(his(occupancy=.55))
    key = (0, "A", "12A", "HIS", "")
    context = prepared_context(s, overrides={key: state})
    r = s[0][0][0]
    hydrogens = [a for a in r if a.is_hydrogen()]
    assert r.name == "HIS" and len(hydrogens) == 1
    assert hydrogens[0].name == h_name and hydrogens[0].flag == "P"
    assert hydrogens[0].occ == pytest.approx(.55)
    assert context.sites[0].protonated_fraction != hydrogens[0].occ
    assert context.is_acceptor(0, "A", r, r.find_atom(acceptor, "\0"))
    assert context.sites[0].source == "user_override"


def test_low_ph_his_contains_two_h_but_is_not_neutral_pi_or_n_donor():
    s = structure(his())
    c = prepared_context(s, ph=4)
    r = s[0][0][0]
    assert c.sites[0].state == "doubly_protonated"
    assert len([a for a in r if a.is_hydrogen()]) == 2
    assert not c.pi_allowed(0, "A", r, "")
    for name in ("ND1", "NE2"):
        assert c.hydrogen_names(0, "A", r, r.find_atom(name, "\0")) == set()
        assert not c.is_acceptor(0, "A", r, r.find_atom(name, "\0"))


def test_his_near_pka_is_explicitly_uncertain():
    s = structure(his())
    c = prepared_context(s, 6.5)
    assert c.sites[0].reason == "mixed_charge_states_near_pKa"
    assert c.sites[0].protonated_fraction == .5
    assert not c.pi_allowed(0, "A", s[0][0][0], "")


def test_his_tautomer_selected_by_local_hbond_not_pi():
    r = his()
    n, h = r.find_atom("ND1", "\0"), r.find_atom("HD1", "\0")
    direction = np.array(h.pos.tolist()) - np.array(n.pos.tolist())
    water = residue("HOH", 20, [atom("O", "O", np.array(h.pos.tolist()) + 1.9 * direction)])
    s = structure(r, water)
    c = prepared_context(s)
    assert c.sites[0].state == "delta"
    assert c.sites[0].candidates[0]["score"] > 0


def test_cys_deprotonation_prevents_cone_resurrection():
    r = thiol()
    r.add_atom(atom("HG", "H", (10, 10, 12.462)))
    s = structure(phe(), r)
    assert core.detect_interactions_in_structure(s, "x")
    c = prepared_context(s, 12)
    assert c.sites[0].state == "deprotonated"
    assert core.detect_interactions_in_structure(s, "x", protonation=c) == []
    assert c.is_acceptor(0, "A", s[0][0][1], s[0][0][1].find_atom("SG", "\0"))


def test_ph_off_keeps_existing_detector_and_ph_fields_are_opt_in():
    r = thiol()
    r.add_atom(atom("HG", "H", (10, 10, 12.462)))
    s = structure(phe(), r)
    baseline = core.detect_interactions_in_structure(s, "x")
    assert "pH" not in baseline[0]
    c = prepared_context(s, 6)
    hits = core.detect_interactions_in_structure(s, "x", protonation=c)
    assert len(hits) == len(baseline)
    for field in ("is_hudson", "is_plevin", "H_xyz_x", "H_xyz_y", "H_xyz_z"):
        assert hits[0][field] == baseline[0][field]
    assert hits[0]["X_protonation_state"] == "protonated"


def test_altlocs_and_insertion_codes_are_separate_states():
    r = his("A")
    for a in his("B"):
        r.add_atom(a)
    s = structure(r)
    overrides = {(0, "A", "12A", "HIS", "A"): "delta", (0, "A", "12A", "HIS", "B"): "epsilon"}
    c = prepared_context(s, overrides=overrides)
    assert [(a.name, a.altloc) for a in s[0][0][0] if a.is_hydrogen()] == [("HD1", "A"), ("HE2", "B")]
    assert len(c.sites) == 2


def test_metal_site_cannot_be_force_protonated():
    r = his()
    n = r.find_atom("ND1", "\0")
    metal = residue("ZN", 40, [atom("ZN", "Zn", np.array(n.pos.tolist()) + [0, 0, 2.1])])
    s = structure(r, metal)
    c = prepared_context(s, overrides={(0, "A", "12A", "HIS", ""): "delta"})
    assert "metal_coordination" in c.sites[0].reason
    assert c.warnings and not any(a.is_hydrogen() for a in s[0][0][0])


def test_backend_failure_is_audited_without_reference_fallback(monkeypatch):
    def fail(*args):
        raise RuntimeError("bad backend")
    monkeypatch.setattr(protonation_pka, "predict_model", fail)
    c = protonation.prepare(structure(his()), "unused", 8, method="propka")
    assert c.sites[0].pka is None and c.sites[0].status == "unresolved"
    assert "no reference-pKa fallback" in c.warnings[0]


def test_override_validation_and_unmatched_site(tmp_path):
    path = tmp_path / "states.json"
    path.write_text(json.dumps([dict(model_index=0, chain="A", residue_id="12A", residue_name="HIS", altloc="", state="delta")]))
    assert protonation.read_overrides(path)[(0, "A", "12A", "HIS", "")] == "delta"
    path.write_text(path.read_text().replace('"delta"', '"HIP"'))
    with pytest.raises(ValueError, match="Unsupported"):
        protonation.read_overrides(path)
    with pytest.raises(ValueError, match="match no site"):
        prepared_context(structure(his()), overrides={(1, "A", "12A", "HIS", ""): "delta"})


@pytest.mark.parametrize("file_type", ["csv", "json", "parquet"])
def test_ph_output_schema_stays_stable_across_null_and_float(tmp_path, file_type):
    if file_type == "parquet":
        pytest.importorskip("pyarrow")
    path = tmp_path / ("x." + file_type)
    with ResultStreamer(path, file_type, False, include_protonation=True) as writer:
        writer.write_chunk([{"pH": 7.0, "X_site_pKa": None}])
        writer.write_chunk([{"pH": 7.0, "X_site_pKa": 6.5}])
    assert path.stat().st_size > 0


def test_cli_options_default_off_and_conflicting_api_mode():
    args = cli._build_parser().parse_args([])
    assert args.ph is None
    args = cli._build_parser().parse_args(["--ph", "auto", "--pka-method", "reference"])
    assert args.ph == "auto"
    with pytest.raises(ValueError, match="rebuilding"):
        api.detect("unused", ph=7, pka_method="reference", h_mode=0)


@pytest.mark.slow
def test_acid_geometry_uses_ccp4_nuclear_restraints():
    lib = hydrogen_prep._get_shared_monlib({"ASP", "GLU"})
    if "ASPprot" not in lib.modifications:
        pytest.skip("CCP4 modifications not installed")
    r = residue("ASP", 3, [atom("CG", "C", (10, 10, 10)), atom("CB", "C", (9, 10, 10)),
                            atom("OD1", "O", (10.65, 11.126, 10)), atom("OD2", "O", (10.65, 8.874, 10))])
    candidates = protonation_geometry.acid_candidates(r, "", lib)
    assert len(candidates) == 4
    for _, group in candidates:
        for parent, hs in group.items():
            assert r.find_atom(parent, "\0").pos.dist(hs[0].pos) == pytest.approx(.966)
            assert hs[0].pos.z == 10
    s = structure(r)
    c = prepared_context(s, 1, overrides={(0, "A", "3", "ASP", ""): "protonated_OD2_trans"})
    assert c.sites[0].status == "resolved"
    assert [a.name for a in s[0][0][0] if a.is_hydrogen()] == ["HD2"]


def facing_ring(parent, hydrogen):
    direction = np.array(hydrogen.pos.tolist()) - np.array(parent.pos.tolist())
    direction /= np.linalg.norm(direction)
    centre = np.array(parent.pos.tolist()) + 3.2 * direction
    v = np.cross(direction, [0, 0, 1])
    if np.linalg.norm(v) < .1:
        v = np.cross(direction, [0, 1, 0])
    v /= np.linalg.norm(v)
    w = np.cross(direction, v)
    atoms = []
    for i, name in enumerate(("CG", "CD1", "CE1", "CZ", "CE2", "CD2")):
        t = np.pi * i / 3
        atoms.append(atom(name, "C", centre + 1.4 * (v * np.cos(t) + w * np.sin(t))))
    return residue("PHE", 50, atoms)


def test_assigned_his_is_detected_and_correctly_labelled():
    r = his()
    s = structure(r, facing_ring(r.find_atom("ND1", "\0"), r.find_atom("HD1", "\0")))
    assert not [h for h in core.detect_interactions_in_structure(s, "x") if h["X_atom"] == "ND1"]
    c = prepared_context(s, overrides={(0, "A", "12A", "HIS", ""): "delta"})
    hits = [h for h in core.detect_interactions_in_structure(s, "x", protonation=c) if h["X_atom"] == "ND1"]
    assert hits and hits[0]["H_source"] == "ph_model"
    assert hits[0]["X_protonation_state"] == "delta"


def test_assigned_acid_oh_is_not_rejected_by_default_acid_filter():
    r = residue("ASP", 3, [atom("CG", "C", (10, 10, 10)), atom("CB", "C", (9, 10, 10)),
                            atom("OD1", "O", (10.65, 11.126, 10)), atom("OD2", "O", (10.65, 8.874, 10))])
    lib = hydrogen_prep._get_shared_monlib({"ASP"})
    candidates = protonation_geometry.acid_candidates(r, "", lib)
    if not candidates:
        pytest.skip("CCP4 modifications not installed")
    h = next(hs["OD2"][0] for name, hs in candidates if name == "protonated_OD2_trans")
    s = structure(r, facing_ring(r.find_atom("OD2", "\0"), h))
    c = prepared_context(s, 1, overrides={(0, "A", "3", "ASP", ""): "protonated_OD2_trans"})
    hits = core.detect_interactions_in_structure(s, "x", protonation=c)
    assert any(h["X_atom"] == "OD2" and h["H_atom"] == "HD2" for h in hits)


def test_selected_model_only_and_coupled_site_is_not_single_henderson_hasselbalch(monkeypatch):
    s = structure(his())
    m = s[0].clone()
    m.num = 2
    s.add_model(m)
    called = []
    def predictor(model, mi, messages=None):
        called.append(mi)
        return {(mi, "A", "12A", "HIS", ""): {"pka": 6.5, "coupled": True}}, "test"
    monkeypatch.setattr(protonation_pka, "predict_model", predictor)
    c = protonation.prepare(s, "unused", 8, model_mode="1")
    assert called == [1]
    assert c.sites[0].reason == "coupled_sites_require_microstate_assignment"
    assert c.sites[0].protonated_fraction is None


def test_terminal_amine_is_not_mislabelled_as_neutral_backbone():
    s = structure(residue("ALA", 1, [atom("N", "N", [10, 10, 10]),
                                     atom("H2", "H", [10, 10, 11])]))
    c = prepared_context(s)
    assert c.sites[0].site_type == "N_terminus"
    assert c.hydrogen_names(0, "A", s[0][0][0], s[0][0][0][0]) == set()


def test_api_and_worker_have_identical_ph_results_and_diagnostics(tmp_path, monkeypatch):
    r = thiol()
    r.add_atom(atom("HG", "H", (10, 10, 12.462)))
    s = structure(phe(), r)
    path = tmp_path / "model.cif"
    s.make_mmcif_document().write_file(str(path))
    def prep(st, h_change_val=4):
        return hydrogen_prep.HydrogenPreparationResult(st, hydrogen_prep.HydrogenPreparationReport("success", 4))
    monkeypatch.setattr(hydrogen_prep, "prepare_hydrogens_memory", prep)
    diagnostics = {}
    api_hits = api.detect(path, ph=6, pka_method="reference", diagnostics=diagnostics)
    task = cli.TaskPacket(path, "model", "direct", "json", 4, str(tmp_path), False,
                          {}, False, "0", "auto", False, False, False, False, True,
                          None, False, 0., False, False, 0., "6", "reference")
    result = cli.process_one_file(task)
    assert result.error is None and result.data == api_hits
    assert result.diagnostics["protonation"] == diagnostics["protonation"]
    json.dumps(result.diagnostics)


def test_propka_adapter_maps_long_chains_and_author_ids(monkeypatch):
    run = pytest.importorskip("propka.run")
    from types import SimpleNamespace
    s = structure(his())
    s[0][0].name = "long_chain"
    def single(filename, stream, write_pka):
        import logging
        logging.getLogger("propka").warning("Missing atoms for HIS 1 A")
        assert write_pka is False
        mapped = gemmi.read_pdb_string(stream.read())
        assert mapped[0][0].name == "A" and mapped[0][0][0].seqid.num == 1
        assert not any(a.is_hydrogen() for a in mapped[0][0][0])
        group = SimpleNamespace(titratable=True, atom=SimpleNamespace(res_name="HIS", chain_id="A", res_num=1),
                                type="HIS", pka_value=6.4, covalently_coupled_groups=[],
                                non_covalently_coupled_groups=[], coupled_titrating_group=None)
        return SimpleNamespace(conformation_names=["1A"], conformations={"1A": SimpleNamespace(groups=[group])})
    monkeypatch.setattr(run, "single", single)
    messages = []
    predictions, _ = protonation_pka.predict_model(s[0], 0, messages)
    assert predictions[(0, "long_chain", "12A", "HIS", "")]["pka"] == 6.4
    assert "long_chain" in messages[0] and "residue 12A HIS" in messages[0]


def test_protonated_cys_context_overrides_charge_without_modifying_heavy_atoms():
    r = thiol(charge=-1)
    r.add_atom(atom("HG", "H", (10, 10, 12.462)))
    s = structure(phe(), r)
    before = heavy_snapshot(s)
    c = prepared_context(s, 6)
    assert s[0][0][1].find_atom("SG", "\0").pos.dist(s[0][0][1].find_atom("HG", "\0").pos) == pytest.approx(1.338)
    assert core.detect_interactions_in_structure(s, "x", protonation=c)
    assert before == heavy_snapshot(s)


def test_shared_x_checks_each_complete_conformer_state():
    r = thiol()
    cb = r.find_atom("CB", "\0")
    cb.altloc = "A"
    cb.occ = .5
    second = cb.clone()
    second.altloc = "B"
    second.pos.y += .2
    r.add_atom(second)
    r.add_atom(atom("HG", "H", (10, 10, 12.462), altloc="A"))
    r.add_atom(atom("HG", "H", (10, 10.1, 12.462), altloc="B"))
    s = structure(phe(), r)
    c = prepared_context(s, 6, overrides={(0, "A", "2", "CYS", "B"): "deprotonated"})
    assert s[0][0][1].find_atom("HG", "A").occ == pytest.approx(.5)
    hits = core.detect_interactions_in_structure(s, "x", protonation=c)
    assert hits and all(h["H_altloc"] == "A" for h in hits)
    assert all(h["X_protonation_state"] == "protonated" for h in hits)


def test_shared_his_ring_retains_consensus_neutral_charge():
    r = his()
    r.add_atom(atom("CB", "C", (13, 10, 10), altloc="A", occ=.5))
    r.add_atom(atom("CB", "C", (13, 10, 10.1), altloc="B", occ=.5))
    s = structure(r)
    c = prepared_context(s, 8)
    assert all(site.charge == 0 for site in c.sites)
    assert c.pi_allowed(0, "A", s[0][0][0], "")


def test_shared_h_is_not_copied_to_different_parent_frameworks():
    r = thiol(altloc="A")
    r.add_atom(atom("HG", "H", (10, 10, 12.462)))
    c = prepared_context(structure(r), 6)
    assert c.sites[0].reason == "missing_dictionary_hydrogen_geometry"


@pytest.mark.parametrize("length,angle", [(float("nan"), 110), (0., 110), (.966, float("nan")), (.966, 180)])
def test_invalid_acid_dictionary_geometry_is_not_guessed(length, angle):
    from types import SimpleNamespace as NS
    bond = NS(id1=NS(atom="OD2"), id2=NS(atom="HD2"), value_nucleus=length)
    restraint = NS(id1=NS(atom="CG"), id2=NS(atom="OD2"), id3=NS(atom="HD2"), value=angle)
    library = NS(modifications={"ASPprot": NS(rt=NS(bonds=[bond], angles=[restraint]))})
    assert protonation_geometry.acid_candidates(residue("ASP", 1, []), "", library) == []


def test_backend_logging_restored_on_failure(monkeypatch):
    run = pytest.importorskip("propka.run")
    import logging
    log = logging.getLogger("propka")
    before = log.level, log.propagate, list(log.handlers)
    def fail(*args, **kwargs):
        log.warning("test chemical warning")
        raise RuntimeError("failed")
    monkeypatch.setattr(run, "single", fail)
    messages = []
    with pytest.raises(RuntimeError, match="failed"):
        protonation_pka.predict_model(structure(his())[0], 0, messages)
    assert messages == ["test chemical warning"]
    assert (log.level, log.propagate, list(log.handlers)) == before
