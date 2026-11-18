"""Default His policy: neutral, with only high-confidence His+ cationic."""
import pytest

from xpid import his_default, protonation_pka
from test_scientific_fixes import structure
from test_protonation import his, heavy_snapshot

KEY = (0, "A", "12A", "HIS", "")


def predicted(monkeypatch, pka, coupled=False):
    monkeypatch.setattr(protonation_pka, "predict_model", lambda model, mi, messages=None:
                        ({KEY: {"pka": pka, "coupled": coupled}}, "test"))


@pytest.mark.parametrize("ph,pka,coupled,charge", [
    (5.5, 7.0, False, 1),     # exactly at the margin
    (5.6, 7.0, False, 0),     # near the pKa: stays neutral, unlike --ph
    (5.5, 7.0, True, 0),      # coupled sites are not called
])
def test_his_plus_needs_the_margin(monkeypatch, ph, pka, coupled, charge):
    predicted(monkeypatch, pka, coupled)
    s = structure(his())
    context = his_default.prepare(s, "unused", "default", ph)
    assert context.sites[0].charge == charge
    assert context.pi_allowed(0, "A", s[0][0][0], "") is (charge == 0)


def test_his_plus_has_no_n_h_donor(monkeypatch):
    predicted(monkeypatch, 8.0)
    s = structure(his())
    context = his_default.prepare(s, "unused", "default", 5)
    r = s[0][0][0]
    assert context.hydrogen_names(0, "A", r, r.find_atom("ND1", "\0")) == set()


def test_undecided_his_stays_pi_without_donor():
    s = structure(his())
    before = heavy_snapshot(s)
    context = his_default.prepare(s, "unused", "neutral")
    r = s[0][0][0]
    assert context.sites[0].state == "neutral_unresolved"
    assert context.pi_allowed(0, "A", r, "")
    assert context.hydrogen_names(0, "A", r, r.find_atom("NE2", "\0")) == set()
    assert not any(a.is_hydrogen() for a in r)
    assert before == heavy_snapshot(s)


def test_missing_ph_is_recorded_and_his_is_neutral(tmp_path):
    path = tmp_path / "x.cif"
    path.write_text("data_x\n_exptl.method 'X-RAY DIFFRACTION'\n")
    s = structure(his())
    context = his_default.prepare(s, path, "default", "auto")
    assert context.ph is None and context.sites[0].charge == 0
    assert any("No pH" in w for w in context.warnings)


def test_backend_failure_leaves_his_neutral(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("boom")
    monkeypatch.setattr(protonation_pka, "predict_model", fail)
    context = his_default.prepare(structure(his()), "unused", "default", 5)
    assert context.sites[0].charge == 0
    assert any("PROPKA failed" in w for w in context.warnings)


@pytest.mark.parametrize("state,h_name", [("delta", "HD1"), ("epsilon", "HE2")])
def test_override_sets_tautomer_with_occupied_h(state, h_name):
    s = structure(his(occupancy=.55))
    his_default.prepare(s, "unused", "neutral", overrides={KEY: state})
    hydrogens = [a for a in s[0][0][0] if a.is_hydrogen()]
    assert [a.name for a in hydrogens] == [h_name] and hydrogens[0].occ == pytest.approx(.55)


def test_override_his_plus_without_ph():
    s = structure(his())
    context = his_default.prepare(s, "unused", "neutral", overrides={KEY: "doubly_protonated"})
    assert not context.pi_allowed(0, "A", s[0][0][0], "")


def test_unmatched_override_rejected():
    with pytest.raises(ValueError, match="match no site"):
        his_default.prepare(structure(his()), "unused", "neutral", overrides={(0, "B", "1", "HIS", ""): "delta"})


def test_rows_keep_their_schema():
    assert his_default.HisContext(None, "", "his_neutral").hit_fields() == {}


def test_invalid_policy_rejected():
    with pytest.raises(ValueError, match="policy"):
        his_default.validate_options("hip")


# --- 2.1.14: legacy stays the default; overrides change only the named sites ----------------------------

def test_api_default_policy_is_legacy():
    import inspect
    from xpid import api
    assert inspect.signature(api.detect).parameters["his_policy"].default == "legacy"


def test_legacy_without_override_leaves_his_alone():
    s = structure(his())
    before = [(a.name, a.occ, a.pos.tolist()) for a in s[0][0][0]]
    context = his_default.prepare(s, "unused", "legacy")
    assert context.sites == []
    assert [(a.name, a.occ, a.pos.tolist()) for a in s[0][0][0]] == before


@pytest.mark.parametrize("state,h_name", [("delta", "HD1"), ("epsilon", "HE2")])
def test_legacy_override_sets_the_named_site(state, h_name):
    s = structure(his())
    context = his_default.prepare(s, "unused", "legacy", overrides={KEY: state})
    assert [a.name for a in s[0][0][0] if a.is_hydrogen()] == [h_name]
    assert context.pi_allowed(0, "A", s[0][0][0], "")


def test_legacy_override_his_plus():
    s = structure(his())
    context = his_default.prepare(s, "unused", "legacy", overrides={KEY: "doubly_protonated"})
    assert not context.pi_allowed(0, "A", s[0][0][0], "")


@pytest.mark.parametrize("kwargs,message", [
    ({"his_overrides": {KEY: "delta"}, "ph": 7.0}, "cannot be combined with ph"),
    ({"his_overrides": {KEY: "delta"}, "h_mode": 0}, "h_mode"),
])
def test_api_rejects_bad_combinations(kwargs, message):
    from xpid import api
    with pytest.raises(ValueError, match=message):
        api.detect("unused.cif", **kwargs)
