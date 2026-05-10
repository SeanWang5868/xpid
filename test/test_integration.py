import gemmi
from xpid import core, config


def _atom(name, element, xyz, b_iso=10.0):
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element(element)
    atom.pos = gemmi.Position(*xyz)
    atom.occ = 1.0
    atom.b_iso = b_iso
    return atom


def _seqid(num):
    return gemmi.SeqId(str(num))


def _structure_with_phe_and_og(ring_b=10.0, include_external_donor=True):
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    phe = gemmi.Residue()
    phe.name = "PHE"
    phe.seqid = _seqid(1)
    ring = {
        "CG": (1.4, 0.0, 0.0),
        "CD1": (0.7, 1.212, 0.0),
        "CE1": (-0.7, 1.212, 0.0),
        "CZ": (-1.4, 0.0, 0.0),
        "CE2": (-0.7, -1.212, 0.0),
        "CD2": (0.7, -1.212, 0.0),
    }
    for name, xyz in ring.items():
        phe.add_atom(_atom(name, "C", xyz, b_iso=ring_b))
    chain.add_residue(phe)

    if include_external_donor:
        ser = gemmi.Residue()
        ser.name = "SER"
        ser.seqid = _seqid(2)
        ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
        ser.add_atom(_atom("HG", "H", (0.0, 0.0, 2.0)))
        chain.add_residue(ser)

    model.add_chain(chain)
    st.add_model(model)
    return st

def test_core_detection_empty():
    st = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    st.add_model(model)
    model.add_chain(chain)
    
    # Empty structure should return empty list
    res = core.detect_interactions_in_structure(st, "test", {}, model_mode=0)
    assert len(res) == 0

def test_model_selection():
    st = gemmi.Structure()
    m1 = gemmi.Model("1")
    m2 = gemmi.Model("2")
    st.add_model(m1)
    st.add_model(m2)
    
    # Dummy atoms to prevent immediate return
    c1 = gemmi.Chain("A")
    m1.add_chain(c1)
    
    # Test valid index
    core.detect_interactions_in_structure(st, "test", {}, model_mode=0)
    # Test 'all'
    core.detect_interactions_in_structure(st, "test", {}, model_mode='all')


def test_donor_atom_filter_accepts_element_symbols():
    st = _structure_with_phe_and_og()

    oxygen_hits = core.detect_interactions_in_structure(
        st, "test", filter_donor_atom=["O"])
    assert any(hit["X_atom"] == "OG" for hit in oxygen_hits)

    atom_name_hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", filter_donor_atom=["OG"])
    assert any(hit["X_atom"] == "OG" for hit in atom_name_hits)

    nitrogen_hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", filter_donor_atom=["N"])
    assert nitrogen_hits == []


def test_hydrogen_atoms_are_not_x_donors():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    asn = gemmi.Residue()
    asn.name = "ASN"
    asn.seqid = _seqid(2)
    asn.add_atom(_atom("HD21", "H", (0.0, 0.0, 3.0)))
    asn.add_atom(_atom("H1", "H", (0.0, 0.0, 2.0)))
    chain.add_residue(asn)

    assert core.detect_interactions_in_structure(st, "test") == []


def test_dictionary_bonded_hydrogens_reject_unbonded_neighbor(monkeypatch):
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    bma = gemmi.Residue()
    bma.name = "BMA"
    bma.seqid = _seqid(2)
    bma.add_atom(_atom("O6", "O", (0.0, 0.0, 3.0)))
    bma.add_atom(_atom("H1", "H", (0.0, 0.0, 2.0)))
    chain.add_residue(bma)

    monkeypatch.setattr(
        config,
        "get_bonded_hydrogens",
        lambda res, atom: {"HO6"} if (res, atom) == ("BMA", "O6") else set(),
    )

    assert core.detect_interactions_in_structure(st, "test") == []


def test_max_b_filters_pi_ring_atoms():
    assert core.detect_interactions_in_structure(
        _structure_with_phe_and_og(ring_b=80.0), "test", max_b=50.0) == []


def test_same_residue_donor_is_excluded():
    st = _structure_with_phe_and_og(include_external_donor=False)
    phe = st[0][0][0]
    phe.add_atom(_atom("N", "N", (0.0, 0.0, 3.0)))
    phe.add_atom(_atom("H", "H", (0.0, 0.0, 2.0)))

    assert core.detect_interactions_in_structure(st, "test") == []
