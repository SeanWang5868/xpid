"""Regression cases for provenance, covalent state and crystallographic context."""
import gemmi
import numpy as np
import pytest

from xpid import (
    cone, cooperativity, core, hbond_acceptors, hbond_annotations,
    hydrogen_prep, monomer_bonds, rotatable_groups,
)


def atom(name, element, xyz, altloc="\0", occ=1.0, charge=0, b=10.0):
    a = gemmi.Atom()
    a.name, a.element, a.pos = name, gemmi.Element(element), gemmi.Position(*xyz)
    a.altloc, a.occ, a.charge, a.b_iso = altloc, occ, charge, b
    return a


def residue(name, seqid, atoms):
    r = gemmi.Residue()
    r.name, r.seqid = name, gemmi.SeqId(str(seqid))
    for a in atoms:
        r.add_atom(a)
    return r


def structure(*residues, spacegroup="P 1"):
    s = gemmi.Structure()
    s.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    s.spacegroup_hm = spacegroup
    m, c = gemmi.Model("1"), gemmi.Chain("A")
    for r in residues:
        c.add_residue(r)
    m.add_chain(c)
    s.add_model(m)
    s.setup_cell_images()
    return s


def phe():
    coordinates = [("CG", (1.4, 0, 0)), ("CD1", (.7, 1.212, 0)),
                   ("CE1", (-.7, 1.212, 0)), ("CZ", (-1.4, 0, 0)),
                   ("CE2", (-.7, -1.212, 0)), ("CD2", (.7, -1.212, 0))]
    return residue("PHE", 1, [atom(n, "C", np.array(p) + 10) for n, p in coordinates])


def ser(h_occ=1.0, parent_occ=1.0, x_alt="\0", h_alt="\0", parent_alt="\0"):
    return residue("SER", 2, [atom("OG", "O", (10, 10, 13), altloc=x_alt),
                               atom("HG", "H", (10, 10, 12), altloc=h_alt, occ=h_occ),
                               atom("CB", "C", (11.42, 10, 13), altloc=parent_alt,
                                    occ=parent_occ)])


@pytest.mark.parametrize("mode,expected", [(0, "input_model"), (1, "adjusted"),
                                          (3, "added"), (4, "added"), (5, "added")])
def test_hydrogen_preparation_records_actual_operation(monkeypatch, mode, expected):
    s = structure(phe(), ser())
    lib = gemmi.MonLib()
    for name in ("SER", "PHE"):
        comp = gemmi.ChemComp()
        comp.name = name
        lib.monomers[name] = comp
    monkeypatch.setattr(hydrogen_prep, "_get_shared_monlib", lambda names: lib)

    def prepare(working, monlib, h_change, model_index):
        r = working[model_index][0][1]
        if mode == 1:
            r.find_atom("HG", "\0").pos.z -= .01
        else:
            r.remove_atom("HG", "\0", gemmi.Element("H"))
            # Gemmi's default new-atom flag is deliberately left unchanged.
            r.add_atom(atom("HG", "H", (10, 10, 12)))
        return hydrogen_prep.TopologyPreparationReport(model_index, "success", 1)

    monkeypatch.setattr(hydrogen_prep, "_prepare_topology_with_retries", prepare)
    hydrogen_prep.prepare_hydrogens_memory(s, mode)
    hits = core.detect_interactions_in_structure(s, "test", cone_mode="none")
    assert len(hits) == 1
    assert hits[0]["H_source"] == expected
    assert hits[0]["H_source"] != "experimental"
    # No-change reuse must not relabel previously prepared coordinates as input.
    hydrogen_prep.prepare_hydrogens_memory(s, 0)
    assert core.detect_interactions_in_structure(s, "test", cone_mode="none")[0]["H_source"] == expected


def test_mode_four_retained_water_and_unknown_h_are_input_model(monkeypatch):
    water = residue("HOH", 3, [atom("O", "O", (10, 10, 13)),
                                atom("H1", "H", (10, 10, 12))])
    unknown = residue("XXX", 4, [atom("H1", "H", (12, 12, 12))])
    s = structure(phe(), water, unknown)
    lib = gemmi.MonLib()
    comp = gemmi.ChemComp()
    comp.name = "HOH"
    lib.monomers["HOH"] = comp
    monkeypatch.setattr(hydrogen_prep, "_get_shared_monlib", lambda names: lib)
    monkeypatch.setattr(hydrogen_prep, "_prepare_topology_with_retries",
                        lambda *args: hydrogen_prep.TopologyPreparationReport(0, "success", 1))
    hydrogen_prep.prepare_hydrogens_memory(s, 4)
    assert s[0][0][1].find_atom("H1", "\0").flag == "I"
    assert s[0][0][2].find_atom("H1", "\0").flag == "I"


def thiol(charge=0, altloc="\0"):
    return residue("CYS", 2, [atom("SG", "S", (10, 10, 13.8), charge=charge, altloc=altloc),
                               atom("CB", "C", (11.812, 10, 13.8), altloc=altloc)])


def crosslink(s, altloc="\0", kind=gemmi.ConnectionType.Covale,
              target_res="LIG", target_atom="C1"):
    link = gemmi.Connection()
    link.name, link.type = "modified_SG", kind
    link.partner1 = gemmi.AtomAddress("A", gemmi.SeqId("2"), "CYS", "SG", altloc)
    link.partner2 = gemmi.AtomAddress("A", gemmi.SeqId("3"), target_res, target_atom, altloc)
    s.connections.append(link)


@pytest.mark.parametrize("kind", [gemmi.ConnectionType.Covale,
                                 gemmi.ConnectionType.Disulf, gemmi.ConnectionType.MetalC])
def test_declared_heavy_link_prevents_virtual_thiol_even_without_short_contact(kind):
    target_res, target_atom, element = {
        gemmi.ConnectionType.Covale: ("LIG", "C1", "C"),
        gemmi.ConnectionType.Disulf: ("CYS", "SG", "S"),
        gemmi.ConnectionType.MetalC: ("ZN", "ZN", "Zn"),
    }[kind]
    s = structure(phe(), thiol(), residue(target_res, 3, [atom(target_atom, element, (20, 20, 20))]))
    assert any(h["X_atom"] == "SG" for h in core.detect_interactions_in_structure(s, "test"))
    crosslink(s, kind=kind, target_res=target_res, target_atom=target_atom)
    diagnostics = {}
    assert core.detect_interactions_in_structure(s, "test", diagnostics=diagnostics) == []
    assert diagnostics["cone_covalent_state_exclusions"]


def test_link_in_alternate_b_does_not_remove_free_thiol_a():
    s = structure(phe(), thiol(altloc="A"),
                  residue("LIG", 3, [atom("C1", "C", (20, 20, 20), altloc="B")]))
    crosslink(s, altloc="B")
    hits = core.detect_interactions_in_structure(s, "test")
    assert len(hits) == 1 and hits[0]["H_altloc"] == "A"


def test_unrecorded_short_heavy_bond_and_thiolate_are_not_sh_donors():
    free = structure(phe(), thiol())
    assert core.detect_interactions_in_structure(free, "test")
    substituted = structure(phe(), thiol(),
                            residue("LIG", 3, [atom("C1", "C", (10, 11.8, 13.8))]))
    assert core.detect_interactions_in_structure(substituted, "test") == []
    assert core.detect_interactions_in_structure(structure(phe(), thiol(charge=-1)), "test") == []


@pytest.mark.parametrize("defect", ["extra_heavy", "no_hydrogen"])
def test_modified_dictionary_prevents_free_cone_donor(monkeypatch, defect):
    s = structure(phe(), ser())
    assert core.detect_interactions_in_structure(s, "test")
    if defect == "extra_heavy":
        monkeypatch.setattr(monomer_bonds, "get_bonded_heavy_atom_names",
                            lambda *args: {"CB", "P"})
    else:
        monkeypatch.setattr(monomer_bonds, "get_bonded_hydrogen_names", lambda *args: set())
    assert core.detect_interactions_in_structure(s, "test") == []


def test_water_exclusion_filters_donors_without_mutating_input_or_hbond_environment():
    donor = ser()
    # Continue the explicit X-H vector; water is 1.8 Å beyond H.
    water = residue("HOH", 3, [atom("O", "O", (10, 10, 10.2)),
                                atom("H1", "H", (10, 10, 11.17))])
    water_donor = residue("HOH", 4, [atom("O", "O", (10, 10, 7)),
                                      atom("H1", "H", (10, 10, 8))])
    s = structure(phe(), donor, water, water_donor)
    before = s.make_mmcif_document().as_string()
    hits = core.detect_interactions_in_structure(s, "test", cone_mode="none")
    assert len(hits) == 1 and hits[0]["X_res"] == "SER"
    assert hits[0]["hbond_acceptor_res"] == "HOH"
    assert s.make_mmcif_document().as_string() == before
    included = core.detect_interactions_in_structure(s, "test", cone_mode="none", include_water=True)
    assert any(h["X_res"] == "HOH" for h in included)
    assert next(h for h in included if h["X_res"] == "SER")["hbond_HA_dist"] == hits[0]["hbond_HA_dist"]


def test_water_environment_retains_alternative_cone_hbond():
    d = ser()
    definition = rotatable_groups.get_rotatable_group("SER", "OG")
    x, parent = np.array([10, 10, 13]), np.array([11.42, 10, 13])
    conformers = cone.generate_conformers(parent, x, definition)
    h = max((c.hydrogen_positions[0] for c in conformers), key=lambda p: p[2])
    oxygen = h + 1.8 * (h - x) / np.linalg.norm(h - x)
    s = structure(phe(), d, residue("HOH", 3, [atom("O", "O", oxygen)]))
    excluded = core.detect_interactions_in_structure(s, "test")
    included = core.detect_interactions_in_structure(s, "test", include_water=True)
    assert len(excluded) == len(included) == 1
    assert excluded[0]["hbond_relation"] == included[0]["hbond_relation"] == "alternative_conformer"


@pytest.mark.parametrize("h_occ,parent_occ", [(0, 1), (1, 0)])
def test_unoccupied_hydrogen_or_cone_parent_cannot_make_a_hit(h_occ, parent_occ):
    s = structure(phe(), ser(h_occ=h_occ, parent_occ=parent_occ))
    mode = "none" if h_occ == 0 else "auto"
    assert core.detect_interactions_in_structure(s, "test", cone_mode=mode) == []


def test_explicit_h_and_shared_x_parent_must_match_pi_altloc():
    ring = phe()
    for a in ring:
        a.altloc = "A"
    explicit = structure(ring, ser(h_alt="B"))
    assert core.detect_interactions_in_structure(explicit, "test", cone_mode="none") == []
    shared = structure(ring, ser(parent_alt="B"))
    assert core.detect_interactions_in_structure(shared, "test") == []
    compatible = structure(ring, ser(parent_alt="A"))
    assert core.detect_interactions_in_structure(compatible, "test")[0]["H_altloc"] == "A"


def context_fixture(altloc="A", occ=1.0, b=10, z=13):
    s = structure(residue("ASP", 3, [atom("OD1", "O", (10, 10, z),
                                          altloc=altloc, occ=occ, b=b)]))
    hit = {"X_chain": "A", "X_res": "SER", "X_id": "2", "X_atom": "OG",
           "X_altloc": "A", "H_altloc": "A", "pi_altloc": "A",
           "X_xyz_x": 10, "X_xyz_y": 10, "X_xyz_z": 10,
           "H_xyz_x": 10, "H_xyz_y": 10, "H_xyz_z": 11,
           "pi_center_x": 10, "pi_center_y": 10, "pi_center_z": 15,
           "symmetry_code": "1_555", "combined_occupancy": .6,
           "H_source": "input_model"}
    return s, hit


@pytest.mark.parametrize("altloc,occ,b,z,min_occ,max_b,expected", [
    ("A", 1, 10, 13, 0, 0, 1), ("\0", 1, 10, 13, 0, 0, 1),
    ("B", 1, 10, 13, 0, 0, 0), ("A", 0, 10, 13, 0, 0, 0),
    ("A", .09, 10, 13, 0, 0, 0), ("A", .15, 10, 13, 0, 0, 1),
    ("A", .4, 10, 13, .5, 0, 0), ("A", .5, 10, 13, .5, 0, 1),
    ("A", 1, 50, 13, 0, 20, 0), ("A", 1, 10, 12, 0, 0, 0),
])
def test_hbond_context_uses_occupied_compatible_conformer(
        altloc, occ, b, z, min_occ, max_b, expected):
    s, hit = context_fixture(altloc, occ, b, z)
    hbond_annotations.annotate_explicit_hbond_context([hit], s, min_occ=min_occ, max_b=max_b)
    hbond_annotations.annotate_cone_hbond_context([hit], s, min_occ=min_occ, max_b=max_b)
    assert hit["hbond_competing"] == hit["hbond_also_hbonded"] == expected


def test_unoccupied_hydrogen_does_not_define_histidine_protonation():
    his = residue("HIS", 3, [atom("ND1", "N", (10, 10, 10), altloc="A"),
                              atom("HD1", "H", (11, 10, 10), altloc="A", occ=0)])
    assert hbond_acceptors.is_hbond_acceptor(his, his.find_atom("ND1", "A"))
    his.find_atom("HD1", "A").occ = 1
    assert not hbond_acceptors.is_hbond_acceptor(his, his.find_atom("ND1", "A"))


def test_hbond_acceptor_image_is_selected_around_h_not_around_ring():
    # H is near a cell face; the acceptor must be the +a copy, while the ring
    # centroid's nearest copy is the ASU atom. Geometry is intentionally isolated.
    s, hit = context_fixture(z=10)
    a = s[0][0][0].find_atom("OD1", "A")
    a.pos = gemmi.Position(1, 10, 10)
    hit.update(X_xyz_x=28, X_xyz_z=10, H_xyz_x=29, H_xyz_z=10,
               pi_center_x=14, pi_center_z=10)
    hbond_annotations.annotate_cone_hbond_context([hit], s)
    assert hit["hbond_competing"] == 1
    assert hit["hbond_HA_dist"] == 2
    assert hit["hbond_acceptor_id"] == "3"
    assert hit["hbond_acceptor_symmetry_code"] == "1_655"


def tshape_structure(translation=0, operation=None):
    ring = phe()
    coordinates = [("CG", (0, 0, 1.4)), ("CD1", (1.212, 0, .7)),
                   ("CE1", (1.212, 0, -.7)), ("CZ", (0, 0, -1.4)),
                   ("CE2", (-1.212, 0, -.7)), ("CD2", (-1.212, 0, .7))]
    donor_atoms = [atom(n, "C", np.array(p) + (10, 10, 15.2)) for n, p in coordinates]
    donor_atoms.append(atom("HZ", "H", (10, 10, 12.708)))
    for a in donor_atoms:
        if operation is not None:
            a.pos = gemmi.Position(*operation.apply_to_xyz([a.pos.x, a.pos.y, a.pos.z]))
        a.pos.x += translation
    return structure(ring, residue("PHE", 2, donor_atoms),
                     spacegroup="P 2" if operation else "P 1")


@pytest.mark.parametrize("translation,rotation", [(0, False), (30, False), (0, True)])
def test_tshape_annotation_uses_whole_donor_ring_image(translation, rotation):
    # The P2 twofold is its own inverse; applying it to ASU coordinates must
    # yield the same detected donor copy and the same perpendicular ring.
    op = gemmi.Op("-x,y,-z") if rotation else None
    s = tshape_structure(translation, op)
    hits = core.detect_interactions_in_structure(s, "test", cone_mode="none", sym_contacts=True)
    selected = [h for h in hits if h["pi_id"] == "1" and h["X_id"] == "2" and h["H_atom"] == "HZ"]
    assert len(selected) == 1
    assert selected[0]["is_pi_pi_tshaped"] == 1
    assert selected[0]["dist_X_centroid"] == pytest.approx(3.8, abs=.002)
    if translation:
        assert selected[0]["symmetry_code"] == "1_455"
    if rotation:
        assert selected[0]["sym_op"] != 0


def test_two_trp_rings_and_ring_alternates_do_not_share_cooperativity():
    base = {"pdb": "test", "model": "1", "pi_chain": "A", "pi_res": "TRP",
            "pi_id": "1", "X_chain": "A", "X_res": "SER", "X_atom": "OG"}
    hits = [dict(base, pi_ring_id="ring1", pi_altloc="A", X_id="2", X_side_of_pi=1),
            dict(base, pi_ring_id="ring2", pi_altloc="A", X_id="3", X_side_of_pi=-1),
            dict(base, pi_ring_id="ring1", pi_altloc="B", X_id="4", X_side_of_pi=-1)]
    cooperativity.annotate_cooperativity(hits)
    assert all(h["coop_total_donors"] == h["coop_same_face_donors"] == 1 for h in hits)
    assert all(h["coop_bi_face"] == h["coop_opp_face_donors"] == 0 for h in hits)


def test_cooperativity_is_independent_of_coordinate_export():
    s = structure(phe(), ser(), residue("SER", 3, [atom("OG", "O", (10, 10, 7)),
                                                    atom("HG", "H", (10, 10, 8))]))
    default = core.detect_interactions_in_structure(s, "test", cone_mode="none")
    coordinates = core.detect_interactions_in_structure(s, "test", cone_mode="none", include_coordinates=True)
    assert len(default) == len(coordinates) == 2
    for h in default + coordinates:
        assert h["coop_bi_face"] == h["coop_same_face_donors"] == h["coop_opp_face_donors"] == 1
        assert h["coop_total_donors"] == 2


@pytest.mark.parametrize("acceptor_alt,occ,b,min_occ,max_b,relation", [
    ("A", 1, 10, 0, 0, "alternative_conformer"),
    ("B", 1, 10, 0, 0, "none"),
    ("A", 0, 10, 0, 0, "none"),
    ("A", .4, 10, .5, 0, "none"),
    ("A", 1, 50, 0, 20, "none"),
])
def test_cone_environment_uses_shared_x_parent_state_and_context_filters(
        acceptor_alt, occ, b, min_occ, max_b, relation):
    x, parent = np.array([10, 10, 13]), np.array([11.42, 10, 13])
    definition = rotatable_groups.get_rotatable_group("SER", "OG")
    conformers = cone.generate_conformers(parent, x, definition)
    h = max((c.hydrogen_positions[0] for c in conformers), key=lambda p: p[2])
    oxygen = h + 1.8 * (h - x) / np.linalg.norm(h - x)
    s = structure(phe(), ser(parent_alt="A"),
                  residue("HOH", 3, [atom("O", "O", oxygen,
                                           altloc=acceptor_alt, occ=occ, b=b)]))
    hits = core.detect_interactions_in_structure(s, "test", min_occ=min_occ, max_b=max_b)
    assert len(hits) == 1
    assert hits[0]["H_altloc"] == "A"
    assert hits[0]["hbond_relation"] == relation


def test_alternate_metal_does_not_block_shared_x_with_parent_a():
    s = structure(phe(), ser(parent_alt="A"),
                  residue("ZN", 3, [atom("ZN", "Zn", (10, 12.3, 13), altloc="B")]))
    assert core.detect_interactions_in_structure(s, "test")
    s[0][0][2].find_atom("ZN", "B").altloc = "A"
    assert core.detect_interactions_in_structure(s, "test") == []


def test_unoccupied_clash_atom_is_ignored_by_scalar_and_vectorized_cone():
    x = np.array([0., 0., 0.])
    h = np.array([0., 0., 1.])
    a = atom("C1", "C", h, occ=0)
    env = [cone.EnvironmentAtom(h, a, residue("LIG", 1, [a]), "A", 0)]
    conformers = [cone.HydrogenConformer(0, (h,))]
    valid, hb = cone._classify_conformer_arrays(conformers, x, env)
    assert bool(valid[0]) is True and not hb.any()
    assert cone._hydrogen_state(x, h, env) == (True, False)


def test_tshape_requires_compatible_donor_ring_state():
    s = tshape_structure()
    for a in s[0][0][0]:
        a.altloc = "A"
    for a in s[0][0][1]:
        a.altloc = "\0" if a.name == "CZ" else "A" if a.name == "HZ" else "B"
    hits = core.detect_interactions_in_structure(s, "test", cone_mode="none")
    assert len(hits) == 1
    assert hits[0]["is_pi_pi_tshaped"] == 0


def test_nonring_sidechain_h_is_not_labelled_as_tshape_donor():
    s = tshape_structure()
    d = s[0][0][1]
    d.remove_atom("HZ", "\0", gemmi.Element("H"))
    d.add_atom(atom("CB", "C", (10, 10, 13.8)))
    d.add_atom(atom("HB2", "H", (10, 10, 12.708)))
    hits = core.detect_interactions_in_structure(s, "test", cone_mode="none")
    assert len(hits) == 1 and hits[0]["X_atom"] == "CB"
    assert hits[0]["is_pi_pi_tshaped"] == 0


def test_covalent_state_reads_single_block_component_dictionary(tmp_path, monkeypatch):
    dictionary = tmp_path / "cys.cif"
    dictionary.write_text('''data_component
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
SG S
CB C
HG H
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
SG CB
SG HG
''')
    monkeypatch.setattr(monomer_bonds.monlib, "find_monomer_cif", lambda code: dictionary)
    assert monomer_bonds.get_bonded_heavy_atom_names("CYS", "SG") == {"CB"}
    assert monomer_bonds.get_bonded_hydrogen_names("CYS", "SG") == {"HG"}
