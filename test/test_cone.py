import math
from types import SimpleNamespace

import gemmi
import numpy as np
import pytest

from xpid import acceptors
from xpid import cone
from xpid import donors
from xpid import xhpi_criteria


def _atom(name, element, xyz=(0.0, 0.0, 0.0), occ=1.0, altloc="\0"):
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element(element)
    atom.pos = gemmi.Position(*xyz)
    atom.occ = occ
    atom.altloc = altloc
    return atom


def _residue(name, atoms):
    residue = gemmi.Residue()
    residue.name = name
    residue.seqid = gemmi.SeqId(1, " ")
    for atom in atoms:
        residue.add_atom(atom)
    return residue


def _angle(a, b, c):
    ba = a - b
    bc = c - b
    cosine = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))


def test_ser_conformer_uses_group_specific_geometry():
    definition = donors.get_definition("SER", "OG")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.42, 0.0, 0.0])

    conformers = cone.generate_conformers(parent, x_pos, definition)

    assert len(conformers) == 360
    assert all(len(item.hydrogen_positions) == 1 for item in conformers)
    h_pos = conformers[0].hydrogen_positions[0]
    assert np.linalg.norm(h_pos - x_pos) == pytest.approx(0.972, abs=1e-8)
    assert _angle(parent, x_pos, h_pos) == pytest.approx(108.539, abs=1e-8)


def test_cys_cone_uses_ccp4_nuclear_geometry():
    definition = donors.get_definition("CYS", "SG")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.812, 0.0, 0.0])

    h_pos = cone.generate_conformers(
        parent, x_pos, definition)[0].hydrogen_positions[0]

    assert np.linalg.norm(h_pos - x_pos) == pytest.approx(1.338, abs=1e-8)
    assert _angle(parent, x_pos, h_pos) == pytest.approx(97.543, abs=1e-8)


def test_methyl_conformer_contains_three_coupled_hydrogens():
    definition = donors.get_definition("ALA", "CB")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.513, 0.0, 0.0])

    conformers = cone.generate_conformers(parent, x_pos, definition)

    assert len(conformers) == 120
    assert all(len(item.hydrogen_positions) == 3 for item in conformers)
    for h_pos in conformers[17].hydrogen_positions:
        assert np.linalg.norm(h_pos - x_pos) == pytest.approx(1.092, abs=1e-8)
        assert _angle(parent, x_pos, h_pos) == pytest.approx(109.742, abs=1e-8)
    h1, h2, h3 = conformers[17].hydrogen_positions
    assert _angle(h1, x_pos, h2) == pytest.approx(109.2, abs=0.1)
    assert _angle(h1, x_pos, h3) == pytest.approx(109.2, abs=0.1)


def test_any_hydrogen_clash_invalidates_complete_methyl_conformer():
    definition = donors.get_definition("ALA", "CB")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.513, 0.0, 0.0])
    conformer = cone.generate_conformers(parent, x_pos, definition)[0]

    blocker = _atom("C1", "C", tuple(conformer.hydrogen_positions[1]))
    ligand = _residue("LIG", [blocker])
    environment = [
        cone.EnvironmentAtom(
            conformer.hydrogen_positions[1], blocker, ligand, "B", 0)
    ]

    valid, constrained = cone.classify_conformers(
        [conformer], x_pos, environment)

    assert valid == []
    assert constrained == []


def test_valid_short_hbond_is_not_rejected_as_steric_clash():
    definition = donors.get_definition("SER", "OG")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.42, 0.0, 0.0])
    conformer = cone.generate_conformers(parent, x_pos, definition)[0]
    h_pos = conformer.hydrogen_positions[0]
    direction = (h_pos - x_pos) / np.linalg.norm(h_pos - x_pos)
    acceptor_pos = h_pos + 1.8 * direction

    oxygen = _atom("OD1", "O", tuple(acceptor_pos))
    asp = _residue("ASP", [oxygen])
    environment = [
        cone.EnvironmentAtom(acceptor_pos, oxygen, asp, "B", 0)
    ]

    valid, constrained = cone.classify_conformers(
        [conformer], x_pos, environment)

    assert valid == [conformer]
    assert constrained == [conformer]


def test_nonlocking_hbond_contact_is_not_a_steric_clash():
    definition = donors.get_definition("SER", "OG")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.42, 0.0, 0.0])
    conformer = cone.generate_conformers(parent, x_pos, definition)[0]
    h_pos = conformer.hydrogen_positions[0]

    # Construct a 1.8 Å H...O contact at a 130° D-H...A angle.  It is a
    # chemically valid contact, but deliberately below the 140° lock cutoff.
    h_to_x = (x_pos - h_pos) / np.linalg.norm(x_pos - h_pos)
    perpendicular = np.cross(h_to_x, np.array([0.0, 0.0, 1.0]))
    if np.linalg.norm(perpendicular) < 1e-8:
        perpendicular = np.cross(h_to_x, np.array([0.0, 1.0, 0.0]))
    perpendicular /= np.linalg.norm(perpendicular)
    angle = np.radians(130.0)
    h_to_a = np.cos(angle) * h_to_x + np.sin(angle) * perpendicular
    acceptor_pos = h_pos + 1.8 * h_to_a

    oxygen = _atom("OD1", "O", tuple(acceptor_pos))
    asp = _residue("ASP", [oxygen])
    environment = [
        cone.EnvironmentAtom(acceptor_pos, oxygen, asp, "B", 0)
    ]

    valid, constrained = cone.classify_conformers(
        [conformer], x_pos, environment)

    assert valid == [conformer]
    assert constrained == []


def test_steric_only_detector_does_not_lock_to_hbond_direction():
    definition = donors.get_definition("SER", "OG")
    parent = np.array([1.42, 0.0, 3.0])
    x_pos = np.array([0.0, 0.0, 3.0])
    conformers = cone.generate_conformers(parent, x_pos, definition)
    away = max(
        conformers,
        key=lambda item: item.hydrogen_positions[0][2],
    )
    away_h = away.hydrogen_positions[0]
    direction = (away_h - x_pos) / np.linalg.norm(away_h - x_pos)
    acceptor_pos = away_h + 1.8 * direction

    oxygen = _atom("OD1", "O", tuple(acceptor_pos))
    asp = _residue("ASP", [oxygen])
    environment = [
        cone.EnvironmentAtom(acceptor_pos, oxygen, asp, "B", 0)
    ]
    ring_context = SimpleNamespace(
        pi_center_arr=np.array([0.0, 0.0, 0.0]),
        pi_normal=np.array([0.0, 0.0, 1.0]),
        p_radius=2.0,
        p_slab_half_thickness=0.5,
    )

    default = cone.evaluate_binary(
        ring_context, definition, parent, x_pos, "O", environment,
        dist_x_plane=3.0, dist_x_centroid=3.0, proj_dist=0.0,
    )
    assert default is not None
    assert default.hbond_relation == "alternative_conformer"


def test_lys_nitrogen_cannot_lock_a_cone_as_acceptor():
    nz = _atom("NZ", "N")
    lys = _residue("LYS", [nz])
    assert acceptors.is_hbond_acceptor(lys, nz) is False


def test_backbone_and_amide_acceptor_typing():
    backbone_o = _atom("O", "O")
    backbone_n = _atom("N", "N")
    amide_o = _atom("OD1", "O")
    amide_n = _atom("ND2", "N")

    ala = _residue("ALA", [backbone_o, backbone_n])
    asn = _residue("ASN", [amide_o, amide_n])

    assert acceptors.is_hbond_acceptor(ala, backbone_o) is True
    assert acceptors.is_hbond_acceptor(ala, backbone_n) is False
    assert acceptors.is_hbond_acceptor(asn, amide_o) is True
    assert acceptors.is_hbond_acceptor(asn, amide_n) is False


def test_hydroxyl_oxygen_remains_acceptor_when_protonated():
    oxygen = _atom("OG", "O")
    hydrogen = _atom("HG", "H", (0.96, 0.0, 0.0))
    ser = _residue("SER", [oxygen, hydrogen])

    assert acceptors.is_hbond_acceptor(ser, oxygen) is True


def test_acceptor_protonation_ignores_incompatible_hydrogen_altloc():
    sulfur = _atom("SD", "S", altloc="B")
    hydrogen_a = _atom("H1", "H", (0.9, 0.0, 0.0), altloc="A")
    met = _residue("MET", [sulfur, hydrogen_a])
    sulfur = next(atom for atom in met if atom.name == "SD")

    assert acceptors.is_hbond_acceptor(met, sulfur) is True

    met.add_atom(_atom("H2", "H", (0.9, 0.0, 0.0), altloc="B"))
    # Gemmi may reallocate the residue's atom vector in add_atom(); reacquire
    # the wrapper instead of retaining a potentially invalid C++ reference.
    sulfur = next(atom for atom in met if atom.name == "SD")
    assert acceptors.is_hbond_acceptor(met, sulfur) is False


def test_altloc_compatibility():
    assert acceptors.altlocs_compatible("", "A")
    assert acceptors.altlocs_compatible("A", "")
    assert acceptors.altlocs_compatible("A", "A")
    assert not acceptors.altlocs_compatible("A", "B")


@pytest.mark.parametrize("parent_order", [("A", "B"), ("B", "A")])
def test_donor_parent_resolution_matches_altloc_independent_of_atom_order(
        parent_order):
    parents = {
        "A": _atom("CG", "C", (0.0, 1.5, 0.0), altloc="A"),
        "B": _atom("CG", "C", (1.5, 0.0, 0.0), occ=0.35, altloc="B"),
    }
    donor = _atom("CD1", "C", (3.0, 0.0, 0.0), occ=0.35, altloc="B")
    residue = _residue(
        "LEU", [parents[alt] for alt in parent_order] + [donor])
    donor = next(atom for atom in residue
                 if atom.name == "CD1" and atom.altloc == "B")

    resolution = donors.resolve_donor_conformers(
        residue, donor, donors.get_definition("LEU", "CD1"),
        model_index=2, chain_name="A")

    assert resolution.issue is None
    assert len(resolution.conformers) == 1
    conformer = resolution.conformers[0]
    assert conformer.model_index == 2
    assert conformer.x_altloc == "B"
    assert conformer.parent_altloc == "B"
    assert conformer.parent_atom.pos.x == pytest.approx(1.5)
    assert conformer.parent_selection == "matching_altloc"
    assert conformer.occupancy == pytest.approx(0.35)


def test_labelled_donor_can_use_one_shared_blank_parent():
    residue = _residue("LEU", [
        _atom("CG", "C", (1.5, 0.0, 0.0), altloc="\0"),
        _atom("CD1", "C", (3.0, 0.0, 0.0), altloc="B"),
    ])
    donor = next(atom for atom in residue if atom.name == "CD1")

    resolution = donors.resolve_donor_conformers(
        residue, donor, donors.get_definition("LEU", "CD1"), 0, "A")

    assert resolution.issue is None
    assert len(resolution.conformers) == 1
    assert resolution.conformers[0].parent_altloc == ""
    assert resolution.conformers[0].parent_selection == "shared_blank_parent"


def test_blank_donor_expands_unique_alternate_parent_conformers():
    residue = _residue("LEU", [
        _atom("CG", "C", (0.0, 1.5, 0.0), occ=0.6, altloc="B"),
        _atom("CG", "C", (1.5, 0.0, 0.0), occ=0.4, altloc="A"),
        _atom("CD1", "C", (3.0, 0.0, 0.0), altloc="\0"),
    ])
    donor = next(atom for atom in residue if atom.name == "CD1")

    resolution = donors.resolve_donor_conformers(
        residue, donor, donors.get_definition("LEU", "CD1"), 0, "A")

    assert resolution.issue is None
    assert [item.parent_altloc for item in resolution.conformers] == ["A", "B"]
    assert all(item.active_altloc in {"A", "B"}
               for item in resolution.conformers)


def test_incompatible_or_duplicate_parent_is_auditable_not_order_dependent():
    incompatible = _residue("LEU", [
        _atom("CG", "C", altloc="A"),
        _atom("CD1", "C", altloc="B"),
    ])
    donor = next(atom for atom in incompatible if atom.name == "CD1")
    resolution = donors.resolve_donor_conformers(
        incompatible, donor, donors.get_definition("LEU", "CD1"), 0, "A")
    assert resolution.conformers == ()
    assert resolution.issue == "incompatible_parent_altloc:A"

    duplicate = _residue("LEU", [
        _atom("CG", "C", (0.0, 0.0, 0.0), altloc="B"),
        _atom("CG", "C", (0.1, 0.0, 0.0), altloc="B"),
        _atom("CD1", "C", altloc="B"),
    ])
    donor = next(atom for atom in duplicate if atom.name == "CD1")
    resolution = donors.resolve_donor_conformers(
        duplicate, donor, donors.get_definition("LEU", "CD1"), 0, "A")
    assert resolution.conformers == ()
    assert resolution.issue == "duplicate_parent_altloc:B"


def test_vectorized_allowed_set_matches_scalar_reference():
    definition = donors.get_definition("SER", "OG")
    parent = np.array([0.0, 0.0, 0.0])
    x_pos = np.array([1.42, 0.0, 0.0])
    conformers = cone.generate_conformers(parent, x_pos, definition)

    oxygen = _atom("OD1", "O", (1.4, 1.8, 0.2))
    carbon = _atom("C1", "C", (1.4, -1.6, 0.1))
    asp = _residue("ASP", [oxygen])
    ligand = _residue("LIG", [carbon])
    environment = [
        cone.EnvironmentAtom(np.array([1.4, 1.8, 0.2]), oxygen, asp, "B", 0),
        cone.EnvironmentAtom(np.array([1.4, -1.6, 0.1]), carbon, ligand, "B", 0),
    ]

    valid_mask, strong_flags = cone._classify_conformer_arrays(
        conformers, x_pos, environment)
    for conformer_index, conformer in enumerate(conformers):
        scalar_states = [
            cone._hydrogen_state(x_pos, h_pos, environment)
            for h_pos in conformer.hydrogen_positions
        ]
        assert valid_mask[conformer_index] == all(
            state[0] for state in scalar_states)
        if valid_mask[conformer_index]:
            assert strong_flags[conformer_index].tolist() == [
                state[1] for state in scalar_states
            ]


def test_vectorized_hudson_plevin_matches_scalar_geometry():
    definition = donors.get_definition("SER", "OG")
    parent = np.array([1.42, 0.0, 3.0])
    x_pos = np.array([0.0, 0.0, 3.0])
    conformers = cone.generate_conformers(parent, x_pos, definition)
    h_positions = np.asarray([
        conformer.hydrogen_positions[0] for conformer in conformers])
    ring_context = SimpleNamespace(
        pi_center_arr=np.array([0.0, 0.0, 0.0]),
        pi_normal=np.array([0.0, 0.0, 1.0]),
        p_radius=2.0,
        p_slab_half_thickness=0.5,
    )
    spatial = xhpi_criteria.prepare_spatial_criteria(
        ring_context, "O", x_pos, 3.0, 0.0)

    theta, xh_pi, hudson, plevin = (
        xhpi_criteria.evaluate_hudson_plevin_batch(
            ring_context, x_pos, h_positions, spatial))

    for index, h_pos in enumerate(h_positions):
        scalar = xhpi_criteria.evaluate_xhpi_geometry(
            ring_context, "O", x_pos, h_pos, 3.0, 3.0, 0.0,
            spatial=spatial, include_direction_metrics=False)
        assert bool(hudson[index]) == bool(scalar["is_hudson"])
        assert bool(plevin[index]) == bool(scalar["is_plevin"])
        if scalar["theta"] is None:
            assert np.isnan(theta[index])
        else:
            assert theta[index] == pytest.approx(
                scalar["theta"], abs=1e-12)
        assert xh_pi[index] == pytest.approx(
            scalar["angle_XH_Pi"], abs=1e-12)
