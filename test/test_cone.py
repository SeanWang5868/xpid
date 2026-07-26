import math
from types import SimpleNamespace

import gemmi
import numpy as np
import pytest

from xpid import acceptors
from xpid import cone
from xpid import donors


def _atom(name, element, xyz=(0.0, 0.0, 0.0), occ=1.0):
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element(element)
    atom.pos = gemmi.Position(*xyz)
    atom.occ = occ
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


def test_cys_cone_uses_non_tetrahedral_sulfur_angle():
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

    valid, constrained = cone.filter_conformers(
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

    valid, constrained = cone.filter_conformers(
        [conformer], x_pos, environment)

    assert valid == [conformer]
    assert constrained == [conformer]


def test_strong_hbond_direction_can_prevent_single_h_xhpi():
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

    unconstrained = cone.evaluate_binary(
        ring_context, definition, parent, x_pos, "O", environment,
        dist_x_plane=3.0, dist_x_centroid=3.0, proj_dist=0.0,
        disable_hbond_constraint=True,
    )
    constrained = cone.evaluate_binary(
        ring_context, definition, parent, x_pos, "O", environment,
        dist_x_plane=3.0, dist_x_centroid=3.0, proj_dist=0.0,
    )

    assert unconstrained is not None
    assert constrained is None


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


def test_altloc_compatibility():
    assert acceptors.altlocs_compatible("", "A")
    assert acceptors.altlocs_compatible("A", "")
    assert acceptors.altlocs_compatible("A", "A")
    assert not acceptors.altlocs_compatible("A", "B")
