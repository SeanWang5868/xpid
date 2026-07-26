from pathlib import Path

import pytest

from xpid import aromatic_rings


def _write_component_cif(
        path: Path,
        component: str,
        bonds: list[tuple[str, str, str]],
        plane_atoms: list[str] | None = None) -> None:
    lines = [
        f"data_comp_{component}",
        "loop_",
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.aromatic",
    ]
    lines.extend(f"{atom_1} {atom_2} {aromatic}"
                 for atom_1, atom_2, aromatic in bonds)
    if plane_atoms:
        lines.extend([
            "loop_",
            "_chem_comp_plane_atom.plane_id",
            "_chem_comp_plane_atom.atom_id",
        ])
        lines.extend(f"plane-1 {atom}" for atom in plane_atoms)
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture(autouse=True)
def _clear_ring_cache():
    aromatic_rings.clear_cache()
    yield
    aromatic_rings.clear_cache()


def test_plane_only_component_is_not_aromatic(tmp_path, monkeypatch):
    cif_path = tmp_path / "PLN.cif"
    atoms = [f"C{index}" for index in range(1, 7)]
    nonaromatic_bonds = [
        (atoms[index], atoms[(index + 1) % 6], "n")
        for index in range(6)
    ]
    _write_component_cif(
        cif_path, "PLN", nonaromatic_bonds, plane_atoms=atoms)
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif",
        lambda component: cif_path if component == "PLN" else None,
    )

    assert aromatic_rings.get_aromatic_rings("PLN") == []


def test_ligand_requires_complete_aromatic_cycle(tmp_path, monkeypatch):
    cif_path = tmp_path / "BRK.cif"
    _write_component_cif(
        cif_path,
        "BRK",
        [
            ("C1", "C2", "y"),
            ("C2", "C3", "y"),
            ("C3", "C4", "y"),
            ("C4", "C5", "y"),
            ("C5", "C6", "y"),
        ],
    )
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: cif_path)

    assert aromatic_rings.get_aromatic_rings("BRK") == []


def test_aromatic_ligand_six_membered_cycle(tmp_path, monkeypatch):
    cif_path = tmp_path / "BEN.cif"
    atoms = [f"C{index}" for index in range(1, 7)]
    _write_component_cif(
        cif_path,
        "BEN",
        [
            (atoms[index], atoms[(index + 1) % 6], "y")
            for index in range(6)
        ],
    )
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: cif_path)

    assert aromatic_rings.get_aromatic_rings("BEN") == [set(atoms)]


def test_fused_aromatic_graph_returns_only_two_minimal_cycles():
    edges = [
        ("A1", "A2"), ("A2", "A3"), ("A3", "A4"),
        ("A4", "A5"), ("A5", "A6"), ("A6", "A1"),
        ("A3", "B1"), ("B1", "B2"), ("B2", "B3"),
        ("B3", "B4"), ("B4", "A4"),
    ]

    rings = aromatic_rings._minimal_aromatic_cycles(edges)

    assert rings == [
        {"A1", "A2", "A3", "A4", "A5", "A6"},
        {"A3", "A4", "B1", "B2", "B3", "B4"},
    ]


def test_trp_aromatic_graph_preserves_six_then_five_ring_order(
        tmp_path, monkeypatch):
    cif_path = tmp_path / "TRP.cif"
    aromatic_edges = [
        ("CG", "CD1", "y"),
        ("CG", "CD2", "y"),
        ("CD1", "NE1", "y"),
        ("CD2", "CE2", "y"),
        ("CD2", "CE3", "y"),
        ("NE1", "CE2", "y"),
        ("CE2", "CZ2", "y"),
        ("CE3", "CZ3", "y"),
        ("CZ2", "CH2", "y"),
        ("CZ3", "CH2", "y"),
    ]
    _write_component_cif(cif_path, "TRP", aromatic_edges)
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: cif_path)

    assert (
        aromatic_rings.get_aromatic_rings("TRP")
        == aromatic_rings.STANDARD_AROMATIC_RINGS["TRP"]
    )


@pytest.mark.parametrize("residue_name", ["PHE", "TYR", "TRP", "HIS"])
def test_standard_residue_falls_back_without_dictionary(
        residue_name, monkeypatch):
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: None)

    assert (
        aromatic_rings.get_aromatic_rings(residue_name)
        == aromatic_rings.STANDARD_AROMATIC_RINGS[residue_name]
    )


def test_standard_residue_falls_back_when_aromatic_graph_is_broken(
        tmp_path, monkeypatch):
    cif_path = tmp_path / "PHE.cif"
    _write_component_cif(
        cif_path,
        "PHE",
        [("CG", "CD1", "y"), ("CD1", "CE1", "y")],
    )
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: cif_path)

    assert (
        aromatic_rings.get_aromatic_rings("PHE")
        == aromatic_rings.STANDARD_AROMATIC_RINGS["PHE"]
    )


def test_standard_fused_residue_falls_back_when_only_one_ring_is_complete(
        tmp_path, monkeypatch):
    cif_path = tmp_path / "TRP.cif"
    _write_component_cif(
        cif_path,
        "TRP",
        [
            ("CG", "CD1", "y"),
            ("CD1", "NE1", "y"),
            ("NE1", "CE2", "y"),
            ("CE2", "CD2", "y"),
            ("CD2", "CG", "y"),
        ],
    )
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: cif_path)

    assert (
        aromatic_rings.get_aromatic_rings("TRP")
        == aromatic_rings.STANDARD_AROMATIC_RINGS["TRP"]
    )


def test_nonstandard_component_has_no_fallback(monkeypatch):
    monkeypatch.setattr(
        aromatic_rings.monlib, "find_monomer_cif", lambda _: None)

    assert aromatic_rings.get_aromatic_rings("PTR") == []
    assert aromatic_rings.get_aromatic_rings("BER") == []
    assert aromatic_rings.get_aromatic_rings("4PO") == []
