import csv
import inspect
from pathlib import Path

import gemmi
from xpid import cli, core, config, detect, output, prep, resolver


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


def test_api_and_core_default_disable_cone_and_p_slab():
    api_sig = inspect.signature(detect)
    core_sig = inspect.signature(core.detect_interactions_in_structure)

    assert api_sig.parameters["use_cone"].default is False
    assert api_sig.parameters["include_p_slab"].default is False
    assert api_sig.parameters["report_xh_candidates"].default is False
    assert api_sig.parameters["include_coordinates"].default is False
    assert api_sig.parameters["residue_pair"].default is None
    assert core_sig.parameters["use_cone"].default is False
    assert core_sig.parameters["include_p_slab"].default is False
    assert core_sig.parameters["report_xh_candidates"].default is False
    assert core_sig.parameters["include_coordinates"].default is False
    assert core_sig.parameters["residue_pair"].default is None


def test_cli_parser_default_disables_cone_and_p_slab():
    parser = cli._build_parser()
    args = parser.parse_args(["dummy.cif"])

    assert args.use_cone is False
    assert args.include_p_slab is False
    assert args.report_xh_candidates is False
    assert args.include_coordinates is False
    assert args.residue_pair is None


def test_cli_parser_accepts_residue_pair_selectors():
    parser = cli._build_parser()
    args = parser.parse_args(["dummy.cif", "--residue-pair", "//A/12", "//A/18"])

    assert args.residue_pair == ["//A/12", "//A/18"]


def test_read_structure_falls_back_to_rcsb_for_corrupt_named_gzip(tmp_path, monkeypatch):
    from xpid import structure_io

    bad_local = tmp_path / "5abc.cif.gz"
    bad_local.write_bytes(b"truncated gzip")
    downloaded = tmp_path / "downloaded" / "5abc.cif.gz"
    sentinel = object()

    def fake_read_structure(path):
        path = Path(path)
        if path == bad_local:
            raise RuntimeError(
                f"Cannot determine uncompressed size of {path}\n"
                "Would it be 12312 -> 1263743296 bytes?"
            )
        if path == downloaded:
            return sentinel
        raise AssertionError(f"unexpected structure path: {path}")

    def fake_download(pdb_code):
        assert pdb_code == "5abc"
        downloaded.parent.mkdir()
        downloaded.write_bytes(b"valid replacement")
        return downloaded

    monkeypatch.setattr(structure_io.gemmi, "read_structure", fake_read_structure)
    monkeypatch.setattr(structure_io, "download_rcsb_mmcif", fake_download)

    assert structure_io.read_structure(bad_local) is sentinel


def test_read_structure_does_not_fallback_when_recovery_disabled(tmp_path, monkeypatch):
    from xpid import structure_io

    bad_local = tmp_path / "5abc.cif.gz"
    bad_local.write_bytes(b"truncated gzip")

    def fake_read_structure(path):
        raise RuntimeError(
            f"Cannot determine uncompressed size of {path}\n"
            "Would it be 12312 -> 1263743296 bytes?"
        )

    def fail_download(pdb_code):
        raise AssertionError("download should not be attempted")

    monkeypatch.setattr(structure_io.gemmi, "read_structure", fake_read_structure)
    monkeypatch.setattr(structure_io, "download_rcsb_mmcif", fail_download)

    try:
        structure_io.read_structure(bad_local, allow_remote_recovery=False)
    except RuntimeError as exc:
        assert "appears corrupted" in str(exc)
        assert "disabled for batch inputs" in str(exc)
    else:
        raise AssertionError("expected corrupt gzip error")


def test_cli_enables_remote_recovery_only_for_single_direct_file(tmp_path):
    single = tmp_path / "5abc.cif.gz"
    single.write_bytes(b"")
    other = tmp_path / "6def.cif.gz"
    other.write_bytes(b"")
    folder = tmp_path / "folder"
    folder.mkdir()
    nested = folder / "7ghi.cif.gz"
    nested.write_bytes(b"")
    pdb_list = tmp_path / "ids.txt"
    pdb_list.write_text("5abc\n", encoding="utf-8")

    parser = cli._build_parser()

    single_args = parser.parse_args([str(single)])
    assert cli._allow_remote_recovery(single_args, [single]) is True

    multi_args = parser.parse_args([str(single), str(other)])
    assert cli._allow_remote_recovery(multi_args, [single, other]) is False

    folder_args = parser.parse_args([str(folder)])
    assert cli._allow_remote_recovery(folder_args, [nested]) is False

    list_args = parser.parse_args([
        "--pdb-list", str(pdb_list), "--pdb-mirror", str(tmp_path)
    ])
    assert cli._allow_remote_recovery(list_args, [single]) is False


def test_directory_input_collects_structure_extensions_recursively(tmp_path):
    root = tmp_path / "structures"
    nested = root / "nested"
    nested.mkdir(parents=True)

    expected = [
        root / "5abc.cif",
        root / "6def.cif.gz",
        root / "model_one.mmcif",
        nested / "model_two.mmcif.gz",
        nested / "custom_name.pdb",
        nested / "custom_name.pdb.gz",
    ]
    ignored = [
        root / "notes.txt",
        nested / "map.ccp4",
    ]
    for candidate in expected + ignored:
        candidate.write_bytes(b"")

    found = resolver.gather_inputs([str(root)], None, None, None)

    assert found == sorted(candidate.resolve() for candidate in expected)


def test_cation_pi_donors_are_excluded_from_core_detection():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    lys = gemmi.Residue()
    lys.name = "LYS"
    lys.seqid = _seqid(2)
    lys.add_atom(_atom("NZ", "N", (0.0, 0.0, 3.0)))
    lys.add_atom(_atom("HZ1", "H", (0.0, 0.0, 2.0)))
    chain.add_residue(lys)

    assert core.detect_interactions_in_structure(st, "test", use_cone=False) == []


def test_residue_pair_filter_keeps_only_selected_inter_residue_contacts():
    st = _structure_with_phe_and_og()
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(3)
    ser.add_atom(_atom("OG", "O", (0.5, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0.5, 0.0, 2.0)))
    chain.add_residue(ser)

    all_hits = core.detect_interactions_in_structure(st, "test", use_cone=False)
    assert {hit["X_id"].strip() for hit in all_hits} == {"2", "3"}

    pair_hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", use_cone=False,
        residue_pair=("//A/1", "//A/2"))
    assert len(pair_hits) == 1
    assert pair_hits[0]["pi_id"].strip() == "1"
    assert pair_hits[0]["X_id"].strip() == "2"


def test_residue_pair_filter_is_direction_agnostic():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", use_cone=False,
        residue_pair=("//A/2", "//A/1"))

    assert len(hits) == 1
    assert hits[0]["pi_id"].strip() == "1"
    assert hits[0]["X_id"].strip() == "2"


def test_residue_pair_filter_returns_empty_for_unmatched_pair():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", use_cone=False,
        residue_pair=("//A/1", "//A/99"))

    assert hits == []


def test_residue_pair_filter_respects_chain_in_selection():
    st = _structure_with_phe_and_og()
    model = st[0]
    chain_b = gemmi.Chain("B")

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
        phe.add_atom(_atom(name, "C", xyz))
    chain_b.add_residue(phe)

    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 2.0)))
    chain_b.add_residue(ser)
    model.add_chain(chain_b)

    hits = core.detect_interactions_in_structure(
        st, "test", use_cone=False, residue_pair=("//A/1", "//A/2"))

    assert len(hits) == 1
    assert hits[0]["pi_chain"] == "A"
    assert hits[0]["X_chain"] == "A"


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


def test_lower_occupancy_altloc_can_contribute_interaction():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 5.2), occ=0.8, altloc="A"))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 4.2), occ=0.8, altloc="A"))
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0), occ=0.2, altloc="B"))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 2.0), occ=0.2, altloc="B"))
    chain.add_residue(ser)

    hits = core.detect_interactions_in_structure(st, "test", use_cone=False)

    assert len(hits) == 1
    assert hits[0]["X_atom"] == "OG"
    assert hits[0]["dist_X_Pi"] == 3.0
    assert hits[0]["is_hudson"] == 1
    assert hits[0]["is_plevin"] == 1


def test_duplicate_altloc_hits_are_collapsed_to_best_occupancy():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0), occ=0.7, altloc="A"))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 2.0), occ=0.7, altloc="A"))
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.1), occ=0.3, altloc="B"))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 2.1), occ=0.3, altloc="B"))
    chain.add_residue(ser)

    hits = core.detect_interactions_in_structure(st, "test", use_cone=False)

    assert len(hits) == 1
    assert hits[0]["dist_X_Pi"] == 3.0


def test_p_model_uses_plane_distance_not_centroid_distance():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (1.9, 0.0, 4.0)))
    ser.add_atom(_atom("HG", "H", (1.9, 0.0, 3.0)))
    chain.add_residue(ser)

    hits = core.detect_interactions_in_structure(
        st, "test", use_cone=False, include_p_slab=True)

    assert len(hits) == 1
    assert hits[0]["dist_X_Pi"] == 4.0
    assert hits[0]["dist_X_centroid"] == 4.428
    assert hits[0]["proj_dist"] == 1.9
    assert hits[0]["h_proj_dist"] == 1.9
    assert hits[0]["is_hudson"] == 0
    assert hits[0]["is_plevin"] == 0
    assert hits[0]["is_p_slab"] == 1


def test_p_model_rejects_h_ray_that_misses_p():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (1.0, 0.0, 3.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(
        st, "test", use_cone=False, include_p_slab=True) == []


def test_p_model_rejects_x_projection_outside_p():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (2.1, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (2.1, 0.0, 2.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(
        st, "test", use_cone=False, include_p_slab=True) == []


def test_p_model_is_disabled_by_default():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (1.9, 0.0, 4.0)))
    ser.add_atom(_atom("HG", "H", (1.9, 0.0, 3.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(st, "test", use_cone=False) == []


def test_xh_candidate_mode_reports_direction_failed_spatial_candidate():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 4.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(st, "test", use_cone=False) == []

    hits = core.detect_interactions_in_structure(
        st, "test", use_cone=False, report_xh_candidates=True)

    assert len(hits) == 1
    hit = hits[0]
    assert hit["is_xh_candidate"] == 1
    assert hit["is_hudson"] == 0
    assert hit["is_plevin"] == 0
    assert hit["is_hudson_spatial"] == 1
    assert hit["is_plevin_spatial"] == 1
    assert hit["hudson_direction_ok"] == 0
    assert hit["plevin_direction_ok"] == 0
    assert hit["xh_centroid_cos"] == -1.0
    assert hit["P_radius"] == 2.0
    assert hit["H_ray_t"] is None
    assert hit["H_plane_t"] is None
    assert "is_p_slab" not in hit


def test_xh_candidate_mode_keeps_positive_labels_and_ray_geometry():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", use_cone=False, report_xh_candidates=True)

    assert len(hits) == 1
    hit = hits[0]
    assert hit["is_hudson"] == 1
    assert hit["is_plevin"] == 1
    assert hit["hudson_direction_ok"] == 1
    assert hit["plevin_direction_ok"] == 1
    assert hit["h_proj_dist"] == 0.0
    assert hit["H_ray_t"] == 2.5
    assert hit["H_ray_entry_dist"] == 2.5
    assert hit["h_plane_proj_dist"] == 0.0
    assert hit["H_plane_t"] == 3.0
    assert hit["H_plane_entry_dist"] == 3.0


def test_xh_candidate_mode_ignores_cone_virtual_hydrogens():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("CB", "C", (1.0, 0.0, 3.0)))
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(
        st, "test", use_cone=True, report_xh_candidates=True) == []


def test_p_slab_accepts_near_edge_directional_ray():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.83, 0.0, 3.97)))
    ser.add_atom(_atom("HG", "H", (1.13, 0.0, 2.9775)))
    chain.add_residue(ser)

    hits = core.detect_interactions_in_structure(
        st, "test", use_cone=False, include_p_slab=True)

    assert len(hits) == 1
    assert hits[0]["proj_dist"] == 0.83
    assert hits[0]["h_proj_dist"] == 1.879
    assert hits[0]["P_slab_half_thickness"] == 0.5
    assert hits[0]["is_p_slab"] == 1


def test_legacy_hudson_hit_is_retained_when_p_slab_fails():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (1.9, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (2.4, 0.0, 2.134)))
    chain.add_residue(ser)

    hits = core.detect_interactions_in_structure(st, "test", use_cone=False)

    assert len(hits) == 1
    assert hits[0]["is_hudson"] == 1
    assert hits[0]["is_plevin"] == 0


def test_coordinate_output_includes_h_xyz_pi_normal_and_x_side():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", use_cone=False,
        include_coordinates=True)

    assert len(hits) == 1
    hit = hits[0]
    assert hit["pi_center_x"] == 0.0
    assert hit["pi_center_y"] == 0.0
    assert hit["pi_center_z"] == 0.0
    assert hit["X_xyz_x"] == 0.0
    assert hit["X_xyz_y"] == 0.0
    assert hit["X_xyz_z"] == 3.0
    assert hit["H_xyz_x"] == 0.0
    assert hit["H_xyz_y"] == 0.0
    assert hit["H_xyz_z"] == 2.0
    assert hit["pi_normal_x"] == 0.0
    assert hit["pi_normal_y"] == 0.0
    assert hit["pi_normal_z"] == 1.0
    assert hit["X_side_of_pi"] == 1


def test_hydrogen_merge_preserves_residues_with_experimental_h(monkeypatch):
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(1)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (9.0, 9.0, 9.0)))
    chain.add_residue(ser)

    thr = gemmi.Residue()
    thr.name = "THR"
    thr.seqid = _seqid(2)
    thr.add_atom(_atom("OG1", "O", (1.0, 0.0, 3.0)))
    chain.add_residue(thr)

    model.add_chain(chain)
    st.add_model(model)

    monkeypatch.setattr(prep, "_get_shared_monlib", lambda codes: gemmi.MonLib())

    def fake_prepare(working, monlib, h_change_val):
        for residue in working[0][0]:
            residue.add_atom(_atom("HX", "H", (1.0, 1.0, 1.0)))

    monkeypatch.setattr(prep, "_prepare_topology_with_retries", fake_prepare)

    prep.add_hydrogens_memory(st)

    ser_h = [atom for atom in st[0][0][0] if atom.element.name.upper() == "H"]
    thr_h = [atom for atom in st[0][0][1] if atom.element.name.upper() == "H"]

    assert len(ser_h) == 1
    assert ser_h[0].name == "HG"
    assert ser_h[0].pos.x == 9.0
    assert len(thr_h) == 1
    assert thr_h[0].name == "HX"


def test_max_b_filters_pi_ring_atoms():
    assert core.detect_interactions_in_structure(
        _structure_with_phe_and_og(ring_b=80.0), "test", max_b=50.0) == []


def test_same_residue_donor_is_excluded():
    st = _structure_with_phe_and_og(include_external_donor=False)
    phe = st[0][0][0]
    phe.add_atom(_atom("N", "N", (0.0, 0.0, 3.0)))
    phe.add_atom(_atom("H", "H", (0.0, 0.0, 2.0)))

    assert core.detect_interactions_in_structure(st, "test") == []


def test_cli_system_summary_defaults_to_hudson_plevin_only():
    rows = [
        {"is_hudson": 1, "is_plevin": 0, "is_p_slab": 0},
        {"is_hudson": 0, "is_plevin": 1, "is_p_slab": 1},
        {"is_hudson": 1, "is_plevin": 1, "is_p_slab": 1},
    ]

    summary = cli.summarize_systems(rows)

    assert summary["hudson"] == 2
    assert summary["plevin"] == 2
    assert summary["hudson_plevin_union"] == 3
    assert summary["hudson_plevin"] == 1
    assert "p_slab" not in summary
    assert "all_three" not in summary


def test_cli_system_summary_can_include_p_slab():
    rows = [
        {"is_hudson": 1, "is_plevin": 0, "is_p_slab": 0},
        {"is_hudson": 0, "is_plevin": 1, "is_p_slab": 1},
        {"is_hudson": 1, "is_plevin": 1, "is_p_slab": 1},
    ]

    summary = cli.summarize_systems(rows, include_p_slab=True)

    assert summary["hudson"] == 2
    assert summary["plevin"] == 2
    assert summary["p_slab"] == 2
    assert summary["hudson_plevin_union"] == 3
    assert summary["hudson_plevin"] == 1
    assert summary["all_three"] == 1


def test_default_csv_output_hides_p_slab_columns(tmp_path):
    rows = [{
        "pdb": "test", "resolution": 1.5,
        "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
        "X_chain": "A", "X_res": "SER", "X_id": "2",
        "X_atom": "OG", "H_atom": "HG", "H_source": "experimental",
        "is_hudson": 1, "is_plevin": 1, "is_p_slab": 1,
        "dist_X_centroid": 3.0, "dist_X_Pi": 3.0,
        "proj_dist": 0.0, "theta": 0.0, "angle_XPCN": 0.0,
        "angle_XH_Pi": 180.0, "P_radius": 2.0,
        "P_slab_half_thickness": 0.5, "h_proj_dist": 0.0, "H_ray_t": 2.5,
        "H_xyz_x": 0.0, "H_xyz_y": 0.0, "H_xyz_z": 2.0,
        "pi_normal_x": 0.0, "pi_normal_y": 0.0, "pi_normal_z": 1.0,
        "X_side_of_pi": 1,
        "is_trp_5ring_acceptor": 0, "is_pi_pi_tshaped": 0, "sym_op": 0,
    }]

    out_path = tmp_path / "hits.csv"
    with output.ResultStreamer(out_path, "csv", verbose=False) as streamer:
        streamer.write_chunk(rows)

    with out_path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))

    assert "is_hudson" in header
    assert "is_plevin" in header
    assert "is_p_slab" not in header
    assert "h_proj_dist" not in header
    assert "H_ray_t" not in header
    assert "H_xyz_x" not in header
    assert "pi_normal_z" not in header
    assert "X_side_of_pi" not in header


def test_csv_output_can_include_coordinate_columns(tmp_path):
    rows = [{
        "pdb": "test", "resolution": 1.5,
        "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
        "X_chain": "A", "X_res": "SER", "X_id": "2",
        "X_atom": "OG", "H_atom": "HG", "H_source": "experimental",
        "is_hudson": 1, "is_plevin": 1,
        "dist_X_centroid": 3.0, "dist_X_Pi": 3.0,
        "proj_dist": 0.0, "theta": 0.0, "angle_XPCN": 0.0,
        "angle_XH_Pi": 180.0,
        "pi_center_x": 0.0, "pi_center_y": 0.0, "pi_center_z": 0.0,
        "pi_normal_x": 0.0, "pi_normal_y": 0.0, "pi_normal_z": 1.0,
        "X_xyz_x": 0.0, "X_xyz_y": 0.0, "X_xyz_z": 3.0,
        "H_xyz_x": 0.0, "H_xyz_y": 0.0, "H_xyz_z": 2.0,
        "X_side_of_pi": 1,
        "is_trp_5ring_acceptor": 0, "is_pi_pi_tshaped": 0, "sym_op": 0,
    }]

    out_path = tmp_path / "coords.csv"
    with output.ResultStreamer(
            out_path, "csv", verbose=False, include_coordinates=True) as streamer:
        streamer.write_chunk(rows)

    with out_path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))

    for col in (
        "pi_center_x", "pi_center_y", "pi_center_z",
        "pi_normal_x", "pi_normal_y", "pi_normal_z",
        "X_xyz_x", "X_xyz_y", "X_xyz_z",
        "H_xyz_x", "H_xyz_y", "H_xyz_z",
        "X_side_of_pi",
    ):
        assert col in header


def test_csv_output_can_include_p_slab_columns(tmp_path):
    rows = [{
        "pdb": "test", "resolution": 1.5,
        "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
        "X_chain": "A", "X_res": "SER", "X_id": "2",
        "X_atom": "OG", "H_atom": "HG", "H_source": "experimental",
        "is_hudson": 1, "is_plevin": 1, "is_p_slab": 1,
        "dist_X_centroid": 3.0, "dist_X_Pi": 3.0,
        "proj_dist": 0.0, "theta": 0.0, "angle_XPCN": 0.0,
        "angle_XH_Pi": 180.0, "P_radius": 2.0,
        "P_slab_half_thickness": 0.5, "h_proj_dist": 0.0, "H_ray_t": 2.5,
        "is_trp_5ring_acceptor": 0, "is_pi_pi_tshaped": 0, "sym_op": 0,
    }]

    out_path = tmp_path / "hits.csv"
    with output.ResultStreamer(out_path, "csv", verbose=False, include_p_slab=True) as streamer:
        streamer.write_chunk(rows)

    with out_path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))

    assert "is_p_slab" in header
    assert "h_proj_dist" in header
    assert "H_ray_t" in header


def test_candidate_csv_output_includes_direction_geometry_without_p_slab_label(tmp_path):
    rows = [{
        "pdb": "test", "resolution": 1.5,
        "pi_chain": "A", "pi_res": "PHE", "pi_id": "1",
        "X_chain": "A", "X_res": "SER", "X_id": "2",
        "X_atom": "OG", "H_atom": "HG", "H_source": "experimental",
        "is_xh_candidate": 1, "is_hudson": 0, "is_plevin": 0, "is_p_slab": 0,
        "is_hudson_spatial": 1, "is_plevin_spatial": 1,
        "hudson_dist_ok": 1, "hudson_proj_ok": 1, "hudson_direction_ok": 0,
        "plevin_dist_ok": 1, "plevin_xpcn_ok": 1, "plevin_direction_ok": 0,
        "dist_X_centroid": 3.0, "dist_X_Pi": 3.0,
        "proj_dist": 0.0, "theta": None, "angle_XPCN": 0.0,
        "angle_XH_Pi": 0.0, "P_radius": 2.0,
        "P_slab_half_thickness": 0.5, "h_proj_dist": None, "H_ray_t": None,
        "H_ray_entry_dist": None, "h_plane_proj_dist": None, "H_plane_t": None,
        "H_plane_entry_dist": None, "xh_centroid_cos": -1.0,
        "xh_lateral_inward_score": None, "delta_h_proj_dist": None,
        "is_trp_5ring_acceptor": 0, "is_pi_pi_tshaped": 0, "sym_op": 0,
    }]

    out_path = tmp_path / "candidates.csv"
    with output.ResultStreamer(
            out_path, "csv", verbose=False, include_xh_candidates=True) as streamer:
        streamer.write_chunk(rows)

    with out_path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))

    assert "is_xh_candidate" in header
    assert "is_hudson_spatial" in header
    assert "hudson_direction_ok" in header
    assert "h_proj_dist" in header
    assert "H_plane_entry_dist" in header
    assert "xh_centroid_cos" in header
    assert "is_p_slab" not in header
