import csv
import inspect
from pathlib import Path

import gemmi
import pytest
from xpid import (
    XPIDError, __version__, api, cli, core, config, detect,
    hydrogen_prep as prep,
    monomer_bonds, output, provenance, resolver,
)


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
        st, "test", cone_mode="none", filter_donor_atom=["O"])
    assert any(hit["X_atom"] == "OG" for hit in oxygen_hits)

    atom_name_hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", cone_mode="none", filter_donor_atom=["OG"])
    assert any(hit["X_atom"] == "OG" for hit in atom_name_hits)

    nitrogen_hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", cone_mode="none", filter_donor_atom=["N"])
    assert nitrogen_hits == []


def test_api_and_core_default_disable_cone_and_p_slab():
    api_sig = inspect.signature(detect)
    core_sig = inspect.signature(core.detect_interactions_in_structure)

    assert api_sig.parameters["cone_mode"].default == "auto"
    assert api_sig.parameters["include_p_slab"].default is False
    assert api_sig.parameters["report_xh_candidates"].default is False
    assert api_sig.parameters["include_coordinates"].default is False
    assert api_sig.parameters["residue_pair"].default is None
    assert core_sig.parameters["cone_mode"].default == "auto"
    assert core_sig.parameters["include_p_slab"].default is False
    assert core_sig.parameters["report_xh_candidates"].default is False
    assert core_sig.parameters["include_coordinates"].default is False
    assert core_sig.parameters["residue_pair"].default is None


def test_api_distinguishes_failure_from_empty_result(monkeypatch, tmp_path):
    monkeypatch.setattr(
        api.structure_io, "read_structure",
        lambda path: (_ for _ in ()).throw(ValueError("bad structure")),
    )

    with pytest.raises(XPIDError, match="bad structure"):
        detect(tmp_path / "bad.cif")
    assert detect(tmp_path / "bad.cif", on_error="empty") == []


def test_provenance_records_real_cone_mode_and_paths(tmp_path):
    args = cli._build_parser().parse_args(["dummy.cif"])
    output_path = tmp_path / "results.csv"

    metadata = provenance.build_metadata(
        args, output_path, "/monomers", file_count=3)

    assert metadata["version"] == __version__
    assert metadata["output"] == str(output_path.resolve())
    assert metadata["monomer_library"] == "/monomers"
    assert metadata["parameters"]["cone_mode"] == "auto"
    args.no_cone = True
    metadata = provenance.build_metadata(
        args, output_path, "/monomers", file_count=3)
    assert metadata["parameters"]["cone_mode"] == "none"


def test_provenance_records_input_resolution_counts(tmp_path):
    args = cli._build_parser().parse_args(["dummy.cif"])
    counts = {
        "pdb_list_entries": 12,
        "pdb_redo_mirror": 5,
        "standard_pdb_mirror": 6,
        "missing": 1,
    }

    metadata = provenance.build_metadata(
        args, tmp_path / "results.csv", "/monomers", file_count=11,
        input_resolution=counts)

    assert metadata["input_resolution"] == counts


def test_cli_parser_default_cone_auto_and_p_slab_off():
    parser = cli._build_parser()
    args = parser.parse_args(["dummy.cif"])

    assert args.no_cone is False
    assert args.include_p_slab is False
    assert args.report_xh_candidates is False
    assert args.include_coordinates is False
    assert args.residue_pair is None


def test_cli_version_reports_installed_package_version(capsys):
    parser = cli._build_parser()
    parser.prog = "xpid"

    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == f"xpid {__version__}\n"


def test_cli_short_v_reports_installed_package_version(capsys):
    parser = cli._build_parser()
    parser.prog = "xpid"

    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["-v"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == f"xpid {__version__}\n"


def test_cli_long_verbose_enables_detailed_output():
    args = cli._build_parser().parse_args(["--verbose", "dummy.cif"])

    assert args.verbose is True


def test_cli_parser_accepts_residue_pair_selectors():
    parser = cli._build_parser()
    args = parser.parse_args(["dummy.cif", "--residue-pair", "//A/12", "//A/18"])

    assert args.residue_pair == ["//A/12", "//A/18"]


def test_cli_sections_use_consistent_plain_text_formatting():
    rendered = cli._format_sections("Xpid initialization", [
        ("Input", [
            ("Unique targets", 3),
            ("Missing", 0),
        ]),
        ("Detection", [
            ("Cone", "auto (rotatable groups)"),
            ("P-slab", cli._on_off(False)),
        ]),
    ])

    assert rendered == (
        "Xpid initialization\n"
        "\n"
        "Input\n"
        "  Unique targets : 3\n"
        "  Missing        : 0\n"
        "\n"
        "Detection\n"
        "  Cone           : auto (rotatable groups)\n"
        "  P-slab         : off"
    )
    assert "[INFO]" not in rendered
    assert "[SUMMARY]" not in rendered


def test_cli_boolean_status_has_only_on_and_off():
    assert cli._on_off(True) == "on"
    assert cli._on_off(False) == "off"


def test_jobs_one_uses_direct_worker_without_pool(monkeypatch):
    tasks = [object(), object()]
    monkeypatch.setattr(
        cli, "process_one_file", lambda task: ("result", task))
    monkeypatch.setattr(
        cli.multiprocessing, "Pool",
        lambda *args, **kwargs: pytest.fail("Pool should not be created"),
    )

    assert list(cli.iter_task_results(tasks, jobs=1)) == [
        ("result", tasks[0]),
        ("result", tasks[1]),
    ]


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


def test_pdb_list_resolution_reports_sources_and_prioritizes_redo(tmp_path):
    pdb_root = tmp_path / "pdb"
    redo_root = tmp_path / "redo"
    pdb_list = tmp_path / "ids.txt"
    pdb_list.write_text("1abc\n2def\n3ghi\n", encoding="utf-8")

    redo_1abc = redo_root / "ab" / "1abc" / "1abc_final.cif"
    pdb_1abc = pdb_root / "ab" / "1abc.cif.gz"
    pdb_2def = pdb_root / "de" / "2def.cif.gz"
    for path in (redo_1abc, pdb_1abc, pdb_2def):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"")

    resolved = resolver.resolve_inputs(
        [], str(pdb_list), str(pdb_root), str(redo_root))

    assert resolved.files == sorted(
        [redo_1abc.resolve(), pdb_2def.resolve()])
    assert resolved.provenance_counts() == {
        "pdb_list_entries": 3,
        "pdb_redo_mirror": 1,
        "standard_pdb_mirror": 1,
        "missing": 1,
    }
    assert resolved.missing_codes == ("3ghi",)
    assert resolved.identity_for(redo_1abc) == ("1abc", "pdb_redo")
    assert resolved.identity_for(pdb_2def) == ("2def", "pdb_mirror")


def test_cation_pi_donors_are_excluded_from_core_detection():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    lys = gemmi.Residue()
    lys.name = "LYS"
    lys.seqid = _seqid(2)
    lys.add_atom(_atom("NZ", "N", (0.0, 0.0, 3.0)))
    lys.add_atom(_atom("HZ1", "H", (0.0, 0.0, 2.0)))
    chain.add_residue(lys)

    assert core.detect_interactions_in_structure(st, "test", cone_mode="none") == []


def test_residue_pair_filter_keeps_only_selected_inter_residue_contacts():
    st = _structure_with_phe_and_og()
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(3)
    ser.add_atom(_atom("OG", "O", (0.5, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0.5, 0.0, 2.0)))
    chain.add_residue(ser)

    all_hits = core.detect_interactions_in_structure(st, "test", cone_mode="none")
    assert {hit["X_id"].strip() for hit in all_hits} == {"2", "3"}

    pair_hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", cone_mode="none",
        residue_pair=("//A/1", "//A/2"))
    assert len(pair_hits) == 1
    assert pair_hits[0]["pi_id"].strip() == "1"
    assert pair_hits[0]["X_id"].strip() == "2"


def test_residue_pair_filter_is_direction_agnostic():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", cone_mode="none",
        residue_pair=("//A/2", "//A/1"))

    assert len(hits) == 1
    assert hits[0]["pi_id"].strip() == "1"
    assert hits[0]["X_id"].strip() == "2"


def test_residue_pair_filter_returns_empty_for_unmatched_pair():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", cone_mode="none",
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
        st, "test", cone_mode="none", residue_pair=("//A/1", "//A/2"))

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
        monomer_bonds,
        "get_bonded_hydrogen_names",
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

    hits = core.detect_interactions_in_structure(st, "test", cone_mode="none")

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

    hits = core.detect_interactions_in_structure(st, "test", cone_mode="none")

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
        st, "test", cone_mode="none", include_p_slab=True)

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
        st, "test", cone_mode="none", include_p_slab=True) == []


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
        st, "test", cone_mode="none", include_p_slab=True) == []


def test_p_model_is_disabled_by_default():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (1.9, 0.0, 4.0)))
    ser.add_atom(_atom("HG", "H", (1.9, 0.0, 3.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(st, "test", cone_mode="none") == []


def test_xh_candidate_mode_reports_direction_failed_spatial_candidate():
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("HG", "H", (0.0, 0.0, 4.0)))
    chain.add_residue(ser)

    assert core.detect_interactions_in_structure(st, "test", cone_mode="none") == []

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", report_xh_candidates=True)

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
        _structure_with_phe_and_og(), "test", cone_mode="none", report_xh_candidates=True)

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
        st, "test", cone_mode="auto", report_xh_candidates=True) == []


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
        st, "test", cone_mode="none", include_p_slab=True)

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

    hits = core.detect_interactions_in_structure(st, "test", cone_mode="none")

    assert len(hits) == 1
    assert hits[0]["is_hudson"] == 1
    assert hits[0]["is_plevin"] == 0


def test_coordinate_output_includes_h_xyz_pi_normal_and_x_side():
    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", cone_mode="none",
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


def test_hydrogen_modes_are_applied_to_every_model_and_synchronized(monkeypatch):
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    for model_number in ("1", "2"):
        model = gemmi.Model(model_number)
        chain = gemmi.Chain("A")
        ser = gemmi.Residue()
        ser.name = "SER"
        ser.seqid = _seqid(1)
        ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
        ser.add_atom(_atom("HG", "H", (9.0, 9.0, 9.0)))
        chain.add_residue(ser)
        model.add_chain(chain)
        st.add_model(model)

    monkeypatch.setattr(prep, "_get_shared_monlib", lambda codes: gemmi.MonLib())
    calls = []

    def fake_prepare(working, monlib, h_change, model_index):
        calls.append((h_change, model_index))
        residue = working[model_index][0][0]
        for index in reversed(range(len(residue))):
            if residue[index].element.name.upper() in {"H", "D"}:
                del residue[index]
        residue.add_atom(_atom("HX", "H", (float(model_index), 1.0, 1.0)))
        return prep.TopologyPreparationReport(
            model_index=model_index, status="success", attempts=1)

    monkeypatch.setattr(prep, "_prepare_topology_with_retries", fake_prepare)

    prep.add_hydrogens_memory(st, h_change_val=3)

    assert calls == [
        (gemmi.HydrogenChange.ReAdd, 0),
        (gemmi.HydrogenChange.ReAdd, 1),
    ]
    for model_index in range(2):
        hydrogens = [
            atom for atom in st[model_index][0][0]
            if atom.element.name.upper() == "H"
        ]
        assert len(hydrogens) == 1
        assert hydrogens[0].name == "HX"
        assert hydrogens[0].pos.x == float(model_index)


def test_hydrogen_topology_failure_is_not_silently_accepted(monkeypatch):
    st = _structure_with_phe_and_og()
    original_h_count = sum(
        atom.element.name.upper() in {"H", "D"}
        for model in st for chain in model for residue in chain
        for atom in residue)
    monkeypatch.setattr(prep, "_get_shared_monlib", lambda codes: gemmi.MonLib())
    monkeypatch.setattr(
        prep, "_prepare_topology_with_retries",
        lambda working, monlib, h_change, model_index:
            prep.TopologyPreparationReport(
                model_index=model_index, status="failed", attempts=1,
                message="synthetic topology failure"),
    )

    with pytest.raises(
            prep.HydrogenPreparationError,
            match="synthetic topology failure") as error:
        prep.prepare_hydrogens_memory(st, h_change_val=4)

    assert error.value.report.status == "failed"
    assert sum(
        atom.element.name.upper() in {"H", "D"}
        for model in st for chain in model for residue in chain
        for atom in residue) == original_h_count


def test_hydrogen_report_marks_missing_monomers_as_partial(monkeypatch):
    st = _structure_with_phe_and_og()
    monkeypatch.setattr(prep, "_get_shared_monlib", lambda codes: gemmi.MonLib())
    monkeypatch.setattr(
        prep, "_prepare_topology_with_retries",
        lambda working, monlib, h_change, model_index:
            prep.TopologyPreparationReport(
                model_index=model_index, status="success", attempts=1),
    )

    result = prep.prepare_hydrogens_memory(st, h_change_val=4)

    assert result.report.status == "partial"
    assert result.report.missing_monomer_components == ("PHE", "SER")


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        (1, gemmi.HydrogenChange.Shift),
        (2, gemmi.HydrogenChange.Remove),
        (3, gemmi.HydrogenChange.ReAdd),
        (4, gemmi.HydrogenChange.ReAddButWater),
        (5, gemmi.HydrogenChange.ReAddKnown),
    ],
)
def test_hydrogen_mode_numbers_match_gemmi(monkeypatch, mode, expected):
    st = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "SER"
    residue.seqid = _seqid(1)
    residue.add_atom(_atom("OG", "O", (0.0, 0.0, 0.0)))
    chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)

    monkeypatch.setattr(prep, "_get_shared_monlib", lambda codes: gemmi.MonLib())
    calls = []
    monkeypatch.setattr(
        prep,
        "_prepare_topology_with_retries",
        lambda working, monlib, h_change, model_index:
            (calls.append((h_change, model_index)) or
             prep.TopologyPreparationReport(
                 model_index=model_index, status="success", attempts=1)),
    )

    prep.add_hydrogens_memory(st, h_change_val=mode)

    assert calls == [(expected, 0)]


def test_dictionary_hydrogen_matching_normalizes_neutron_d_names(monkeypatch):
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ser = gemmi.Residue()
    ser.name = "SER"
    ser.seqid = _seqid(2)
    ser.add_atom(_atom("OG", "O", (0.0, 0.0, 3.0)))
    ser.add_atom(_atom("DG", "D", (0.0, 0.0, 2.0)))
    chain.add_residue(ser)
    monkeypatch.setattr(
        monomer_bonds, "get_bonded_hydrogen_names",
        lambda residue, atom: {"HG"} if (residue, atom) == ("SER", "OG") else set(),
    )

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", annotate_cooperativity=False)

    assert len(hits) == 1
    assert hits[0]["H_atom"] == "DG"


def test_explicit_cys_deuterium_uses_sulfur_specific_bond_cutoff(monkeypatch):
    st = _structure_with_phe_and_og(include_external_donor=False)
    cys = gemmi.Residue()
    cys.name = "CYS"
    cys.seqid = _seqid(2)
    cys.add_atom(_atom("SG", "S", (0.0, 0.0, 3.5)))
    cys.add_atom(_atom("DG", "D", (0.0, 0.0, 2.162)))
    st[0][0].add_residue(cys)
    monkeypatch.setattr(
        monomer_bonds, "get_bonded_hydrogen_names",
        lambda residue, atom:
            {"HG"} if (residue, atom) == ("CYS", "SG") else set(),
    )

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", annotate_cooperativity=False)

    assert len(hits) == 1
    assert hits[0]["X_atom"] == "SG"
    assert hits[0]["H_atom"] == "DG"
    assert hits[0]["H_source"] == "experimental"


def test_unknown_ligand_requires_unambiguous_covalent_hydrogen(monkeypatch):
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ligand = gemmi.Residue()
    ligand.name = "ZZZ"
    ligand.seqid = _seqid(2)
    ligand.add_atom(_atom("O1", "O", (0.0, 0.0, 3.0)))
    ligand.add_atom(_atom("N1", "N", (0.2, 0.0, 3.0)))
    ligand.add_atom(_atom("H1", "H", (0.1, 0.0, 2.0)))
    chain.add_residue(ligand)
    monkeypatch.setattr(
        monomer_bonds, "get_bonded_hydrogen_names",
        lambda residue, atom: None)

    assert core.detect_interactions_in_structure(
        st, "test", cone_mode="none") == []


def test_unknown_ligand_accepts_unique_element_appropriate_xh_bond(monkeypatch):
    st = _structure_with_phe_and_og(include_external_donor=False)
    chain = st[0][0]
    ligand = gemmi.Residue()
    ligand.name = "ZZZ"
    ligand.seqid = _seqid(2)
    ligand.add_atom(_atom("O1", "O", (0.0, 0.0, 3.0)))
    ligand.add_atom(_atom("C1", "C", (2.0, 0.0, 3.0)))
    ligand.add_atom(_atom("H1", "H", (0.0, 0.0, 2.0)))
    chain.add_residue(ligand)
    monkeypatch.setattr(
        monomer_bonds, "get_bonded_hydrogen_names",
        lambda residue, atom: None)

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none", annotate_cooperativity=False)

    assert len(hits) == 1
    assert hits[0]["X_atom"] == "O1"


def test_max_b_filters_pi_ring_atoms():
    assert core.detect_interactions_in_structure(
        _structure_with_phe_and_og(ring_b=80.0), "test", max_b=50.0) == []


def test_same_residue_donor_is_excluded():
    st = _structure_with_phe_and_og(include_external_donor=False)
    phe = st[0][0][0]
    phe.add_atom(_atom("N", "N", (0.0, 0.0, 3.0)))
    phe.add_atom(_atom("H", "H", (0.0, 0.0, 2.0)))

    assert core.detect_interactions_in_structure(st, "test") == []


def test_same_residue_symmetry_copy_is_a_valid_contact(monkeypatch):
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(10, 10, 10, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    phe = gemmi.Residue()
    phe.name = "PHE"
    phe.seqid = _seqid(1)
    ring = {
        "CG": (6.4, 5.0, 0.5), "CD1": (5.7, 6.212, 0.5),
        "CE1": (4.3, 6.212, 0.5), "CZ": (3.6, 5.0, 0.5),
        "CE2": (4.3, 3.788, 0.5), "CD2": (5.7, 3.788, 0.5),
    }
    for name, xyz in ring.items():
        phe.add_atom(_atom(name, "C", xyz))
    phe.add_atom(_atom("N", "N", (5.0, 5.0, 8.0)))
    phe.add_atom(_atom("H", "H", (5.0, 5.0, 9.0)))
    chain.add_residue(phe)
    model.add_chain(chain)
    st.add_model(model)
    monkeypatch.setattr(
        monomer_bonds, "get_bonded_hydrogen_names",
        lambda residue, atom:
            {"H"} if (residue, atom) == ("PHE", "N") else set(),
    )

    assert core.detect_interactions_in_structure(
        st.clone(), "test", cone_mode="none",
        sym_contacts=False, annotate_cooperativity=False) == []

    hits = core.detect_interactions_in_structure(
        st, "test", cone_mode="none",
        sym_contacts=True, annotate_cooperativity=False)

    assert len(hits) == 1
    assert hits[0]["sym_op"] == 0  # P1 translation has operation index 0.
    assert hits[0]["symmetry_code"] == "1_554"
    assert hits[0]["H_xyz_z"] == -1.0


def test_cone_missing_parent_is_reported():
    diagnostics = {}

    hits = core.detect_interactions_in_structure(
        _structure_with_phe_and_og(), "test", diagnostics=diagnostics,
        annotate_cooperativity=False)

    assert hits == []
    assert diagnostics["cone_missing_parent_groups"] == {
        "1:A:2:SER:OG:CB"
    }


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


def test_parquet_streaming_uses_declared_types_across_chunks(tmp_path):
    pyarrow = pytest.importorskip("pyarrow")
    pyarrow_parquet = pytest.importorskip("pyarrow.parquet")
    base = {
        "pdb": "first",
        "model": "1",
        "resolution": 1.5,
        "hbond_acceptor_atom": None,
        "hbond_acceptor_res": None,
        "hbond_acceptor_chain": None,
        "hbond_HA_dist": None,
        "hbond_DHA_angle": None,
        "hbond_vs_xhpi_score": None,
    }
    out_path = tmp_path / "streamed.parquet"

    with output.ResultStreamer(
            out_path, "parquet", verbose=True) as streamer:
        streamer.write_chunk([base])
        streamer.write_chunk([{
            **base,
            "pdb": "second",
            "hbond_acceptor_atom": "O",
            "hbond_acceptor_res": "ASP",
            "hbond_acceptor_chain": "A",
            "hbond_HA_dist": 1.9,
        }])

    table = pyarrow_parquet.read_table(out_path)
    assert table.column("hbond_acceptor_atom").to_pylist() == [None, "O"]
    assert table.column("hbond_acceptor_res").to_pylist() == [None, "ASP"]
    assert table.schema.field(
        "hbond_acceptor_atom").type == pyarrow.large_string()


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
