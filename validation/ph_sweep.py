#!/usr/bin/env python3
"""Real-structure pH sensitivity audit, NOT a protonation-accuracy benchmark.

Run with XPID and its [ph] extra installed. Inputs are fetched from RCSB;
every unmodified API result and diagnostic is saved alongside its input hash.
The audit hook observes preparation without replacing its predictions/geometry.
No production algorithm is changed, cached, or mocked by this script.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import logging
from pathlib import Path
import platform
import re
import time
from unittest.mock import patch
from urllib.request import urlopen

import gemmi

from xpid import api, monlib, protonation, protonation_pka
from xpid._version import __version__


IDENTITY = ("model", "X_chain", "X_res", "X_id", "X_atom", "X_altloc",
            "H_atom", "H_altloc", "pi_chain", "pi_res", "pi_id",
            "pi_ring_id", "pi_altloc", "symmetry_code")
GEOMETRY = ("is_hudson", "is_plevin", "dist_X_Pi", "theta", "angle_XPCN",
            "angle_XH_Pi", "H_xyz_x", "H_xyz_y", "H_xyz_z")


def save(path, value):
    # Detector diagnostics contain sets (e.g. quarantined residue keys).
    def encode(obj):
        if isinstance(obj, set):
            return sorted(obj)
        raise TypeError(f"Unsupported diagnostic type: {type(obj)}")
    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=encode) + "\n")


def heavy_snapshot(structure):
    return [(mi, c.name, str(r.seqid), r.name, a.name, a.altloc,
             a.pos.tolist(), a.occ, a.b_iso, a.charge)
            for mi, m in enumerate(structure) for c in m for r in c for a in r
            if not a.is_hydrogen()]


def audit_context(structure, context):
    checks = Counter()
    for site in context.sites:
        mi, chain_name, resid, resname, label = site.key
        residue = next(r for c in structure[mi] if c.name == chain_name
                       for r in c if str(r.seqid).strip() == resid and r.name == resname)
        hs = {a.name: a for a in residue if a.is_hydrogen()
              and protonation_pka.alt(a) == label}
        for parent, names in site.donor_h.items():
            for name in names:
                assert name in hs and hs[name].occ >= 0.1, (site.key, name)
                checks["assigned_h_present_and_occupied"] += 1
        if site.site_type != "sidechain":
            continue
        all_names = {h for names in protonation.SITE_H[resname].values() for h in names}
        assigned = {h for names in site.donor_h.values() for h in names}
        assert all_names.intersection(hs) == assigned, (site.key, "unexpected site H")
        checks["site_h_matches_assignment"] += 1
        if site.status == "resolved" and resname == "HIS":
            expected = 2 if site.state == "doubly_protonated" else 1
            assert len(assigned) == expected, site.key
            checks["his_h_count"] += 1
        if resname in {"ASP", "GLU", "CYS", "TYR"} and site.state == "deprotonated":
            assert not assigned, site.key
            checks["deprotonated_acid_has_no_site_h"] += 1
        if site.reason in {"covalent_link_or_coordination", "possible_covalent_link",
                           "metal_coordination_requires_specialist_assignment"}:
            assert site.status == "unresolved" and not assigned, site.key
            checks["linked_or_metal_site_not_forced"] += 1
    return dict(checks)


def run(path, ph):
    diagnostics, checks = {}, {}
    original = protonation.prepare

    def observed_prepare(structure, *args, **kwargs):
        before = heavy_snapshot(structure)
        context = original(structure, *args, **kwargs)
        assert before == heavy_snapshot(structure), "pH preparation moved a heavy atom"
        checks.update(audit_context(structure, context))
        checks["heavy_atoms_unchanged"] = True
        return context

    started = time.monotonic()
    with patch.object(protonation, "prepare", observed_prepare):
        rows = api.detect(path, ph=ph, include_coordinates=True, diagnostics=diagnostics)
    for row in rows:
        assert not (row["X_res"] == "LYS" and row["X_atom"] == "NZ")
        assert not (row["X_res"] == "ARG" and row["X_atom"] in {"NE", "NH1", "NH2"})
        if ph is not None:
            assert row.get("pi_protonation_state") != "doubly_protonated"
            if row["X_atom"] in protonation.SITE_H.get(row["X_res"], {}):
                assert row["X_protonation_state"] not in {"unresolved", "deprotonated", "doubly_protonated"}
    checks["reported_state_consistency"] = True
    return rows, diagnostics, checks, round(time.monotonic() - started, 3)


def identity(row):
    return tuple(row.get(field) for field in IDENTITY)


def compare(before, after):
    left = {identity(r): r for r in before}
    right = {identity(r): r for r in after}
    assert len(left) == len(before) and len(right) == len(after), "duplicate interaction identity"
    return {"added": [right[k] for k in sorted(right.keys() - left.keys())],
            "removed": [left[k] for k in sorted(left.keys() - right.keys())],
            "shared_geometry_changed": [
                {"before": left[k], "after": right[k]} for k in sorted(left.keys() & right.keys())
                if any(left[k].get(f) != right[k].get(f) for f in GEOMETRY)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="New output directory (must not exist)")
    parser.add_argument("--ids", nargs="+", default=["1crn", "1ubq", "1lyz", "1ca2", "2zoi", "3fhp", "7bbc", "5fjj"])
    parser.add_argument("--ph", type=float, nargs="+", default=[4.0, 7.4, 10.0])
    args = parser.parse_args()
    if any(not re.fullmatch(r"[A-Za-z0-9]{4}", code) for code in args.ids):
        parser.error("Expected four-character PDB IDs")
    if len(set(args.ph)) != len(args.ph) or len(set(args.ids)) != len(args.ids):
        parser.error("Duplicate pH values or IDs")
    for ph in args.ph:
        protonation_pka.numeric_ph(ph)
    args.out.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(filename=args.out / "audit.log", level=logging.WARNING)
    root = Path(monlib.get_monomer_library_path())
    summary = {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "xpid": __version__,
               "python": platform.python_version(), "gemmi": gemmi.__version__,
               "propka": version("propka"), "monomer_library": str(root),
               "monomer_hashes": {}, "structures": {},
               "interpretation": "Fixed-coordinate sensitivity and software consistency, not experimental accuracy"}
    for code in protonation.SITE_H:
        cif = monlib.find_monomer_cif(code)
        if cif:
            summary["monomer_hashes"][code] = hashlib.sha256(cif.read_bytes()).hexdigest()
    save(args.out / "summary.json", summary)
    for code in args.ids:
        code = code.lower()
        print(f"{code}: downloading and testing", flush=True)
        folder = args.out / code
        folder.mkdir()
        url = f"https://files.rcsb.org/download/{code.upper()}.cif"
        with urlopen(url, timeout=120) as response:
            data = response.read()
        path = folder / f"{code}.cif"
        path.write_bytes(data)
        block = gemmi.cif.read(str(path)).sole_block()
        record = {"url": url, "sha256": hashlib.sha256(data).hexdigest(),
                  "title": block.find_value("_struct.title"),
                  "methods": list(block.find_values("_exptl.method")), "runs": {}}
        try:
            record["experimental_ph"], record["experimental_ph_source"] = protonation_pka.resolve_ph(path, "auto")
        except ValueError as exc:
            record["experimental_ph_error"] = str(exc)
        summary["structures"][code] = record
        results, sites_by_ph = {}, {}
        for ph in [None] + sorted(args.ph) + (["auto"] if "experimental_ph" in record else []):
            label = "default" if ph is None else f"ph_{ph}"
            rows, diagnostics, checks, seconds = run(path, ph)
            results[label] = rows
            save(folder / f"{label}_hits.json", rows)
            save(folder / f"{label}_diagnostics.json", diagnostics)
            sites = diagnostics.get("protonation", {}).get("sites", [])
            sites_by_ph[label] = sites
            record["runs"][label] = {
                "hits": len(rows), "seconds": seconds, "checks": checks,
                "resolved": sum(s["status"] == "resolved" for s in sites),
                "unresolved": sum(s["status"] != "resolved" for s in sites),
                "his_states": dict(Counter(s["state"] for s in sites if s["key"][3] == "HIS" and s["site_type"] == "sidechain")),
                "unresolved_reasons": dict(Counter(s["reason"] for s in sites if s["status"] != "resolved")),
                "warnings": len(diagnostics.get("protonation", {}).get("warnings", []))}
            if ph is not None:
                delta = compare(results["default"], rows)
                save(folder / f"{label}_vs_default.json", delta)
                record["runs"][label]["vs_default"] = {key: len(value) for key, value in delta.items()}
            print(f"  {label}: {len(rows)} hits ({seconds}s)", flush=True)
            save(args.out / "summary.json", summary)
        # Independent runs of the same fixed structure must yield the same pKas,
        # and the model's protonated fraction must not increase with pH.
        previous = {}
        for ph in sorted(args.ph):
            for site in sites_by_ph[f"ph_{ph}"]:
                key = (tuple(site["key"]), site["site_type"])
                if key in previous:
                    old = previous[key]
                    assert old["pka"] == site["pka"], (code, key, "unstable pKa")
                    if site["protonated_fraction"] is not None:
                        assert site["protonated_fraction"] <= old["protonated_fraction"], (code, key)
                previous[key] = site
        record["ph_fraction_monotonic"] = True
        for first, second in zip(sorted(args.ph), sorted(args.ph)[1:]):
            delta = compare(results[f"ph_{first}"], results[f"ph_{second}"])
            save(folder / f"ph_{first}_to_{second}.json", delta)
        # A fresh complete API call checks reproducibility, not cached predictions.
        repeated, _, _, _ = run(path, args.ph[0])
        assert repeated == results[f"ph_{args.ph[0]}"], (code, "non-deterministic")
        record["repeat_identical"] = True
        save(args.out / "summary.json", summary)
    print(f"Complete: {args.out / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
