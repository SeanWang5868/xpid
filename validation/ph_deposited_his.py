#!/usr/bin/env python3
"""Compare a saved pH sweep with deposited neutron HIS H/D records.

This does not inspect density, refine occupancy, correct pH to pD, or treat
missing H/D atoms as evidence of deprotonation. It is a discrepancy audit,
not an accuracy calculation. Run after ph_sweep.py on its output directory.
"""
import argparse
import json
from pathlib import Path

import gemmi

from xpid.monomer_bonds import canonical_hydrogen_name
from xpid.protonation_pka import alt


def inspect(root):
    summary = json.loads((root / "summary.json").read_text())
    observations = []
    for code, entry in summary["structures"].items():
        if not any("NEUTRON" in method.upper() for method in entry["methods"]):
            continue
        structure = gemmi.read_structure(str(root / code / f"{code}.cif"))
        runs = {}
        for label in entry["runs"]:
            if label == "default":
                continue
            diagnostic = json.loads((root / code / f"{label}_diagnostics.json").read_text())["protonation"]
            for site in diagnostic["sites"]:
                if site["site_type"] == "sidechain" and site["key"][3] == "HIS":
                    runs.setdefault(tuple(site["key"]), {})[label] = {
                        "pH": diagnostic["pH"], "pH_source": diagnostic["pH_source"],
                        **{k: site[k] for k in ("state", "status", "charge", "pka", "reason")}}
        for key, predictions in runs.items():
            mi, chain, resid, resname, label = key
            residue = next(r for c in structure[mi] if c.name == chain for r in c
                           if str(r.seqid).strip() == resid and r.name == resname)
            hydrogens = []
            for parent_name, hydrogen_name in (("ND1", "HD1"), ("NE2", "HE2")):
                parents = [a for a in residue if a.name == parent_name and alt(a) in {"", label}]
                if len(parents) != 1:
                    continue
                parent = parents[0]
                for atom in residue:
                    if (not atom.is_hydrogen() or canonical_hydrogen_name(atom.name) != hydrogen_name
                            or atom.occ < 0.1):
                        continue
                    # H/D alternate labels can encode isotope occupancy; retain
                    # the raw label rather than pretending it is a heavy-atom conformer.
                    distance = parent.pos.dist(atom.pos)
                    if 0.5 < distance <= 1.25:
                        hydrogens.append({"parent": parent_name, "name": atom.name,
                                          "element": atom.element.name, "altloc": alt(atom),
                                          "occupancy": atom.occ, "parent_occupancy": parent.occ,
                                          "distance": distance})
            observations.append({"pdb": code, "site": key, "deposited_hydrogens": hydrogens,
                                 "predictions": predictions})
    result = {"interpretation": "Deposited coordinate evidence only; no density validation or accuracy percentage",
              "observations": observations}
    output = root / "deposited_his_audit.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"{len(observations)} HIS sites: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    inspect(parser.parse_args().directory)
