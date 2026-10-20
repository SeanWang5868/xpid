#!/usr/bin/env python3
"""Compare the audit with the unmodified detection path on local structures."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from xpid import api
from xpid._version import __version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, required=True, help="New directory")
    parser.add_argument("--ph", type=float)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    summary = {"xpid": __version__, "pH": args.ph, "structures": [],
               "interpretation": "Software invariance, not evidence of physical H-bond occupation"}
    for i, path in enumerate(args.inputs):
        baseline_start = time.monotonic()
        baseline = api.detect(path, ph=args.ph, include_coordinates=True)
        baseline_seconds = time.monotonic() - baseline_start
        d = {}
        audit_start = time.monotonic()
        audited = api.detect(path, ph=args.ph, include_coordinates=True, hbond_audit=True, diagnostics=d)
        seconds = time.monotonic() - audit_start
        assert baseline == audited, f"Audit changed hits: {path}"
        # Earlier pH-validation files provide an independently saved pre-audit
        # baseline with the same API defaults and coordinate export.
        prior = path.parent / "default_hits.json"
        old_match = None
        if args.ph is None and prior.exists():
            old_match = json.loads(prior.read_text()) == audited
            assert old_match, f"Pre-audit baseline changed: {path}"
        report = d["hbond_audit"]
        record = {"input": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "hits": len(audited), "identical": True, "pre_audit_baseline_identical": old_match,
                  "baseline_seconds": round(baseline_seconds, 3), "audit_seconds": round(seconds, 3),
                  "pairs": len(report["pairs"]), "pairs_omitted": report["pairs_omitted"],
                  "negative_pairs": sum(not p["detected"] for p in report["pairs"]),
                  "polar_disjoint_pairs": sum(p["typed_contacts"]["disjoint_only"]
                                               and p["donor_class"] == "polar_XH" for p in report["pairs"]),
                  "uncertain_contact_pairs": sum(p["uncertain_contacts"]["hbond_conformers"] > 0
                                                 for p in report["pairs"])}
        summary["structures"].append(record)
        (args.out / f"{i}_{path.stem}_audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        (args.out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        print(f"{path.name}: {len(audited)} unchanged hits; {len(report['pairs'])} audit pairs", flush=True)


if __name__ == "__main__":
    main()
