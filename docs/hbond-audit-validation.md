# H-bond audit implementation validation

Local check, 20 October 2026, XPID 2.1.12. This is a software-invariance and
chemical-guard audit, not an experimental validation of H-bond occupation.

Final suite: **268 passed in 68.96 s**, including 38 audit-specific cases.
The optional PROPKA 3.5.1 tests were enabled using an isolated dependency
directory; no base-environment installation was required. The 2.1.12 wheel
was built successfully and checked to contain the new audit module.

Eight previously downloaded RCSB structures were processed through the
coordinate-inclusive API with audit off/on. Every returned row was identical
in both runs and also matched its independently saved pre-audit baseline from
the pH validation. No tolerance or ignored-column comparison was used.

| Input | Unchanged hits | Audited Cone pairs |
| --- | ---: | ---: |
| 1UBQ | 2 | 1 |
| 1LYZ | 21 | 5 |
| 1CA2 | 54 | 20 |
| 2ZOI | 23 | 11 |
| 3FHP | 14 | 6 |
| 7BBC | 129 | 57 |
| 1CRN | 2 | 1 |
| 5FJJ | 565 | 190 |

The 291 audited donor-conformer/ring pairs include 14 negative pairs. This is
not a count of all donors or all interactions: rigid donors and pairs removed
by earlier filters are outside the audit. No diagnostic truncation occurred.

Artifacts and input hashes are in the ignored local directory
`_test_hbond_audit/final_eight/`; the reproducible runner is
`validation/hbond_audit_check.py`. Runtime measurements are single ordered
runs with shared caches and must not be interpreted as a controlled benchmark.

Additional checks include pH 7.4 with PROPKA (2ZOI: 23 unchanged hits;
7BBC: 126 unchanged hits), CLI JSON/Parquet diagnostics,
multiple models, P-slab, nearest periodic images, altloc/occupancy/B filters,
metal-associated water, the S–H search-radius bound, and explicit diagnostic
truncation without detection truncation. The audit's contact counts are model
descriptions, not measured occupancy or probabilities.
The two-input, two-worker Parquet test produced 152 identical rows and
identical schemas with audit off/on; both worker audit records were saved.

Next scientific step: inspect selected neutron H/D sites and their density
support against the recorded compatible/disjoint candidate directions.
No new H-bond-based filtering policy is justified by these regression tests.
