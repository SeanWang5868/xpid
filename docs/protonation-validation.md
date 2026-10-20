# pH extension: implementation validation

Local audit, 20 October 2026. This records software checks, not independent
validation of protonation accuracy. The feature remains opt-in.

Environment: Python 3.13.5, Gemmi 0.7.4, NumPy 2.3.4, PROPKA 3.5.1. PROPKA
was installed in an isolated test directory; the base environment was not changed.

## Default-path regression

The previous code was run from a detached worktree at `74a9041`. Both versions
processed the same local `5fjj.cif.gz` with default Cone, no symmetry contacts,
verbose JSON and coordinate export. Both produced **565 rows**, with identical
row contents and byte-identical output files (including H coordinates and all
annotations).

- Input SHA-256: `fd5ab0ca0109557b4fd47ddc77dfd3338afc482c74a74b50b6ea7216cfeae214`
- Both output SHA-256: `5a502afb2ee052070393c25d378056f185c36d429ea3a63bf5306a2e9252333f`

## Real-structure pH smoke test

```bash
xpid 5fjj.cif.gz --ph 7.4 --jobs 2 --file-type parquet \
  --include-coordinates --out-dir ph_smoke --provenance
```

This completed with **539 rows**, 772 resolved sites and 122 unresolved sites
(including terminal sites). Among 60 HIS sites: 11 delta, 7 epsilon, 8 doubly
protonated, and 34 unresolved. These assignments are model outputs, **not
measured accuracy or demonstrated corrections**. The conditional pH result
should not be treated as a recall comparison with the default result.

PROPKA emitted 266 warning messages, including incomplete local protonation
of some ASN groups. The warnings were retained in diagnostics rather than
silenced; recognisable remapped labels also carry the original input identity.
A direct native-PDB PROPKA check on chain A also produced ASN warnings, so
they cannot simply be attributed to XPID's residue-number mapping. Inspect
these warnings before interpreting predicted site pKas quantitatively.

Additional smoke runs covered reference-pKa CSV output, no-Cone output and
Cone with symmetry contacts. `--ph auto` correctly rejected this particular
5FJJ file because its experimental pH is `?`; it did not silently use 7.
Unit tests cover successful automatic reading of diffraction, EM and NMR pH.

## Test and build checks

Final result: **230 tests passed**, including **47 pH-specific cases**.

The automated suite includes dedicated checks for HIS H counts and occupancy,
charge/tautomer separation, donor and H-bond acceptor consistency, coupled-site
handling, zero-occupancy filtering, acidic OH detection, nuclear X–H distances,
deprotonated Cone suppression, covalent/metal protection, altloc and shared-X
identity, invalid dictionary geometry, overrides, output schemas and API/CLI
agreement. Run it with the optional pH dependency installed:

```bash
python -m pytest -q
python -m pip wheel --no-deps --no-build-isolation . --wheel-dir wheels
```

The package version has been corrected to **2.1.12** at the user's request.
`xpid --version` reports 2.1.12; a 2.1.12 wheel was built successfully. The full
suite was rerun: **230 passed in 67.51 s**. No GitHub/PyPI/Zenodo release was
published. Published citation metadata is therefore left at the existing
release, pending the user's release workflow.

## Eight-structure sensitivity check (2.1.12)

Eight distinct RCSB inputs were tested through the complete Python API with
default Cone, h-mode 4, coordinates included and no symmetry contacts. Each
run independently recalculated PROPKA predictions; no predicted pKas or
prepared structures were reused between runs. An observational preparation
hook checked invariants without changing the production calculation.

| PDB | No pH option | pH 4.0 | pH 7.4 | pH 10.0 |
| --- | ---: | ---: | ---: | ---: |
| 1CRN | 2 | 2 | 2 | 2 |
| 1UBQ | 2 | 2 | 2 | 2 |
| 1LYZ | 21 | 21 | 21 | 21 |
| 1CA2 | 54 | 46 | 51 | 53 |
| 2ZOI | 23 | 22 | 23 | 23 |
| 3FHP | 14 | 13 | 13 | 14 |
| 7BBC | 129 | 119 | 126 | 129 |
| 5FJJ | 565 | 532 | 539 | 544 |

These are reported interaction rows, **not precision or recall**. All rows
removed relative to the corresponding no-pH run in this panel had HIS as the
pi partner. No additional rows appeared relative to that baseline, and no
shared-row geometry or Hudson/Plevin decisions changed. At higher pH, some
HIS rings returned to the eligible neutral-ring set. For example, the
2ZOI MET A91 CE → HIS A108 pair is excluded at pH 4 and included at 7.4;
this is a model consequence, not evidence that the higher-pH call is correct.
No newly detected acidic OH or HIS NH interactions were demonstrated by this
panel; their implementation checks remain synthetic unit tests.

The no-HIS control 1CRN remained unchanged. 1UBQ and 1LYZ contain HIS but their
reported pairs also remained unchanged. All tested preparations preserved
heavy coordinates/occupancies/charges exactly, matched assigned H names and
counts to their site states, and did not create site H on excluded covalent/
metal-linked sites. Repeated runs were identical, pKas were independent of
the requested pH, and the computed protonated fraction decreased monotonically
with pH. These checks test consistency of the implementation, not the pKa model.

`--ph auto` read 9.0 for 2ZOI (23 rows) and 6.3 for 3FHP (13 rows).
The other six inputs lacked a unique usable pH in the supported metadata;
the supplied grid must not be presented as their experimental conditions.
Additional pH-12 sensitivity checks returned 2, 53 and 554 rows for 1CRN,
1CA2 and 5FJJ respectively. Fixed coordinates at pH 12 do **not** establish
protein stability or correct structure at that pH.

There were **50 complete API calls** across the two batches, including fresh
repeat runs and a repeated 1CA2 control. Full inputs, input/dictionary SHA-256
hashes, results, pair-level differences, warnings and site diagnostics are
retained locally under `_test_ph_2.1.12/six_structures/` and
`_test_ph_2.1.12/additional/` (ignored by Git). Reproduce the main grid with:

```bash
python validation/ph_sweep.py --out /path/to/new/ph_sweep
python validation/ph_deposited_his.py /path/to/new/ph_sweep
# Optional extreme-pH sensitivity check, not a physiological benchmark:
python validation/ph_sweep.py --out /path/to/new/ph_extreme \
  --ids 1crn 1ca2 5fjj --ph 4 7.4 10 12
```

### Deposited neutron H/D spot checks: agreement and a counterexample

The follow-up script records 17 HIS sites in 2ZOI, 3FHP and 7BBC, preserving
H/D atom names, occupancies and N–H/D distances. Missing atoms are not treated
as proof of an unprotonated N. No density maps were inspected in this test.
Only 2ZOI and 3FHP had a usable automatic pH for condition-matched spot checks.

- **2ZOI HIS A3:** at metadata pH 9, the predicted epsilon state agrees with
  deposited DE2 (occupancy 1.0; N–D 0.980 Å).
- **3FHP HIS B5:** at metadata pH 6.3, the model predicts a doubly protonated
  site (pKa 7.500); the deposited model contains DD1 and DE2, each occupancy 1.
  HIS D5 remains unresolved near its predicted pKa, despite deposited D on
  both ring N atoms. The Zn-linked HIS B10/D10 sites are conservatively
  unresolved rather than assigned from the generic pKa calculation.
- **Important counterexample — 2ZOI HIS A108:** predicted pKa 5.555 gives a
  neutral charge at metadata pH 9 (tautomer unresolved). The deposited model
  instead has DD1 occupancy 0.90 and DE2 occupancy 0.85, each at an N–D
  distance of approximately 0.98 Å. This is inconsistent with treating the
  deposited site as a single neutral tautomer and warns against accepting
  the resulting neutral-pi assignment as experimentally established.

The A108 disagreement is **not resolved by this test**. Possible contributions
include pKa-model limitations, crystal/environment effects, interpretation of
pH versus pD, and deposited-model uncertainty. No cause was established, no
occupancy was refined, and no parameters were fitted to make this site agree.
The 90% population policy is not 90% predictive accuracy. The pH option should
remain opt-in, not replace default preparation on the strength of this panel.

Warnings were retained: e.g. 1LYZ emitted 271 PROPKA warnings per run. Its
GLU35 C–OE1 distance is approximately 0.96 Å; the acid geometry guard rejects
this input rather than constructing a protonated carboxyl from implausible
heavy-atom geometry. Such stress cases are useful for software behaviour but
are not suitable as unqualified protonation ground truth.

The next scientific check is a held-out neutron validation with independently
supported protonation/tautomer assignments and known pH/pD conditions. This has
**not** been performed for the new pH option; the spot checks above and existing Cone validation do
not establish the accuracy of these new chemical-state assignments.
