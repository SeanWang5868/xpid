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

The 2.2.0 wheel was built and checked separately. No GitHub/PyPI/Zenodo release
was published. Published citation metadata is therefore left at the existing
release, pending the user's release workflow.

The next scientific check is a held-out neutron validation with independently
supported protonation/tautomer assignments and known pH/pD conditions. This has
**not** been performed for the new pH option; the existing Cone validation does
not establish the accuracy of these new chemical-state assignments.
