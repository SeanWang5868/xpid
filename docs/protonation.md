# Optional pH-dependent protonation

This is a **conditional chemical model**, not an experimental determination of
proton positions, a refinement program, or a constant-pH simulation. It changes
which protein side-chain X–H groups can participate in detection, but does not
change Hudson/Plevin cutoffs, the Cone sampling or representative ranking.
The default (no `--ph`) retains the existing Gemmi preparation and zero-H-
occupancy filtering.

## Use

```bash
pip install "xpid[ph]"
xpid model.cif --ph 7.4 --file-type csv --out-dir ph_results --provenance
xpid model.cif --ph auto --file-type csv --out-dir ph_results --provenance
```

For a local checkout, install with `pip install -e '.[ph]'`.
H modes 3, 4 and 5 are supported. Modes 0–2 are rejected with `--ph` because
their preserve/shift/remove semantics conflict with state rebuilding.
Numeric pH must be finite and in 0–14; this range is an input safeguard, not a
claim that the structure or every unmodelled group remains valid at extreme pH.

`auto` uses `_exptl_crystal_grow.pH` for diffraction, `_em_buffer.pH` for EM,
`_pdbx_nmr_exptl_sample_conditions.pH` for NMR, or explicit PDB REMARK 200 pH.
Missing/conflicting values and ranges are errors. No value of 7 or range
midpoint is silently substituted. Experimental buffer pH is not necessarily
the effective pH of every site in a cryocooled crystal.
No heavy-water pD correction or isotope-dependent pKa shift is inferred.

PROPKA is optional, but required for the default `--ph` method. A deliberately
simpler sensitivity analysis can use `--pka-method reference`. Its model-
compound values (from PROPKA 3.5) are ASP 3.8, GLU 4.5, HIS 6.5, CYS 9.0,
TYR 10.0, LYS 10.5 and ARG 12.5. These are **not structure-specific pKas**.
Backend failures are diagnosed; they never trigger a reference-pKa fallback.

## Chemical decisions and limits

For an uncoupled site, the estimated protonated fraction is
`1 / (1 + 10**(pH - pKa))`. Only a charge state with an estimated population
of at least 0.90 is assigned automatically. This is a conservative reporting
policy, **not a calibrated confidence interval**: pKa prediction errors remain.
Near-pKa sites and PROPKA-flagged coupled sites are unresolved; no independent
fraction is reported for coupled sites. Results from different possible
microstates are never unioned into an unconditional positive call.

| Site | Treatment |
| --- | --- |
| HIS | Separate charge assignment from neutral delta/epsilon tautomer choice; doubly protonated imidazolium is not a neutral π acceptor or N–H donor in this analysis. A neutral ring can remain a π acceptor when its tautomer is unresolved. |
| ASP/GLU | Deprotonated carboxylates have no O–H donor. Protonated candidates cover both oxygens and both planar O–H directions using CCP4 `ASPprot`/`GLUprot` nuclear-distance and angle restraints. Only a resolved candidate is used. Carboxyl O–H is not a freely rotating hydroxyl Cone. |
| CYS/TYR | A deprotonated site cannot be recreated by Cone. Resolved neutral SH/OH uses the existing chemical Cone. Covalent/coordination links, missing site atoms and suspected short extra bonds prevent automatic assignment. |
| LYS/ARG | Charge states are assessed, but their N donors remain excluded under XPID's existing cation–π policy. Neutral amine/guanidine donor geometry is not generated. Neutral LYS NZ can be annotated as an H-bond acceptor; neutral ARG tautomer/acceptors remain unresolved. |
| Termini | No pKa assignment in this release. Recognised terminal/chain-break amines and terminal carboxyls are explicitly unresolved and excluded at affected donor/acceptor atoms. A backbone NH elsewhere is not titrated. |
| Ligands, nonstandard residues, other groups | Existing dictionary preparation is retained; their protonation is **not inferred from pH**. The feature is not a general ligand protonator. |

Heavy coordinates, occupancies, names and deposited formal-charge fields are
not changed. State-dependent charge and donor/acceptor typing live in a separate
context shared by detection, Cone and H-bond annotation. HIS is never renamed
to HIP/HID/HIE: these are not portable CCP4 chemical-component identities.

HIS candidate directions come from Gemmi's dictionary-generated H placeholders.
Assigned pH-state Hs use dictionary **nuclear** X–H lengths. Other prepared Hs
retain the existing Gemmi convention; this option is not a global conversion
of all H distances. Acid candidates use the named CCP4 modifications directly;
missing restraints yield an unresolved site, not guessed geometry.

Tautomer/orientation selection uses a deliberately limited, deterministic local
geometric discriminator: reject H–heavy distances below 1.45 Å (except the
parent and its bonded neighbours), then score H-bond-like contacts with
1.5–2.6 Å H···A and angles ≥140°. A unique surviving candidate or a score
advantage ≥0.25 resolves the choice. Competing titratable side chains do not
vote in this local score; no global H-bond-network optimization is claimed.
These selection thresholds are heuristics requiring independent validation,
not energy parameters. The score never uses the aromatic-ring location or an
XH–π result. Conventional H-bond annotation still does **not** veto Cone hits.

Site identities include model index, full chain name, residue sequence ID
(including case-sensitive insertion code), residue name and altloc. No A/B
coordinates are combined to form a hydrogen. For PROPKA's PDB interface, chain
names and residue numbers are mapped reversibly; one highest-occupancy
atom-name conformer is passed per residue. Other alternate states receive no
copied pKa and remain unresolved unless explicitly assigned. This is not an
alternate-conformer ensemble pKa calculation. The adapter is limited to 62
chains and 9999 residues per chain; unsupported input is diagnosed. Only models
selected by `--model` are submitted to PROPKA, and no `.pka` files are created.
PROPKA sees the finite input coordinates, not a periodic crystal or an explicit
solvent ensemble; its predictions do not include all crystal-packing effects.

## Explicit microstates

Use only for an independently justified assignment, not to improve the number
of detected interactions. For one input structure:

```json
[
  {"model_index": 0, "chain": "A", "residue_id": "57", "residue_name": "HIS", "altloc": "", "state": "delta"}
]
```

```bash
xpid model.cif --ph 7.4 --protonation-overrides states.json --out-dir ph_results
```

HIS states are `delta`, `epsilon`, `doubly_protonated`. CYS/TYR/LYS/ARG accept
`protonated` or `deprotonated` (neutral LYS/ARG geometry remains unsupported).
ASP/GLU accept `deprotonated` or e.g. `protonated_OD2_trans` /
`protonated_OE1_cis`: `cis` points in-plane toward the other oxygen's side;
`trans` points away. Overrides do not bypass missing geometry, coordination,
covalent links, quarantine or the candidate clash check. An unmatched identity
is an error; an unsafe matching override is reported as not applied.

## Output and validation

Every pH run adds `pH`, `protonation_method`, donor state/source/pKa/protonated
fraction, and HIS π-state columns. `H_source=ph_model` means a selected model,
**not an observed proton**. Assigned H occupancy inherits its heavy parent,
capped by the active labelled heavy-atom conformer occupancy:
the estimated chemical-state population is separate and is never multiplied
into crystallographic occupancy or the interaction `combined_occupancy`.

The existing `*_diagnostics.json` contains every assessed site, including sites
with no hits, reasons for unresolved/excluded states, candidate scores, backend
version and pH source. Always inspect it before comparing pH runs. A zero hit
count at an unresolved site does not establish absence of an interaction.

```python
from xpid import detect

audit = {}
hits = detect("model.cif", ph=7.4, pka_method="propka", diagnostics=audit)
sites = audit["protonation"]["sites"]
```

Run the new chemical-state tests with
`python -m pytest test/test_protonation.py -q`, and the full suite with
`python -m pytest -q`. They check H counts/occupancies, charge vs tautomer,
altloc identities, metal protection, deprotonated Cone suppression, acidic
O–H detection, default-path compatibility, API/CLI agreement and output schemas.
These tests and a real-structure smoke run establish implementation behaviour,
**not the accuracy of automatic protonation assignments**. Independent neutron
validation across pH, metal sites and tautomer classes is still required before
using this option as a production replacement for the default preparation.

Method references: [PROPKA documentation and citation guidance](https://propka.readthedocs.io/en/latest/),
[Gemmi monomer-library chemistry](https://gemmi.readthedocs.io/en/latest/chemistry.html).
