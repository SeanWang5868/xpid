# Hydrogen-bond compatibility audit

```bash
xpid model.cif --hbond-audit --out-dir audited --provenance
```

This optional analysis **does not change detected interactions, Cone geometry,
steric thresholds, selected representative conformers, or existing result
columns**. It does not lock H atoms to acceptors or apply a hydrogen-bond veto.
Version remains 2.1.12. It requires Cone auto and cannot be combined with
`--no-cone` or `--xh-candidates`.

Read `structures[*].hbond_audit` in `xpid_results_diagnostics.json` (or the
corresponding custom output stem). This file is written independently of
`--provenance`. The ordinary JSON/CSV/Parquet interaction output is unchanged.

## What is recorded

For every eligible rotatable donor-conformer/ring pair reaching Cone after
the existing chemical and spatial prefilters, the audit records:

- Model, donor/parent altloc identities, donor periodic image, and target ring.
- Number of generated, sterically allowed, and XH–pi-positive conformers.
- Nearby acceptors with chemical typing, state source, occupancy, image and
  coordinates in the same frame as the donor.
- Independent counts of same-H coexistence, different-H coexistence within
  a complete conformer, H-bond-only conformers and XH–pi-only conformers.
- `disjoint_only`: both kinds of contact exist, but never in the same allowed
  conformer. This means **potential directional competition**, not occupation
  of the H or evidence that the XH–pi interaction should be deleted.

`typed_contacts` combines contacts to acceptors supported by the typing rules;
`uncertain_contacts` separately records possible contacts whose chemical or
conformational state is unresolved. Each acceptor also has its own relationship
counts. Thus an aggregate compatible contact does not hide competition with
another acceptor. Negative pairs are retained too, including pairs with no
allowed conformers. `detected` is the pair-level Cone result before final hit
deduplication, not a count of rows in the results table.

Each acceptor includes one geometric witness (`example`: rotation angle,
zero-based hydrogen index, H···A distance and X–H···A angle). It is the first
sampled witness, **not the strongest or most likely interaction**.

## Chemical scope and limitations

The audit uses the existing contact geometry: 1.30–2.80 Å H···A and at least
120° X–H···A. These are screening conditions, not an energy test, proof of a
hydrogen bond, or a universal calibrated constraint for every donor class.
CH contacts are explicitly labelled `weak_CH_contact_only`; SH geometry is
not independently calibrated by this feature. No population is inferred from
the number of sampled orientations.

Known backbone/amide carbonyls, hydroxyl O, water O and standard thioether
acceptors are typed conservatively. Unknown ligand O/S/N, uncertain terminal
states and unassigned His/carboxyl/thiol states are not promoted to definite
acceptors merely because H is absent. Occupied His N–H/D is not an acceptor;
zero-occupancy H is not proof of deprotonation. An available resolved pH state
is used with its original source recorded, not presented as experimental truth.
Unresolved states, declared covalent/metal links, nearby metal or unexpected
short heavy-atom bonds remain uncertain. Water orientation is not optimized.

The audit uses occupied, filter-compatible atoms and the existing alternate
selection convention. Matching nonblank labels on different residues do not
establish their correlation: such acceptors are marked uncertain when the
donor also has a labelled alternate. Other cross-residue alternate combinations
are not exhaustively enumerated. It is not a joint conformational ensemble.

The audit retains same-residue noncovalent contacts, excluding the actual X,
parent and dictionary-connected 1–3 neighbours. Image positions are resolved
around the donor. Crystallographic rotation images follow the detector's cell
image setup; nearest lattice translations remain possible even without
`--sym-contacts`, consistent with the existing environment logic.
`--include-water` does not control whether water can be an environmental acceptor.

The production steric mask is reused exactly, including its existing H-bond
short-contact exemptions and acceptor rules. **The stricter audit typing does
not retroactively change that mask.** Its separate search radius is X–H length
plus 2.8 Å, covering the longer S–H geometry without enlarging the production
steric search. Sites rejected before Cone are outside this audit's scope.

There is no local-network optimization, pH feedback, rigid-H reorientation,
receiver lone-pair directional test, or experimental density validation here.

## Reproducible details and Python API

Masks are hex-encoded `numpy.packbits` with big bit order. The valid mask has
one bit per conformer; `pi_mask` and each acceptor's `contact_mask` use the
reported `[conformers, hydrogens]` shape in row-major order. Remove padding bits
before reshaping. Rotation angles are the existing 1° sampling starting at 0°,
with the group's existing period and hydrogen ordering.

To bound diagnostics memory, at most 10,000 pair records are stored per
structure. `pairs_seen` and `pairs_omitted` explicitly report truncation;
the CLI warns when it occurs. Detection itself is never truncated. These
records are not an unbiased sample when the cap is reached.

```python
from xpid import detect

diagnostics = {}
hits = detect("model.cif", hbond_audit=True, diagnostics=diagnostics)
audit = diagnostics["hbond_audit"]
```

The Python API requires a diagnostics dictionary when the audit is requested.
It can be combined with `ph=7.4` if the optional pH dependencies are installed.

Run `python -m pytest test/test_hbond_audit.py -q` for targeted regressions.
For local real structures:

```bash
python validation/hbond_audit_check.py model1.cif model2.cif --out new_audit_check
```

The script compares complete coordinate-inclusive API results with the audit
off/on and retains input hashes and pair diagnostics. A passing comparison
establishes software invariance, **not the physical accuracy of competition
assignments**. See [implementation validation](hbond-audit-validation.md).
