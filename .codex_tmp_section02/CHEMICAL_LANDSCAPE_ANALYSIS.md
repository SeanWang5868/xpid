# Section-02 chemical-landscape analysis

This workflow closes the evidence boundary for the manuscript's high-resolution
PDB characterization section.  It does not alter XPID positives and does not
introduce COD, DFT, refinement restraints or external energetic models.

## Why the background must be regenerated

Schema v8 adds two candidate-level annotations that cannot be recovered from a
positive-only table:

1. complete donor-ring/acceptor-ring geometry for Aromatic C-H candidates,
   independent of the X-H direction;
2. fixed library-H Hudson and Plevin angles for Backbone N-H
   candidates, including candidates that fail the primary directional cutoff;
3. a separate 3.0--5.5 A aromatic ring-image background, independent of any
   donor atom or X-H direction, for the T-shape enrichment denominator.

For symmetry-related contacts, the neighbour-search crystallographic operation
and the marked donor's lattice translation are applied to every atom in the
donor ring or N-H group.  The QC fails rather than silently dropping a candidate
if this complete-group transformation cannot be reconstructed.

## Step 1: regenerate background schema v8

Use the same command and inputs as the accepted background run.  The updated
generator writes a new timestamped output directory and leaves the accepted v6
run untouched.

```bash
python generate_background.py \
  --xpid-source ~/scratch/sofeware/xpid/src/ \
  --positive ~/xhpi_paper/paper_section_02_xhpi_characterization/xpid_result/Xray_1.6A_30July.parquet \
  --metadata ~/xhpi_paper/paper_section_02_xhpi_characterization/xpid_result/Xray_1.6A_30July_metadata.json \
  --pdb-list ~/xhpi_paper/paper_section_02_xhpi_characterization/pdb_id/Xray_1.6A_29July_pdb_ids.txt \
  --pisces-cull ~/xhpi_paper/paper_section_02_xhpi_characterization/pdb_id/cullpdb_pc25.0_res0.0-1.6_len40-10000_R0.25_Xray_d2026_07_29_chains4791 \
  --pdb-mirror /scratch/bql506/pdb_mirrors/pdb/mmCIF/ \
  --redo-mirror /scratch/bql506/pdb_mirrors/pdb-redo/ \
  --monomer-lib /run/user/198171/cache/xpid/ccp4-monomers/ \
  --jobs 64
```

The final QC must report 4,664 processed structures, no processing failures,
no unmatched positives and schema version 8.

## Step 2: run the manuscript analysis

Replace `BACKGROUND_V8` with the new output directory.

```bash
python analyze_chemical_landscape.py \
  --background BACKGROUND_V8 \
  --bootstrap 2000 \
  --seed 20260801

python analyze_residue_donor.py \
  --background BACKGROUND_V8 \
  --bootstrap 2000 \
  --seed 20260802 \
  --force
```

For a fast three-structure smoke test:

```bash
python analyze_chemical_landscape.py \
  --background BACKGROUND_V8 \
  --max-structures 3 \
  --bootstrap 100 \
  --output BACKGROUND_V8/analysis/chemical_landscape_test \
  --force
```

## Outputs

- `global_summary.csv`: full-model unique-pair count, structure coverage and
  per-structure contact-density summary; raw XPID rows are reported separately.
- `global_structure_summary.csv`: the auditable PDB-level numerator,
  standard-residue denominator and pairs per 1,000 residues.
- `main_contact_composition.csv`: exact 173,732-pair master composition.
- `acceptor_composition.csv`: acceptor composition in the master set.
- `acceptor_propensity.csv`: selected-acceptor eligible, participating and
  contact-multiplicity estimands.
- `donor_acceptor_matrix.csv`: distance, H-independent ring-face and positive
  matrices. Acceptor enrichment is calculated within each donor class against
  that donor class's ring-face acceptor distribution, with PDB-cluster 95% CI;
  it is a structural realization ratio, not an affinity estimate.
- `tshape_summary.csv`: the exact Aromatic C-H positive-contact T-shape
  fraction, including symmetry contacts, and unique-ring-pair enrichment over
  the fully matched H-independent non-symmetry nearby-ring background.
- `gly_backbone_nh_funnel.csv`: eligible -> distance -> ring-face -> direction
  funnel for Gly, non-Gly and all backbone N-H donors.
- `backbone_nh_angle_sensitivity.csv`: primary and prespecified +/-5 degree
  directional-threshold sensitivity in the blank-altloc donor/ring subset;
  the four-stage Gly funnel itself retains all contacts.
- `distance_by_donor_class.csv`: medians and interquartile ranges.
- `bootstrap_estimates.csv`: PDB-cluster 95% confidence intervals for major
  proportions, distance contrasts, T-shape enrichment and Gly/non-Gly rate ratio.
- `qc_report.json` and `analysis_metadata.json`: scope, provenance and audit.

`analysis/residue_donor/` additionally contains the latest schema-v8
`residue_contribution.csv`, `residue_atom_propensity.csv`,
`donor_class_propensity.csv`, `donor_class_geometry.csv` and
`backbone_nh_by_residue.csv`. Residue and donor-class rate intervals resample
PDB structures; detailed atom rows remain descriptive.

The primary statistical unit for uncertainty is the PDB structure.  Contacts
within a structure are descriptive subsamples and are not treated as independent
replicates.  The analysis reports effect sizes and confidence intervals, not
contact-level significance tests.

The analysis QC fails if the 280,843 full-model unique-pair count, 4,644
positive-structure count or 173,732-pair master count changes, if any Aromatic
C-H candidate lacks complete ring geometry, if a positive nearby ring pair is
absent from the H-independent background, or if the recomputed primary
Backbone N-H angular classification differs from the recorded classification.
