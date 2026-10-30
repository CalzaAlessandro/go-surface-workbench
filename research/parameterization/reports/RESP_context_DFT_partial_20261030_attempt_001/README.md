# Interim DFT–RESP controls: four completed structures

This is a dated partial calculation snapshot, not a transferable parameter set. The primary new derivation remains RESP; these explicit DFT–RESP controls use ωB97X-D/6-311++G(d,p) MK potentials after accepted geometry, frequency and wavefunction-stability gates. They are separate from the unstable HF/6-31G* diagnostic references.

Four structures were captured. Three rigid MK orientations per structure enter one atomic RESP fit: there are four physical contexts, not twelve independent chemical samples. No new bonded/frcmod coefficients or production charge libraries were written.

## Verified evidence

The downloaded archive SHA-256 is `6ec622c478c71d7abee9fa8a581fa60d540f7ef4ac7fa4a6adf2623e77bf4b16`. All 100 collected files matched their recorded hashes. Before collection, every original fit artifact was verified on the cluster. The corrected MK gate independently verified actual finite SCF energies and matching electronic methods in all 12 logs. All 120 candidate/orientation score rows were recalculated locally on the original espgen grids; supported scores reproduce the recorded results. Charge identity and total-charge checks passed at the documented RESP output-rounding tolerance.

## Training ESP comparison

The table reports the arithmetic mean relative ESP RMSE across the three numerical orientations of each one fitted geometry. This aggregation describes numerical orientations; it is not a sampling error estimate. The new atomic fit was trained on these same grids, so its smaller training residual is not held-out evidence of superior MD or adsorption predictions.

| Finite structure | New atomic RESP fit | Current fixed blocks | HierGO OPLS charge table |
|---|---:|---:|---:|
| pyrene_pristine | 0.1100 | 0.4514 | 0.2362 |
| pyrene_diol_cis | 0.1076 | 0.5487 | 0.4708 |
| pyrene_diol_trans | 0.1208 | 0.3191 | 0.2395 |
| pyrene_dione | 0.0655 | 0.3775 | 0.3596 |

The controls compare charge representations only. They do not execute the complete original OPLS or GAFF2–AIM force fields, and they do not replace their published Lennard–Jones coefficients. GOPY entries with a charge mismatch are excluded rather than normalized. OLD whole-sheet charges and GAFF2–AIM do not have an established exact mapping onto these new finite fragments and remain unsupported here; missing coverage is not zero error.

## What the raw charges do and do not show

For the trans diol control, the two fitted OH groups share C = +0.429174 e, O = −0.546598 e and H = +0.328281 e. Each C–O–H partition therefore has net charge +0.210857 e; the rest of the neutral fragment compensates it. The current fixed OH block uses +0.250/−0.640/+0.390 e and is exactly neutral. This difference quantifies a representation constraint in this finite control; it does not establish an optimum transferable OH block.

For the cis control, the two hydroxylated-carbon charges are +0.190507 and +0.790822 e. Their OH partition sums are −0.007683 and +0.506003 e. Do not average these directly into a library: geometry/environment asymmetry and fitting sensitivity have not yet been separated by the broader panel and additional reference checks.

The pristine finite pyrene control already assigns nonzero charges to interior carbons (−0.029809 to +0.222699 e). An oxidized-minus-matching-pristine comparison is therefore required before attributing background charges to oxidation or deriving a distance correction. These finite-parent values must not be copied as the charge of macroscopic pristine graphene.

## Preserved outputs

- `per_atom_charges.csv`: exact atom identities, roles and fitted charges.
- `functional_group_charges.csv`: group partitions, without enforcing local neutrality.
- `functional_group_pairs.csv`: graph/centroid distances and raw paired-group charges.
- `background_carbon_distance_charges.csv`: interior scaffold and hydrogen-terminated edge separated.
- `role_charge_ranges.csv`: raw atom counts, means and ranges, with the current fixed role value.
- `training_ESP_comparisons.csv`: all 120 orientation/candidate rows, including exclusions.
- `PARTIAL_REVIEW.json`: gate and numerical-reproduction evidence.
- `SOURCE_SNAPSHOT_MANIFEST.json`: source paths, file hashes and collection time.

`dataset.json` in the local study records absolute paths on the analysis machine. To reproduce the context tables elsewhere, extract the matching archive and build a new input list of `{model_path, charges_path, fit_evidence_path}` entries pointing to each `<model>/fit/{model.json,charges.json,FIT_COMPLETE.json}`; then run `PYTHONPATH=src python -m go_param.context dataset.json --output-dir NEW_OUTPUT`. The original calculation archive remains outside the Git repository.

## Remaining scientific decisions

The remaining accepted-parent controls, matched C₂₄/C₅₄ isolated controls and distance/face panel are still required. The data do not yet choose between fixed modular RESP, an environment lookup, a compact correction rule or a RESP-calibrated published redistribution method. General transfer, protonation, conformation/fitting sensitivity and coherent charge/LJ performance remain unvalidated.
