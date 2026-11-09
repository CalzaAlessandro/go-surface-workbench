# Sixteen completed DFT–RESP comparison controls (7 pyrene C16 + 9 coronene C24)

All 16 accepted-parent controls (7 pyrene $\text{C}_{16}$ + 9 coronene $\text{C}_{24}$, including the recovered `C24_OH_OH_near_trans` job `1408584`) completed gas-phase $\omega\text{B97X-D/6-311++G(d,p)}$ geometry optimization, harmonic frequency verification, external wavefunction stability verification, three-orientation rigid Merz–Kollman ESP evaluation, and Amber `resp` fitting.

- **Pyrene ($\text{C}_{16}$) source snapshot**: `../../calculations/dft_resp_final_context_snapshot_attempt_001` (SHA-256 `84a8a6024b0e23dbec09bb34a4b5cff6aeb18b0e7b85200b8531e0f1c06c5b29`)
- **Coronene ($\text{C}_{24}$) source snapshot**: `../../calculations/c24_completed_context_snapshot_20261101_attempt_001` (SHA-256 `21c8e0288509b2ca22b6672451d87531e0f98c35baeac0c3529ba90275acea0d`, 72 verified files)
- **Pending cluster controls**: `C54_EPX_EPX_near_trans`, `C54_pristine_parent_reference`, `C54_single_basal_EPX_control` (running on Ciamician with 5-day walltime).

## Files in this report
- `context_analysis.json`: full structured output from `go_param.context` across 16 contexts and 15 functional-group pairs.
- `per_atom_charges.csv`: atomic DFT–RESP charges, roles, coordinates, and scaffold flags.
- `functional_group_charges.csv`: unaveraged net group charges and scaffold-anchor charges.
- `functional_group_pairs.csv`: graph bond distances, 3D centroid distances, and paired net charges.
- `background_carbon_distance_charges.csv`: interior scaffold vs H-terminated edge charges by bond distance to nearest functional group.
- `training_ESP_comparisons.csv`: 480 candidate/orientation rows (16 models $\times$ 3 orientations $\times$ 10 candidate presets).

> **Scope note**: These are unconstrained atomic DFT–RESP comparison controls, not the final Amber-compatible modular library preset. Training-grid ESP residuals compare electrostatic reproduction on the QM parent grid and do not by themselves validate condensed-phase adsorption or full force-field accuracy.
