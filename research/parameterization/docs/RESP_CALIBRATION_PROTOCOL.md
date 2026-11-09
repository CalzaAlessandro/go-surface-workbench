# Proposed RESP calibration of optional static redistribution

This is a calibration design, not a fitted SQE-GO model. Update 2026-11-03:
HF–RESP is the primary reference; DFT–RESP and AM1-BCC remain controls. Fixed
modular, context-library and static redistribution candidates remain under
comparison. The final integer charge follows the product's composition and
protonation and is not forced to zero. Charges stay fixed throughout MD and
are inherited by MM/GBSA subsets. See the current
[strategy benchmark](CHARGE_STRATEGY_BENCHMARK_20261103.md).

## Baseline and training references

Retain the new jointly constrained modular RESP library as an L0 candidate,
with each source/state/fit scope explicit. Sixteen conventional restricted-HF
RESP references on accepted DFT geometries are now available; they do not
claim wavefunction-stability validation. The earlier instability diagnostic
and the conventional reference definition remain separate scopes under
[AMBER_RESP_ACCEPTANCE.md](AMBER_RESP_ACCEPTANCE.md). No production parameters
are selected merely because a fit completed.

Build target data on mixed, terminal and differently sized fragments. Use
whole-fragment RESP targets with correct total charge and justified equivalences,
without pre-imposing zero charge on every background atom or identical group
charges across all chemical environments. Training only against the L0 vector
would erase the very redistribution being investigated. Retain the original QM
ESP grids and molecular moments as independent descriptions of fit quality.

The prospective panel varies OH/OH, EPX/EPX and OH/EPX distance and relative
face, includes two PAH sizes and EPH/EPT/CAB/COOH/COO controls. Larger or charged
edge/basal combinations still need expansion; the panel is not full coverage.

## SQE-family calibration proposal

Use charge-conserving bond transfers added to a reference per-atom vector.
Chemical atom/bond classes distinguish aromatic, hydroxylated and epoxidized
carbon, ring/terminal environments, carbonyl and protonation. Avoid unnecessary
environment subclasses before the data support them.

Initially fix the reference vector to the selected RESP block library, then
fit a small set of response parameters: electronegativity differences, atomic
hardness and bond hardness. Keep the distance-smearing model explicit and test
its sensitivity; do not freely optimize every parameter family on a small
dataset. Fix parameter gauges, check identifiability, and require stable solves
on all supported input graphs. The reference-charge choice and any extra
constraints make this a GO calibration/variant, not the published parameter set.

Use molecule-normalized, atom-type-balanced charge residuals alongside common-
grid ESP and moment errors, plus regularization against unsupported large
corrections. Select relative weights by held-out performance. Do not let the
many nearly uncharged carbons dominate a metric and hide errors at functional
groups. Compare predictions against current charges, constant L0 RESP and
available native/derived controls with their scopes kept distinct.

Always preserve the molecular charge. Compare an unrestricted global-transfer
candidate with deliberately local/block-conserving candidates if data warrant
them: neutral blocks and frozen-zero scaffold remove response channels, and
their effect must be measured rather than assumed harmless. Native global
conservation does not enforce each block's net charge.

Hold out entire physical fragments, separations and sizes. Three rigid rotations
of one QM geometry are one chemical context, not three independent validation
examples. Test size dependence, distorted/folded geometries, charged termini and
unseen functional-group combinations. Unsupported contexts remain explicit.

Finally evaluate coherent charge/LJ/bonded models using relevant QM interaction
curves, hydration and adsorption observables. Lower training RESP error alone
does not establish better MD or MM/GBSA performance.

## Published method families to compare

- [SQE+qp](https://doi.org/10.1186/s13321-021-00528-w): published reference-charge
  plus bond-response architecture; existing example parameters target NPA.
- [QTPIE](https://doi.org/10.1016/j.cplett.2007.02.065): alternative charge-flow
  model with improved dissociation behavior; GO parameters still need testing.
- [IPolQ](https://doi.org/10.1021/jp311851r): fixed charges incorporating an
  explicit-water electrostatic environment through iterative MD/QM. Solvent
  and force-field coupling are integral to the protocol.
- [RESP2](https://doi.org/10.1038/s42004-020-0291-4): combines gas/aqueous QM
  information to tune fixed-charge polarity, alongside nonbonded optimization.
- [EspalomaCharge](https://doi.org/10.1021/acs.jpca.4c01287): a learned
  electronegativity/hardness model approximating AM1-BCC; extrapolation to
  extended GO is not established by its molecular benchmarks.

These are alternatives at different levels, not interchangeable correction
coefficients. IPolQ/RESP2 are possible expensive reference branches; a small
SQE-family or local-increment model is a candidate for routine preparation.
