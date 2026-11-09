# Environment RESP pilot implementation plan

> Execute inline with executing-plans and test-driven-development. User approved
> the environment-library and RESP-assumption work on 2026-11-05. No Git needed.

**Goal:** prepare an invariant first-neighbor charge classification, run an actual
Amber RESP pilot on existing grids, and quantify the localization assumptions.

**Spec:** the approved work is items 1–2 in the conversation, together with
CONTEXT_DEPENDENT_CHARGES.md and RESP_CALIBRATION_PROTOCOL.md. Canonical research
sources remain PARAMETERIZATION; mirror reviewed source/docs/reports in WORKBENCH.

## Scientific design and constraints

- New first-neighbor classes are inspired by GAFF2–AIM's organization, not its
  DDEC6 charges or LJ parameters. Published presets stay intact.
- Functional carbon environments use the multiset of primary chemical roles of
  directly bonded atoms. A group's carbon environments define its shared group
  context; individual carbon environments may distinguish asymmetric carbons.
- Group O/H charges inherit group context. Equivalent COO oxygen labels are
  insensitive to the localized resonance drawing. Graph signatures never use
  atom IDs, model IDs, coordinates, iteration order or arbitrary Kekulé bonds.
- For this first controlled L1 comparison, retain L0's zero basal background,
  neutral blocks, COO block −1 including its anchor, and shared finite-H caps.
  Separately audit those assumptions; L1 does not establish their correctness.
- Use actual external Amber RESP with the existing hyperbolic restraint,
  equivalences and orientation weighting. No regression surrogate is called RESP.
- Train on the seven predeclared C16 contexts; test on original C24. Fit all
  available contexts only as an explicitly in-sample diagnostic if engine
  capacity permits. New classes remain unsupported; no hidden L0 fallback.
- Report class support by physical model, not by atom or rigid orientation.
  Split entire physical models; first-neighbor classes ignore relative face,
  which remains a limitation to evaluate rather than a hidden geometry rule.
- Keep raw QM inputs/results and active jobs unchanged. New jobs use isolated
  calculation folders, bounded CPUs and a five-day maximum.
- No final reference, production library, bonded or LJ parameter changes.

## Tasks and verification

1. Add tests for permutation/rotation invariance, distinct chemical contexts,
   COO resonance, input immutability, strict unseen-class rejection, and exact
   block/total conservation. Run and see the missing-feature failure.
2. Implement `go_param.environments.classify_model`, `build_problem`, `coverage`
   and `predict`; reuse existing RESP card generation and result validation.
   Run the new tests and the full parameterization suite.
3. Freeze and execute a new actual-RESP pilot using verified archived grids.
   Check convergence, shared labels, every group/total charge, raw-grid hashes,
   held-out support and matched electrostatic scores. Archive input/output hashes.
4. Integrate an independent scalar audit of existing molecular HF/DFT RESP
   charges, including basal/edge separation and COO anchor decomposition.
   Produce Italian HTML tables and CSV/JSON, explicitly retain reference caveats.
5. Review all new code/results, verify links/hashes, mirror only necessary
   code/docs/reports and update campaign status. No production promotion.

## Review focus

- Hidden atom-number or geometry dependence in class keys.
- Treating COO resonance oxygen labels as distinct chemistry.
- Training atom counts mistaken for independent class support.
- Comparing different supported subsets or reporting in-sample fit as transfer.
- Dropping dependent constraints without retaining every group's output check.

## Execution ledger

- Started: inspected existing modular RESP and charge-context analyzers; reuse
  their validated partition and external RESP conventions.
- Ruling: preserve zero-background/block constraints in the first L1 pilot to
  isolate the effect of context classes; a separate audit assesses their limits.
  Cost: L1 cannot yet describe inter-block charge transfer or basal polarization.
- Ruling: no Git operations, in accordance with the user's existing instruction.
- RED: seven tests failed because the new environment module did not exist.
  GREEN: the seven invariant/conservation tests pass; full suite 245/245 passed.
- First-neighbor coverage before fitting: training on seven C16 models covers
  2/9 original C24 contexts. Adding two new loading references covers 4/9.
  These counts are physical models, not atoms or rigid orientations.
- Four actual Amber RESP fits submitted as job 1469243: L0/L1 on C16 only,
  then L0/L1 on C16 plus the two explicitly declared new loading controls.
- Ruling: all16 joint fit omitted because 1641 atom copies exceed the existing
  1024-center executable. The nine-model expanded training has 837 copies;
  no engine source or installed software was modified.
- Added user requirements: transfer across approximately 4–10 nm and low/high
  oxidation are separate acceptance axes; see SCALABILITY_ACCEPTANCE_20261105.md.
- Added scope: numeric frcmod audit and literature/QM-needs comparison in
  isolated reports; new ether/lactone/lactol reference fragments remain unfitted.

- Complete: four native fits passed all class/group/total gates; 81 raw grids verified; old L0 exactly reproduced. Matched L1 does not improve L0 on this limited panel.
- Complete: independent core review found no important defect; extended-fragment review found MMFF aromaticity drift. Added a failing coumarin SDF round-trip regression, restored standard aromaticity after MMFF, and exported attempt_002 without changing historical attempt_001.
- Verification: complete suite 265/265 passed; numerical frcmod reconstruction and source hashes verified. Root HTML tables/plots inspected; dense numeric tables remain horizontally scrollable.
- Live QM refresh: 4 expanded HF fits completed, 2 rejected (one stability gate, one incomplete optimization); reference acetate requires review of small negative frequencies. These newer cases did not enter the frozen L1 pilot. Final QM reference and production transfer remain open.
