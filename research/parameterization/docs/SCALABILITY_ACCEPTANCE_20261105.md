# Charge-model transfer across sheet size and oxidation

User requirements confirmed 2026-11-05: preparation methods must work on about
4 nm and 10 nm GO sheets, with both low and high oxidation, then retain fixed
charges during MD with small molecules, proteins or DNA. A new whole-sheet QM
fit for every production geometry is not an acceptable routine dependency.

## Two independent validation axes

- Compare different sheet sizes at matched functional-group fractions,
  protonation, edge chemistry and spatial organization as closely as feasible.
- Compare oxidation levels at fixed sheet size. Record O/C and C/O with clear
  conventions, plus actual group counts; equal C/O does not mean equal chemistry.
- Include isolated and clustered OH/epoxide, mixed groups and edge distances;
  defect/ether/ester groups require explicit validated coverage.
- Report different edge-to-area ratios as a physical size effect. Do not make
  small-flake H caps stand in for basal carbon or mix them in training averages.
- Test distorted/folded geometries separately. Graph classes do not change with
  rotation/folding, but the electrostatic validity of unchanged charges still
  needs checking. SQE's 3D response may change during preparation; q then stays fixed.

Current fixed-count C16/C24/C54 fragment comparisons conflate size and coverage
in some series. They are exploratory evidence, not acceptance of the 4–10 nm
or low-to-high-oxidation range. No numerical validity range is claimed yet.

## Acceptance gates for each candidate

1. Same class definitions/reference parameters across sizes; exact atom mapping,
   finite charges and integer total from protonation and chemistry, including
   negative sheets and grafted groups. No size-dependent global charge rescaling.
2. Explicit coverage for low/high oxidation and mixed dense neighborhoods;
   unseen classes are unsupported rather than silently replaced.
3. Entire structures, loadings and combinations excluded from parameter fitting;
   rotations, atoms or duplicate conformers are not independent test structures.
4. Report ESP/moments and region/group diagnostics, not only atom-charge RMSE.
   Assess total redistribution and response far from a local chemical change.
5. Record time/peak memory vs actual atom and bond counts. Lightweight assignment
   must not depend on a new QM solve for each sheet.
6. Before MD use, validate electrostatics with the chosen LJ/bonded/water/ion
   model and representative hydration/adsorption tests. Preserve protein/DNA
   force-field and ligand preparation provenance separately.

## Implementation status and limits

- L0 and the first-neighbor classifier assign local classes on a supplied,
  validated molecular graph. Their production-scale scientific transfer remains
  unverified; no claim of passing a 10 nm oxidized sheet is made here.
- `go_param.sqe._operators` currently builds dense atom interactions and a dense
  incidence/response solve. This research implementation must be benchmarked
  before use on a large sheet; a sparse/iterative or otherwise bounded solver
  needs separate numerical and physical validation if the current cost is too high.
- A graph-only class model cannot reproduce every long-range electrostatic or
  conjugation effect automatically. Finite-fragment fits and first-neighbor
  classification do not themselves prove size consistency.
- New ether/lactone/lactol fragments are chemistry references for later fitting,
  not large-sheet validation or assigned charge/force-field presets.
- The current research partitioner also inherits the atomic RESP input's
  conservative formal-charge range −98 to +98. Routine assignment on highly
  charged large sheets must use a graph/chemistry validator separated from
  this fit-input restriction; it must not change or neutralize the chemistry
  merely to pass a RESP header limit. This gate is not implemented here.

No existing model, charge library or running QM job is changed by these gates.
