# GO QM reference validation implementation plan

Goal: compare electronic reference methods without changing the conventional
Amber HF–RESP branch or any production force-field parameter.

Architecture: a separate fixed-geometry Gaussian validation module, frozen
input campaigns, and gated postanalysis of direct potentials and dipoles.
System Gaussian/Amber only; Python/NumPy for analysis. Original structures,
published parameter bundles, active jobs and fixed MD charges stay unchanged.

Spec: DFT_REFERENCE_FIXED_CHARGES_20261104.md and the user's explicit request
to validate the QM reference while preserving Amber/GAFF/protein consistency.

## Scientific design

- Existing GO controls: pristine, OH, OH/EPX, carbonyl, COOH/COO on C16, plus
  two mixed C24 contexts. Freeze accepted coordinates, Q and trial spin.
- DFT comparison: current ωB97X-D/Pople control, and ωB97X-D, PBE0 and PW6B95
  with aug-cc-pVTZ. PW6B95/aug-cc-pVDZ tests a cheaper RESP2-inspired basis.
- Small chemical anchors: benzene, phenol, benzene oxide, p-benzoquinone,
  acetic acid and acetate. Accepted geometry/frequency/stability comes first.
- On anchors, compare MP2 and CCSD generalized densities with explicit
  Density=Current. Do not interpret an SCF density as a correlated result.
  Screen the underlying HF determinant and record single-reference diagnostics.
- For benzene oxide and acetate, test CCSD(T) finite-field dipoles using two
  field strengths. Calibrate Gaussian's field sign against an analytic CCSD
  dipole on water. Do not claim a CCSD(T) density from a CCSD density.
- A water runtime probe establishes syntax/property extraction, not GO accuracy.
- Compare potentials on matched coordinates/points, vector dipoles with a
  common origin, absolute errors and regularized relative errors. Separate
  method, basis and size effects; no winner selected by agreement alone.

## Implementation and checks

- [x] Add failing tests for input parity, fixed geometry, correlated density,
  explicit field units/sign and rejection of unsupported methods.
- [x] Implement go_param.reference_qm; check output gates against actual probe
  logs, including negative controls for missing stability/density evidence.
- [x] Pass the runtime probe, freeze manifests, submit isolated GO and anchor
  arrays with five-day limits and bounded concurrency.
- [x] Record jobs, failures and acceptance separately; produce tables and mirror
  source/docs. Submission is not reference validation.
- [ ] Retrieve accepted method comparisons, diagnose failed/unstable cases and
  select the independent QM reference. This remains scientifically pending.

## Review focus

RHF instability; wrong density selected by default; inconsistent dipole frame
or field sign; unmatched grids; incomplete correlated convergence; finite edges
and neutral/ionic states; confusing electronic accuracy with a consistent
Amber charge/LJ/water model. No new charge set is installed by this validation.

Literature/database audit is an additional input: distinguish Amber parameter
libraries, structural/DFTB databases, and DFT energy/force datasets from actual
electrostatic reference data. No large dataset or code is adopted unexamined.

## Executed status

- Gaussian runtime job 1469071: 21 intended positive cases pass, one deliberate
  SCF-density negative control rejected. CCSD analytic vs finite-field dipole
  differs by 0.000034 D; field-step difference 0.0000064 D. Water is a runtime
  control, not validation of GO chemistry or a CCSD(T) ground-state guarantee.
- End-to-end job 1469072: three accepted density/property outputs, native Amber
  espgen, 17,719 matching grid points per method and all three pairwise comparisons.
- GO controls array 1469073: 8 existing fixed geometries × 6 levels (HF/6-31G*,
  current DFT/Pople plus four DFT/basis controls). First two tasks running.
- Chemical anchors array 1469076: 6 molecules × 9 property levels after geometry
  gates; two anchors also have gated CCSD(T) two-step finite-field controls.
- All real jobs use 5-day limits, 16 CPU/64 GB per task, maximum 2 tasks per
  array. Original whole-sheet/expanded RESP jobs continue unchanged.
- Early RHF→UHF instabilities are explicitly rejected as stable QM references;
  their raw conventional RHF outputs remain available as labelled Amber
  diagnostics, not silently replaced with broken-symmetry RESP charges.
- 238 local tests passed, including actual Gaussian output excerpts and negative
  stability, density, incomplete-log and atomic-geometry controls.

## Published structural pools

See `../references/go_prior_art_20261104/README.md` and its hashed catalogs.
Four small finite CSIRO candidates and 75 periodic GO-MACE frames are available
locally for curation, plus one large published GO-MACE sheet for structural
coverage review. None is yet a fitted or accepted electrostatic reference.
Parent-sheet and near-duplicate grouping must precede train/test splitting.
The reference-method controls already launched are not additional random GO
generation: they reuse fixed, previously accepted geometries and chemical anchors.
