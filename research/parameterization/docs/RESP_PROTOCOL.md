# Atomic-reference RESP protocol

`go_param.resp.fit_atomic(model_path, gaussian_logs, output_dir, *, espgen='espgen', resp='resp')` fits the **whole finite model**. It uses installed AmberTools executables; no Python substitute for RESP or fitted-charge library is included. The CLI is:

```sh
PYTHONPATH=src python -m go_param.resp \
  --model models/pyrene_pristine/model.json \
  --logs orientation_001.log orientation_002.log orientation_003.log \
  --output-dir fits/pyrene_pristine
```

These filenames denote actual completed Gaussian logs supplied by the caller. Three orientations are the intended campaign protocol; the helper also accepts a single orientation for pilot QA. Changing the number is recorded explicitly.

## QM and atom identity

The caller must first approve the parent geometry and frequency evidence and
retain the electronic-state assessment. A standalone ESP calculation cannot
prove these parent properties. The initial strict research branch additionally
requires a stable wavefunction; conventional restricted-HF controls retain
their explicit instability status rather than being silently replaced with
UHF or DFT. See [AMBER_RESP_ACCEPTANCE.md](AMBER_RESP_ACCEPTANCE.md).
Each supplied log must end normally, request the declared ESP method
(HF/6-31G* by default), explicit symmetry disabling (`NoSymm` or `Symmetry=None`),
and MK ESP, and agree with model charge, multiplicity, element sequence, and
stable 1-based atom IDs. HF ESP on an accepted DFT geometry is supported through
separate ESP/geometry method arguments and explicit provenance. The last MK
atomic centers must match the last preceding Gaussian atom identity table.
Declared bonds and possible new short contacts are checked with explicit
C/H/O covalent-distance bounds; these checks do not constitute a reaction
or electronic-structure analysis.

All supplied geometries must have the same ordered pair distances and superpose by a proper rotation and translation within 0.00005 Å. They represent **one physical conformer**, not independent conformers. Duplicate log paths/content are rejected. Atom order remains the input contract; swaps of graph-symmetric, chemically equivalent atoms cannot be distinguished by geometry alone.

`espgen` headers retain a fixed I5 atom field and a minimum-width point-count
field: the latter expands beyond 99,999 points. The Python parser preserves
the complete count, checks exact atom/grid record counts and rejects unsupported
multipole records. Both native 64-column numeric records and the manual's
65-column records are supported. Native RESP builds with fixed `(2i5)` input
also need separately verified input-decoding support for expanded counts;
a Python parser change alone does not fix the external executable. The C54
calculation-only engine was checked against the installed pilot and retained
every original point. Atomic ESP coordinates must reproduce Gaussian MK centers
after conversion from bohr to Å. Grid coordinates and potential remain in atomic
units throughout the fit and independent metric calculation.

Freshly generated optimization, stability-check, and rigid MK inputs use `Symmetry=None`. The [Gaussian Symmetry documentation provided by CONFLEX](https://www.conflex.co.jp/gaussian_support/symmetry.php) distinguishes this from `NoSymm`: `NoSymm` preserves the input orientation but still attempts point-group identification, while `Symmetry=None` bypasses symmetry handling and point-group assignment entirely. This change addresses the observed Gaussian `Logic error in ASyTop` in a rounded optimized pilot geometry. It changes no coordinates, atom identities, rotations, or electronic-structure settings. Previously prepared attempt inputs and logs remain preserved; retries require new output directories. The route gates accept either explicit disabling form so earlier valid evidence remains readable. Successful numerical execution of fresh inputs must still be verified independently.

## RESP constraints

Stage 1 uses hyperbolic restraints, `qwt=0.0005`, and unrestrained hydrogen charges (`ihfree=1`). Declared symmetry classes are tied, including the carboxylate oxygen pair despite localized single/double input bond orders. Atoms reserved for stage 2 are independently fitted in stage 1. These settings and card formats follow the [RESP documentation](https://q4md-forcefieldtools.org/REDServer-Development/Documentation/resp/) and [Amber respgen implementation](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/respgen.c).

Stage 2 refits only tetrahedral carbon atoms bearing two or three hydrogens, plus those hydrogens, using `qwt=0.001`; every other atom is frozen to its stage-1 value. Stage-2 symmetry follows the declared classes. The current pyrene models have no such CH2/CH3 groups. They therefore retain stage-1 charges and report `noop_no_methyl_methylene`. Their all-frozen `stage2.in` is an audit artifact and is not executed. No stage-2 charge change is fabricated.

Independent charge representatives are shared across all orientations. One total-charge constraint is applied to the first orientation; subsequent charge headers use Amber's `-99` sentinel to omit redundant constraints. Equivalencing imposes the same physical integer total in each orientation. The exact [Amber RESP source](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/resp.F) applies the squared molecular weight to ESP residuals and a restraint to every atomic copy. Consequently, molecular weights are `sqrt(mean_grid_count / this_grid_count)` and the listed qwt values remain unchanged. Dividing the combined objective by the orientation count yields equally weighted mean-squared ESP errors and one physical set of atomic restraints. Weights are recorded at the F10.5 precision actually read by RESP. Repeating equally sampled orientations therefore does not strengthen the ESP relative to its restraints. This is a fit of one atomic reference, separate from future joint modular fitting. Orientation sensitivity motivates retaining multiple orientations, as discussed by [Dupradeau et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2918240/).

## Evidence and interpretation

Every invocation creates a new `attempt_NNN`; source model/log snapshots, executable/code digests, command lines, subprocess output, ESP files and RESP decks/results are retained. A failed attempt receives `FIT_FAILED.json`; only all passing gates produce `FIT_COMPLETE.json`. Gates require convergence evidence, finite charge files of the exact expected size, the whole-model integer charge, symmetry/orientation equivalence, and preservation of frozen stage-2 charges. Charge-sum tolerance accounts explicitly for RESP's six-decimal output.

`charges.csv` and `charges.json` preserve atom IDs, elements, names, roles and symmetry classes. Independent per-orientation metrics include ESP RMSE in a.u., RRMS, and the fitted point-charge dipole in Debye. Its origin is explicitly `(0, 0, 0)` in that Gaussian MK coordinate frame; ionic dipoles depend on the origin. These are point-charge dipoles, not an independent validation of a QM dipole. No arbitrary RRMS acceptance cutoff is imposed: metric assessment and out-of-sample scientific validation remain necessary before parameter publication or transfer to larger flakes.

The unit suite uses explicitly identified synthetic ESP data and executables that only copy predeclared fixture data. It verifies parsing, constraints, failure gates, preservation and metrics; it does **not** establish real Gaussian or AmberTools numerical validation. Real `module load amber/25` fitting QA on the compute host is a separate required campaign check. There are no topology, bonded-parameter, modular-transfer, or bulk-GO validation claims in a helper success marker.

Gaussian's combined `Opt Freq` job can print a separate internal frequency route such as `Geom=AllCheck Guess=TCheck SCRF=Check Test GenChk RHF/6-31G(d) Freq`, omitting the symmetry keyword inherited from the explicit parent route. The reference gate accepts this specific checkpoint-based internal signature only after an explicit HF/6-31G* Opt+Freq route with symmetry disabled. The internal route must independently retain HF/6-31G*; explicit optimization and stability routes still require their symmetry keyword. This parser compatibility preserves an unstable reference's rejection and does not turn it into an accepted reference.
