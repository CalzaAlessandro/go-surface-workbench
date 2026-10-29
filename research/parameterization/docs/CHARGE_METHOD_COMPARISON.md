# Comparing computed charge methods

The campaign compares actual per-atom charges from distinct methods, their motif sums and their sensitivity to finite-flake environment. A method must first complete its own calculation and pass atom mapping, chemistry, convergence and integer-charge checks. Published motif tables remain separate controls. None of these comparisons alone validates a complete force field.

| Method | What is calculated | Required evidence and limits |
| --- | --- | --- |
| HF/6-31G* RESP | Atom-centered charges fitted to molecular electrostatic potentials with restraints and declared equivalences | Validated reference geometry/electronic state and actual Amber RESP fit; orientation dependence and buried-atom identifiability require attention. The pristine restricted-HF instability remains a diagnostic scope. |
| AM1-BCC | Semiempirical AM1 populations followed by parameterized bond-charge corrections | Actual AmberTools `antechamber -c bcc` and converged SQM calculation. These charges are not RESP charges. The AM1 optimization can change the starting HF geometry and is recorded separately. |
| CM5 | A mapped Hirshfeld population analysis of an electronic density | Requires an appropriate density calculation and explicit CM5 output; unavailable CM5 charges are not inferred from RESP or AM1-BCC. |
| DDEC6 | Electron-density partitioning designed for molecules and materials | Requires a suitable converged density and actual DDEC6 analysis, including source-model context; periodic sheet charges cannot simply be copied onto finite PAHs. |

RESP's original formulation is described by [Bayly et al.](https://doi.org/10.1021/j100142a004). AM1-BCC combines AM1 atomic charges with bond-charge corrections parameterized against HF/6-31G* ESP data: [Jakalian et al., method](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291096-987X%2820000130%2921%3A2%3C132%3A%3AAID-JCC5%3E3.0.CO%3B2-P) and [parameterization/validation](https://onlinelibrary.wiley.com/doi/10.1002/jcc.10128/abstract). CM5 is a separate population-based charge model, described by [Marenich et al.](https://pubs.acs.org/doi/10.1021/ct200866d). For DDEC6 see [Manz and Limas, part 1](https://pubs.rsc.org/lv/content/articlehtml/2016/ra/c6ra04656h), its [correction](https://pubs.rsc.org/en/content/articlelanding/2022/ra/d2ra90050e), and [part 2](https://doi.org/10.1039/C6RA05507A). These sources motivate distinct controls; they do not establish a winner for the present graphene-oxide fragments.

## Actual AM1-BCC worker

`go_param.alternatives` depends only on Python's standard library and installed AmberTools. It imports no RDKit or NumPy on the compute node. The same worker runs locally or under Slurm:

```sh
PYTHONPATH=src python -m go_param.alternatives \
  --model path/to/optimized/model.json --output-dir calculations/am1_bcc/model_id
# Add --slurm for the system amber/25 worker; this submits parameterization only.
```

The Python API is `run_am1_bcc(model_path, output_dir, *, antechamber='antechamber')`.

Each call creates a new `attempt_NNN`, snapshots the complete model, and writes a V2000 SDF with explicit graph, atom maps and formal charges. SDF coordinates use its standard four-decimal representation. The computational input is an equivalent named MOL2 with 12 coordinate decimals and element-compatible four-character names such as `C001`, `O017`, and `H027` (element plus the global atom ID, limited to 999 atoms) This deliberate format choice makes atom names an independent identity check: the [Amber MDL reader](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/mdl.c) reads atom symbols and does not preserve arbitrary SDF atom-map names. Five-character names such as `A0001` are rejected: actual AmberTools runs rewrote and truncated those names, causing duplicate output names. Exact names, uniqueness and mapping remain mandatory; previous failed attempts are preserved. No coordinate optimization or charge fitting occurs during this Python conversion.

The worker invokes standard `antechamber -fi mol2 -fo mol2 -c bcc -at gaff2 -nc <integer> -rn GOF`. Strict mode passes no SQM convergence overrides, charge rescaling, `-dr no`, or alternative charge procedure. Verbose native status (`-s 2`) exposes the actual AM1-BCC subcommands. The actual generated `sqm.in`, `sqm.out`, MOL2, subprocess log, executable digests and all scratch artifacts remain in the attempt. AM1 Hamiltonian, SQM geometry convergence and completion must be explicit. Final SQM coordinates must preserve atom ordering, elements and bond-distance chemistry. Exported MOL2 additionally must preserve names, atom IDs, connectivity, nonaromatic bond orders, finite coordinates/charges and integer total within recorded output precision. Aromatic `ar` versus Kekulé single/double representations are recorded explicitly.

A passing attempt exports `candidate.json` in the comparison module's candidate shape, plus `charges.json`, `charges.csv`, and `AM1_BCC_COMPLETE.json`. The candidate name is `new_am1_bcc_gaff2`. The original model's diagnostic scope remains in both candidate and completion evidence. A failure creates `AM1_BCC_FAILED.json`; previous attempts are never overwritten. Completion establishes the AM1-BCC charge control, not final bonded/Lennard-Jones parameters or transferability. Stub-based unit tests check transport and failure gates and are not real AM1-BCC validation.

## What to compare

Compare mapped atomic charges, symmetry-related charge spread, motif charge sums, total charge and dipoles with an explicit origin. Evaluate all charge sets on the same HF geometry and ESP grids as a common diagnostic; retain each method's actual charge-generation geometry separately. AM1 geometry changes are recorded rather than relabeled as an HF/RESP calculation. For charged models the dipole depends on the chosen origin.

Training-ESP residuals favor methods trained on that target and cannot alone select the best force field. Independent hydration/adsorption observables, solvent structure, energetics and held-out sizes or functional-group environments are necessary. Bonded terms, Lennard-Jones terms, combining rules and their compatibility with the charge model must also be assessed. An `frcmod` generated by analogy is a parameter-completeness artifact, not proof of those observables or of bulk-GO accuracy.


## Native precision and an explicitly derived comparison candidate

Actual AmberTools calculations exposed a precision boundary before the MOL2 writer: [SQM prints Mulliken populations using three decimal places](https://github.com/Amber-MD/AmberClassic/blob/main/src/sqm/qm2_print_charges.F90), and [antechamber reads these printed values before applying BCC](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/charge.c). Pairwise bond-charge corrections conserve their input sum. Consequently, six-decimal MOL2 output can still have a raw integer-charge discrepancy from the earlier three-decimal populations. The worker's default `--charge-policy strict` continues to reject that discrepancy beyond ordinary six-decimal output tolerance.

For an explicitly authorized *derived* comparison, `--charge-policy constrained-rounding` adds a distinct candidate called `new_am1_bcc_charge_constrained_vsraw`. The full unmodified native charges and sum remain in `raw_am1_bcc.json` and `am1_bcc.mol2`. Before deriving anything, the worker verifies the SQM population atom map and three-decimal precision, preservation of its total in the pre-BCC AC file, and agreement of the final MOL2 with post-BCC AC charges. It snapshots the active system `BCCPARM.DAT` and independently verifies every native BCC increment against the actual typed bonds and parameters; absent BCC evidence or parameter coverage is a hard failure. Zero increments are allowed only when explained by the active parameters. The pristine model's native C16–H91 and C16–C16 correction entries are zero, so unchanged pre/post charges alone do not imply a missing BCC stage.

Only a residual bounded by N × 0.00051 e is eligible for this policy. The derived values subtract the same residual/N from every atom: the minimum sum-of-squared-changes correction under the integer-total constraint. This preserves existing equalities. Raw and corrected sums, each atomic correction, maximum change and precision bound are recorded. It is not called unmodified standard AM1-BCC, and the policy is never applied to published GOPY/OLD/reference tables. On the eight inspected native attempts, the resulting maximum atomic changes would range from zero to approximately 0.000077 e; this numerical audit does not by itself establish final parameter validity.

## Resonance representations and fixed chemistry

Aromatic metadata can be prepared locally with RDKit before node execution:

```python
from go_param.alternatives import prepare_resonance_model
prepare_resonance_model("original/model.json", "NEW/model_with_resonance.json")
```

The new model retains all IDs, coordinates, formal charges and bonds and adds RDKit aromatic-edge metadata tied to an exact graph digest. Nodes remain standard-library-only. The worker permits `ar`/Kekulé representations only on those verified aromatic edges or edges already explicitly aromatic in the input. Numeric single/double swaps elsewhere must form closed, even, alternating carbon cycles of at least six atoms and conserve the bond valence at every affected atom. This handles the observed 14-carbon conjugated perimeter without allowing isolated sigma or functional-group changes.

The separate carboxylate exception requires two terminal oxygens on the same three-coordinate carbon, their declared resonance equivalence, the expected localized single/double starting graph, a −1 formal charge across the two oxygens, and compatible oxygen types. Amber's delocalized representation can write both C–O bonds as single while native AC files use delocalized bond code 9; these changes are recorded. COOH, carbonyl, hydroxyl, epoxide, connectivity and atom-identity changes do not receive a general bond-order exemption. All geometry and native-BCC evidence gates still apply. A single reproducible native pilot must pass before launching another full set of cases.
