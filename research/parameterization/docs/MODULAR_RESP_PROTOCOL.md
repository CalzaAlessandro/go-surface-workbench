# L0 modular GO charge protocol

This protocol derives **shared functional-group library charges**, with background scaffold carbon fixed at zero. It implements the user's intentional localization/scalability model: current block charges were averaged over finite fragments; neutral scaffold carbon is a modeling choice, not missing chemistry to be fitted across an entire sheet. L0 preserves that framework when comparing current charges, newly constrained RESP blocks and a separately labeled AM1-BCC projection.

The seven existing usable fragments support a **PILOT_LIMITED_CONTEXT** fit, not final transferability. The standalone `pyrene_epoxide` calculation failed an upstream QM bond gate and is excluded from training. EPX currently occurs only in the mixed diol/epoxide model, so it has one physical chemical context despite three rigid orientations. An RHF-instability diagnostic reference stays diagnostic; a successful constrained fit does not validate its electronic state or final parameters.

## Exact atom partition and indexed library scope

The helper uses stable model IDs, primary atom roles and intended graph connectivity. Secondary `scaffold_C` annotations on functional carbons never freeze those carbons. Every atom must belong to exactly one recognized block or the fixed background; uncovered, overlapping or malformed chemistry fails before fitting. Model and atomic formal charges must agree with the block decomposition.

| Indexed unit | L0 block definition | Target charge |
|---|---|---:|
| GGG | Primary `scaffold_C`, neutral three-carbon junction; fixed q = 0 | 0 e |
| GRT | Terminal edge C–H, shared fitted C/H nuisance charges | 0 e |
| HYD | Basal C–O–H, roles `hydroxyl_C/O/H` | 0 e |
| EPX | Basal C–C–O epoxide triangle, two equivalent `epoxide_C` and one `epoxide_O` | 0 e |
| EPH | Epoxide with one unsubstituted C and one CH; `EPH_C`, `EPH_CH`, `EPH_H`, `EPH_O`; its two carbon charges remain distinct | 0 e |
| EPT | Epoxide with two CH groups; shared `EPT_C`, shared `EPT_H`, and `EPT_O` | 0 e |
| CAB | Carbonyl C=O, roles `carbonyl_C/O` | 0 e |
| CBX | COOH **including its scaffold anchor C**, carbonyl C/O, hydroxyl O/H | 0 e |
| CBD | COO− **including its scaffold anchor C**, carboxyl C and equivalent two O atoms | −1 e |

All nine units are registered, but only observed roles are fitted. `unfitted_terminal_blocks` explicitly lists EPH/EPT for the present seven-model pilot, and `unobserved_library_units` also identifies other units absent from a smaller pilot. Registry values from the current libraries are reference data, never reported as newly fitted values. Current CAB is +0.581/−0.581 e; EPH uses C +0.053, CH +0.137, H +0.120 and O −0.310 e; EPT uses C +0.020, H +0.120 and O −0.280 e. EPH's unindexed orphan EPX raw subtable is not an additional fitted library unit. OFF connectivity flags are not interpreted as chemical bond orders.

COOH and COO anchors remain inside their respective blocks, as in the user's current libraries. GRT charges are fitted shared nuisance parameters; they are not zeroed or silently merged with substituted anchors. COO uses the intended explicit one-single/one-double bond graph while sharing both partial oxygen charges; equivalence is a charge constraint rather than an invented fractional QM bond.

## System RESP constraint construction

`build_problem` collapses atom charges to shared role variables and constructs one stoichiometric sum equation for each observed block type. The complete seven-model pilot has 17 fitted role variables, six independent group equations and 11 remaining dimensions. Primary background carbon contributes no fitted variable. Each group occurrence must have the same role stoichiometry; full coverage and block formal totals are checked before system RESP. Rank and augmented-rank checks reject dependent or infeasible constraints after sharing.

The primary optimizer is the installed Amber **RESP** executable, with the existing espgen-format datasets. No local least-squares or alternative optimizer replaces its hyperbolic fit. Generated cards use `iqopt=2` and an explicit all-zero `qin.zero`; scaffold atoms have `ivary=-1`. All molecule headers use Amber's −99 sentinel, omitting redundant whole-model constraints. One global atom-group charge equation per observed block is emitted, and within-/between-model role equivalences propagate it to every occurrence and orientation. The resulting model totals follow from neutral blocks plus each −1 COO block and are independently checked after fitting.

The input-card conventions are verified against the [Amber RESP source](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/resp.F). Molecular weights are squared by RESP, while hyperbolic restraints are applied to heavy-atom copies with `ihfree=1`. Exactly three rigid orientations per physical model are required. The input weight is √(mean grid count / orientation grid count), with qwt = 0.0005 unchanged. Dividing the complete objective by three gives equal model contributions through orientation-mean ESP mean-squared errors. It also leaves a physical heavy-atom-occurrence prior: a role present many times receives repeated restraints. This multiplicity is recorded explicitly; it is not one restraint per library variable. Literal 1/3 input weights would give the wrong balance.

The present campaign has no aliphatic CH₂/CH₃ stage-2 group. Stage 2 is therefore an explicit no-op; no second optimization is run or claimed. Adding such chemistry requires a reviewed stage-2 extension.

## Inputs, execution and evidence

The dataset JSON contains already validated atomic-reference ESP files, in bohr and atomic units, bound to the matching optimized model and stable atom order:

```json
{"datasets": [{
  "model_path": "/path/to/accepted_or_diagnostic_fit/model.json",
  "esp_paths": ["/path/orientation_001.esp", "/path/orientation_002.esp", "/path/orientation_003.esp"],
  "scope": "rhf_instability_diagnostic_control"
}]}
```

Existing atomic-reference processing uses system `espgen` to produce these files. The modular helper consumes that established format; it does not infer element identities from an element-free ESP file. Upstream Gaussian/espgen atom-identity validation remains required. Here atom count, finite coordinates, intended topology, full pair-distance preservation and proper rigid rotation against the supplied optimized model are checked. Reflections or changed physical geometries fail. Use the optimized atomic-fit `model.json`, not initial force-field coordinates.

```sh
PYTHONPATH=src python3 -m go_param.modular prepare --dataset dataset.json --output-dir modular_attempt_001
PYTHONPATH=src python3 -m go_param.modular run --output-dir modular_attempt_001 --resp resp
# Or perform both sequentially:
PYTHONPATH=src python3 -m go_param.modular fit --dataset dataset.json --output-dir modular_attempt_002 --resp resp
```

Preparation refuses to overwrite an existing directory. It retains input ESP snapshots/hashes, problem/constraint/weighting plans, cards and initializer. Execution records the actual system executable hash and command, requires finite complete charges and explicit convergence, verifies frozen background, all role equivalences, each group sum and every model total. It exports shared role JSON/CSV, per-occurrence neutrality/model checks, and current-versus-modular ESP metrics on the same snapshots. Six-decimal RESP output tolerances are propagated without renormalizing fitted charges.

Synthetic executable fixtures test only input/output plumbing and constraint validation. They do not establish scientific fits or successful real Amber execution; compute-node numerical QA is a separate gate. A helper success status remains `PILOT_LIMITED_CONTEXT`, with `final_parameter_validation=false`.

### Actual execution and external engine capacity

The two-model system-RESP pilot and the full seven-model external-RESP pilot
completed on 2026-10-30. The complete cards contain 600 atom copies and six
charge constraints. The installed RESP default capacity of 500 centers cannot
run those cards, despite having only 17 shared fitted roles after equivalencing.
Use an explicitly identified external executable with sufficient capacity via
`--resp`; no executable or Amber source bundle is part of this repository.

A calculation-only AmberClassic build at upstream commit
`76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002` increased maxq to 1024 and matched
the installed two-model pilot exactly at output precision. Inspection also
identified a multi-molecule statistics inconsistency: weighted ESP norm was
combined with unweighted cross/quadratic terms, yielding a negative reported
SSE despite finite converged charges. A reporting-only call correction uses
weighted matrices consistently; the full charge vector is unchanged. Independently
recomputed weighted SSE is 30.43361928 versus 30.434 in the corrected rounded
report. The framework's NaN rejection remains strict. Global metrics use the
original grids and actual weights rather than the program's last-molecule
SIGMA denominator. See `reports/RESP_external_*_QA_20261030.json` and
`reports/modular_RESP_seven_model_diagnostic_pilot.json`.

These execution checks do not remove the RHF diagnostic scope or establish
held-out transferability. The six-model AM1-BCC projection is also complete as
a separate comparator; CAB is absent because its native BCC evidence is incomplete.

## AM1-BCC comparator: explicit constrained block projection

`project_am1_bcc(problem, raw_by_model)` retains raw atom-resolved AM1-BCC evidence and solves a distinct constrained L2 projection into the same shared-role namespace, fixed-zero scaffold and block-sum equations. Each model contributes its mean squared atomic-charge displacement; fixed background changes contribute constant loss. The role-diagonal normal equations and the independent block equations have an analytic equality-constrained solution.

The method is labeled `AM1_BCC_constrained_L2_block_projection_NOT_native_AM1_BCC`. It is neither a native AM1-BCC result nor an ESP/RESP fit, and raw AM1-BCC values are never overwritten. Differences from the native atomic model, including background removal and group localization, are part of this explicitly derived comparator. JSON input maps each model ID to its original `atoms` list with `atom_id` and `charge_e`; source evidence can carry additional provenance unchanged.

```sh
PYTHONPATH=src python3 -m go_param.modular project-bcc --dataset dataset.json \
  --raw-bcc raw_by_model.json --output-dir bcc_L0_projection_001
```

## Distance-sensitive training and held-out validation

The existing generic mixed fragment is insufficient to establish transferable blocks. The next panel must vary **distance and environment** without immediately creating a new parameter class for each count:

- OH–OH, EPX–EPX and OH–EPX separation ladders, sampled by graph distance and 3D motif-center distance.
- Edge COO−–basal O separation ladders, separating electrostatic distance effects from finite-edge/cap effects.
- At least two PAH sizes, same-face and opposite-face basal arrangements, and comparable boundary conditions.
- Dedicated EPH/EPT terminal cases before claiming all indexed units are fitted.

Fit one global L0 block model on a declared training subset. Hold out entire physical fragments/separations/sizes, never merely a rotated copy of a training geometry. Compare held-out ESP errors, molecular moments and unconstrained atomic/motif charge responses versus both graph and 3D distance. Keep group identity, local motif geometry, protonation, total charge and method/state labels explicit. Charge variation within a motif need not imply varying its constrained net charge.

Only if held-out behavior demonstrates a reproducible local-environment failure should a small number of chemically interpretable nearest-neighbor classes be considered. Such an extension must retain group neutrality/COO− charge, test rank again, and show improvement on genuinely held-out contexts. L1/environment classes are not implemented or fitted in this L0 pilot.
