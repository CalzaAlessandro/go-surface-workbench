# GO parameterization study

Independent study of GO functional-group electrostatics and force fields.
Original GO_BUILDER, user libraries, OLD parameters and existing MD campaigns
are preserved. Gaussian and AmberTools are system installations.

## Latest verified results: 2026-11-04

The [QM reference validation and published-structure audit](reports/QM_REFERENCE_VALIDATION_20261104_attempt_001/REPORT.html)
records the completed Gaussian/density/finite-field probes, new method-comparison
arrays, and externally sourced GO/nanoflake candidates. HF/6-31G*–RESP remains
the conventional Amber baseline; no production charge or LJ parameter changed.
Arrays 1469073 (eight GO controls) and 1469076 (six chemical anchors) are running;
their submission is not a validated choice of reference. External structures are
catalogued with provenance but are not yet curated charge-fitting labels.
See the [validation protocol](docs/QM_REFERENCE_VALIDATION_PLAN_20261104.md) and
[source-pool audit](references/go_prior_art_20261104/README.md).

The [table-based QM/literature report](reports/DIRECT_QM_ESP_LITERATURE_20261104_attempt_001/REPORT.html)
contains the direct HF/DFT potential comparison, RESP/BCC mixing diagnostic and
native GAFF2–AIM parameter inventory. The user-facing report now uses tables.
The [reference recommendation](docs/DFT_REFERENCE_FIXED_CHARGES_20261104.md)
separates electronic validation from the effective fixed-charge model for water.
It does not introduce new calculations or production parameters.

The [current HTML report](reports/QM_SQE_UPDATE_20261104_attempt_001/REPORT.html)
records the completed first SQE calibration to molecular HF–RESP, the
matched C16→C24 transfer comparison, whole-sheet QM array 1468880 and the
55-context expanded RESP array 1468993. Both whole-sheet low-spin trials
have entered SCF; the first four expanded references are running Gaussian.
These are active calculations, not completed new charge references. Every
new production trial has a five-day limit.

The recalibrated modular L0 remains better than this first SQE on the matched
held-out HF potential metric. No final charge method is selected. See the
[SQE protocol and limitations](docs/SQE_HF_RESP_PILOT_20261104.md).

## Earlier verified controls, updated through 2026-11-03

The [Espaloma control](reports/ESPALOMA_CONTROL_20261103_attempt_001/REPORT.html)
completed on 16 graphs, with 12 fully covered native BCC comparisons. It uses
the published model with an explicit formal-charge adapter; the report
records the release API charge issue and source stereo-metadata warnings.
Charge-vector similarity is not force-field validation. See the
[charge/LJ/bonded coupling note](docs/CHARGES_BONDED_LJ_COUPLING.md).

The [HTML charge report](reports/CHARGE_COMPARISON_20261103_attempt_001/RESP_AMBER_REPORT.html)
contains the actual RESP charges, comparison with the current library,
matched HF/DFT controls and held-out-context scores. Sixteen HF/6-31G* RESP
controls on accepted DFT geometries and sixteen matching DFT–RESP controls
are complete. Three recovered C54 references bring the DFT panel to nineteen;
these are not 34 independent chemical contexts. Twelve native AM1-BCC
controls have full parameter coverage; four exclusions are explicit.
Four shared-block RESP L0 fits are complete, including C16-only training
evaluated on seven supported held-out C24 contexts. No production library
or force-field winner has been selected.

The C54 single-epoxide recovery and gated postprocessing are complete in the
[19-context update](reports/CHARGE_CONTROLS_UPDATE_20261103_attempt_001/REPORT.html).
The [campaign status](reports/CAMPAIGN_STATUS.json) distinguishes the current
observation from historical snapshots.

The mixed-level workflow selects the DFT parent with
`--method 'wB97XD/6-311++G(d,p)'` and the conventional HF ESP target with
`--esp-method 'HF/6-31G*'` in `02_fit_and_compare.py`. The subsequent fit uses
Amber RESP. An ESP error metric evaluates the fitted RESP charges, not a
separate Gaussian charge set.

## Original study scope (historical snapshot through 2026-10-30)

Eight explicit finite pyrene-based models cover pristine scaffold, epoxide,
paired OH (cis/trans), paired carbonyl, COOH, COO− and mixed diol/epoxide.
Initial coordinates are force-field geometries, not QM results. The initial
HF/6-31G* branch remains diagnostic because of electronic instability. An
explicit ωB97X-D/6-311++G(d,p) branch supplies accepted-parent DFT–RESP controls,
with geometry, frequency and stability gates before three rigid MK
orientations, atomic RESP and matched-grid charge comparisons. The separate
limited-context shared-block pilot is still diagnostic. Held-out validation
and full force-field comparison remain pending.

The pristine pilot has a restricted-surface minimum with positive frequencies
but an RHF→UHF instability. Strict acceptance remains false. Its RHF-RESP path
is an explicit diagnostic control; a separate fixed-geometry UHF test is used
for interpretation. Lower broken-symmetry energy alone will not select final
charges or establish an exact singlet state. See docs/QM_GATES.md.

## Reproduce the tools

Use `uv sync --group dev`, then `uv run pytest -q`. Gaussian, espgen and resp
must be available externally. The isolated cluster study is recorded in
campaign.json; actual receipts and artifact gates determine execution status.

Inspect a reference with:

```bash
PYTHONPATH=src uv run python -m go_param.qm inspect models/pyrene_pristine/model.json /path/to/gaussian.log
```

Run the accepted-reference postprocessing locally, in an initialized molecular
software environment, using a fresh output directory:

```bash
uv run python 02_fit_and_compare.py --method 'wB97XD/6-311++G(d,p)' --model models/pyrene_diol_trans/model.json --log /path/to/accepted_DFT_parent.log --output-dir results/diol_trans/attempt_001
```

Append `--slurm --partition ciamician_cpu --modules gaussian/16 amber/25`
to submit the same scientific worker through Slurm. Module names and partition
are site settings. `--dependency afterany:JOBID` adds a scheduler dependency;
the reference artifact gate still runs. `--rhf-diagnostic` is an explicit
control for the exact documented RHF→UHF instability only; it never makes the
reference accepted or creates a final-parameter validation claim.

The selected Slurm Python must have NumPy; the site-specific campaign workers
probe Amber's Python before other runtimes. Always inspect completion/failure
JSON and scientific evidence, not just queue disappearance or zero exit.

## Outputs and limits

- `reference_gate.json`: convergence, frequency, identity, topology-distance and stability evidence.
- `orientations/`: immutable prepared geometries, inputs, actual logs and hashes.
- `fits/attempt_*/`: Amber RESP inputs, convergence logs, charges and ESP metrics.
- `charge_comparison.csv/json`: shared-grid charge-only controls with coverage/status.
- `POSTPROCESS_COMPLETE.json` or `POSTPROCESS_DIAGNOSTIC_COMPLETE.json`: distinct scopes.

No complete scientific force-field winner has been selected. OLD whole-sheet
and published periodic DDEC6 models are retained as references where exact
fragment mappings do not exist. Original GOPY charge sums are reported without
renormalization. OPLS charge-only controls do not constitute the full OPLS model.

See docs/PROTOCOL.md, docs/RESP_PROTOCOL.md, docs/COMPARISON_PROTOCOL.md and
docs/FORCE_FIELD_CANDIDATES.md for the scientific choices and primary sources.

## Transferable blocks and chemical context

The [assumption assessment](docs/ASSUMPTIONS_LITERATURE_ASSESSMENT.md) separates
physical total-charge requirements from block-neutrality and zero-background
approximations. The [L0 protocol](docs/MODULAR_RESP_PROTOCOL.md) fits shared roles
using external Amber RESP, retaining current charges as a matched control.
The [distance panel](docs/MIXED_DISTANCE_PANEL.md) contains 34 initial paired
configurations and five edge controls on two PAH sizes. Its eight-model first
QM subset was submitted as job 1401678; submission is not scientific completion.
The [literature group inventory](docs/MISSING_GROUPS_LITERATURE.md) distinguishes
ethers, phenolic OH, lactols, lactones and other conditional additions. These
groups have not been silently added to current libraries or material recipes.
The [context-dependent charge proposal](docs/CONTEXT_DEPENDENT_CHARGES.md)
describes tabulated and fitted local corrections; it is not implemented.

The eight original DFT geometry/Hessian calculations have finished. Seven pass
the preliminary identity, connectivity-distance, positive-frequency and
wavefunction-stability checks. The standalone pyrene epoxide changes its
intended C–C connectivity and remains excluded. These are references for later
bonded fitting, not fitted force constants.

The seven initial DFT–RESP controls and their independent postanalysis have
completed: see reports/RESP_context_DFT_seven_controls_20261030_attempt_001.
Seven HF probes on their fixed DFT geometries completed Gaussian normally
but reported external RHF-to-UHF instability; that additional strict gate
does not imply SCF failure or automatic Amber incompatibility. See
docs/DFT_LEVEL_AND_RHF_STABILITY.md and docs/AMBER_RESP_ACCEPTANCE.md.
The near-OH mixed reference timed out at the original eight-hour limit and
is under isolated checkpoint-based recovery as job 1408584, with five days.
