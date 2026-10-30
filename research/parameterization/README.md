# GO parameterization study

Independent study of GO functional-group electrostatics and force fields.
Original GO_BUILDER, user libraries, OLD parameters and existing MD campaigns
are preserved. Gaussian and AmberTools are system installations.

## Current scope

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

RESP and AM1-BCC method definitions, raw/constrained rounding controls and
resonance checks are documented in [CHARGE_METHOD_COMPARISON.md](docs/CHARGE_METHOD_COMPARISON.md).
Three provisional GAFF/GAFF2/AMBER99SB-compatible bonded/LJ candidates and
their unresolved gaps are in [NEW_FRCMOD_CANDIDATES.md](docs/NEW_FRCMOD_CANDIDATES.md).

## Current experiment status

The real LEaP coverage test constructed 11 of 24 fragment/candidate topologies;
the remaining cases have unresolved terms. These source-selected baselines
are not complete production force fields. Seven HF atomic RESP controls remain
explicitly diagnostic because of RHF electronic instability. Seven AM1-BCC
models pass the current evidence gates; the dione has one bond without an explicit BCC parameter and remains a separately documented native-zero-fallback diagnostic.
All eight original DFT geometry/Hessian calculations finished; seven pass the
preliminary geometry, frequency and electronic-state gates. The standalone
pyrene epoxide changes its intended connectivity. No new bonded constants have
yet been fitted. An eight-model mixed-distance/edge subset is now running.

The machine-readable [campaign status](reports/CAMPAIGN_STATUS.json) records
the verification timestamp and exact scope. See [bonded QM protocol](docs/BONDED_QM_PROTOCOL.md)
and the input-only examples under `examples/dft_bonded_reference/`.

The new [assumption assessment](docs/ASSUMPTIONS_LITERATURE_ASSESSMENT.md),
[L0 shared-block protocol](docs/MODULAR_RESP_PROTOCOL.md),
[distance/face/size panel](docs/MIXED_DISTANCE_PANEL.md) and
[missing-group inventory](docs/MISSING_GROUPS_LITERATURE.md) separate current
evidence from further validation. Seven-model modular RESP remains an
RHF-instability diagnostic pilot; the AM1-BCC block projection is separately
labeled and is not native AM1-BCC.

The proposed distribution architecture has independent material recipes,
charge choices and compatible parameter packages, with an external-GO entry
path: see [preset/external-input design](docs/PRESETS_AND_EXTERNAL_PREPARATION_DESIGN.md).
An optional fixed-charge redistribution mode is a proposal, not implemented.
[SQE+qp method selection](docs/STATIC_CHARGE_METHOD_SELECTION.md) records the
published solver candidate, existing local implementation and NPA/RESP/GO
validation limits.

Current priority is RESP context data and coherent charge/Lennard-Jones
comparisons; deriving new bonded/frcmod coefficients is deferred. The
[experimental charge constraints](docs/EXPERIMENTAL_CHARGE_CONSTRAINTS.md) and
[nonbonded benchmarks](docs/EXPERIMENTAL_NONBONDED_BENCHMARKS.md) distinguish
measured titration/adsorption behavior from non-unique atomic partial charges.

Explicit [DFT-RESP support](docs/EXPLICIT_RESP_METHODS.md) preserves HF defaults
and labels the selected QM method throughout the gates and fit. The actual
seven-control DFT-RESP panel and independent postanalysis have completed.
The [complete context snapshot](reports/RESP_context_DFT_seven_controls_20261030_attempt_001/README.md)
contains seven physical contexts and 210 training ESP comparison rows.
The [context analyzer](docs/RESP_CONTEXT_ANALYSIS.md) exports individual atom
and group values versus scaffold graph distances, separately identifying
interior carbons and H-terminated edges. The 39-model starting panel now has
[four supplementary isolated controls](docs/ISOLATED_CONTROLS.md). Their
QM→RESP sequence is now submitted as array 1406363: the two C₂₄ controls are
running and the two C₅₄ controls are queued. They are not completed references.

The [interim four-structure report](reports/RESP_context_DFT_partial_20261030_attempt_001/README.md)
contains individually resolved atom/group charges and 120 locally reproduced
training ESP comparisons. It is a fixed snapshot of four completed structures,
so its count differs from the subsequent live completion count. The source
files, actual MK SCF evidence and charge totals were checked; this is not
transferable force-field validation.

Array 1406371 schedules mixed-panel RESP continuations with
[`aftercorr`](https://slurm.schedmd.com/sbatch.html#OPT_dependency), followed by
an independent scientific reference gate in each worker. New calculation
arrays and the previously pending mixed references have a five-day wall-time
limit. One far-OH parent and its RESP control completed. The near-OH parent
timed out at its original eight-hour limit; recovery job 1408584 uses a copied
checkpoint in a new folder, recomputes frequencies and then runs gated
DFT-RESP. Its obsolete RESP dependency task was canceled. The current
[live snapshot](reports/CHARGE_SCOPE_LIVE_STATUS_20261030_attempt_002.json)
records the separate original, recovered and downstream job states.

The additive `go_param.hf_fixed` helper tests fixed-geometry HF evidence after
an independently accepted DFT parent. All seven native probes completed
Gaussian normally, with external RHF-to-UHF instability; the additional
strict research gate therefore remains false. This is distinct from SCF
failure or automatic Amber incompatibility. The
[HF probe report](reports/HF_on_DFT_stability_20261030_attempt_001/README.md),
[DFT level and stability interpretation](docs/DFT_LEVEL_AND_RHF_STABILITY.md)
and [Amber RESP acceptance policy](docs/AMBER_RESP_ACCEPTANCE.md) preserve
that distinction. These probes did not fit new RESP charges.
