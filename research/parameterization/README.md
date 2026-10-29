# GO parameterization study

Independent study of GO functional-group electrostatics and force fields.
Original GO_BUILDER, user libraries, OLD parameters and existing MD campaigns
are preserved. Gaussian and AmberTools are system installations.

## Current scope

Eight explicit finite pyrene-based models cover pristine scaffold, epoxide,
paired OH (cis/trans), paired carbonyl, COOH, COO− and mixed diol/epoxide.
Initial coordinates are force-field geometries, not QM results. The campaign
runs HF/6-31G* optimization, frequencies and wavefunction stability before
three rigid-orientation MK single points, atomic RESP and matched-grid charge
controls. These are atomic reference fits; transferable block fitting,
held-out validation and full force-field comparison remain pending.

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
uv run python 02_fit_and_compare.py --model models/pyrene_epoxide/model.json --log /path/to/gaussian.log --output-dir results/epoxide/attempt_001
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
DFT geometry/Hessian references are running independently.

The machine-readable [campaign status](reports/CAMPAIGN_STATUS.json) records
the verification timestamp and exact scope. See [bonded QM protocol](docs/BONDED_QM_PROTOCOL.md)
and the input-only examples under `examples/dft_bonded_reference/`.
