# Explicit finite-model RESP methods

`go_param.qm`, `go_param.resp`, and `go_param.pipeline` accept exactly two explicit
method strings. `HF/6-31G*` remains the default. `wB97XD/6-311++G(d,p)` selects a
separately labeled DFT-RESP control. Selecting the DFT method does not reinterpret
HF evidence or establish force-field transferability.

The selected method applies to parent optimization/frequency/stability evidence,
the optimized model metadata, each of the three rigid MK single points, atomic
RESP input validation, and fit provenance. DFT-generated inputs additionally use
`Integral=UltraFine`. Both methods retain `SCF=(Tight,XQC)`, disabled symmetry, and
the existing MK grid IOp settings.

Parent acceptance still requires successful optimization, complete nonlinear
3N−6 frequency evidence, no negative frequencies (small negative modes require
manual review), final finite energy, explicitly stable wavefunction, correct
charge/multiplicity and atom order, and the original topology distance screen.
All Gaussian routes must match the selected method/basis, and all observed
`SCF Done` energy labels must match its electronic method. DFT energies are
reported as `final_scf_energy_hartree` with `final_scf_energy_method`; they are
never relabeled as HF energies. The legacy `final_hf_energy_hartree` field and
HF geometry status remain available.

Gaussian's internal checkpoint frequency route may inherit symmetry only from
a preceding explicit Opt/Freq route that independently passes method and
symmetry checks. It must contain `Geom=AllCheck Guess=TCheck SCRF=Check GenChk`
and the matching method/basis. The observed HF signature additionally requires
`Test`; the verified wB97XD internal route omits that keyword. This exception
does not apply to an explicit stability job or an arbitrary frequency route.

The `--rhf-diagnostic` option and `allow_rhf_instability_diagnostic=True` remain
HF-only. Combining either with the DFT method raises an error, including for a
stable DFT parent. DFT parent rejection cannot use the RHF instability bypass.

## Command-line usage

Run the strict evidence gate before any calculation:

```sh
python -m go_param.qm inspect models/pyrene_pristine/model.json parent.log \
  --method 'wB97XD/6-311++G(d,p)'
```

Prepare three MK inputs without running Gaussian:

```sh
python -m go_param.qm prepare-orientations \
  models/pyrene_pristine/model.json parent.log prepared_dft_resp \
  --method 'wB97XD/6-311++G(d,p)' --cpus 8 --mem-gb 6
```

The controller runs the parent gate, generates/runs the three MK single points,
fits with external AmberTools `espgen`/`resp`, and compares on the training ESP:

```sh
python -m go_param.pipeline \
  --model models/pyrene_pristine/model.json --log parent.log \
  --output-dir dft_resp_control_attempt_001 \
  --method 'wB97XD/6-311++G(d,p)' --cpus 8 --mem-gb 6
```

Adding `--slurm` submits the same controller with the method preserved in the
worker command. `--partition`, `--dependency`, `--gaussian`, and `--modules`
retain their existing meanings. Output directories must be new.

For already generated MK logs, the atomic fitting entry point is:

```sh
python -m go_param.resp --model prepared_dft_resp/model.json \
  --logs orientation_01.log orientation_02.log orientation_03.log \
  --output-dir fits --method 'wB97XD/6-311++G(d,p)'
```

This standalone fitter does not approve parent optimization, frequencies, or
stability: upstream acceptance remains mandatory. Every supplied MK route must
match the selected method. An optimized model's `metadata.qm_method`, when
present, must also match; legacy models without that field remain supported.

## Python API

All existing calls preserve their HF default. The new keyword is `method`:

```python
from go_param.qm import gaussian_input, inspect_reference, prepare_orientations
from go_param.resp import fit_atomic
from go_param.pipeline import run

method = 'wB97XD/6-311++G(d,p)'
text = gaussian_input(model, method=method)
gate = inspect_reference(model, parent_log, method=method)
manifest = prepare_orientations(model, parent_log, fresh_output, method=method)
# Invoke only with completed MK logs and an upstream-approved parent:
fit = fit_atomic(optimized_model_path, mk_logs, fit_output, method=method)
# Alternatively, run all gated downstream stages:
result = run(model_path, parent_log_path, fresh_pipeline_output, method=method)
```

RESP weighting and priors remain unchanged: one physical conformer, equalized
orientation ESP weight using the original point-count rule, shared charges,
stage-1 qwt 0.0005, stage-2 qwt 0.001, ihfree 1, and the original stage-2 CH₂/CH₃
policy. No NPA fallback or alternative numerical fitter is introduced.

DFT orientation/controller scopes are `dft_resp_control`; completed controller
status is `DFT_RESP_CONTROL_AND_TRAINING_ESP_CONTROLS_COMPLETE`. Fit protocol,
provenance, each snapshotted MK log entry, and training-comparison JSON/CSV
identify the QM method. `final_parameter_validation` remains false. A completed
atomic fit and training-grid comparison do not validate modular transfer,
held-out adsorption, sampling, or a final force field.

The automated tests use synthetic Gaussian text and explicit executable
transport fixtures. They verify gates and data transport, not QM or RESP
scientific results.
