# Gaussian reference gates and MK orientation inputs

`go_param.qm` uses only the Python standard library. It prepares Gaussian text and evaluates explicit evidence; it does not run Gaussian, submit jobs, fit charges or establish that text came from an authentic calculation. The parser tests use labeled synthetic fixtures, never reported as QM results.

## Reference input

`gaussian_input(model_record, cpus=8, mem_gb=6)` retains the ordered atoms, explicit coordinates, formal charge and multiplicity. The first link requests:

```text
#P HF/6-31G* Opt=(Tight,CalcFC,MaxCycles=200) Freq Pop=MK SCF=(Tight,XQC) NoSymm IOp(6/33=2,6/41=10,6/42=17)
```

The second link reads the same `reference.chk`:

```text
#P HF/6-31G* Geom=AllCheck Guess=Read Stable NoSymm SCF=(Tight,XQC)
```

`Stable` checks the optimized reference wavefunction; it does not request `Stable=Opt` or silently change the reference state. An unstable wavefunction requires a separately reviewed protocol. Resource values must be positive integers. Element support is explicit, coordinates must be finite, atom IDs unique, and electron/multiplicity parity consistent.

## Inspection and acceptance

The controller API is `inspect_reference(model_path_or_record, log_path_or_text)`. Its result includes `accepted`, `issues`, `manual_review_required`, `optimized_atoms`, frequency and HF energy evidence, final orientation and topology screens. `optimized_model` is a complete model record only when accepted; otherwise it is null. `validate_gaussian_log` exposes the evidence gate without constructing that model; `parse_gaussian_log` extracts raw evidence without claiming acceptance.

All of these gates must pass:

- Explicit HF/6-31G* Opt/Freq and HF/6-31G* stability-check route evidence (the printed RHF/6-31G(d) spelling is equivalent). All relevant Opt/Freq/Stable routes must use that level; a DFT optimization followed by an HF stability check cannot pass. `Stable=Opt` routes are excluded from this reference gate.
- At least two explicit `Normal termination of Gaussian` records; no recognized Gaussian error termination, convergence failure, stopped optimization or write error.
- Explicit `Optimization completed.` evidence.
- A finite final `SCF Done` HF energy. RHF, UHF and ROHF labels are recognized; numerical D exponents are accepted.
- Exactly 3N − 6 finite frequencies from the last harmonic-frequency section, appropriate to the nonlinear fragment catalog. A frequency below −20 cm⁻¹ rejects the reference. Any frequency from −20 to below 0 cm⁻¹ requires manual review and prevents automatic acceptance. The script offers no silent manual override.
- Explicit stable-wavefunction evidence. A later instability message invalidates an earlier stable message.
- Explicit charge and multiplicity records, all matching the intended model.
- A final Input or Standard orientation with exactly N atoms, sequential Gaussian center numbers, and the expected ordered atomic numbers. Missing, malformed or nonfinite coordinates reject the reference. The model's stable atom IDs, names and roles are preserved when replacing coordinates.
- A distance screen against the intended model bonds. Intended bond lengths must fall between 0.65 and 1.30 times the sum of element covalent radii. Nonbonded pairs closer than 1.20 times that sum are flagged as novel close pairs. Lost bonds, excessively short bonds and novel close pairs prevent automatic acceptance. The radii and factors are included in the result.

The distance screen **does not validate bond order, aromaticity or stereochemistry**. It is deliberately a structural anomaly gate. Passing it does not establish a global minimum or validate the initial closed-shell ansatz beyond the specified stability evidence. The route/basis guard checks the printed route sections; the archived generated input and job records remain necessary provenance. The parser does not reconstruct every Gaussian keyword from an output archive. Do not apply this 3N − 6 frequency rule to linear molecules without an explicit protocol change.

A CLI inspection prints portable JSON and exits with status 0 only if accepted, otherwise 2:

```sh
PYTHONPATH=src python3 -m go_param.qm inspect models/pyrene_epoxide/model.json reference.log
```

Rejected nonfinite numerical evidence is printed as JSON null, with the relevant issue retained. The program does not infer missing frequencies, stability, geometry or charge from a successful process exit.

## Orientation preparation

`prepare_orientations(model_path, log_path, output_dir, cpus=8, mem_gb=6)` checks acceptance before creating outputs and refuses to overwrite an existing output directory. It writes:

- `model.json`: the accepted, unrotated optimized geometry and source-log hash.
- `validation.json`: evidence and gate results.
- `orientation_01` through `orientation_03`: each contains its rotated `model.json` and `mk.gjf` input.
- `manifest.json`: preparation status, source model/log and output hashes, optimized model path, and each orientation's input/model paths and SHA-256 values.

The three deterministic rotations use Euler rotations about laboratory x, then y, then z, with angle triples (17°, 31°, 47°), (−29°, 13°, 71°), and (43°, −37°, 19°), about the atom-coordinate centroid. They preserve all internal distances. Element-resolved coordinate-cloud checks reject a rotation that exactly reproduces the original or a previous orientation within 10⁻⁷ Å through a molecular symmetry. These are distinct grid-orientation controls for one rigid geometry, not independent conformers.

Each `mk.gjf` requests HF/6-31G* MK single-point evaluation with tight SCF/XQC, the same MK output options, and `NoSymm`. It requests no geometry optimization or frequency calculation. Source log hashes are computed from original file bytes when a path is supplied; model/input hashes likewise identify written content. Geometry/model records preserve atom ordering for downstream ESP/RESP mapping.

```sh
PYTHONPATH=src python3 -m go_param.qm prepare-orientations \
  models/pyrene_epoxide/model.json reference.log calculations/accepted_orientations \
  --cpus 8 --mem-gb 6
```

Manifest `orientations` entries include absolute `input_path` and `model_path`; the top-level `optimized_model_path` identifies the unrotated accepted model. The RESP helper can receive that optimized model and the three actual MK logs: it must validate ordered identities and rigid geometry equivalence across the rotations. Comparisons against a particular grid must use the matching orientation geometry or a verified coordinate transformation.

The manifest status `prepared_inputs_only_not_executed` describes prepared inputs. Completion of orientation jobs and verified ESP extraction are separate gates. No partial charge, RESP fit or electrostatic ranking follows from input preparation alone.

## Explicit constrained-RHF diagnostic branch

The default gate remains strict. `prepare_orientations(..., allow_rhf_instability_diagnostic=True)` (CLI `--rhf-diagnostic`) permits a separately labeled control only when the sole validation issue is `stability_not_confirmed` **and** the log explicitly contains the exact Gaussian statement `The wavefunction has an RHF -> UHF instability.` All optimization, frequency, identity, method/route, charge, geometry and topology gates must otherwise pass. An absent stability calculation or an arbitrary instability message cannot qualify.

The reference report retains `accepted: false`; `inspect_reference` continues to provide no accepted optimized model for it. The diagnostic preparation constructs geometry from the already inspected `optimized_atoms`, labels the model `wavefunction_stability_validated: false`, and gives the manifest/model scope `rhf_instability_diagnostic_control`. Manifest `reference_accepted` and `final_parameter_validation` remain false. The branch allows a reproducible conventional restricted-HF ESP/RESP sensitivity control; it does not validate the unstable reference or authorize final parameters.

An RHF→UHF instability is an electronic orbital-space instability when restricted spin constraints are relaxed. Positive nuclear vibrational frequencies on the restricted surface do not establish unrestricted electronic stability. Conversely, the triplet stability-vector ⟨S²⟩ printed beside an eigenvalue is not the final spin contamination of a separately optimized UHF determinant. A separate UHF singlet wavefunction stabilization test at the same geometry should report final energy, actual determinant ⟨S²⟩ and a final stability check before interpreting the electrostatic sensitivity. It must remain a separate protocol branch, with reoptimization/frequencies if its nuclear geometry is subsequently changed. A lower broken-symmetry UHF determinant alone does not establish the exact spin-pure ground state.

The Gaussian [Stable manual](https://www.chem.kuleuven.be/peeters/manuals/g03/g_ur/k_stable.htm) distinguishes testing from `Stable=Opt` wavefunction reoptimization and distinguishes internal instabilities from allowing RHF to become UHF. This primary manual copy is hosted by KU Leuven; the current Gaussian website was unavailable during this check. The restricted diagnostic and the UHF sensitivity branch must retain their separate method/state labels in fitted-charge and ESP reports.
