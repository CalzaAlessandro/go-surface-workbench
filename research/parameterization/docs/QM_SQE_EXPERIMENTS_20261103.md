# Whole-sheet QM and RESP-calibrated SQE experiments

The user authorized these experiments on 2026-11-03, including independently
reconstructed sheets if electronic problems require them. Original structures,
libraries, calculations and published parameters remain immutable.

## Accepted direction and implementation choices

- Primary charge target: the 16 conventional HF/6-31G* Amber RESP references.
  Electronic-stability limitations remain recorded. DFT is a separate control.
- Whole-sheet calculations use original GOX935 (1007 atoms, Q=-16, 5859
  electrons) and REGOX980 (913 atoms, Q=-6, 4994 electrons). A non-Kekulizable
  neutral-carbon graph does not prohibit a coordinate-based unrestricted QM
  calculation. Test spin explicitly before changing chemical composition.
- System ORCA 6.1.0 is available on compute nodes; CP2K is not in the module
  catalog. Use unrestricted PBE/ma-def2-SVP, RI-J/def2/J, TightSCF and a CPCM
  water environment for the highly charged aqueous-material controls. This is
  not a replacement for the canonical gas-phase HF-RESP protocol.
- Start a doublet/quartet pair for GOX and singlet/triplet pair for REGOX;
  the latter includes an alpha/beta-symmetry-breaking initial rotation. These
  are trial determinants, not an exhaustive determination of the ground spin.
- Each production trial has a five-day Slurm limit. Verify the installed
  executable and syntax with a tiny system before submitting the full sheets.
- Never fix parity by arbitrarily changing Q. A future rebuilt sheet must
  report C/O, all functional-group counts, protonation, H/C and defects relative
  to its source, not just claim equivalence from oxygen percentage.

## SQE research model

Implement the published SQE initial-charge equations with a fixed RESP-derived
reference vector. Bond transfers conserve each connected component's charge.
Atom roles, not arbitrary Kekule bond orders or the localized O formal-charge
label, define chemistry. Both COO oxygens share coefficients; unequal local
fields may still produce unequal final charges.

First fit role electronegativities with fixed positive hardness, bond hardness
and Gaussian widths; test a small, declared sensitivity grid. The first model
uses the current block library as its explicitly labeled baseline. No NPA
coefficients are relabeled as RESP coefficients. Fitted charges are not obtained
by fitting to the constrained L0 library itself.

Use molecule-normalized and role-balanced reference residuals with ridge
regularization. Fit C16, evaluate the seven supported C24 contexts without
using them to select the default hyperparameters; exclude the unsupported EPH
and EPT roles from that transfer test. Then fit all16 as a separate candidate,
without claiming independent validation for unique COOH/COO contexts.

## Meaningful checks

- Numerical agreement with the native published SQE implementation on a small
  molecule, including the initial-charge term.
- Conservation of Q and component charges, bond orientation and atom-order
  invariance, rigid-coordinate invariance, chemical resonance invariance.
- Synthetic recovery of observable charges from a known response model.
- Positive-definite solves, finite outputs and explicit unsupported-role errors.
- Whole-sheet QM success requires converged SCF, finite energy, correct electron
  count and trial spin, and population data; process exit alone is insufficient.

## Presets and application modes

Published `gaff2-aim` and `opls-hiergo` are atomic force-field bundles: charges,
bonded/LJ terms, mixing rules and 1-4 scaling travel together. Missing chemistry
is an explicit unsupported case. GAFF2-AIM's current basal-only source cannot
silently acquire carboxylate/edge terms from GAFF or OPLS.

Keep the application modes distinct: fixed library, minimum-distance chemical
graph context, and static redistribution (SQE/Espaloma). Unimplemented or
uncalibrated modes are not exposed as production-ready. The material recipe
and dimensions remain independent of the selected complete force-field bundle.

## Execution plan

- [x] Check native QM runtime and trial input; submit traceable full-sheet jobs.
- [x] Add failing numerical/chemical invariance tests, implement the SQE core,
      verify against native controls, and run the bounded HF-RESP calibration.
- [x] Record complete-bundle preset policy and application-mode readiness,
      without inventing missing parameters or enabling unfitted methods.
- [x] Verify outputs, mirror research sources/docs, update HTML/status reports.

## Execution update, 2026-11-04

- Whole-sheet array 1468875 exhausted memory during MPI startup, before SCF.
  It provides no electronic conclusion. Inputs and failure logs are preserved.
- Recovery array 1468880 changes resources only: 32 MPI ranks, 1024 GB/node,
  24000 MB/rank, explicit `ciamician_cpu_caps` QoS, five days, two trials at a
  time. Both low-spin sheets have entered actual SCF iterations. Convergence,
  spin comparison and population analysis remain pending.
- SQE fit 1468879 completed. Seven C16 training contexts transfer to seven
  supported C24 contexts; EPH/EPT remain explicitly unsupported in that split.
  Default settings were declared before observing C24. SQE improves the current
  library but does not outperform the recalibrated fixed-block L0 on matched
  held-out HF ESP. All16 fitting is separately labelled in-sample.
- Expanded array 1468993 contains 55 unique starting contexts (27 C24, 28 C54),
  four simultaneous tasks, five days each. Functional-group counts include
  OH/EPX loading series, CAB, COOH/COO, multiple acids, mixed protonation and
  EPH/EPT controls. Near/far placement uses minimum C–C graph distance.
  One requested C24 far EPX2 case failed the bounded placement/matching search;
  no smaller count was substituted and no mathematical impossibility is claimed.
- Each new context first requires an accepted ωB97X-D/6-311++G(d,p) geometry,
  frequency and wavefunction-stability check; only then is conventional
  HF/6-31G* MK ESP fitted by external Amber RESP. Some references may fail this
  gate. Prepared inputs and scheduler success are not completed RESP results.

No Git operations or new bonded/frcmod fitting are required for this experiment.
