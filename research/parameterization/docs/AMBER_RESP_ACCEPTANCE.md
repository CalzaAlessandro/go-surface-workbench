# Acceptance policy for the new Amber-oriented RESP preset

The user requires scientific consistency with the intended Amber force field,
in addition to readable Amber topology/library files. This policy applies to
the new RESP-derived preset. Published GAFF2–AIM and OPLS source presets retain
their own charge/parameter definitions unchanged.

## Canonical reference branch

- The reference charge target for the new conventional RESP/GAFF-family branch
  is a gas-phase HF/6-31G* Merz–Kollman electrostatic potential followed by
  Amber RESP hyperbolic fitting. See the [AmberTools manual, GAFF and Antechamber
  sections](https://ics.uci.edu/~sources/amber/AmberTools.pdf).
- Geometry optimization and ESP generation are separate protocol choices.
  An accepted DFT-optimized geometry may be used for a subsequent HF/6-31G*
  ESP calculation, provided both levels and the checks at each level are
  recorded. DFT geometry acceptance does not prove HF wavefunction stability.
- Require the declared charge, multiplicity, intended chemistry, atom mapping,
  geometry checks and electronically acceptable HF reference. An unstable
  restricted HF solution remains diagnostic. Do not silently change spin,
  protonation, functional-group chemistry or ESP level to obtain acceptance.
- Use external Amber RESP, retaining stage-1 qwt 0.0005 and unrestrained H,
  justified equivalences and the physical integer total. Apply the documented
  second-stage CH2/CH3 treatment where applicable; explicitly record a no-op
  for fragments without those groups. Keep orientation/conformer weighting
  and any additional modular constraints explicit.
- Neutral functional blocks and zero background-carbon charges are candidate
  transfer approximations, not requirements of the atomic reference fit.
  Measure their effect before choosing the library or redistribution model.

## Scope of the current calculations

The ongoing omegaB97X-D/6-311++G(d,p) ESP fits are DFT–RESP comparison controls.
They do not satisfy the HF ESP target above and cannot automatically become
the canonical new RESP preset. Preserve their completed outputs and accepted
geometries as independent evidence.

The current controller couples the parent reference method to the MK method.
A mixed-level HF-ESP-on-DFT-geometry branch still needs explicit implementation
and validation with separate geometry and ESP provenance. This policy does
not claim that branch has already run, or that a production selector enforces
these requirements.

## Production acceptance

Whole-model charge fits, executable regression tests and low training ESP
residuals are insufficient for production acceptance. Validate transfer across
sizes/environments and the coherent charge/LJ/bonded combination for the
intended solvent and protonation conditions. Preserve atom order and charge
totals in Amber libraries, topologies and MM/GBSA subsets. If a motif lacks an
acceptable canonical reference, retain an explicit unresolved status; adopting
an alternative charge model requires a separately identified and validated
charge/parameter branch.
