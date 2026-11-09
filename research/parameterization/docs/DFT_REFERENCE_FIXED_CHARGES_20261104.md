# A defensible reference for fixed-charge GO models

Recommendation for the research protocol, not a newly validated charge preset.
No new electronic calculations or MD parameters are introduced by this note.
Charges are assigned once during preparation and remain fixed in MD and GBSA.
The total is the integer dictated by chemistry/protonation, including nonzero Q.

## Separate the questions

| Question | Reference | What it establishes |
|---|---|---|
| Electronic fidelity | Direct QM electrostatic potential and dipole, with a checked method/basis/state | Whether the charge model reproduces chosen electronic properties |
| Point-charge representability | Individual RESP fit to that QM potential | A flexible fitting control; not independent physical validation |
| Transfer to new GO | Whole fragments, arrangements and sizes excluded from calibration | Generalization of a library, context rule or static SQE |
| Fixed-charge MD in water | Controlled prepolarization, coherent LJ/water model, independent interaction and condensed-phase tests | Whether the complete effective model serves the intended environment |
| Published presets | Original GAFF2–AIM or OPLS bundles in their supported chemistry and boundary conditions | Comparison with an established model, not a unique atomic-charge truth |

RESP is a fitting procedure. HF/6-31G* is the electronic target of the
conventional Amber-oriented branch; neither is a universal ground truth for GO.
The existing HF branch stays unchanged as a named baseline. A DFT-based branch
can also produce fixed charges for Amber, with separate parameter validation.

## Working reference and acceptance gate

The current gas-phase ωB97X-D/6-311++G(d,p) potential is a working DFT control.
Its ranking for GO electrostatics has not been established. Before selecting
a final reference, use a representative subset with pristine scaffold, OH,
epoxide, COOH/COO and mixed environments to check:

1. Declared chemistry, protonation, electron count, trial spin and wavefunction
   stability. Converged SCF alone does not identify the ground state.
2. Sensitivity to a second functional and an augmented basis. Functional names
   or a large basis alone do not establish accuracy. On suitably small,
   single-reference cases, use a higher-level correlated-property calculation
   where feasible; open-shell/static-correlation cases need separate scrutiny.
3. Direct potential and vector dipole agreement. Keep a common origin for
   charged-system dipoles. Report absolute potential errors alongside relative
   errors: monopoles and small reference potentials can distort relative scores.
4. Finite-edge, size and geometry effects. A finite cluster is not a periodic
   sheet. Hold out entire physical contexts; rigid rotations are not new ones.

The dependence of partial-charge definitions on basis and correlation is
documented by [Mehta and Martin](https://doi.org/10.1002/jcc.27294)
([author manuscript](https://arxiv.org/abs/2309.12184)). Their small-molecule
benchmark is evidence for method checks, not proof of an optimal GO functional.

## Fixed charges for water

I recommend comparing gas- and solvent-polarized QM references before choosing
the effective fixed-charge target. Full aqueous polarization is not automatically
the best fixed-charge approximation because its self-polarization cost is absent
from a simple Coulomb model. The published
[RESP2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8204736/)
combines gas/water-derived charges and validates them with LJ parameters in a
fixed-charge model. Its calibrated balance is not assumed transferable to GO.

This is a methodological precedent, not implementation of RESP2 here. Fitting
should target the solute electrostatics consistently; an already screened PB
solution is not an interchangeable target for vacuum-Coulomb atomic charges.
After derivation, the resulting library/SQE assignment is frozen in the topology.
There is no frame-dependent charge update or polarizable MD in this proposal.

## What the completed controls currently mean

We have evaluated direct HF and DFT potentials from existing calculations.
The small DFT residual of an individual DFT–RESP fit is expected: it is fitted
on those data. It cannot prove that the chosen functional is correct or that
its charges give the best aqueous MD. The current panel also does not establish
the optimal prepolarization or final charge/LJ combination. These remain explicit
gates before selecting a transferable production preset.
