# Proposed context-dependent GO charge library

Status: scientific proposal, not implemented or fitted. The user's subsequent
clarification sets standard library charges as the default, with optional static
redistribution on the completed structure; see
[published method selection](STATIC_CHARGE_METHOD_SELECTION.md) and
[preset/external-input design](PRESETS_AND_EXTERNAL_PREPARATION_DESIGN.md).
Preserve the current
constant-block L0 model as a baseline. First establish the useful range of its
zero-background and block-neutrality assumptions using independent contexts.

## Alternatives and recommendation

1. A small table of charge vectors for explicitly defined graph environments
   is transparent but needs adequate examples in every class. Unseen classes
   remain unsupported rather than receiving invented parameters.
2. A regularized rule predicts corrections from graph environment and, only
   if justified, electrostatic descriptors. It shares evidence across classes
   but introduces assumptions about interpolation and response.
3. Updating charges throughout MD requires a coherent polarizable/fluctuating-
   charge energy and force model. It is a different force-field project.

Recommendation: begin with a few interpretable environment classes, then test
whether a low-complexity correction rule improves genuinely held-out cases.
The user has proposed tables or a fitted rule; neither alternative has yet
been selected as a production implementation.

## Chemical context and electrostatic embedding are distinct

The covalent neighborhood changes hybridization, conjugation and substitution.
Use explicit nearby groups in the QM region, varying separation, relative face,
edge/terminal context and PAH size. Relative face must be local to the sheet,
so folding a GO receptor does not change a label merely through global z.

Distant electrostatics can be investigated with controlled external point
charges or fields. Such an embedded-fragment reference includes environmental
polarization at the selected QM/MM boundary. A graph-only lookup table is a
surrogate for sampled chemical contexts; it is not a self-consistent embedding
calculation on every new sheet.

Any embedding study needs converged QM-region size and boundary/cap treatment.
Do not place nearby classical charges so that artificial boundary overpolarization
dominates the fit. If the user's current charges supply the outer environment,
vary the starting scheme and region size: reproducing that imposed environment
alone would be circular validation. An embedded RESP target must contain the
polarized QM contribution with external-charge potential removed/handled
consistently; otherwise fitting that potential again and retaining the MM charges
would double count it. This separation must be verified for the actual Gaussian
output before an embedded fit is accepted.

## Conserving charge without one arbitrary scaling factor

Use base role charges plus learned increments. Paired increments of +δ and −δ
on a selected bond conserve net charge by construction. Restricting them to a
neutral EPX block preserves C1 + C2 + O = 0 while allowing unequal carbon
corrections in asymmetric environments. A single multiplier of all EPX charges
cannot reproduce such asymmetric changes. A COO block should retain −1 e if
that block constraint is chosen; multiplying its entire vector changes that sum.

Initially constrain the total correction of each block to zero and leave generic
background carbons at zero, matching L0. These constraints exclude net charge
transfer between blocks and background polarization. Quantify the lost accuracy;
if evidence requires broader redistribution, test enlarged blocks or explicitly
charge-conserving inter-block increments as a separate model. Scaling is not a
way around those physical limitations.

Electrostatic-response terms should use fields or potential differences rather
than an arbitrary absolute potential reference. Neither their coefficients nor
a distance-decay formula are known yet. They must be calibrated from accepted
QM references, regularized and checked outside the training fragments.

## Amber implementation boundary

Assign the final per-atom vector during system preparation, record its chemical
context/source and retain it throughout ordinary fixed-charge MD and MM/GBSA.
Do not overwrite a global EPX.lib once per site: a single shared template cannot
represent simultaneous different environments. Use derived templates or assign
mapped per-atom topology charges in an auditable preparation step. Derive
complex/receptor/ligand topologies from that same assignment for MM/GBSA.

This gives environment-dependent *preparation*, not adaptive polarization during
the trajectory. Charges alone do not validate the selected LJ/bonded terms or
the treatment of solvation/adsorption. No new class is production-ready solely
because its charge sum passes.

## Primary precedents and limits

- [Pinto et al., GAFF2–AIM](https://doi.org/10.1021/acs.jpcb.2c06037): the
  inspected SI, printed S5, averages charge/LJ parameters in nearest-neighbor
  environment classes. It supports contextual, reusable parameters. Its DDEC6
  derivation and coupled LJ model are distinct from RESP and from the proposed
  local neutral-block increments.
- [Jakalian, Jack and Bayly, AM1-BCC II](https://doi.org/10.1002/jcc.10128):
  additive bond-charge corrections supplement molecular AM1 population charges.
  This is a precedent for conserving increments, not evidence that a correction
  to the current fixed GO blocks is native AM1-BCC or an embedded QM calculation.

The existing L0 constrained-HF pilot remains electronically diagnostic. The
34 distance/face configurations plus five edge controls are initial graphs;
only a first eight-model DFT subset is currently running. No context-response
coefficient, validated charge table or fluctuating-charge MD method is claimed.
