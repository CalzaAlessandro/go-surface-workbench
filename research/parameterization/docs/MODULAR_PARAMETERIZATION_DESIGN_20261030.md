# Transferable GO blocks: agreed target and implementation plan

User intent: rederive and test reusable GO charge/parameter blocks, including
mixed groups at different separations, terminal epoxides, carbonyls and edge
chemistry. The current background-zero approximation is deliberate. Preserve
old/current libraries and simulations; regenerate definitive material models
after parameter selection.

## Compared approaches

1. Primary L0: joint multi-model RESP with background scaffold C fixed to zero,
   neutral HYD/EPX/EPH/EPT/CAB/CBX/GRT blocks and CBD charge −1. Keep carboxyl
   anchors and edge caps inside their own blocks. Shared role charges across
   models and orientations enforce portability; do not average atom charges
   then silently redistribute charge afterwards.
2. Same L0 feasible space projected against the independently computed AM1-BCC
   controls, labelled a transferable projection rather than native AM1-BCC.
3. Optional later context refinement: a few chemically justified neighbor
   classes or neutral local bond-charge increments, only if held-out size,
   distance and face comparisons show L0 inadequate. Finite local rules must
   preserve assembled-system charge for arbitrary sheet size.

GAFF2-AIM is a pertinent comparator: its SI averages DDEC6 charges and volume-
derived LJ parameters within first-neighbor environment classes. It retains
GAFF2 bonded terms; it does not impose all unoxygenated carbons at zero or
supply the full finite edge/COOH/COO/terminal catalog required here.

## Execution checklist

- [x] Inspect indexed source library units and distinguish EPX/EPH/EPT chemistry.
- [ ] Run actual Amber constrained RESP pilot on existing HF diagnostic ESPs.
- [ ] Export L0 coefficients and compare current blocks in the same constraints.
- [ ] Build mixed-distance/face/size panel plus terminal and edge controls.
- [ ] Validate panel chemistry before staged QM submission.
- [ ] Obtain larger-patch/QM references and held-out estimates of transfer error.
- [ ] Fit/refine ring-context-specific bonded terms using accepted Hessians and
      targeted scans; verify 1–4 and LJ conventions separately.
- [ ] Select charges and full parameters together using matched interactions
      and aqueous/adsorption benchmarks; then rebuild the four materials.

The current seven-model pilot lacks standalone accepted EPX and terminal
EPH/EPT contexts. It cannot establish final transferability. All raw atomic
fits remain diagnostic controls, not residue libraries.
