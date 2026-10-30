# Transferable GO blocks: agreed target and implementation plan

User intent: rederive and test reusable GO charge/parameter blocks, including
mixed groups at different separations, terminal epoxides, carbonyls and edge
chemistry. The current background-zero approximation is deliberate. Preserve
old/current libraries and simulations; regenerate definitive material models
after parameter selection.

## Compared approaches

1. Primary reference program: atomic/multiconformer RESP on isolated or
   appropriately balanced groups, mixed combinations and separation ladders.
   Retain the atomic context response before selecting a transferable model.
   L0 candidate: joint multi-model RESP with background scaffold C fixed to zero,
   neutral HYD/EPX/EPH/EPT/CAB/CBX/GRT blocks and CBD charge −1. Keep carboxyl
   anchors and edge caps inside their own blocks. Shared role charges across
   models and orientations enforce portability; do not average atom charges
   then silently redistribute charge afterwards.
2. After observing the RESP context data, compare (a) distance/environment
   lookup tables, (b) a fitted interpretable correction rule, and (c) an
   appropriate published redistribution method recalibrated to RESP. Do not
   select SQE or any other response formula before those data support it.
3. AM1-BCC and its explicitly named L0 projection remain secondary controls,
   not the primary derivation or production charge preset. Published GAFF2–AIM
   and OPLS remain distinct exceptions to the new-RESP policy. Any local rule
   must preserve assembled-system charge for arbitrary sheet size.

GAFF2-AIM is a pertinent comparator: its SI averages DDEC6 charges and volume-
derived LJ parameters within first-neighbor environment classes. It retains
GAFF2 bonded terms; it does not impose all unoxygenated carbons at zero or
supply the full finite edge/COOH/COO/terminal catalog required here.

## Execution checklist

- [x] Inspect indexed source library units and distinguish EPX/EPH/EPT chemistry.
- [x] Run actual Amber constrained RESP pilot on existing HF diagnostic ESPs.
- [x] Export L0 coefficients and compare current blocks in the same constraints.
- [x] Build mixed-distance/face/size panel plus terminal and edge controls.
- [x] Validate panel chemistry before staged QM submission.
- [x] Add strict, explicitly labeled DFT-RESP support without changing HF defaults.
- [x] Register published GAFF2–AIM charge/LJ preset with unchanged source files.
- [x] Implement descriptive atomic/group RESP context analysis with scaffold
      graph distances and separate interior/edge-carbon reporting.
- [ ] Complete the staged DFT-RESP controls and the expanded mixed panel.
- [ ] Expand distance ladders to COOH/COO pairs and their combinations with
      basal OH/EPX and carbonyl environments; do not treat the present basal
      panel as complete coverage of the source library.
- [ ] Compare each functionalized parent with its corresponding pristine
      parent to distinguish oxidation-induced charge response from finite-edge
      effects. A single neutral basal OH needs an explicit open-shell reference;
      a remote OH partner is a balanced proxy and must be labeled as such.
- [ ] Obtain larger-patch/QM references and held-out estimates of transfer error.
- [ ] Deferred per the user's current scope: fit/refine ring-context-specific
      bonded terms using accepted Hessians and
      targeted scans; verify 1–4 and LJ conventions separately.
- [ ] Select charges and full parameters together using matched interactions
      and aqueous/adsorption benchmarks; then rebuild the four materials.

The current seven-model pilot lacks standalone accepted EPX and terminal
EPH/EPT contexts. It cannot establish final transferability. All raw atomic
fits remain diagnostic controls, not residue libraries.
