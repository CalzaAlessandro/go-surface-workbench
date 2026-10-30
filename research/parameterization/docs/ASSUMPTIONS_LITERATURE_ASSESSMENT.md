# Assessment of the GO modeling assumptions

**The constraints make sense as a deliberately local, transferable fixed-charge model. They are hypotheses to validate, rather than consequences that the cited literature establishes exactly.** They simplify large-sheet assembly while suppressing some environment-dependent polarization.

| Assumption | Assessment | What establishes its useful range |
|---|---|---|
| **Exact whole-model charge** | Physically required for the specified composition/state: 0 for neutral models; −1 e for the selected carboxylate. | Check identity, protonation and total; internal partial-charge partitioning remains model-dependent. |
| **Each OH/epoxide/carbonyl/COOH block neutral; COO block −1** | A sensible assembly convention, not an independently conserved quantum charge for every subfragment. Stauffer's OH/epoxide examples support such a convention in a GO model [1]; they do not establish every block boundary or environment. | Test common boundary definitions, retained anchors and held-out ESP/moments. Neighboring blocks can redistribute density even when their imposed net charges remain fixed. |
| **Every background scaffold C exactly zero** | A strong, useful sparsity/localization approximation. “Most redistribution is local” does not imply that every other atom has exactly zero partial charge. | Compare constrained L0 predictions with appropriate less-restricted references across distances and sizes. No universal cutoff or acceptable error is asserted here. |
| **Shared role charges versus averaging fragment fits** | Both can compress many environments into a reusable library. A joint shared-role fit directly optimizes the stated common objective; post-hoc averaging need not produce that optimum. This does not establish that the earlier averaging was wrong. | Compare both on held-out fragments with matched background, caps and state labels. |
| **One charge per role across distances/faces/sizes** | A reasonable L0 starting point, with limited transferability until tested. Nearest-neighbor dependence has direct precedent in GAFF2–AIM [4]. | Hold out OH–OH, EPX–EPX, OH–EPX and COO–basal-O distance ladders, both faces and multiple PAH sizes. Rotations of one geometry are not independent chemical contexts. |
| **Finite termini and caps** | Necessary boundary choices, but their fitted charges can absorb finite-fragment effects. EPH/EPT are chemically distinct from basal EPX. | Include terminal cases and size tests; retain GRT and substituted-anchor roles separately. Missing terminal data remain explicitly unfitted. |
| **Fixed protonation at a chosen pH** | A conditional state model, not a prediction of pH-dependent ionization or charge regulation. A chosen COO fraction is an input assumption. | State solution conditions and compare plausible protonation scenarios; do not infer a universal deprotonation fraction from a successful charge fit. |
| **Charges combined with bonded/LJ parameters** | Must be assessed as a coherent model. Good ESP agreement alone does not validate ring stiffness, torsions, dispersion or adsorption thermodynamics. | Check actual topology/matching, geometry/Hessians, torsional scans and relevant hydration/adsorption observables with consistent 1–4 conventions. |

## What the thesis citations actually support

**[1] Stauffer et al. (2014), DOI [10.1063/1.4890503](https://doi.org/10.1063/1.4890503):** directly concerns GO atomic charges and amino-acid adsorption in explicit water. The archived primary Table I contains neutral OH/epoxide charge examples. This is a model precedent, not universal proof of exact-zero background charges. The full PDF timed out; Table I retains its indexed-primary provenance.

**[2] Smith et al. (2011), DOI [10.1021/jp202757p](https://pubs.acs.org/doi/abs/10.1021/jp202757p):** the inspected primary abstract studies bicyclo[1.1.1]pentane-carboxylic-acid induction using QTAIM/dipoles. This is broader transferability context, not GO-specific proof of vanishing scaffold charges.

**[3] Liu and Wilcox (2012), DOI [10.1021/es204071g](https://pubs.acs.org/doi/abs/10.1021/es204071g):** the inspected primary abstract concerns functionalized microporous carbons/CO₂ and Bader analysis. Its application and charge definition differ from aqueous GO RESP; it establishes no universal exact-zero rule.

**[4] Pinto et al., DOI [10.1021/acs.jpcb.2c06037](https://doi.org/10.1021/acs.jpcb.2c06037):** the archived primary SI explicitly averages parameters within nearest-neighbor atom-type classes, with different hydroxyl/epoxide environments (printed p. S5; `references/papers/Pinto2023_SI.txt`, lines 111–133). This supports environment-aware modularity, not required strict zeros.

## Current evidence and conclusion

The two-model and full seven-model constrained-HF pilots have succeeded technically, but remain **RHF-instability diagnostic controls**. The installed engine's 500-center limit was handled with a separately checked external build; it was an implementation limit, not a scientific rejection or validation of L0. See the dated capacity/reporting QA and modular pilot reports.

Retain L0 as a baseline; assess held-out predictions and physical observables before introducing justified environment classes. Neither convergence nor exact satisfaction of imposed constraints establishes the physics.
