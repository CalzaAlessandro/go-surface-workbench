# GO parameterization and force-field comparison

User authorization, 29 October 2026: derive new GO functional-group parameters,
compare current modular HF/6-31G* RESP, historical embedded QM/MM, pertinent
published GO charge models and AMBER-family/GAFF/GAFF2/OPLS models. Preserve
existing campaigns and GO_BUILDER. The user delegated method choices.

## Molecular design

Use finite, explicitly H-terminated PAH patches. The initial training family
uses pyrene C16H10: pristine, basal epoxide, cis and trans paired basal OH,
paired edge carbonyl (quinone), COOH, COO− and mixed diol/epoxide. Each has
explicit formal charge, even electron count, no radicals and fixed atom IDs.
A single OH added to neutral pristine PAH is an odd-electron species; paired
OH avoids silently fitting that structure as a closed-shell singlet. Paired
carbonyls similarly provide a quinone reference. Compare larger patches as
held-out transferability tests, rather than assuming one size is sufficient.

## New electrostatic models

1. Optimize reference conformers with Gaussian HF/6-31G*, check termination,
   molecular identity, charge/spin and frequencies. Compute MK ESP with NoSymm.
2. Use additional rigid orientations of each accepted geometry to reduce grid
   orientation bias. Orientations are not additional physical conformers.
3. Derive a two-stage RESP molecular reference with explicit graph symmetry and
   integer total charge, equal total weight per distinct conformer.
4. Fit a separately labelled transferable block model jointly across molecules:
   constrain neutral/charged building-block totals and chemically equivalent
   roles, and account for removed caps during fitting. Compare a local zero-
   scaffold approximation with a model retaining first-shell polarization.
   Never zero fitted atoms afterwards without refitting/quantifying ESP loss.
5. Assess on held-out patterns/sizes and molecule–water/adsorbate interactions.
   A fit to the training ESP alone does not establish improved adsorption.

## Force-field separation

Compare electrostatics at common geometries first. Then generate distinct full
parameter sets: GAFF and GAFF2 share the new RESP candidate where applicable,
while their bonded/LJ parameters remain distinct; the user's current GAFF2,
GOPY-AMBER analogy, HierGO-OPLS and GO-specific published models retain their
own documented definitions. OPLS charges alone do not constitute an OPLS model.
GAFF2-AIM is a published alternative with DDEC6-derived charges and modified
nonbonded terms; do not label arbitrary charge averages as the original model.

Protein force fields and grafted amino acids/PEG require their own consistent
covalent junction templates. No unvalidated parameter set replaces active MD.

## Evidence and execution

Gaussian and AmberTools are external system software. Run QM on Slurm CPU
nodes, in an isolated campaign directory. Preserve every attempt, input/log,
hash and job receipt. Separate prepared, submitted, running, failed and verified
results. Quantitative comparisons require matching chemistry/protonation;
historical whole-sheet charges are descriptive where no exact mapping exists.

Primary methodological reference:
[R.E.D. multiple-molecule/conformation fitting](https://pmc.ncbi.nlm.nih.gov/articles/PMC2918240/).
GO-specific references and source hashes are catalogued under references/.
