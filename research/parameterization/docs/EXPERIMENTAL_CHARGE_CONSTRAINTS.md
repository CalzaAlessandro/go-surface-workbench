# Experimental constraints on GO charge models

The project keeps RESP as its primary new charge derivation. Published
GAFF2–AIM and HierGO/OPLS remain separate, unchanged reference presets. New
bonded/frcmod fitting is deferred; the current priority is electrostatics and
coherent charge/Lennard-Jones comparisons.

## What is measured

Partial atomic charges are model-dependent partitions, not uniquely measured
numbers. The focused primary literature checked here provides acid/base,
electrokinetic and adsorption data, not a universal experimental EPX/OH/C
partial-charge library.

| Evidence | Proper role in this study | What it does not determine |
|---|---|---|
| Potentiometric titration | Dissociable-site populations and sample-specific protonation response versus pH | A unique partial charge on each carbon or oxygen |
| Zeta potential | Electrokinetic response under matched pH, ionic strength and dispersion conditions | Direct atomic charge or a unique microscopic surface-charge density |
| XPS/complementary spectroscopy | Composition and chemical-environment constraints, with assignment uncertainty | A unique atomic partial-charge vector from peak fitting |
| Adsorption ITC/equilibria | Combined charge/LJ/solvation/protonation behavior under matched experimental conditions | Independent values for a single charge or LJ coefficient |

## A directly relevant mixed-group experiment

[Tao et al. (2023)](https://doi.org/10.1021/acsomega.3c01596) combine four GO
sample titrations with DFT fragment hypotheses. Their Table 1 titration fits
give four apparent acidity ranges: 3.76–3.95, 5.66–6.36, 8.42–8.61 and
10.38–10.68. Proposed paired-group assignments are model interpretations;
the fragment pKa predictions in Figure 4 are not the same measurements.

For our design this motivates testing coupled neighboring groups and keeping
COOH/COO and other protonation states explicit. It does not establish those
constants or structures for LayerOne, nor determine our RESP coefficients.
Fit targets and independent validation observations must remain separate.

[Konkena and Vasudevan (2012)](https://doi.org/10.1021/jz300236w) combine
titration, electrokinetics and spectroscopy for GO/rGO; their GO response
includes apparent pK values 4.3, 6.6 and 9.0. These are sample- and
condition-dependent evidence of heterogeneous acidity, not a universal
deprotonation fraction for any constructed sheet.

## Linking these data to fixed-charge simulations

Choose and report a protonation model before comparing adsorption. A
fixed-protonation MD simulation does not itself predict a titration curve;
testing pKa numerically needs an appropriate protonation/free-energy protocol.
Changing LJ or intragroup partial charges cannot substitute for selecting the
right chemical state. Keep morphology, accessible surface, salt, buffer,
temperature and standard-state conventions explicit.

For charge/LJ validation start with several chemically distinct probes rather
than fitting one attraction. The primary ITC and water-uptake references in
[the nonbonded benchmark note](EXPERIMENTAL_NONBONDED_BENCHMARKS.md) supply
candidate observations. Before using numerical targets, verify the relevant
GO composition, site definition and preparation against the modeled material.
