<a id="cariche-lennard-jones-e-bonded-accoppiamento-e-protocollo"></a>
# Charges, Lennard-Jones and bonded terms: coupling and protocol

Updated on 3 November 2026. This note describes the choice of model;
it introduces no new GAFF, LJ or frcmod coefficients.

<a id="perché-sono-accoppiati"></a>
## Why the terms are coupled

The MM energy contains bonded, Coulomb and Lennard-Jones terms. Observables
depend on their sum. Charges and LJ jointly affect hydration, water
orientation, contact distances and adsorption. A change in one term can be
partly compensated by another: a single adsorption value does not uniquely
identify the charges and LJ parameters.

Intramolecular electrostatic and LJ contributions, including 1–4 interactions
and their conventions, also enter angle and torsion energy profiles. Bonded
coefficients derived in the context of a particular charge model can compensate
for part of its energy surface. This does not mean that every charge change
requires refitting all bonded terms; it requires checking their consistency,
especially for new covalent junctions.

The study [Links between the charge model and bonded parameter force constants
in biomolecular force fields](https://pmc.ncbi.nlm.nih.gov/articles/PMC5669220/)
examines this relationship in protein force fields. It provides a methodological
precedent, rather than a measurement of the effect in our GO models.

<a id="cosa-hanno-fatto-in-gaff2aim"></a>
## What was done in GAFF2–AIM

The paper and Supporting Information by
[Pinto et al.](https://doi.org/10.1021/acs.jpcb.2c06037) describe:

| Component | Reported source |
|---|---|
| Bonded terms | GAFF2 parameterization of GO fragments through ACPYPE, reused for the sheets |
| Charges | DFT electron density of the sheets, partitioned with DDEC6/Chargemol |
| LJ | Atomic volumes from the same density, using the Tkatchenko–Scheffler procedure described in the SI |
| Environment classes | Parameters averaged over classes defined by nearest neighbors |

Sections S1 and S4–S5 of the SI do not describe a new AIM bonded fit to QM scans.
The environment labels assigned to atom types are distinct from the source of
the bonded coefficients. The published GROMACS files contain explicit bonded
terms; the authors describe the amber99sb-ildn.ff-modified folder as containing
GO parameters obtained with ACPYPE.

The coordinated derivation of charges and LJ differs from fitting all
coefficients simultaneously to experiments. The paper tests the complete model
through water properties and amino-acid adsorption. It does not establish that
pairing DDEC6 with arbitrary LJ parameters preserves those results.

Preserved local sources: references/papers/Pinto2023_SI.txt, S1 and S4–S5;
references/gaff2_aim/README.md and GO_30/GO_sheet_3_0.top.
This inspection is not a numerical comparison of every published coefficient
against every version of gaff2.dat.

<a id="conseguenza-per-il-nostro-procedimento"></a>
## Implications for our procedure

Studying charges first is useful for defining the HF–RESP reference, checking
the effects of constraints and selecting a transferable assignment method.
Potential and charge comparisons test electrostatics; they do not select the
final force field for MD.

Subsequent comparisons of observables must identify the charge–LJ pair, water
model and 1–4 conventions. A RESP candidate with GAFF2 and a GAFF2–AIM preset
must each be evaluated in their own context; combinations from different
sources are explicitly identified as hybrid controls.

Proposed sequence:

1. Derive and verify charges with HF–RESP, retaining the controls.
2. Retain the stated reference bonded/LJ terms during diagnostic comparisons,
   without treating them as finally validated.
3. Check and, if necessary, calibrate the charge–LJ pair using multiple probes,
   interaction geometries and relevant observables.
4. Subsequently check bonded geometries and energy profiles, with particular
   attention to epoxides, edges and new GO–graft junctions.

Selecting the charge model and selecting the final parameters remain separate
research stages, with joint validation required before production.
