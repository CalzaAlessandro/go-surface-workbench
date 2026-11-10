<a id="scelta-delle-cariche-go-riproducibilità-trasferibilità-e-scala"></a>
# Selecting GO charges: reproducibility, transferability and scale

Decisions confirmed by the user on 3 November 2026. This is the protocol for
the next comparison, rather than a newly implemented or validated method.

Update on 4 November: the first SQE fit is complete as a pilot study;
QM calculations on the two whole sheets and fragment expansion are running.
The [SQE protocol](SQE_HF_RESP_PILOT_20261104.md) and
[current report](../reports/QM_SQE_UPDATE_20261104_attempt_001/REPORT.html)
describe the results and limitations. This page remains the plan for choosing
among A/B/C; it does not approve a final method.

<a id="obiettivo-e-perimetro"></a>
## Objective and scope

Derive a reproducible charge model for GO of a user-selected size, suitable
for Amber with small molecules, proteins and subsequent functionalization.
HF/6-31G* MK followed by RESP remains the main reference. DFT–RESP and AM1-BCC
are separate controls. GAFF/GAFF2, Lennard-Jones and frcmod fitting remain
deferred: the current decision concerns how to assign charges.

The classical groups already covered remain in scope: GGG/GRT, HYD, EPX,
EPH/EPT, CAB, CBX and CBD. Ethers, lactones/lactols and other chemistries
present in models such as HierGO, together with PEG, amino acids and peptides,
are extensions requiring parameterization. They must not silently receive
charges from superficially similar groups.

<a id="vincolo-obbligatorio-sulla-carica-totale"></a>
## Mandatory total-charge constraint

The sum of partial charges must equal the integer formal charge of the final
product. GO is not required to be neutral. Protonation, charged groups and
reaction products determine the total.

Charge-accounting example: 16 carboxylates and 4 ammonium centers at +1 e,
with no other ionic groups, give Q = −12 e. Neutral PEG does not change the
total if the newly formed junction and proton balance are also neutral. A new
graft requires the product graph with explicit H atoms and formal charges:
simply summing precursor charges is insufficient when atoms or protons may
be lost.

Covalently functionalized GO, adsorbed molecules and counterions have separate
charge accounts. Each independent chemical component retains its own charge
state; any neutrality requirement for the solvated box is separate. The total
does not depend on the choice of L0, environment classes or redistribution.
It is imposed during fitting or by construction, rather than by arbitrarily
rescaling all atoms at the end.

A formal COO⁻ charge of −1 e does not require the partial charges on that
particular atom subset to sum to −1 e in a molecular fit. Neutral individual
blocks and an exactly zero basal background are additional constraints of
the modular model; they must be tested separately from mandatory conservation
of the total charge.

<a id="tre-famiglie-da-confrontare"></a>
## Three families to compare

| Family | Assignment during preparation | What to test |
|---|---|---|
| A · Fixed L0 library | One charge per role, shared by every occurrence; baseline already calculated | Transferability without environment classes, effects of a zero background and block boundaries |
| B · Environment-based library | A small number of chemical classes with shared RESP vectors; neighbor type, bond distance, edge and junction | Benefits on entire fragments excluded from fitting; coverage and parameter count |
| C · Static redistribution | Baseline values plus charge-conserving corrections assigned to the complete GO structure | RESP calibration, identifiability, size/geometry dependence and assignment cost |

All three produce fixed charges during MD and MM/GBSA. Family B does not
describe variable-charge dynamics. Dependence on instantaneous Euclidean
distance during MD would require a different energy model.

For B, graph distance is the first descriptor to test: it does not change
when the sheet bends. Group type, basal/edge position and local face are
distinct descriptors. Any 3D distance enters only as an additional hypothesis
tested on representative geometries; it does not automatically become a
universal 1/r correction.

For C, +δ/−δ corrections on bonds conserve the total by construction. The
baseline must already have the correct total; disconnected components do not
exchange charge through nonexistent bonds. A local increment rule and an SQE
solver are different candidates, not synonyms. SQE+qp is a published precedent,
but its published parameters are not a GO-specific RESP calibration.

<a id="primo-confronto-isolare-leffetto-dei-vincoli"></a>
## First comparison: isolate the effects of constraints

Before multiplying classes by distance, compare:

1. Current L0: zero internal background and blocks with imposed charges.
2. A variant with an extended basal environment, keeping charge boundaries
   controlled. Test environments one and two bonds away as benchmark
   hypotheses, rather than as an established physical radius.
3. A variant with charge constrained over a larger functional region rather
   than only the original block. If two regions overlap, merge the domains
   or use charge-conserving increments: one atom cannot belong to two
   independent sums constructed as though they were disjoint.

Change one assumption at a time. Atomwise RESP charges for each fragment
remain a more flexible reference, distinct from a library applicable to any
new sheet. Internal-carbon values from finite fragments are not automatically
transferable to a macroscopic surface.

<a id="dato-diagnostico-sul-carbossilato"></a>
### Diagnostic carboxylate data

For the already calculated pyrene_COO model, the molecular HF–RESP fit gives
COO_C = +0.807216 e and each COO_O = −0.786679 e. The current library uses
+0.760000 e and −0.790000 e; the shared HF L0 gives +1.731213 e and −1.049356 e.
The anchoring carbon differs among the three models.

The molecular and L0 fits have different constraints: this comparison does
not establish which constraint alone causes the difference. It requires the
controlled comparison described above. In particular, it does not show that
the HF–RESP protocol requires +1.73 e on the carboxyl carbon, or that the
original library is wrong. COOH and COO⁻ each have only one physical context:
three rotations do not constitute three independent examples.

<a id="dati-generalizzazione-e-criterio-di-scelta"></a>
## Data, generalization and selection criterion

Use the same HF references, geometries, grids, weights and RESP conventions
for all candidates. Reuse available real data before launching new QM
calculations. Evaluate the residual of the potential generated by RESP
charges, parameter sensitivity, domain sums and the response to changes in
size. Analyze internal and edge carbons separately.

Exclude entire fragments, separations and sizes from fitting, rather than
one of their rotations. A chemical class absent from training remains
unsupported. Subsequent targeted expansion must include COOH/COO⁻ at multiple
sizes and near/far from OH/EPX, followed by representative new covalent
junctions. For conformers, distinguish physical variations from numerical
rotations.

Prefer the simplest model that gives stable results outside training and
requires suitable time and memory for the intended sheets. A local-environment
rule does not require a new whole-sheet HF calculation for every generated
structure. For a global solver, cost and size dependence must be measured;
linear scaling is not assumed.

<a id="riproducibilità-e-scalabilità-della-costruzione"></a>
## Reproducibility and scalability of construction

Separate the chemical recipe, geometric realization and charge preset.
Record the seed, dimensions, actual group counts, protonation, connectivity
and IDs, library version/hash and the charge of each atom. The same recipe
with the same seed must reproduce the same realization; different seeds
produce an ensemble of surfaces for robustness checks.

Basal groups scale with area; edge sites scale with the perimeter and defect
boundaries. With carboxylates confined to the outer edge, preserving the same
COO/C fraction or charge density across sizes is not guaranteed. The recipe
must state which densities it preserves and reject unrealizable counts,
without inventing sites or changing the charge to compensate. This is a
specification to verify, rather than a newly added feature.

For PEG or amino acids, characterize the repeating chain, termini and
GO–graft junction separately. The new bond changes the chemical environment:
a free PEG/amino-acid model does not validate charges in the covalent product.
GO, complex and MM/GBSA subset topologies must inherit the same final
assignment, without recalculating it separately.

<a id="appoggio-sperimentale"></a>
## Experimental support

Atomic partial charges have no unique, directly measurable experimental value.
Experiments constrain composition, protonation and collective properties.
There is therefore no universal experimental table that can independently
select among the three assignment schemes.

- XPS and complementary analyses: the recipe and group proportions, with
  assignment uncertainty; not a unique partial charge for each atom.
- Potentiometric titration: the number of dissociable sites and the sample's
  pH response, useful for defining model protonation and total charge.
  [Tao et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10339393/) combine
  titrations with DFT structural hypotheses; the latter remain interpretations.
- Functionalization: [Orth et al.](https://doi.org/10.1016/j.jcis.2016.01.013)
  use titration and spectroscopy to characterize introduced groups.
  These do not directly measure junction RESP coefficients.
- Adsorption and calorimetry: [Pandit and De](https://doi.org/10.1021/acs.jpcc.6b11571)
  measure GO–amino-acid interactions. These provide a subsequent test of the
  complete model, with consistent water, protonation and nonelectrostatic
  parameters. Charges alone are not currently being fitted to a single
  adsorption energy.

LayerOne materials require data from the corresponding sample: results for
GO produced by other synthesis routes provide precedents and transferability
controls, rather than measurements of our sample's properties. Proteins and
peptides provide a later, more complex test; small probes first help separate
the effects of the model more clearly.

<a id="precedenti-metodologici"></a>
## Methodological precedents

[R.E.D., Dupradeau et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2918240/)
describe RESP fits over multiple molecules/conformations with fragment
constraints. This provides a precedent for building consistent libraries;
it does not by itself validate GO. Section S5 of the
[GAFF2–AIM SI, Pinto et al.](https://doi.org/10.1021/acs.jpcb.2c06037)
uses nearest-neighbor classes: a chemical-environment precedent distinct
from our RESP target. [SQE+qp](https://link.springer.com/article/10.1186/s13321-021-00528-w)
provides a family of charge-conserving redistribution models with a charge
reference; GO/RESP calibration still needs to be performed and validated.
