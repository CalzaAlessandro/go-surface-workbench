> Subsequent update: the third C54 has also been completed. The DFT panel now has 19 contexts. See the [new report with Espaloma and SQE](../CHARGE_CONTROLS_UPDATE_20261103_attempt_001/REPORT.html). This folder preserves the previous snapshot with 18 contexts.

<a id="cariche-resp-per-amber-confronto-del-3-novembre-2026"></a>
# RESP charges for Amber: comparison of 3 November 2026

16 atom-specific HF/6-31G* RESP references on DFT geometries,
16 corresponding DFT–RESP controls, four L0 shared-block RESP fits,
and 12 AM1-BCC controls with complete BCC parameter coverage have been completed.
Two additional C₅₄ structures bring the DFT–RESP panel to 18 structures.
No new production library has been selected.

HF and DFT share the same first 16 geometries: these are not 34 independent
chemical contexts. COOH and COO⁻ each have only one physical context in the panel.

The version to read is the [complete HTML page](RESP_AMBER_REPORT.html),
with charges, methods, selectable plots, and considerations for the next steps.

<a id="resp-ed-esp-cosa-viene-calcolato"></a>
## RESP and ESP: what is calculated

**The new charges are RESP charges, obtained with the Amber RESP executable.**
ESP denotes the quantum electrostatic potential used as input.
The ESP residual measures how well the RESP charges reproduce that potential:
it does not mean that Gaussian's unrestrained ESP charges were adopted.
See also the [RESP documentation](https://q4md-forcefieldtools.org/RED/resp/).

The primary reference is the gas-phase MK potential at HF/6-31G*,
calculated on geometries optimized at ωB97X-D/6-311++G(d,p), followed by the
hyperbolic RESP fit. The DFT geometry and the HF level of the potential are recorded
separately. The DFT–RESP controls instead use the ωB97X-D/6-311++G(d,p) potential.
AM1-BCC is an alternative method, explicitly kept separate from RESP.

The atom-specific fit uses qwt = 0.0005, ihfree = 1, integer-charge constraints,
and the declared equivalence constraints. The three rigid orientations have recorded weights
that compensate for the different numbers of grid points. They represent one physical geometry,
not three independent conformers. Since there are no CH₂/CH₃ groups to refine,
RESP stage 2 is declared inactive: a second fit is not simulated.
Details are in [RESP_PROTOCOL.md](../../docs/RESP_PROTOCOL.md).

The external RHF→UHF instability observed in previous checks is distinct
from an SCF error or an automatic incompatibility with Amber. This panel
retains the restricted HF reference without claiming that its electronic stability
has been validated. The DFT geometries passed their respective checks.

<a id="l0-e-gruppi-analizzati"></a>
## L0 and the groups analyzed

L0 is a **RESP fit with one shared charge per chemical role**, independent
of the distance from other groups. Unfunctionalized interior carbon is
fixed at 0 e; each neutral block sums to 0 e and the carboxylate block to −1 e.
The COOH/COO⁻ anchor carbons belong to their respective blocks.
Edge carbons and hydrogens have roles distinct from the interior background.
This is a test of the transferability of the modular hypothesis, not an approved preset.

The panel covers:

| Chemistry | Library unit | Available contexts |
|---|---|---|
| Interior skeleton and C–H edge | GGG, GRT | C₁₆, C₂₄, and C₅₄ DFT controls |
| Basal hydroxyls | HYD | Two cis/trans OH groups, near/far, OH/EPX mixtures |
| Basal epoxides | EPX | Single, pair, and mixtures; C₅₄ DFT pair recovered |
| Edge epoxides | EPH, EPT | Two dedicated C₂₄ controls |
| Carbonyls | CAB | Fragment with two carbonyls |
| Protonated carboxyl | CBX | COOH with anchor carbon |
| Carboxylate | CBD | COO⁻ with anchor carbon |

C₁₆/C₂₄/C₅₄ indicate the size of the starting aromatic skeleton;
the carboxyl carbon is additional. The cases named OH_EPX contain
two OH groups and one epoxide. This is not yet an exhaustive panel of all
GO group distances, densities, faces, protonation states, and combinations.
Lactones and other proposed functionalities have not been parameterized here.

<a id="risultati-descrittivi"></a>
## Descriptive results

The three orientations of each structure are averaged first, followed by the
median across structures. RRMS is the relative root-mean-square error in the potential;
lower values indicate better reproduction of that QM reference.

| Charges evaluated against the same HF ESP | Median RRMS, 16 structures |
|---|---:|
| Atom-specific HF RESP, specific to each fragment | 0.0797 |
| Atom-specific DFT RESP on the same fragments | 0.1123 |
| HF L0 RESP, shared fit across the 16 cases | 0.1883 |
| Complete current library | 0.3505 |

The atom-specific fits are also evaluated against the data used to determine them. The advantage
of these fits does not demonstrate better dynamics or adsorption affinity.
The L0 fit across all 16 cases is also a training result.

A separate check excludes all C₂₄ structures from the fit: L0 is estimated
on the seven contexts with a C₁₆ skeleton and applied to the seven C₂₄ cases with covered roles.
EPH and EPT remain unavailable because they are absent from the training set.

| Comparison on the same 7 C₂₄ cases | Median RRMS, HF ESP |
|---|---:|
| HF L0 RESP trained only on C₁₆ | 0.1792 |
| HF L0 RESP trained on the 16 cases, thus including C₂₄ | 0.1777 |
| Complete current library | 0.3501 |

This is an initial check outside the training set, favoring the L0 candidate
on the available panel. It does not demonstrate transferability to macroscopic sheets,
new chemistries, or all oxidation conditions.

However, the charge table shows strong compensation for COO⁻: in HF L0
across 16 cases, carboxyl C = +1.731213 e, each O = −1.049356 e, and anchor
carbon = −0.632501 e. Before adoption, sensitivity to constraints,
correlation between variables, and additional carboxylate contexts must be checked.
A RESP decomposition does not uniquely measure charge transfer.

<a id="am1-bcc-e-modelli-pubblicati"></a>
## AM1-BCC and published models

Four cases have an uncovered native BCC term and are excluded from the
main comparison: pyrene_dione, C24_single_basal_EPX_control,
C24_EPH_edge_control, and C24_OH_OH_near_trans. The missing types are
16–12–2 or 12–16–2. They are replaced neither with zero nor with invented parameters.
The native charges are preserved; the candidate labeled AM1-BCC* in the plots
contains a small constrained correction for total-charge rounding.
The L0 projection of AM1-BCC is yet another candidate, distinct from the native method.

The complete CSV also includes controls from HierGO/OPLS, GOPY, and Stauffer,
and the actual coverage of each candidate. Comparisons with the OLD whole-sheet
libraries and GAFF2-AIM do not have a complete atom mapping onto the
fragments here: **they have not been assigned a fictitious score**. Hybrid controls that
combine different sources are labeled as such. Ranking only the
potential residuals does not select a charge/Lennard-Jones pair.

<a id="file-da-leggere"></a>
## Files to read

- [Selectable HF/DFT plot](comparison.html), [HF PNG](HF_ESP.png), [DFT PNG](DFT_ESP.png); the SVGs can be exported.
- [L0 RESP charges by role](L0_shared_role_charges.csv), with four separate fits for HF/DFT and full/C₁₆ training.
- [Summary numerical results](SUMMARY_RESULTS.json).
- [Coverage of all candidates](CANDIDATE_COVERAGE.csv) and [comparisons on identical subsets](MATCHED_SUBSET_SCORES.csv).
- [Residuals by structure](RRMS_by_physical_context.csv), including unavailable cases and their respective reasons.
- [HF–DFT RESP differences by role](matched_HF_DFT_role_differences.csv).
- Individual charges, group sums, pairs, and distances from groups: [HF, 16 contexts](context_HF_16/context_analysis.json) and [DFT, 18 contexts](context_DFT_18/context_analysis.json), with CSVs in the same folders.

<a id="recupero-e-provenienza"></a>
## Recovery and provenance

Job 1452975 completed the complementary synthesis. Actual fits that had already
finished were recovered, checking artifact hashes and MK logs;
Gaussian, SQM, and Antechamber were not rerun. The new shared fit
uses RESP from an identified AmberClassic source, with increased capacity for
atom copies across different molecules and orientations. The solver is unchanged;
comparison with the earlier pilot using the installed RESP passes within 10⁻⁶ e.

Job 1453077 recovered the DFT–RESP fits for C54_pristine_parent_reference and
C54_EPX_EPX_near_trans. The espgen header expands the point count beyond
99,999: both the Python parser and the Fortran reader required a read operation
that would not truncate the additional digit. The change to the separate executable
concerns header reading and reporting, not the fit equations.
The pilot regression agrees at output precision.
All points were retained: 120,100/120,005/120,107 in the first case and
122,272/122,477/122,602 in the second. The original inputs and outputs remain
preserved, and no system installation was modified.

Local result snapshots in `PARAMETERIZATION/calculations/`:

| Archive | SHA-256 |
|---|---|
| complementary_results_20261103_snapshot_001.tar.gz | 32ea75cec824cd06e740687215ec4d59267461f462ddca88e1509927c1d46552 |
| c54_recovered_results_snapshot_20261103_attempt_001.tar.gz | 471798f22cd9d6ae5306d9b301c878a466cc11ac03ba32a3b1d04a786cca4759 |

The remote directories are isolated under
`/export/scratch/alessandro.calza3/GO_PARAMETERIZATION_20261029_fqOFQv/calculations/`.
The synthesis and report code is preserved in the respective attempt
directories, together with the input manifests and executable checks.

At the last check on 3 November, C54_single_basal_EPX_control
(1406363_3) was still running, after approximately 3 days and 11 hours, against
the five-day limit. Job 1453111 waits for its completion and will attempt
postprocessing only after checking the QM artifacts. The dependency is
afterany: even if the parent fails, the scientific checks remain mandatory.

Charge selection, sensitivity/transferability checks,
a consistent charge–Lennard-Jones comparison, and solvent and adsorption
checks remain outstanding. The new bonded/frcmod fit and the final reconstruction
of the four materials remain deferred.

<a id="decisioni-successive-confermate-il-3-novembre"></a>
## Subsequent decisions confirmed on 3 November

HF–RESP is the primary method; DFT–RESP and AM1-BCC are controls. The total charge follows
from the product chemistry and protonation and is an integer, which may be nonzero.
The choice between a fixed library, local classes, and static redistribution remains
open. GAFF/frcmod work is deferred. The [protocol for the next comparison](../../docs/CHARGE_STRATEGY_BENCHMARK_20261103.md)
includes transferability, scale, new junctions, and experimental references.

The [COO diagnostic](COO_charge_diagnostic.csv) compares the current library,
molecular HF RESP, and L0: the high C value appears in L0, whereas the
molecular fit gives C = +0.807216 e and O = −0.786679 e. Controlled comparisons
are needed to attribute the difference to the various constraints. The [interior carbons](basal_C_summary_by_context.csv)
are summarized by context and distance; these are not measurements of physical
charge transfer. No new fits were performed in this update.
