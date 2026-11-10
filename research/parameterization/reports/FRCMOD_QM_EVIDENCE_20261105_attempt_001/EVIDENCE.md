<a id="evidenza-e-protocollo-qm-per-i-frcmod-go"></a>
# Evidence and QM protocol for GO frcmod files

Date: 5 November 2026. Status: scientific proposal, no fit and no new QM job. Only existing documents/sources were read. This dossier does not replace the independent numerical file audit or a comparison using resolved topologies.

<a id="decisione"></a>
## Decision

Targeted physical checks are needed before promoting a frcmod to production; it has not yet been demonstrated that all terms need refitting. Exact term attribution, chemical context, and reuse of acceptable QM references already available come first. The final charge choice remains a prerequisite for torsional fitting and joint q–LJ validation. QM Hessians and scans are references independent of the MM model: they can be examined now, but they do not justify a new frcmod while q remains unresolved.

<a id="fonti-primarie-e-limiti-di-accesso"></a>
## Primary sources and access limitations

| ID | Exact source | Verification and locator | Support; limitation |
|---|---|---|---|
| E01 | Pinto AV, Ferreira P, Fernandes PA, Magalhães AL, Ramos MJ. *Development of Nanoscale Graphene Oxide Models for the Adsorption of Biological Molecules*. J Phys Chem B 127(2), 557–566 (2023; online 25 October 2022). DOI https://doi.org/10.1021/acs.jpcb.2c06037 ; https://pubmed.ncbi.nlm.nih.gov/36282235/ | Original abstract and bibliographic metadata read online on 5 November 2026; direct ACS access unavailable through the tool. | The work changes the nonbonded parameters, inherits the GAFF2 bonded parameters, and compares water/amino-acid adsorption; it does not certify every local GAFF2 variant or new edges. |
| E02 | Supporting Information for E01; https://api.figshare.com/v2/articles/21395006 | Local copy `references/papers/Pinto2023_SI.txt`, lines 12–18, 49–89, 111–132, 142–155; local PDF preserved with JSON metadata. The figshare endpoint was not reopened in the browser in this session. | ACPYPE fragments → GAFF2 library; PBE/DZVP-MOLOPT-SR-GTH; DDEC6/Chargemol and atomic volumes → LJ; nearest-neighbor classes. No new bonded fit against Hessians/scans is described here. |
| E03 | GAFF2-AIM authors' repository https://github.com/arvpinto/GO_sheets ; snapshot https://github.com/arvpinto/GO_sheets/tree/1e59f43bd6c4d8c19d0a2a208edf53ce11e00bc5 | Online README read; local snapshot `references/gaff2_aim/` and `references/gaff2_aim_source_parameters.json`. | Complete topologies and periodic coordinates; the modified ff folder incorporates GO parameters from ACPYPE. Do not treat the folder name amber99sb-ildn as an independent origin of a bonded fit. |
| E04 | Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. *Development and testing of a general amber force field*. J Comput Chem 25, 1157–1174 (2004). DOI https://doi.org/10.1002/jcc.20035 | Publisher abstract opened online; metadata and general scope. | GAFF covers organic compounds and also uses empirical/heuristic models; this is neither a dense GO benchmark nor a versioned source for every GAFF2 coefficient. |
| E05 | Cerutti DS, Debiec KT, Case DA, Chong LT. *Links between the charge model and bonded parameter force constants in biomolecular force fields*. J Chem Phys 147, 161730 (2017). DOI https://doi.org/10.1063/1.4985866 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC5669220/ | Full article read online, Abstract and Introduction; particularly the discussion of gas-/aqueous-phase surfaces. | Methodological evidence of q–bonded and 1–4 coupling; protein system, with no effect size automatically transferable to GO. |
| E06 | Allen AEA, Payne MC, Cole DJ. *Harmonic Force Constants for Molecular Mechanics Force Fields via Hessian Matrix Projection*. J Chem Theory Comput 14(1), 274–281 (2018). DOI https://doi.org/10.1021/acs.jctc.7b00785 ; https://www.repository.cam.ac.uk/items/aa3d9bf5-0a37-4322-9906-929dadd11ab8 | Institutional record and authors' abstract consulted through search; ACS 403; indexed accepted-manuscript PDF https://api.repository.cam.ac.uk/server/api/core/bitstreams/06646b1c-4830-47ae-8eba-85e7b8c3549a/content | Modified Seminario derives bond/angle parameters from the QM Hessian; it does not derive LJ or a reaction surface, nor does it automatically make a fragment transferable to a sheet. |
| E07 | Garcia NA, Awuah JB, Zhao C, Vuković F, Walsh TR. *Simulation-ready graphene oxide structures with hierarchical complexity: a modular tiling strategy*. 2D Mater 10, 025007 (2023). DOI https://doi.org/10.1088/2053-1583/acb0e1 | Publisher Crossmark metadata verified online: https://crossmark.crossref.org/dialog-content?doi=10.1088%2F2053-1583%2Facb0e1&domain=pdf ; IOP full text blocked by robots. | Bibliographic attribution verified; no claim about lactones inferred from the title/abstract or secondary citations. |
| E08 | HierGO code/authors https://github.com/IFM-molecular-simulation-group/HierGO ; commit 60382a1a1807ef0003ffcc3903392e913b1e5e9a | Online README and local snapshot `references/local_snapshot/opls/PROVENANCE.json`, `Oplsatom2`, `atom-typing.py`, `OPLS-top.f90`. | Actual inventory of the preserved subset: OH, epoxide, non-epoxide ether, carbonyl, COOH, terminal C–H, and sp² C; no separate lactone or COO⁻ table. This does not prove absence from every possible upstream resource. |
| E09 | Primary Amber-MD parmchk2 implementation https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/parmchk2.c | Source and CORR/penalty definitions consulted online. | A similarity score and term coverage are not equivalent to a measure of accuracy on GO. |

<a id="fonti-locali-lette-e-affermazioni-legate"></a>
## Local sources read and associated claims

- `docs/CHARGES_BONDED_LJ_COUPLING.md`: distinction of the q–LJ pair, inherited GAFF2-AIM bonded parameters, final validation not performed.
- `docs/BONDED_QM_PROTOCOL.md`: finite ωB97X-D/6-311++G(d,p) reference, separate Opt/Freq/Stability checks, QA to be completed. Do not suggest a new undifferentiated campaign.
- `docs/NEW_FRCMOD_CANDIDATES.md`: unfitted candidates and the `cx-cx-cx` collision outside a C–C–C triangle; context takes precedence over a simple type match.
- `docs/FRCMOD_AUDIT.md:9,15–21,25–39,55–57`: parmchk provenance + manual analogies reported by the user; override ≠ complete FF; context-dependent zeros, continuations, and values not physically validated.
- `docs/FORCE_FIELD_CANDIDATES.md`: complete conventions and preset scope; GAFF2-AIM for basal OH/EPX, not an already parameterized model of COOH/COO⁻ edges.
- `parameters/comparison_sources/HOSHIN_OLD_CURRENT_20261105/README.md` and `MANIFEST.json`: byte-identical sources and Hoshin/OLD/current origin, without evidence of a new QM fit.
- `references/go_prior_art_20261104/AmberGO_Tutorial.txt:743–772`: reduced GAFF selection for practical management/coverage; the tutorial's adjective “validated” does not replace an independent physical comparison.
- `reports/dft_bonded_reference_QA_20261030.json`: report dated 30 October, 7/8 preliminary passes and `final_parameter_validation: false`. In `pyrene_epoxide`, expected bond distance 15–16 = 2.110471977 Å and topology screen failed, despite positive frequencies and electronic stability. The raw Gaussian/fchk files were not recovered/checked here; no conclusion about the physical reaction or current job status.

<a id="inventario-hiergo-direttamente-verificato"></a>
## Directly verified HierGO inventory

`Oplsatom2` contains 17 rows, not 17 groups: HYD 3, EPX 2, ETH 2, CAB 2, CARB 5, terminal GRAT 2, and GRA 1. These are atom roles; the names do not automatically define all graphs, protonation, or the final topology.

`atom-typing.py:204–238` uses geometric neighbors: O with two C neighbors, with those two C atoms mutually neighboring → OE/epoxide; otherwise OT/ether. It does not provide a separate class for basal versus edge ether or for lactones. A lactone contains C(=O)–O–C closed into a ring and requires a specific ester environment: the mere presence of OT + OC does not prove correct assignment. CARB includes the acidic H and is not equivalent to lactone or COO⁻.

Discrepancy to keep visible: `Oplsatom2:6` reports ETH O −0.1700 e, whereas the comments in `atom-typing.py:174,344` report −0.2300 e. `OPLS-top.f90:40–68` reads Oplsatom2. The preserved operational table is therefore the source to follow for the comparison, without silently correcting the comments. Do not import these charges into GAFF2 as a complete HierGO preset.

<a id="riferimenti-aggiuntivi-scelta-metodologica-non-conclusioni-bibliografiche"></a>
## Additional references: methodological choice, not bibliographic conclusions

All scan counts, distances, and thresholds reported in REPORT.html are experimental-design proposals, not published results or new calculations. The final numerical criteria must be set with respect to the desired observable and QM convergence; no errors or CPU times have been invented.

This proposal reuses the already documented QM method for consistency of comparison, without declaring it the definitive level for reactions, periodic sheets, or anions in solution. For interaction energies, dispersion, BSSE/convergence checks, and declared geometries are proposed; the actual decision on the level and solvent must be made in the separate protocol before calculations.

Memory used only for routing: `MEMORY.md:47–60`, session 01a07be4-18ac-7111-bcc4-ab20caec3b9d; the operational claims above were reread from current sources. No memory changes.
