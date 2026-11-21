# Force-field presets, custom combinations and references

*Status verified as of 2026-11-21* (Block 2 AMBER Core validated across `step2a`–`step2d`;
`resources/reference/gaff2_v3/PROVENANCE.json`).

The CLI now ships versioned AMBER preset selection (`03_prepare_amber.py --list-presets`,
`--preset gaff2_v3 | gaff2_current`, `--charge-model {auto,l0_v2,current_lib,custom}`,
`--charges-json`) plus multi-engine exporters (`07_export_gromacs.py`,
`08_export_openmm.py`). **`gaff2_v3` (`resources/reference/gaff2_v3`:
`GO_gaff2_v3.frcmod`, 18 `*_gaff2_v3.lib` with 4-color DSATUR epoxides,
`L0_v2_charges.json`) is the validated and only recommended AMBER preset** (pairing cleanly with default `L0_v2` charges or external custom charges via `--charge-model custom --charges-json`), while **`gaff2_current` (`resources/reference/current`: `GO_gaff2.frcmod` + legacy libraries)** is **deprecated** (known bonded/charge defects) and kept untouched only for archival backwards compatibility (selected with `--preset gaff2_current` or implicit `--reference-dir` without `--preset`).

## Candidates and verified campaign status (2026-11-21)

The [research catalogue](../research/benchmarks/2026-11-13/preset_catalogue.json)
and [common-QM benchmark](../research/benchmarks/2026-11-13/README.md) were
expanded to **58** models (`52/58` in-domain Priority 1+2 GO; `6/58` Priority-3
rejected with `ValueError`) and evaluated across the Block 2 6-suite cross-matrix
with **unscaled `1.0` harmonic frequencies**.

| Use | Candidate | Verified evidence & readiness (2026-11-21) |
|---|---|---|
| **Validated CLI default (Block 2)** | **`gaff2_v3 + L0_v2` (`resources/reference/gaff2_v3`)** | **11/11 Go/No-Go gates PASS**: `ALL_52 F_med = 14.78` (`F_max = 25.39`), `H_med = 17.64%`, `G_med = 0.100 Å`, unscaled `Freq_med = 35.0 cm⁻¹`, `0` imag modes (`52/52`); size-transfer `C54` (excluded from parameter fits; unseen sheets in Block 3): `C54_neutral_9 F_med = 17.32`, `C54_anionic_5 F_med = 12.61`; torsions `pyrene_COOH = 0.045`, `pyrene_COO = 0.012 kcal/mol`; π–π stacking `MAE = 1.231 kcal/mol` (`z_min = 3.4 Å`, `-5.806 kcal/mol`); `0/129` Amber namespace collisions; `WORKBENCH` topology parity `< 1e-4 kJ/mol` (`OpenMM 8.6.1` native vs pure-Py XML `NoCutoff`/`CutoffNonPeriodic`/`PME`, and GROMACS `.top/.itp` via `OpenMM GromacsTopFile`) / `< 0.03 kJ/mol` (`sander` Version 24.0, `amber24u3+at25`) |
| Legacy CLI preset (deprecated) | `gaff2_current` (`resources/reference/current` + `GO_gaff2.frcmod`) | **Deprecated archival baseline** (kept only for regression compatibility; known defects): `ALL_52 F_med = 25.46` (`F_max = 102.27`), `H_med = 21.36%`, unscaled `Freq_med = 46.5 cm⁻¹`, `1` imag mode (`pyrene_COOH`), `C54-14` ESP RRMS `28.79%`; `GO_gaff2_v3.frcmod` fixes `COOH`/`COO-` typing (`0` imag) |
| Conventional Amber controls | `GAFF2 / AM1-BCC`, `GAFF2 / HF-RESP`, `GAFF2 / ABCG2` | Native parameter generation; `GAFF2_ABCG2` audited across **58/58** models (`ALL_52 F_med = 18.18`, `F_max = 98.39`, `H_med = 20.10%`, unscaled `Freq = 45.7 cm⁻¹`; `C54_neutral_9` ESP `55.32%`, `C54_anionic_5` `8.69%`) |
| Learned bonded terms | `Grappa 1.4.0 / OpenFF 2.0 / AM1-BCC` | Best harmonic curvature (`17.77%` Hessian RRMS on 49 common models; force `11.46 kcal/(mol·Å)`, unscaled freq `28.56 cm⁻¹`) |
| Alternative learned bonded terms | `Espaloma 0.3.2 / OpenFF 2.0` (`AM1-BCC` / `EspalomaCharge`) | Best median force error (`9.28` / `9.56 kcal/(mol·Å)`) and unscaled freq MAE (`27.75` / `27.56 cm⁻¹`) on 49 common models; charge variants tested with bonded/LJ held fixed |
| Complete OpenFF combination | `Sage 2.3.0 / AshGC 1.0` | Native coupled combination (`17.19 kcal/(mol·Å)`, `18.65%` Hessian, `40.59 cm⁻¹` unscaled on 49 common models); `Sage 2.3 / EspalomaCharge` shows only a shallow π–π minimum (`−0.531 kcal/mol` at `3.8 Å`, MAE `6.406 kcal/mol`) |
| Blocked pre-`L0_v2` charge candidates | `L0_strict40_closed`, `Local_BCC_ESP_strict40`, `SQE_ESP_direct_strict40`, `SQE_L0_strict40` | **Blocked archival drafts**: superseded by `L0_v2` |
| Blocked pre-`v3` bonded candidate | `Local_GAFF2_v2_qm_ref` (`GO_gaff2_v2_qm_ref_20261108.frcmod`) | **Blocked archival draft** (`46`-model overlap, overloaded `c2`, `cx-cx-cx` ambiguity): superseded by clean-room `GO_gaff2_v3.frcmod` |
| Published OPLS baseline | `HierGO_OPLS_native` (`hiergo_opls`); `HierGO_OPLS_docx` (`Generating_OPLS_top.docx`) | Unpatched `6/58`, `docx`-analog `36/58` (unscaled freq `56.69` vs `64.58 cm⁻¹`); `COO⁻` fails neutrality (`-0.36 e` instead of `-1.00 e`) and lacks `EPH`/`EPT` |
| Molecule-specific QM Hessian control | `Q-Force 0.6.14` (Hessian-only fit, GAFF2 LJ/AM1-BCC) | In-sample Hessian fit (`3.928 kcal/(mol·Å)`, `13.75%`); no torsion scans and not transferable |

GAFF 1.81 remains a reproducible baseline (with `GO_gaff2_old.frcmod` separately failing `tleap` on 3 neutral `C24` models due to missing `c2-ca-c6-c6`, `c6-c6-cx-c6`/`c6-cx-c6-c6`, and `c6-c6-c6-c6` torsions). Published HierGO/OPLS and
GAFF2–AIM models can be retained as complete reference combinations, with
their own chemistry and conventions. On the finite-fragment C54 HF/6-31G* ESP
comparison, HierGO/OPLS outperforms `current_gaff2` on neutral fragments
(22.13% vs 32.48%), while raw unclosed GOPY motif charges fail charge
conservation (`−0.45` to `−1.48 e`) and Stauffer ESP/SESP cover only basal
OH/epoxide without scaffold or edge caps. GAFF2–AIM couples periodic DDEC6
charges and AIM LJ with GAFF2 bonded terms; it is distinct from ordinary GAFF2
and does not transfer automatically onto finite fragments. No new RESP2, MBIS
or DDEC6 charge derivations are part of this work.

## Complete preset records

A preset must record exact sources/checkpoints and hashes; charge protocol,
charge/protonation domain; bonded terms; LJ representation and mixing;
exclusions and 1–4 coefficients; water/ion model; tested chemistry, boundary
conditions and observables; and custom graft-junction terms. Gas-phase
benchmark presets have no water model assigned.

The actual 54 pilot systems were checked against their graphs: all 1–2 and
1–3 nonbonded pairs are excluded, 1–4 electrostatics are multiplied by 5/6,
1–4 LJ by 1/2, and LJ uses arithmetic σ/geometric ε mixing
(Lorentz–Berthelot). Amber denominators are SCEE = 1.2 and SCNB = 2 for these
scaled pairs. Individual torsion records may contain other defaults;
the decisive evidence is the resolved pair interaction. See the
[convention audit](../research/benchmarks/2026-11-13/data/verification.json).

These observed conventions must not become forced rules for future presets.
An OPLS-style combination can use different mixing/scaling: preserve its
declared rules or name an intentional hybrid.

## Selection within a force-field family

Evaluate each complete preset for a declared host force field and application.
First check atom typing, parameter coverage, charge and protonation, the
functional form, the effective exclusions and 1–4 terms, and the treatment of
cross interactions with water, ions, and any biomolecular partner. For grafts,
check the bonded and electrostatic terms across the junction as well as the
unchanged parts of the partner topology. Engine parity verifies the resolved
parameters. Physical accuracy requires separate reference comparisons.

Compare candidates on common cases. Report errors by observable, chemical
environment, and charge state. Inspect the largest regressions alongside
median performance. Give greater weight to the reference conventions and
observables that matter for the intended use:

- For an Amber RESP workflow, use HF/6-31G* RESP as the electrostatic reference.
  Check agreement with the atomic charges and the exterior ESP, together with
  torsions, hydration, and cross interactions. AM1-BCC is a separate Amber
  baseline. Compare published coupled Amber-family presets with their own
  charge and Lennard-Jones terms.
- For CHARMM/CGenFF, assess native atom typing, charges and water interactions
  using the CHARMM conventions. Preserve the CHARMM functional form and pair
  rules. See [the CGenFF parameterization study](https://doi.org/10.1002/jcc.21367).
- For OPLS-AA, GROMOS, OpenFF, and other families, retain the documented charge
  model, site representation, bonded terms, mixing and 1–4 conventions of each
  tested preset. Treat learned bonded terms and charge substitutions as named
  combinations with their own evidence.

Freeze the fragment-derived candidates and selection criteria before evaluating
realistic whole sheets. Minimize the generated sheets, verify their chemistry,
then test whether the frozen candidates reproduce the relevant structure and
interactions at that scale. Whole-sheet QM can provide an independent check
when the electronic state and calculation pass their acceptance criteria.

## Free custom combinations

Users can choose charges, bonded parameters, LJ/nonbonded terms, mixing,
exclusions and 1–4 settings separately. Scientific validation status is
metadata; **lack of validation does not block expert use**. Exported provenance
must identify the exact choice. Changing Sage's native charge model creates
a custom Sage-based combination, rather than an unchanged Sage 2.3.0 preset.

For example, Grappa bonded terms, GAFF2 LJ and the local GO charge correction
would be a custom experiment. It has no score in this panel and must not
inherit the score of Grappa/OpenFF/AM1-BCC. Component selection does not
silently refit LJ or compensate charges to reproduce a named preset.

Technical errors still require a fix: missing atom mappings or parameters,
non-finite coefficients, inconsistent units or an engine unable to represent
the chosen functional form/pair rules exactly. Report such limitations
without silently changing the selection. A charge swap updates both particle
charges and affected 1–4 products. Unknown chemical environments need explicit
assignments, without hidden fallbacks.

Amber can execute other charge sets; its file format does not require AM1-BCC
or RESP. Physical consistency is separate: valence, LJ and water models may
have been fitted with a particular charge model.
[Sage 2.3.0](https://doi.org/10.1021/acs.jctc.6c00169), for example, was developed
with AshGC. A better isolated ESP fit does not establish a better liquid or
GO–protein interaction model.

Protein/DNA charges and internal parameters stay with their chosen molecular
force field. Check atom-type namespaces, wildcard terms, load precedence and
resolved energies against a baseline topology when adding GO libraries.
Noncovalent contacts need cross-nonbonded validation. Covalent grafts also
need junction charges, angles, torsions and QM fragments; geometry construction
or topology export alone is insufficient.

## QM and experimental references

| Reference | Useful comparisons | Conditions to match |
|---|---|---|
| Finite-fragment QM | ESP, forces/Hessians, geometries, torsions, graft junctions | Charge, spin, protonation, capping, size sensitivity, geometry and QM protocol |
| Published periodic QM | Basal-plane geometry, corrugation, mechanics, interlayers | Composition, pattern, strain, boundaries and environment |
| Elemental analysis | Atomic O/C and composition feasibility | Mass-to-atom conversion, dry/hydrated basis, impurities and sample identity |
| XPS / solid-state NMR | Group populations and local environments | Assignments, deconvolution and overlapping signals |
| Titration | Accessible ionizable groups / protonation | Same product, pH, ionic strength and chemical interpretation |
| AFM / TEM | Apparent height, lateral size and domains | Substrate/environment and resolution |
| XRD | Interlayer spacing in stacked materials | Hydration, stacking, salt and preparation |

Finite QM is suitable for local tests; a new periodic VASP campaign is not
required to make them meaningful. Published periodic results add references
for their own observables. Do not pool finite and periodic absolute energies
into one accuracy score, or compare single-flake thickness with XRD spacing.

The P935 demonstration uses the
[vendor TDS](https://filgen.jp/Product/Bioscience4/Abalonyx/TDS_GO_Pastes_and_powders.pdf)
dry elemental intervals only. Four examples satisfy atomic O/C 0.381–0.399.
Their basal OH:epoxide ratio and zero edge oxygen population are assumptions,
despite FTIR evidence for edge groups in the product. They are not full P935
reconstructions. Existing LayerOne metadata keeps vendor measurements,
author-derived integer recipes and thesis assumptions distinct.

Measured intervals generate families of feasible compositions and patterns.
Preserve uncertainty and report infeasible constraints instead of relaxing
them silently. Do not equate XPS peak fractions with exact group counts
without a measurement model, or zeta potential with atomic charges.
