# Force-field presets, custom combinations and references

*Status verified as of 2026-11-19* (from the companion expanded-validation archive:
`FORCE_FIELD_SELECTION.md` and `ERRATA_AND_UPDATES_20261119.md`).

The intended interface offers named combinations and freely selectable
components. The [research catalogue](../research/benchmarks/2026-11-13/preset_catalogue.json)
records ten initial exercised combinations (nine transferable presets plus one
molecule-specific Q-Force Hessian-only in-sample control), subsequently expanded
across the verified **58/58** HF/6-31G\* RESP reference campaign (`strict40`
calibration + `14 C54` size holdout) and **49** common MM models (`14`
transferable paired methods evaluated with **unscaled `1.0` harmonic
frequencies**). The general preset bank is **not yet integrated into the
production CLI**: **`resp-current` (`resources/reference/current` + `GO_gaff2.frcmod`)
remains the only operational CLI assignment path and default**. Existing
commands and defaults are documented in [CHARGES.md](CHARGES.md),
[CHARGE_PRESETS.md](CHARGE_PRESETS.md), and the [README](../README.md).

## Candidates and verified campaign status (2026-11-19)

The [common-QM benchmark](../research/benchmarks/2026-11-13/README.md) supplies
the original pilot conventions, updated by the 2026-11-19 expanded audit and
errata. A preferred candidate for future validation is **not yet a validated
production default**.

| Use | Candidate | Verified evidence & readiness (2026-11-19) |
|---|---|---|
| Operational CLI default | `resp-current` (`resources/reference/current` + `GO_gaff2.frcmod`) | **Only operational CLI default** (`28.79%` C54-14 median ESP RRMS; unscaled freq MAE `46.52 cm⁻¹`); retained until `L0_v2` + `GO_gaff2_v3` pass Block 2 gates |
| Conventional Amber controls | `GAFF2 / AM1-BCC`, `GAFF2 / HF-RESP`, `GAFF2 / ABCG2` | Native parameter generation; `GAFF2_ABCG2` audited across **58/58** models (unscaled freq MAE `46.67 cm⁻¹` on `n=58`, `45.73 cm⁻¹` on `n=52`) |
| Learned bonded terms | `Grappa 1.4.0 / OpenFF 2.0 / AM1-BCC` | Best harmonic curvature (`17.77%` Hessian RRMS on 49 common models; force `11.46 kcal/(mol·Å)`, unscaled freq `28.56 cm⁻¹`) |
| Alternative learned bonded terms | `Espaloma 0.3.2 / OpenFF 2.0` (`AM1-BCC` / `EspalomaCharge`) | Best median force error (`9.28` / `9.56 kcal/(mol·Å)`) and unscaled freq MAE (`27.75` / `27.56 cm⁻¹`) on 49 common models; charge variants tested with bonded/LJ held fixed |
| Complete OpenFF combination | `Sage 2.3.0 / AshGC 1.0` | Native coupled combination (`17.19 kcal/(mol·Å)`, `18.65%` Hessian, `40.59 cm⁻¹` unscaled on 49 common models); `Sage 2.3 / EspalomaCharge` shows only a shallow π–π minimum (`−0.531 kcal/mol` at `3.8 Å`, MAE `6.406 kcal/mol`) |
| Blocked pre-`L0_v2` charge candidates | `L0_strict40_closed`, `Local_BCC_ESP_strict40`, `SQE_ESP_direct_strict40`, `SQE_L0_strict40` | **Blocked from promotion**: `L0_strict40` has severe `-COO-` charge pathology (`COO_C = +1.716076`, `COO_O = -1.030563`, `COO_anchor = -0.654951 e`); `Local_BCC_ESP_strict40` is rank-deficient (`35/36`, `cond = 4.26e+18`, `COO_anchor = +0.839..+0.923 e`); `SQE_ESP_direct_strict40` has unphysical `role_theta_map` weights (`augmented_rank = 24`) |
| Blocked pre-`v3` bonded candidate | `Local_GAFF2_v2_qm_ref` (`GO_gaff2_v2_qm_ref_20261108.frcmod`) | **Blocked from promotion**: `46`-model training overlap, overloaded `c2` (`CBX`/`CBD`), inter-epoxide `cx-cx-cx` ambiguity (`117.95°` vs `60.0°`), and force regressions on `16/52` (`16/37` neutral/mixed) cases despite improving `15/15` anions |
| Published OPLS baseline | `HierGO_OPLS_native` (`hiergo_opls`); `HierGO_OPLS_docx` (`Generating_OPLS_top.docx`) | Unpatched `6/58`, `docx`-analog `36/58` (unscaled freq `56.69` vs `64.58 cm⁻¹`); `COO⁻` fails neutrality (`-0.36 e` instead of `-1.00 e`) and lacks `EPH`/`EPT` |
| Molecule-specific QM Hessian control | `Q-Force 0.6.14` (Hessian-only fit, GAFF2 LJ/AM1-BCC) | In-sample Hessian fit (`3.928 kcal/(mol·Å)`, `13.75%`); no torsion scans and not transferable |
| Future transferable refits (Block 2) | `L0_v2` / `L1_reg` charges + `GO_gaff2_v3` bonded (`CBX`/`CBD` split) | Planned clean-room refit targets before any production default change |

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
