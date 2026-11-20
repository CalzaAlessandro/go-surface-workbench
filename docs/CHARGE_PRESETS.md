# Charge presets and the RESP-first study

*Status verified as of 2026-11-19* (from the companion expanded-validation archive:
`FORCE_FIELD_SELECTION.md` and `ERRATA_AND_UPDATES_20261119.md`).

The source catalog is `resources/charge_presets.json`. It records charge family,
coupled parameter expectations, provenance, and readiness. It remains a research
catalog: **`resp-current` (`resources/reference/current`) is the only operational
CLI assignment path and default**, and no generic production preset selector is
yet integrated.

| ID | Charge model | Status (2026-11-19) |
|---|---|---|
| `resp-current` | Existing user RESP-derived fragment libraries (`resources/reference/current`) | **Only operational CLI default**; legacy reference (`28.79%` median C54-14 ESP RRMS) |
| `resp-go-new` | New RESP-derived GO blocks/context response (`L0_strict40`, `Local_BCC_ESP`, `SQE_*`) | **Blocked from production promotion** pending `L0_v2` / `L1_reg` refits; known pathologies documented below |
| `gaff2-aim` | Published environment-averaged DDEC6 plus AIM LJ | Audited on native periodic `GO_18`/`GO_30` sheets; generic finite-GO mapping unsupported |
| `opls-hiergo` | Published HierGO/OPLS model | Audited (`6/58` unpatched, `36/58` with `Generating_OPLS_top.docx` analogs); lacks `COO-` (`-0.36 e`), `EPH`/`EPT` |

GAFF2–AIM source topologies/coordinates at nominal C/O 1.8 and 3.0 are retained
with their associated force-field files and license under
`resources/reference/gaff2_aim/`. The source manifest verifies the original
download hashes. These fixed source sheets are parameter references, not
material-size presets: users always choose generated sheet dimensions.
The published basal OH/epoxide data do not establish finite-edge, carbonyl,
COOH/COO⁻, or graft-junction parameters. No missing group charge is invented and
no source coefficient is silently renormalized.

The new derivation remains RESP-based. First measure group and surrounding-carbon
responses on isolated/balanced and mixed fragments at controlled graph distances,
faces, protonation, and sizes. Then compare a regularized fixed-block library
(`L0_v2`), a full-rank context/proximity rule (`L1_reg`), and a physically
constrained redistribution method calibrated to RESP.
Existing NPA/AM1-BCC/ABCG2 tools remain comparison controls rather than default new
charge derivations. Source presets are not refitted.

Published presets are complete, indivisible force-field bundles. Their charges,
bonded terms, LJ parameters, mixing rules, and 1–4 conventions must be preserved
together. Selecting `gaff2-aim` or `opls-hiergo` does not permit a charges-only
override or an SQE/Espaloma correction under the same published name. A new
derivative would require a separate identity, provenance, and validation.
For HierGO/OPLS, `COO−` is not covered by the archived 17-row `Oplsatom2` charge table
(`-0.36 e` instead of `-1.00 e`). No missing chemistry is silently borrowed from
another force field. Amber export also requires a validated conversion and an
explicit compatible GB-radius policy.

## Verified 58/58 Reference Panel and Blocked Pre-`L0_v2` Candidates

For the Amber-oriented preset, the canonical reference target is gas-phase
HF/6-31G\* MK ESP (three rigid orientations per model) with Amber RESP fitting.
The completed reference campaign covers **58/58 models** (`51` expanded C16/C24/C54
+ `7` extended functional-group controls; `C54_EPX_EPX_near_trans` completed in
job `1479457`), with **16** separate ωB97X-D/6-311++G(d,p) DFT-RESP controls
retained strictly as comparison benchmarks. No `C54` model entered `strict40`
calibration or hyperparameter selection (`all43` remains excluded as a diagnostic
set). External RHF-to-UHF instability is retained as a scientific flag and is not
automatic Amber incompatibility; see the
[Amber RESP acceptance policy](../research/parameterization/docs/AMBER_RESP_ACCEPTANCE.md).

Across the **14 `C54` size-holdout models** (evaluated on the same HF/6-31G\* MK
ESP grids), the pre-`L0_v2` research candidates exhibit specific pathologies that
**block them from production promotion** prior to the Block 2 `L0_v2` / `GO_gaff2_v3`
refits:

| Method | Category | C54-14 Median ESP RRMS | Blocking Defect / Status (2026-11-19) |
|---|---|---:|---|
| `Atomic_HF_RESP_closed` | Per-model in-sample QM reference | `6.87%` | In-sample reference only; not a transferable predictor |
| `Local_BCC_ESP_strict40` | Graph bond-increment (`strict40`) | `14.68%` | **Blocked**: rank-deficient normal matrix (`35/36`, `cond = 4.26e+18`); unphysical `COO_anchor = +0.839..+0.923 e` |
| `SQE_ESP_direct_strict40` | Static SQE fitted to ESP (`strict40`) | `16.02%` | **Blocked**: ill-conditioned (`augmented_rank = 24`, singular values `0.010..0.058`) and unphysical `role_theta_map` weights (`-7.69..+3.32 eV`) compensating for `L0_strict40` |
| `L0_strict40_closed` | Fixed role library (`strict40`) | `16.55%` | **Blocked**: severe `-COO-` charge pathology (`COO_C = +1.716076`, `COO_O = -1.030563`, `COO_anchor = -0.654951 e`) |
| `resp-current` (`current_fixed_role_library`) | Legacy user library (`resources/reference/current`) | `28.79%` | **Operational CLI default** until `L0_v2` + `GO_gaff2_v3` pass Block 2 gates |
| `SQE_L0_strict40` | Static SQE fitted to `L0`+RESP atomic charges | `46.46%` | **Blocked**: improves charge RMSD (`0.125 e`) but degrades C54 ESP RRMS on `8/13`–`9/14` models (`53.98%` on neutral C54) |

Note on `GAFF2_ABCG2` (AmberTools 25): audited across **58/58** models as an
independent coupled comparator (`-pl 8` and single-point AM1 for `N > 55`), it
remains a reference control alongside `AM1-BCC` rather than a replacement for
the RESP-based `L0_v2` derivation.

All final charges remain fixed during MD. Structure-specific response is a
preparation step and must preserve the chemically determined integer total
charge, which can be nonzero. Whole-sheet QM convergence diagnostics were
stopped on 2026-11-19 (`8` Slurm jobs cancelled) and remain open for Block 3
with closed-shell constrained `xTB` (`S1`/`S2`) pre-relaxation. Experimental
constraints and coupled LJ effects are discussed in
`research/parameterization/docs/EXPERIMENTAL_CHARGE_CONSTRAINTS.md` and
`EXPERIMENTAL_NONBONDED_BENCHMARKS.md`.
