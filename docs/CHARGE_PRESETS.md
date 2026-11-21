# Charge presets and the RESP-first study

*Status verified as of 2026-11-21* (Block 2 validated: `gaff2_v3` / `L0_v2c_a030`
bundle in `resources/reference/gaff2_v3/PROVENANCE.json` and the `step2a`–`step2d`
cross-validation logs).

The source catalog is `resources/charge_presets.json` (`schema_version: 3`). It
records charge family, coupled parameter expectations, provenance, and readiness.
**`resp-l0-v2` (`gaff2_v3`, `L0_v2c_a030` under `resources/reference/gaff2_v3`) is
now the validated and only recommended AMBER CLI default** (`03_prepare_amber.py`
/ `--preset gaff2_v3`), while **`resp-current` (`gaff2_current` under
`resources/reference/current`)** is **deprecated** (known bonded/charge defects)
and kept untouched only for archival backwards compatibility (`--preset gaff2_current`
or implicit `--reference-dir` without `--preset`). Expert per-atom external overrides
are supported via `--charge-model custom --charges-json` and logged with file
hashes in `PROVENANCE.json`.

| ID | Charge model | Status (2026-11-21) |
|---|---|---|
| `resp-l0-v2` | Block 2 dual-regime `L0_v2` (`L0_v2c_a030`: fixed roles when `N_COO == 0`, `exp(-0.30 * d_graph)` sp² delocalization `d_max = 10`, `delta_sheet = -0.345574 e` when `N_COO > 0`) | **Validated & recommended AMBER CLI default** (`18.85%` `strict40`, `15.35%` `C54_all_14`, `17.05%` `C54_neutral_9`, `8.85%` `C54_anionic_5` HF/6-31G\* MK ESP RRMS; physical `CBD`: `COO_C = +0.835809`, `COO_O = -0.799783`, `COO_anchor = +0.109331 e`) |
| `resp-current` | Legacy user RESP-derived fragment libraries (`resources/reference/current`) | **Deprecated archival baseline** (kept only for `--preset gaff2_current` regression compatibility; `32.15%` `strict40`, `28.79%` `C54_all_14`, `32.48%` `C54_neutral_9`, `18.43%` `C54_anionic_5` ESP RRMS) |
| `resp-go-new` | Early pre-`L0_v2` Block 1/2 drafts (`L0_strict40`, `Local_BCC_ESP`, `SQE_*`) | **Blocked archival baselines** superseded by `L0_v2`; known pathologies preserved below |
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

## Verified 58/58 Reference Panel: Promoted `L0_v2` vs Controls & Superseded Drafts

For the Amber-oriented preset, the canonical reference target is gas-phase
HF/6-31G\* MK ESP (three rigid orientations per model) with Amber RESP fitting.
The completed reference campaign covers **58/58 models** (`51` expanded C16/C24/C54
+ `7` extended functional-group controls; `52/58` in-domain Priority 1+2 GO and
`6/58` Priority-3 explicitly rejected with `ValueError`), with **16** separate
ωB97X-D/6-311++G(d,p) DFT-RESP controls retained as comparison benchmarks. No `C54` model entered `strict40` parameter fitting (`all43` remains excluded
as a diagnostic set). Because the `14 C54` models were also inspected when
selecting among the seven `L0_v2` regularization/cutoff candidates in `step2a`,
they serve here as a **size-transfer evaluation set** rather than an unseen final
blind test (completely unseen sheets are deferred to Block 3). External RHF-to-UHF
instability is retained as a scientific flag and is not automatic Amber
incompatibility; see the
[Amber RESP acceptance policy](../research/parameterization/docs/AMBER_RESP_ACCEPTANCE.md).

Across those **14 `C54` size-transfer models** (`9` neutral + `5` anionic on the
same HF/6-31G\* MK ESP grids), Block 2 evaluated `L0_v2` (`L0_v2c_a030`) alongside
`current_lib`, `GAFF2_ABCG2`, and the blocked pre-`L0_v2` candidates:

| Method | Category | C54-14 Med (Neu-9 / Ani-5) ESP RRMS | Status / Notes (2026-11-21) |
|---|---|---:|---|
| `Atomic_HF_RESP_closed` | Per-model in-sample QM reference | `6.87%` (`7.68%` / `1.20%`) | In-sample reference only; not transferable |
| **`resp-l0-v2` (`L0_v2c_a030`)** | **Block 2 dual-regime (`strict40`-fit)** | **`15.35%`** (`17.05%` / `8.85%`) | **Promoted CLI default** (`strict40 = 18.85%`); passes all physical charge gates (`COO_C = +0.835809`, `COO_O = -0.799783`, `COO_anchor = +0.109331 e`), exact permutation & Kekulé invariance (`< 1e-12 e`), and integer charge conservation |
| `Local_BCC_ESP_strict40` | Graph bond-increment (`strict40`) | `14.68%` | **Blocked**: rank-deficient normal matrix (`35/36`, `cond = 4.26e+18`); unphysical `COO_anchor = +0.839..+0.923 e` |
| `SQE_ESP_direct_strict40` | Static SQE fitted to ESP (`strict40`) | `16.02%` | **Blocked**: ill-conditioned (`augmented_rank = 24`, singular values `0.010..0.058`) and unphysical `role_theta_map` weights (`-7.69..+3.32 eV`) |
| `L0_strict40_closed` | Early fixed role library (`strict40`) | `16.55%` (`18.06%` / `16.31%`) | **Blocked**: severe `-COO-` charge pathology (`COO_C = +1.716076`, `COO_O = -1.030563`, `COO_anchor = -0.654951 e`) |
| `resp-current` (`current_lib`) | Legacy library (`resources/reference/current`) | `28.79%` (`32.48%` / `18.43%`) | **Retained backwards-compatible legacy preset** (`strict40 = 32.15%`) |
| `GAFF2_ABCG2` | AmberTools 25 ABCG2 control | `39.53%` (`55.32%` / `8.69%`) | Comparison control (`strict40 n=37 = 37.54%`); good on anionic C54 (`8.69%`) and dense-OH forces (`C54_OH6 12.70`), but poor on neutral C54 ESP (`55.32%`) and worse on `C54_OH4` forces (`24.37` vs `v3_L0v2 18.37`) |
| `SQE_L0_strict40` | Static SQE fitted to `L0`+RESP charges | `46.46%` | **Blocked**: degrades C54 ESP RRMS on `8/13`–`9/14` models (`53.98%` on neutral C54) |

All final charges remain fixed during MD. Structure-specific response is a
preparation step and must preserve the chemically determined integer total
charge, which can be nonzero. Whole-sheet QM convergence diagnostics were
stopped on 2026-11-19 (`8` Slurm jobs cancelled) and remain open for Block 3
with closed-shell constrained `xTB` (`S1`/`S2`) pre-relaxation. Experimental
constraints and coupled LJ effects are discussed in
`research/parameterization/docs/EXPERIMENTAL_CHARGE_CONSTRAINTS.md` and
`EXPERIMENTAL_NONBONDED_BENCHMARKS.md`.
