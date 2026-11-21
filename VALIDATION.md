# Software and Force-Field Validation Status

## Local Software Regression Suite
- Verified with `uv sync --frozen` + `uv run pytest -q`: **`407 passed`** (covers Blocco 1 graph/geometry construction, placement failure handling, anti-grafting, covalent functionalization, external-engine preparation, Blocco 2 `gaff2_v3` / `L0_v2` graph perception, 4-color DSATUR epoxide coloring, `PROVENANCE.json` SHA-256 verification, precision GROMACS `.top`/`.itp` export, and pure-Python OpenMM XML export).
- Quickstart end-to-end smoke test (`uv run python examples/quickstart/run_quickstart.py`) passes and verifies `GO_gaff2_v3.frcmod`, 18 `v3_libraries/*_gaff2_v3.lib`, `PROVENANCE.json`, and `preset == "gaff2_v3"`.

## Block 2 Core AMBER Validation (`gaff2_v3` + `L0_v2`, 2026-11-20/21)
- **Go/No-Go Gates (`step2b`)**: **11/11 PASS** (`all_gates_passed: true`; bonded terms fitted strictly on `strict37` + 2 C16 torsion scans):
  - `0/129` AMBER namespace collisions across 18 lowercase GO-only types (`ga, gb, gc, gd, ge, gh, gk, gm, go, gp, gq, gu, gw, gx, gy, gz, g1, g6`);
  - `0` imaginary harmonic modes across `52/52` in-domain models (`max min F_RMS = 1.88e-4 kJ/(mol·nm) < 0.01`);
  - `ALL_52`: `F_med = 14.78`, `F_max = 25.39 kcal/(mol·Å)`, `H_med = 17.64%`, `G_med = 0.100 Å`, unscaled `Freq_med = 35.0 cm⁻¹`;
  - Size-transfer `C54` panel (excluded from all parameter fits; used for `L0_v2` candidate selection in `step2a`; unseen large sheets scheduled in Block 3): `C54_neutral_9` (`F_med = 17.32`, `F_max = 21.41`) and `C54_anionic_5` (`F_med = 12.61`, `F_max = 17.78`);
  - `pyrene_COOH` torsion MAE `0.0445 kcal/mol` (`pyrene_COO`: `0.0123 kcal/mol`);
  - `52/52` AMBER `sander` vs OpenMM single-point parity (`max|ΔE| = 0.0141 kJ/mol`, `max|ΔF| = 0.0317 kJ/(mol·nm)`);
  - `4/4` synthetic + `ff19SB+TIP3P+Na` stress tests PASS.
- **HF/6-31G\* MK ESP (`L0_v2c_a030`, `step2a`/`step2c`)**: `strict40` median `18.85%`, `C54_all_14` `15.35%`, `C54_neutral_9` `17.05%`, `C54_anionic_5` `8.85%`; permutation & Kekulé invariant (`< 1e-12 e`); integer net charge conserved.
- **6-Suite Coupled Cross-Matrix (`step2c`)**: exact `1e-17..1e-9` parity checks; `52/58` in-domain supported and `6/58` Priority-3 explicitly `UNSUPPORTED`; `v3_L0v2` improves π–π stacking (`MAE = 1.231 kcal/mol` vs `v3_current 2.154`) and outperforms `v3_ABCG2` across dense-OH models.
- **`WORKBENCH` Topology Parity & 58-Model Cross-Check (`step2d`, `ciamician`)**:
  - **Evaluators**: native **`sander` Version 24.0 (`amber24u3+at25`, module `amber/25`)** & `tleap`, **OpenMM 8.6.1** (`AmberPrmtopFile` native vs pure-Python `export_openmm` XML across `NoCutoff`, `CutoffNonPeriodic`, and boxed `PME`), and **OpenMM `GromacsTopFile`** reading `.top`/`.itp` on exact `.inpcrd` coordinates (native `gmx` binary is not installed on the `ciamician` cluster).
  - `OpenMM pure-Py XML` vs `OpenMM native`: `max|ΔE| < 1e-6 kJ/mol`, `max|ΔF| < 5e-5 kJ/(mol·nm)`;
  - `GROMACS .top/.itp` (exact `.inpcrd` coords) vs `OpenMM native`: `max|ΔE| = 5.43e-6 kJ/mol`, `max|ΔF| = 3.36e-5 kJ/(mol·nm)`;
  - **`sander` Version 24.0 (`amber24u3+at25`)** vs `OpenMM native`: `max|ΔE| = 0.0267 kJ/mol`, `max|ΔF| = 0.0212 kJ/(mol·nm)`;
  - `WORKBENCH` vs `step2c` across all 58 benchmark models: `52/52` in-domain accepted, `6/6` Priority-3 rejected, `0.0 e` maximum charge discrepancy.

## Remaining Unvalidated Domains (Block 3+)
- Whole-sheet closed-shell `xTB`/QM benchmarks (`Block 3`), proximity/redistribution charge models (`L1_reg`/`SQE`, `Block 4`), periodic covalent bonded exports (`Block 5`), and non-AMBER force fields (`OPLS-AA`, `OpenFF`, `CHARMM`, `Block 6`) remain outside current production validation claims. Historical external logs remain outside the git tree (`research_trials/`).
