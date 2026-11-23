# Software and Force-Field Validation Status

## Automated Regression Suite
- Verified with `uv sync --frozen` + `uv run pytest -q`: **420 passing tests** covering molecular graph and geometry construction, placement failure handling, anti-grafting stereochemistry, covalent epoxide functionalization, external-engine preparation, `gaff2_v3` / `L0_v2` graph perception, 4-color DSATUR epoxide coloring, `PROVENANCE.json` SHA-256 verification, precision GROMACS `.top`/`.itp` export, pure-Python OpenMM XML export, `irregular_island` Lerf–Klinowski placement, and `--clar` closed-shell Kekulé / Hückel $\pi$-gap screening.
- End-to-end quickstart execution (`uv run python examples/quickstart/run_quickstart.py`) passes and verifies `GO_gaff2_v3.frcmod`, 18 `v3_libraries/*_gaff2_v3.lib`, `PROVENANCE.json`, and `preset == "gaff2_v3"`.

## Core Amber Force-Field Validation (`gaff2_v3` + `L0_v2`)
- **Parameterization acceptance criteria:** **11/11 benchmarks passed** (bonded parameters fitted on `37` core finite models plus 2 pyrene-carboxyl torsion scans):
  - `0/129` Amber atom-type namespace collisions across 18 lowercase GO-only types (`ga, gb, gc, gd, ge, gh, gk, gm, go, gp, gq, gu, gw, gx, gy, gz, g1, g6`);
  - `0` imaginary harmonic modes across `52/52` in-domain finite models (`max min F_RMS = 1.88e-4 kJ/(mol·nm) < 0.01`);
  - Complete `52`-model in-domain benchmark: `F_med = 14.78`, `F_max = 25.39 kcal/(mol·Å)`, `H_med = 17.64%`, `G_med = 0.100 Å`, unscaled `Freq_med = 35.0 cm⁻¹`;
  - Out-of-fit size-transfer `C54` test panel (excluded from bonded parameter fitting): `C54_neutral_9` (`F_med = 17.32`, `F_max = 21.41`) and `C54_anionic_5` (`F_med = 12.61`, `F_max = 17.78`);
  - `pyrene_COOH` torsion MAE `0.0445 kcal/mol` (`pyrene_COO`: `0.0123 kcal/mol`);
  - `52/52` Amber `sander` vs OpenMM single-point energy/force parity (`max|ΔE| = 0.0141 kJ/mol`, `max|ΔF| = 0.0317 kJ/(mol·nm)`);
  - `4/4` synthetic dense multi-group + `ff19SB+TIP3P+Na` co-loading stress tests passing.
- **Electrostatic potential (`HF/6-31G*` MK ESP, `L0_v2`):** `40`-model training median RRMS `18.85%`; unseen `14`-model `C54` panel `15.35%` (`C54_neutral_9` `17.05%`, `C54_anionic_5` `8.85%`); permutation- and Kekulé-invariant (`< 1e-12 e`); exact integer formal-charge conservation.
- **Multi-suite cross-comparison:** Exact `1e-17..1e-9` parity checks; `52/58` finite reference models supported in-domain and `6/58` unsupported edge motifs cleanly rejected; `gaff2_v3` + `L0_v2` improves $\pi$–$\pi$ stacking energetics (`MAE = 1.231 kcal/mol` vs legacy `2.154`) and improves dense-OH electrostatic accuracy over AM1-BCC/ABCG2.
- **Cross-engine topology parity (`58`-model verification):**
  - **Evaluated engines:** native Amber **`sander` 24.0 (`AmberTools 25`)** and `tleap`, **OpenMM 8.6.1** (`AmberPrmtopFile` native vs pure-Python `export_openmm` XML across `NoCutoff`, `CutoffNonPeriodic`, and boxed `PME`), and **OpenMM `GromacsTopFile`** reading exported `.top`/`.itp` on exact `.inpcrd` coordinates;
  - `OpenMM pure-Py XML` vs `OpenMM native`: `max|ΔE| < 1e-6 kJ/mol`, `max|ΔF| < 5e-5 kJ/(mol·nm)`;
  - `GROMACS .top/.itp` (exact `.inpcrd` coords) vs `OpenMM native`: `max|ΔE| = 5.43e-6 kJ/mol`, `max|ΔF| = 3.36e-5 kJ/(mol·nm)`;
  - **`sander` 24.0** vs `OpenMM native`: `max|ΔE| = 0.0267 kJ/mol`, `max|ΔF| = 0.0212 kJ/(mol·nm)`;
  - Across all 58 benchmark models: `52/52` in-domain accepted, `6/6` unsupported edge-motif models rejected, `0.0 e` maximum charge discrepancy.

## Domains Under Active Development
- Larger-sheet (`80–190` atoms) closed-shell QM electrostatic/vibrational benchmarks and charge regularization refinements (`L1_reg`/`SQE`) remain in progress.
- Extended edge-oxygen motifs (phenolic `sp2`-OH, planar furanic ethers, and edge lactones) are not yet included in `gaff2_v3` and are explicitly rejected until validated parameters are finalized.
- Periodic covalent MD benchmarks across external engines and non-Amber force fields (`OPLS-AA`, `OpenFF`, `CHARMM36m`) remain outside current production validation claims.
