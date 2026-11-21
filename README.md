<h1 align="center">GO Workbench</h1>

<p align="center">
  <em>Linear-scaling graphene oxide builder, validated AMBER parameterization, and multi-engine topology export for docking and molecular dynamics.</em>
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.24119849"><img src="https://img.shields.io/badge/Archive-10.5281%2Fzenodo.24119849-2457A7?style=flat" alt="Archive DOI: 10.5281/zenodo.24119849"></a>
  <a href="https://github.com/CalzaAlessandro/go-surface-workbench/actions/workflows/tests.yml"><img src="https://img.shields.io/github/actions/workflow/status/CalzaAlessandro/go-surface-workbench/tests.yml?branch=main&amp;style=flat" alt="Tests"></a>
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-%E2%89%A53.11-58616B?style=flat" alt="Python: 3.11 or later"></a>
  <img src="https://img.shields.io/badge/Status-Milestone%20block--2--amber--core-2E8B57?style=flat" alt="Status: Milestone block-2-amber-core">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache--2.0-607D8B?style=flat" alt="License: Apache-2.0"></a>
</p>

<p align="center">
  <img src="docs/images/banner.webp" width="100%" alt="Illustration of a graphene oxide blueprint becoming an assembled molecular model on blue drafting paper, with a ruler, pencil and eraser">
</p>

**GO Workbench** constructs chemically valid finite and periodic graphene oxide (GO) models in $O(N)$ time up to $24{,}500+$ atoms, assigns validated AMBER **`gaff2_v3` + `L0_v2`** parameters with `0/129` namespace collisions, and exports verified inputs and topologies for **molecular docking** (`PDB`, `MOL2`, `SDF`) and **vacuum or solvated MD** (**AMBER**, **GROMACS**, and **OpenMM**).

---

## Development Roadmap & Block Status

| Block | Scope | Status | Milestone Tag |
|---|---|---|---|
| **Block 1 — Structure Builder** | Finite & periodic GO (`rectangular`, `hexagonal`, `circular`), nm sizing, HierGO-style pores, `pH`-aware `C/O` & `O/C` solver, `random`/`island`/`dispersed` patterns, Priority 1+2 groups, $O(N)$ motif-aware gate, `JSON`/`PDB`/`MOL2`/`XYZ`/`EXTXYZ`/`SDF V3000`/`CIF P1` exports | **Done** | [`block-1-builder`](https://github.com/CalzaAlessandro/go-surface-workbench/releases/tag/block-1-builder) |
| **Block 2 — Core AMBER (`gaff2_v3` + `L0_v2`) & Multi-Engine Parity** | `GO_gaff2_v3.frcmod`, 18 lowercase GO types (`0/129` Amber collisions), 18 `*_gaff2_v3.lib`, 4-color DSATUR epoxide typing, `L0_v2c_a030` + custom charge JSON, SHA-256 `PROVENANCE.json`, precision GROMACS `.top`/`.itp` & pure-Python OpenMM `.xml` exporters (`11/11` gates, `407 passed`) | **Done** | [`block-2-amber-core`](https://github.com/CalzaAlessandro/go-surface-workbench/releases/tag/block-2-amber-core) |
| **Block 3 — Realistic Sheet QM/RESP Validation** | Closed-shell Clar/Hückel blind mid-size (`80–87` atoms) & realistic (`185–190` atoms) GO sheets on HPC (gas & PCM DFT + HF/6-31G* RESP) | **In progress** | — |
| **Block 4 — Advanced Charges, `gaff2_v3.1` & Irregular Islands** | Proximity charges (`L1_reg`/`SQE`), additive Priority-3 motifs (phenols, furanic ethers, lactones), organic irregular `island` mode & optional `clar`/Hückel filter | **Planned** | — |
| **Block 5 — Periodic Topologies** | Direct periodic topology writers from `model.json` bond images for GROMACS, OpenMM, AMBER, and OpenFF | **Planned** | — |
| **Block 6 — Multi-Force-Field & Literature Catalog** | `OPLS-AA`, `OpenFF`, and `CHARMM36m`/CGenFF presets alongside a tiered literature catalog | **Planned** | — |
| **Block 7 — Condensed-Phase MD, MM-GBSA/PBSA & Biomolecules** | Explicit-water MD, 3-level MM-GBSA/PBSA calibration + `idecomp`, docking (`Vina`/`PatchDock`) & plug-and-play downstream solvation (`gmx solvate`, OpenMM `Modeller.addSolvent`), Linux + macOS CI | **Planned** | — |

---

## Installation & Quickstart

Requires **Python ≥ 3.11** (standard-library `tomllib`; dependencies managed via [pyproject.toml](pyproject.toml) and `uv.lock`: NumPy, SciPy, ParmEd, RDKit). From the repository root:

```bash
uv sync --frozen
uv run python examples/quickstart/run_quickstart.py --output-dir work/quickstart
```

When installed from a wheel or sdist, all stages are also exposed through the module entry point:

```bash
python -m go_builder --help
```

External scientific engines (`AmberTools`/`tleap`, `sander`, `GROMACS`, `OpenMM`, `Gaussian`/`xTB`, `AutoDock Vina`/`GNINA`) remain external to the Python package; see [Engine setup](docs/ENGINE_SETUP.md) and [Quickstart details](docs/QUICKSTART.md).

---

## Core Workflows

### 1. Build, Check & Export GO Structures (`01_build_GO.py`)

Construct models directly from CLI flags, commented `.toml`/`.json` configs ([`examples/configs/`](examples/configs/)), or packaged presets:

```bash
# Direct CLI: 4x4 nm GO at pH 7.4, 2 islands, 5% mixed pores, multi-format export
uv run python 01_build_GO.py \
  --size 4 4 --oc 0.35 --ph 7.4 --place island --islands 2 --holes mixed:5 --seed 42 \
  --out work/go.json --pdb work/go.pdb --mol2 work/go.mol2 --sdf work/go.sdf --cif work/go.cif \
  --log work/go.log

# Hexagonal zigzag nanoflake (2.5 nm diameter) for QM / xTB / ASE workflows
uv run python 01_build_GO.py \
  --shape hexagonal --edge zigzag --diam 2.5 --oc 0.25 --place dispersed --seed 17 \
  --out work/hex.json --extxyz work/hex.extxyz --sdf work/hex.sdf

# Config templates, shipped TOML recipes, packaged presets, and strict gate check
uv run python 01_build_GO.py --init-config work/my_go.toml
uv run python 01_build_GO.py --config examples/configs/go_island_porous.toml --out work/island.json
uv run python 01_build_GO.py --list-presets
uv run python check_model.py --model work/go.json --strict
```

Every built model preserves exact constructed connectivity (`geometry_status = "constructed_unrelaxed"`) and passes the $O(N)$ motif-aware integrity gate (`0` clashes, bond lengths, skew-cell minimum image, valence & charge checks) with single-core CPU scaling:

| System | Atoms | Build time | Gate time | Clashes |
|---|---:|---:|---:|---:|
| Finite `12×12` (`O/C=0.35`, pH 7.4) | `913` | `0.16 s` | `0.03 s` | `0` |
| Finite `24×24` (`O/C=0.35`, pH 7.4) | `3,666` | `0.67 s` | `0.12 s` | `0` |
| Finite `36×36` (`O/C=0.35`, pH 7.4) | `8,294` | `1.60 s` | `0.26 s` | `0` |
| Periodic `64×64` (`O/C=0.30`) | `24,576` | `6.30 s` | `1.33 s` | `0` |

### 2. Parameterize for AMBER (`gaff2_v3` + `L0_v2`) & Export to GROMACS / OpenMM

Run the core pipeline `01_build_GO.py` → [`02_open_epoxides.py`](02_open_epoxides.py) (optional) → [`03_prepare_amber.py`](03_prepare_amber.py) → [`07_export_gromacs.py`](07_export_gromacs.py) / [`08_export_openmm.py`](08_export_openmm.py):

```bash
uv run python 01_build_GO.py --recipe examples/finite_GO.json --output work/finite_GO.json --pdb work/finite_GO.pdb --json
uv run python 03_prepare_amber.py --list-presets
uv run python 03_prepare_amber.py --model work/finite_GO.json --output-dir work/amber_run --radii mbondi2 --run

# Verified cross-engine exports from Amber topology + coordinates:
uv run python 07_export_gromacs.py --topology work/amber_run/go.prmtop --coordinates work/amber_run/go.inpcrd --output-prefix work/gmx/go
uv run python 08_export_openmm.py --topology work/amber_run/go.prmtop --coordinates work/amber_run/go.inpcrd --output work/openmm/go.xml
```

- **In-domain coverage (`gaff2_v3`)**: finite Priority 1+2 GO chemistry (`basal_C`, `edge_CH`, `OH`, basal `epoxide`, `epoxide_eph`, `epoxide_ept`, `carbonyl`, `COOH`, `COO−`) using `GO_gaff2_v3.frcmod`, 18 lowercase GO `.lib` units (`0/129` AMBER namespace collisions), 4-color DSATUR epoxide typing, and dual-regime `L0_v2c_a030` charges. Priority-3 motifs and unparameterized junctions are explicitly rejected.
- **Custom charges & archival fallback**: pass `--charge-model custom --charges-json <file>` to inject external per-atom charges; `--preset gaff2_current` is deprecated and kept solely for archival backwards compatibility.
- **Docking & MD regimes**: `PDB` (`CONECT`), `MOL2` (Tripos + charges), and `SDF V3000` support downstream docking (`AutoDock Vina`, `PatchDock`, `assemble_complex`); AMBER/GROMACS/OpenMM exporters support both **vacuum/gas-phase** and **explicitly solvated** systems (prepared in `tleap` via `solvateBox`/`addIonsRand` alongside standard biomolecular force fields). Direct `gmx solvate` split `.itp` helpers and OpenMM `ForceField` XML templates are scheduled in Block 7.

---

## Validation Summary

- **Software regression suite**: **`407 passed`** (`uv run pytest -q`), covering graph construction, motif gate, covalent functionalization, wheel/sdist packaging, `gaff2_v3` / `L0_v2` perception, `PROVENANCE.json` SHA-256 (`20/20` files), OpenMM XML serialization, and GROMACS `.top`/`.itp` + `SETTLE` round-trip verification.
- **AMBER core parameterization (`gaff2_v3` + `L0_v2`)**: **`11/11` Go/No-Go gates (`step2b`) PASS**; `0` imaginary modes across `52/52` in-domain models (`6/58` Priority-3 explicitly `UNSUPPORTED`); `0/129` AMBER namespace collisions; `sander` 24.0 vs OpenMM 8.6.1 vs GROMACS parity **`max|ΔE| = 0.0267 kJ/mol`** and **`max|ΔF| = 0.0317 kJ/(mol·nm)`** (`4/4` synthetic + `ff19SB+TIP3P+Na` stress tests PASS).

Full tables, evaluator details, and current scope limits are recorded in [VALIDATION.md](VALIDATION.md).

---

## Documentation

| Topic | Document |
|---|---|
| Installation & quickstart outputs | [docs/QUICKSTART.md](docs/QUICKSTART.md) |
| CLI flags, recipes & full workflow reference | [docs/REFERENCE.md](docs/REFERENCE.md) |
| Validation gates & parity results | [VALIDATION.md](VALIDATION.md) |
| Charge models (`L0_v2` & custom JSON) | [docs/CHARGE_PRESETS.md](docs/CHARGE_PRESETS.md) |
| Force-field presets & reference provenance | [docs/REFERENCE_AND_PRESETS.md](docs/REFERENCE_AND_PRESETS.md) |
| Feature comparison vs literature | [docs/COMPARISON.md](docs/COMPARISON.md) |
| Linear-scaling benchmarks | [docs/PERFORMANCE.md](docs/PERFORMANCE.md) |
| External simulation & QM engines | [docs/ENGINE_SETUP.md](docs/ENGINE_SETUP.md) |
| Covalent functionalization fragments | [docs/FRAGMENTS.md](docs/FRAGMENTS.md) |
| Trajectory alignment utilities | [docs/ALIGNMENT.md](docs/ALIGNMENT.md) |
| External dataset notices & licenses | [THIRD_PARTY.md](THIRD_PARTY.md) |

---

## Citation, Contributing & License

If you use **GO Workbench** or its parameter libraries, please cite:

> Alessandro Calza (2026). *GO Workbench*, version 0.1.0. [DOI: 10.5281/zenodo.24119849](https://doi.org/10.5281/zenodo.24119849) · [CITATION.cff](CITATION.cff)

<details>
<summary>BibTeX</summary>

```bibtex
@software{Calza2026GOWorkbench,
  author = {Alessandro Calza},
  title = {GO Workbench},
  version = {0.1.0},
  year = {2026},
  doi = {10.5281/zenodo.24119849},
  url = {https://github.com/CalzaAlessandro/go-surface-workbench}
}
```

</details>

- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md) and open issues on [GitHub Issues](https://github.com/CalzaAlessandro/go-surface-workbench/issues).
- **License**: Original code and documentation in this repository are licensed under the [Apache 2.0 License](LICENSE). Third-party tools and reference datasets under `resources/reference/` retain their upstream terms documented in [THIRD_PARTY.md](THIRD_PARTY.md).
