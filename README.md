<h1 align="center">GO Workbench</h1>

Generate graphene oxide sheets with explicit connectivity, assign Amber parameters and export topologies for Amber, OpenMM or GROMACS.

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.24119849"><img src="https://img.shields.io/badge/Archive-10.5281%2Fzenodo.24119849-2457A7?style=flat" alt="Archive DOI: 10.5281/zenodo.24119849"></a>
  <a href="https://github.com/CalzaAlessandro/go-surface-workbench/actions/workflows/tests.yml"><img src="https://img.shields.io/github/actions/workflow/status/CalzaAlessandro/go-surface-workbench/tests.yml?branch=main&amp;style=flat" alt="Tests"></a>
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-%E2%89%A53.11-58616B?style=flat" alt="Python: 3.11 or later"></a>
  <img src="https://img.shields.io/badge/Status-Active%20development-9A6700?style=flat" alt="Status: active development">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache--2.0-607D8B?style=flat" alt="License: Apache-2.0"></a>
</p>

<p align="center">
  <img src="docs/images/banner.webp" width="100%" alt="Illustration of a graphene oxide blueprint becoming an assembled molecular model on blue drafting paper, with a ruler, pencil and eraser">
</p>

## Contents

- [Getting started](#getting-started)
- [Usage and examples](#usage-and-examples)
- [Data and files](#data-and-files)
- [Methods and validation](#methods-and-validation)
- [Documentation](#documentation)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

> [!NOTE]
> Structure generation and the core Amber parameterization (`gaff2_v3` with `L0_v2` charges) are validated for finite graphene oxide sheets containing Priority 1 and 2 oxygen groups. Larger-sheet charge validation, periodic topology writers and non-Amber force fields are in progress; see [VALIDATION.md](VALIDATION.md) before running production simulations.

## Getting started

Requires Python 3.11 or later. Dependencies are listed in [pyproject.toml](pyproject.toml) and `uv.lock` (NumPy, SciPy, ParmEd and RDKit). From the repository root:

```bash
uv sync --frozen
uv run python examples/quickstart/run_quickstart.py --output-dir work/quickstart
```

When installed from a wheel or source distribution, all stages are available through `python -m go_builder --help`. External programs such as AmberTools (`tleap`, `sander`), GROMACS, OpenMM, Gaussian, xTB and AutoDock Vina are optional at install time and must be provided by the host environment when executing simulations or quantum-chemical jobs.

[Quickstart and expected outputs](docs/QUICKSTART.md) ·
[Installation details](docs/REFERENCE.md#installation--environment-setup) ·
[Engine setup](docs/ENGINE_SETUP.md)

## Usage and examples

Sheets can be constructed with [01_build_GO.py](01_build_GO.py) using command-line flags, TOML or JSON configuration files (examples in [`examples/configs/`](examples/configs/)), or packaged presets:

```bash
# Build a 4 x 4 nm GO sheet at pH 7.4 with two oxidized islands and 5% pores
uv run python 01_build_GO.py \
  --size 4 4 --oc 0.35 --ph 7.4 --place island --islands 2 --holes mixed:5 --seed 42 \
  --out work/go.json --pdb work/go.pdb --mol2 work/go.mol2 --sdf work/go.sdf --cif work/go.cif

# Inspect configuration recipes and verify a constructed model
uv run python 01_build_GO.py --init-config work/my_go.toml
uv run python 01_build_GO.py --config examples/configs/go_island_porous.toml --out work/island.json
uv run python 01_build_GO.py --list-presets
uv run python check_model.py --model work/go.json --strict
```

Models export to `JSON` (authoritative graph), `PDB` (with `CONECT` records), `MOL2`, `XYZ`, `EXTXYZ`, `SDF` (V3000) and `CIF`.

Parameterization and multi-engine conversion follow [03_prepare_amber.py](03_prepare_amber.py) → [07_export_gromacs.py](07_export_gromacs.py) or [08_export_openmm.py](08_export_openmm.py), with optional covalent epoxide opening via [02_open_epoxides.py](02_open_epoxides.py):

```bash
# Build, assign gaff2_v3 parameters + L0_v2 charges, and compile with tleap
uv run python 01_build_GO.py --recipe examples/finite_GO.json --output work/finite_GO.json --pdb work/finite_GO.pdb
uv run python 03_prepare_amber.py --model work/finite_GO.json --output-dir work/amber --radii mbondi2 --run

# Export the completed Amber topology and coordinates to GROMACS and OpenMM
uv run python 07_export_gromacs.py --topology work/amber/go.prmtop --coordinates work/amber/go.inpcrd --output-prefix work/gmx/go
uv run python 08_export_openmm.py --topology work/amber/go.prmtop --coordinates work/amber/go.inpcrd --output work/openmm/go.xml
```

User-supplied per-atom charges can be passed with `--charge-model custom --charges-json <file>`. The archival baseline `--preset gaff2_current` is deprecated and kept for backwards compatibility. Exported `PDB` and `MOL2` files support rigid or flexible receptor docking (`AutoDock Vina`, `PatchDock`), while Amber, GROMACS and OpenMM outputs cover both gas-phase runs and explicitly solvated systems prepared in `tleap`.

## Data and files

| What you need | Where to find it |
|---|---|
| Constructed structure and graph | `work/finite_GO.json`, `finite_GO.pdb` (or `.mol2`, `.sdf`, `.xyz`, `.cif`) |
| Ready Amber sheet | `work/amber/go.prmtop` and `go.inpcrd` (after `--run`) |
| Amber input scripts and residue units | `work/amber/tleap.in`, `GO_gaff2_v3.frcmod`, `v3_libraries/` |
| Charge assignment and checksum records | `work/amber/assigned_model.json`, `library_audit.json`, `PROVENANCE.json` |
| Ready GROMACS topology | `work/gmx/go.top`, `go.itp` and `go.gro` |
| Ready OpenMM system | `work/openmm/go.xml` |

All CLI stages and command-line flags are listed in the [command reference](docs/REFERENCE.md#scripts--usage-workflow).

## Methods and validation

Finite and periodic sheets (`rectangular`, `hexagonal`, `circular`) are generated with specified lateral dimensions, pH-aware `C/O` or `O/C` compositions, vacancy distributions and spatial oxidation modes (`random`, `island`, `dispersed`). Every built structure keeps explicit connectivity (`geometry_status = "constructed_unrelaxed"`) and passes pre-export checks for steric clashes, bond lengths, minimum-image distances, valences and net charge. Single-core construction runs in $O(N)$ time, requiring `0.16 s` for `913` atoms (`12×12`, `O/C = 0.35`, pH `7.4`), `0.67 s` for `3,666` atoms (`24×24`), `1.60 s` for `8,294` atoms (`36×36`), and `6.30 s` for a `24,576`-atom periodic `64×64` sheet, with zero steric clashes. Geometric checks verify graph and local coordination validity; they do not establish electronic structure or conformational equilibration.

The default `gaff2_v3` preset covers finite GO models bearing basal carbons, aromatic CH edges, hydroxyls, basal and edge epoxides, carbonyls, carboxylic acids and carboxylates, using 18 dedicated lowercase Amber atom types (`0/129` collisions with `ff19SB`, `GAFF2`, water and ions), 4-color epoxide typing, and dual-regime `L0_v2` charges. Unsupported Priority-3 motifs (phenols, planar ethers, lactones) and unparameterized junctions are rejected cleanly. Across 52 in-domain benchmark models (`6/58` Priority-3 models excluded by design), harmonic frequency analysis gives zero imaginary modes, Amber `sander` 24.0, OpenMM 8.6.1 and GROMACS energies agree within `0.0267 kJ/mol` (`max|ΔF| = 0.0317 kJ/(mol·nm)`), all 20 library files verify against SHA-256 checksums in `PROVENANCE.json`, and the regression suite passes `407` automated tests. Experimental parameterization studies in `research/` are kept separate from the reference libraries.

### Development status

| Block | Scope | Status | Milestone |
|---|---|---|---|
| Block 1 | Structure builder, placement modes, $O(N)$ integrity checks and multi-format exports | Done | [`block-1-builder`](https://github.com/CalzaAlessandro/go-surface-workbench/releases/tag/block-1-builder) |
| Block 2 | Core Amber `gaff2_v3` + `L0_v2` libraries, custom charges, GROMACS and OpenMM exporters (`407` tests) | Done | [`block-2-amber-core`](https://github.com/CalzaAlessandro/go-surface-workbench/releases/tag/block-2-amber-core) |
| Block 3 | Closed-shell QM optimization and HF/DFT RESP benchmarks on `80–190`-atom sheets | In progress | — |
| Block 4 | Proximity charge models (`L1_reg`/`SQE`), additive Priority-3 motifs (`gaff2_v3.1`), irregular island builder | Planned | — |
| Block 5 | Direct periodic topology writers for GROMACS, OpenMM, Amber and OpenFF | Planned | — |
| Block 6 | `OPLS-AA`, `OpenFF` and `CHARMM36m` presets and tiered literature catalog | Planned | — |
| Block 7 | Explicit-water MD, MM-GBSA/PBSA calibration, docking and direct solvation helpers, Linux/macOS CI | Planned | — |

## Documentation

| Topic | Guide |
|---|---|
| First example | [Quickstart](docs/QUICKSTART.md) |
| Commands and complete workflows | [Command reference](docs/REFERENCE.md) |
| Validation and limitations | [Validation status](VALIDATION.md) |
| Charge models and presets | [Charge presets](docs/CHARGE_PRESETS.md) |
| Force-field presets and references | [Force-field presets & references](docs/REFERENCE_AND_PRESETS.md) |
| Feature comparison vs literature | [Comparison](docs/COMPARISON.md) |
| Scaling timings | [Benchmarks](docs/PERFORMANCE.md) |
| Covalent functionalization | [Fragments](docs/FRAGMENTS.md) |
| External programs | [Engine setup](docs/ENGINE_SETUP.md) |
| Trajectory alignment | [Alignment](docs/ALIGNMENT.md) |
| Third-party notices and datasets | [Third-party licenses](THIRD_PARTY.md) |

## Citation

If you use this software or its parameter libraries, cite the archived release:

Alessandro Calza (2026). *GO Workbench*, version 0.1.0.
[DOI: 10.5281/zenodo.24119849](https://doi.org/10.5281/zenodo.24119849).

Citation metadata: [CITATION.cff](CITATION.cff).

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

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for setup and checks. Report problems or propose changes through [GitHub Issues](https://github.com/CalzaAlessandro/go-surface-workbench/issues). Include the affected file or command and a short description.

## License

Original source code, parameter libraries generated for this project and documentation are distributed under the [Apache License, Version 2.0 (Apache-2.0)](LICENSE). Third-party software and external reference datasets in `resources/reference/` retain their original terms and attribution; see [THIRD_PARTY.md](THIRD_PARTY.md).
