<a id="go-builder-go-surface-workbench"></a>
<h1 align="center">GO Workbench</h1>

Generate graphene oxide sheets with explicit connectivity and optional covalent functionalization.

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.24119849"><img src="https://img.shields.io/badge/Archive-10.5281%2Fzenodo.24119849-2457A7?style=flat" alt="Archive DOI: 10.5281/zenodo.24119849"></a>
  <a href="https://github.com/CalzaAlessandro/go-surface-workbench/actions/workflows/tests.yml"><img src="https://img.shields.io/github/actions/workflow/status/CalzaAlessandro/go-surface-workbench/tests.yml?branch=main&amp;style=flat&amp;label=Tests" alt="Software tests"></a>
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-%E2%89%A53.11-58616B?style=flat" alt="Python: 3.11 or later"></a>
  <img src="https://img.shields.io/badge/Status-Active%20development-9A6700?style=flat" alt="Status: active development">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache--2.0-607D8B?style=flat" alt="License: Apache-2.0"></a>
</p>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/banner.png">
  <img src="docs/images/banner.png" width="100%" alt="Illustration of a graphene oxide blueprint becoming an assembled molecular model on blue drafting paper, with a ruler, pencil and eraser">
</picture>

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
> New parameters, charges and topologies are preliminary.
> See the [validation status](VALIDATION.md) before using them in simulations.

## Getting started

Requires Python 3.11 or later; the bundled quickstart uses Python 3.12.
From the repository root:

```bash
uv sync --frozen
uv run python examples/quickstart/run_quickstart.py --output-dir work/quickstart
```

The example builds a small sheet, checks its connectivity and prepares Amber
inputs using Python dependencies. Install AmberTools to compile a topology;
later QM, docking and simulation steps use their respective external engines.

[Quickstart and expected outputs](docs/QUICKSTART.md) ·
[Installation details](docs/REFERENCE.md#installation--environment-setup) ·
[Engine setup](docs/ENGINE_SETUP.md)

## Usage and examples

Build and functionalize a GO sheet with:
`01_build_GO.py` → `02_open_epoxides.py` (optional) → `03_prepare_amber.py`.
Use `prepare_fragment.py` to prepare grafting reagents for reuse across loadings
and random seeds.

To run the quickstart steps separately:

```bash
uv run python 01_build_GO.py --recipe examples/finite_GO.json --output work/finite_GO.json --pdb work/finite_GO.pdb --json
uv run python check_model.py --model work/finite_GO.json --json
uv run python 03_prepare_amber.py --model work/finite_GO.json --output-dir work/amber --radii mbondi2
```

## Data and files

| Output under `work/quickstart/` | Contents |
|---|---|
| `finite_GO.json`, `finite_GO.pdb` | Sheet connectivity and coordinates |
| `amber/tleap.in`, `amber/GO_gaff2.frcmod` | LEaP input and parameter modifications |
| `amber/assigned_model.json`, `amber/library_audit.json` | Atom-type and library assignments |
| `amber/manifest.json`, `amber/normalized_libraries/` | Input manifest and unit libraries |

The [command reference](docs/REFERENCE.md#scripts--usage-workflow) lists the
build, functionalization, parameterization, docking and export scripts.

## Methods and validation

The [charge presets](docs/CHARGE_PRESETS.md) include the published GAFF2–AIM
parameters, the current libraries and the HierGO/OPLS reference.
[Parameterization studies](research/parameterization/README.md) cover RESP
development; experimental parameters are kept separate from the reference libraries.

[VALIDATION.md](VALIDATION.md) lists completed regression benchmarks and
parameterization work still in progress. [COMPARISON.md](docs/COMPARISON.md) compares
methodologies against HierGO and GOPY. Original reference files are indexed in
`resources/reference/SOURCES.json`; HierGO/OPLS tracks commit SHAs, file hashes,
and upstream licenses in `resources/reference/opls/PROVENANCE.json`.

Graph checks establish connectivity and coordination. They do not establish
electronic structure, parameter accuracy or equilibration.
See the [charge-derivation reference](docs/REFERENCE.md#charge-derivation-audits--validation)
for the QM and ESP workflow.

## Documentation

<a id="installation--environment-setup"></a>
<a id="scripts--usage-workflow"></a>
<a id="generating-a-graphene-oxide-sheet"></a>
<a id="construction-details"></a>
<a id="covalent-functionalization-amino-acids-peg-peptides-cyclodextrins"></a>
<a id="amber-topologies-periodicity-docking--mmgbsa"></a>
<a id="receptor-preparation-and-3d-bounding-boxes"></a>
<a id="autodock-vina--empirical-gnina-docking"></a>
<a id="complex-assembly-for-md--mmgbsa"></a>
<a id="charge-derivation-audits--validation"></a>
<a id="local-execution--slurm-integration"></a>
<a id="validation-status--original-baselines"></a>
<a id="citation--license"></a>

<details>
<summary>Topic index</summary>

- [Installation & Environment Setup](docs/REFERENCE.md#installation--environment-setup)
- [Scripts & Usage Workflow](docs/REFERENCE.md#scripts--usage-workflow)
- [Generating a Graphene Oxide Sheet](docs/REFERENCE.md#generating-a-graphene-oxide-sheet)
- [Construction Details](docs/REFERENCE.md#construction-details)
- [Covalent Functionalization (Amino Acids, PEG, Peptides, Cyclodextrins)](docs/REFERENCE.md#covalent-functionalization-amino-acids-peg-peptides-cyclodextrins)
- [Amber Topologies, Periodicity, Docking & MM/GBSA](docs/REFERENCE.md#amber-topologies-periodicity-docking--mmgbsa)
- [Receptor Preparation and 3D Bounding Boxes](docs/REFERENCE.md#receptor-preparation-and-3d-bounding-boxes)
- [AutoDock Vina & Empirical GNINA Docking](docs/REFERENCE.md#autodock-vina--empirical-gnina-docking)
- [Complex Assembly for MD & MM/GBSA](docs/REFERENCE.md#complex-assembly-for-md--mmgbsa)
- [Charge Derivation: Audits & Validation](docs/REFERENCE.md#charge-derivation-audits--validation)
- [Local Execution & Slurm Integration](docs/REFERENCE.md#local-execution--slurm-integration)
- [Validation Status & Original Baselines](docs/REFERENCE.md#validation-status--original-baselines)
- [Citation & License](#citation)

</details>

| Topic | Guide |
|---|---|
| First example | [Quickstart](docs/QUICKSTART.md) |
| Commands and complete workflows | [Command reference](docs/REFERENCE.md) |
| Covalent functionalization | [Fragments](docs/FRAGMENTS.md) |
| Charge models and presets | [Charge presets](docs/CHARGE_PRESETS.md) |
| External programs | [Engine setup](docs/ENGINE_SETUP.md) |
| Validation and limitations | [Validation](VALIDATION.md) |
| Development plans and timings | [Roadmap](docs/ROADMAP.md) · [Benchmarks](docs/PERFORMANCE.md) |

[Trajectory alignment](docs/ALIGNMENT.md).

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

This project is licensed under the [Apache 2.0 License](LICENSE).
