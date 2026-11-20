<a id="go-builder-go-surface-workbench"></a>
<h1 align="center">GO Workbench</h1>

Generate graphene oxide sheets with explicit connectivity and optional covalent functionalization.

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.24119849"><img src="https://img.shields.io/badge/Archive-10.5281%2Fzenodo.24119849-2457A7?style=flat" alt="Archive DOI: 10.5281/zenodo.24119849"></a>
  <a href="https://github.com/CalzaAlessandro/go-surface-workbench/actions/workflows/tests.yml"><img src="https://img.shields.io/github/actions/workflow/status/CalzaAlessandro/go-surface-workbench/tests.yml?branch=main&amp;style=flat&amp;label=Tests" alt="Software tests"></a>
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-%E2%89%A53.11-58616B?style=flat" alt="Python: 3.11 or later"></a>
  <img src="https://img.shields.io/badge/Status-Milestone%20block--1--builder-2E8B57?style=flat" alt="Status: Milestone block-1-builder">
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

> [!IMPORTANT]
> **Project status (`block-1-builder`)**
> - **Ready to use — Blocco 1 (Structure Builder):** finite & periodic GO construction (`rectangular`, `hexagonal`, `circular`), nm sizing, HierGO-style pore distributions, `pH`-aware `C/O` & `O/C` composition solving, `random`/`island`/`dispersed` patterns, Priority 1 & 2 oxygen groups (`OH`, `epoxide`, `COOH`, `COO−`, `carbonyl`, `epoxide_eph`, `epoxide_ept`), $O(N)$ pre-export motif-aware integrity gating, and multi-format exports (`JSON`, `PDB`, `MOL2`, `XYZ`, `EXTXYZ`, `SDF V3000`, `CIF P1`) via direct CLI or `.toml`/`.json` configs.
> - **Active research / preliminary — Blocco 2+:** new force-field parameters, charges, and topologies remain preliminary pending peer-reviewed validation. See the [validation status](VALIDATION.md) and [charge presets](docs/CHARGE_PRESETS.md) before using derived parameters in production simulations.

## Getting started

Requires Python 3.11 or later (`tomllib` stdlib support; bundled quickstart also tested on Python 3.12+). Dependencies (`pyproject.toml` and `uv.lock`) include NumPy, SciPy, ParmEd, and RDKit. From the repository root:

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

### Blocco 1 Builder: Direct CLI, `.toml`/`.json` Configs & Presets

`01_build_GO.py` (equivalent to `go.py build`) works directly from CLI flags without
a config file, from commented `.toml`/`.json` configs, or from packaged presets
(`LayerOne` + illustrative templates):

```bash
# 1) Direct CLI: 4x4 nm oxidized GO at pH 7.4, 2 islands, 5% mixed pores, multi-format export
uv run python 01_build_GO.py \
  --size 4 4 --oc 0.35 --ph 7.4 --place island --islands 2 --holes mixed:5 --seed 42 \
  --out work/go.json --pdb work/go.pdb --mol2 work/go.mol2 --sdf work/go.sdf --cif work/go.cif \
  --log work/go.log

# 2) Hexagonal zigzag nanoflake (2.5 nm diameter) for QM / ASE / xTB workflows
uv run python 01_build_GO.py \
  --shape hexagonal --edge zigzag --diam 2.5 --oc 0.25 --place dispersed --seed 17 \
  --out work/hex.json --extxyz work/hex.extxyz --sdf work/hex.sdf

# 3) Annotated TOML/JSON config template & shipped TOML examples
uv run python 01_build_GO.py --init-config work/my_go.toml
uv run python 01_build_GO.py --config examples/configs/go_island_porous.toml --out work/island.json

# 4) Packaged presets & strict post-build geometry check (exit 1 on failure)
uv run python 01_build_GO.py --list-presets
uv run python 01_build_GO.py --preset template-go-ph7 --out work/go_ph7.json --pdb work/go_ph7.pdb
uv run python check_model.py --model work/go.json --strict
```

Three ready-to-run commented TOML examples are shipped in `examples/configs/`:
[`go_standard_ph7.toml`](examples/configs/go_standard_ph7.toml),
[`go_island_porous.toml`](examples/configs/go_island_porous.toml), and
[`go_hex_qm.toml`](examples/configs/go_hex_qm.toml) (also compatible with
[`prepare_composition.py`](prepare_composition.py)).

### Export Formats & Pre-Export Integrity Gate

Every constructed model retains explicit connectivity (`geometry_status = "constructed_unrelaxed"`,
no forced minimization) and passes the motif-aware $O(N)$ gate (`1-2`/`1-3`/`1-4+` clashes,
motif bond lengths, skew-cell minimum image, valence & charge checks) before export:
`--out` (authoritative JSON graph), `--pdb` (wwPDB `CONECT` + formal charges),
`--mol2` (contiguous-residue Tripos MOL2, ParmEd-verified), `--xyz`, `--extxyz`
(`Lattice`, `pbc`, escaped metadata), `--sdf` (Molfile V3000 with Kekulé orders,
formal charges & `RAD=2` on unpaired $\pi$-sites, RDKit `sanitize=True` verified),
`--cif` (`P 1` with `>= 20 Å` non-periodic clearance), and `--log`.

### Scalability & Literature Comparison

Workbench scales linearly ($O(N)$) from `<100` to `24,500+` atoms on a single CPU core
with `0` clashes (`Gate: PASS`):
- `Finite 12×12` (`913` atoms, `O/C = 0.35`, `pH 7.4`): **`0.16 s`** build, `0.03 s` gate
- `Finite 24×24` (`3,666` atoms, `O/C = 0.35`, `pH 7.4`): **`0.67 s`** build, `0.12 s` gate
- `Finite 36×36` (`8,294` atoms, `O/C = 0.35`, `pH 7.4`): **`1.60 s`** build, `0.26 s` gate
- `Periodic 64×64` (`24,576` atoms, `4,915` groups, `O/C = 0.30`): **`6.30 s`** build, `1.33 s` gate

See [Feature comparison & scalability](docs/COMPARISON.md) for the updated matrix vs
`GOPY`, `HierGO`, `MakeGraphitics`, and `GrapheneGUI`, and [Benchmarks](docs/PERFORMANCE.md)
for the reproducible 10-case scaling table (`scripts/benchmark_build.py --blocco1`).

### Full Pipeline & Amber Quickstart Steps

Build and functionalize a GO sheet with:
`01_build_GO.py` → `02_open_epoxides.py` (optional) → `03_prepare_amber.py`.
Use `prepare_fragment.py` to prepare grafting reagents for reuse across loadings
and random seeds.

To run the minimal quickstart steps separately:

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

The [charge presets](docs/CHARGE_PRESETS.md) and [force-field presets & references](docs/REFERENCE_AND_PRESETS.md)
document the published GAFF2–AIM parameters, current libraries, HierGO/OPLS reference,
and the 2026-11-19 expanded 58-fragment audit (extending the historical
[13 Nov 2026 pilot](research/benchmarks/2026-11-13/README.md)).
[Parameterization studies](research/parameterization/README.md) cover RESP
development; experimental parameters are kept separate from the reference libraries.

[VALIDATION.md](VALIDATION.md) lists completed regression benchmarks and
parameterization work still in progress. [Feature comparison & scalability](docs/COMPARISON.md)
compares capabilities against HierGO, GOPY, MakeGraphitics, and GrapheneGUI. Original
reference files are indexed in `resources/reference/SOURCES.json`; HierGO/OPLS tracks
commit SHAs, file hashes, and upstream licenses in `resources/reference/opls/PROVENANCE.json`.

Graph checks establish connectivity, coordination, and constructed steric/bond-length
validity. They do not establish electronic structure, force-field accuracy, or equilibration.
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
| Feature comparison & scalability | [Comparison vs literature](docs/COMPARISON.md) |
| Covalent functionalization | [Fragments](docs/FRAGMENTS.md) |
| Charge models and presets | [Charge presets](docs/CHARGE_PRESETS.md) · [Force-field presets & references](docs/REFERENCE_AND_PRESETS.md) |
| 13 Nov 2026 pilot benchmark | [Pilot report](research/benchmarks/2026-11-13/README.md) · [Per-case tables](research/benchmarks/2026-11-13/PER_CASE_RESULTS.md) |
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
