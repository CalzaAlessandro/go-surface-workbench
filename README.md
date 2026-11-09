# GO Builder (`go-surface-workbench`)

[![Status: Under Active Development](https://img.shields.io/badge/status-active%20development-orange.svg)]()
[![Python: >=3.11](https://img.shields.io/badge/python-3.11+-blue.svg)]()
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> [!IMPORTANT]
> **Work in Progress / Active Research Repository**
> This repository contains the source code, parameterization protocols, and validation toolchains for the GO Builder workbench.
> The project is currently under active scientific development. Force-field parameters, atomic charges, and model topologies are preliminary and subject to revision prior to peer-reviewed publication.
> 
> Molecular engine outputs (Gaussian/ORCA/Amber logs and trajectories) are tracked separately. External dependencies and molecular dynamics engines must be configured as described in [Engine setup](docs/ENGINE_SETUP.md).

Provisional repository name: **go-surface-workbench**.
See [scientific roadmap](docs/ROADMAP.md) and [measured performance benchmarks](docs/PERFORMANCE.md).

Studies on charge derivation methods and candidate force fields are documented in
[research/parameterization](research/parameterization/README.md); newly derived
parameters remain isolated from the libraries used in legacy simulation campaigns.

The [charge presets catalog](docs/CHARGE_PRESETS.md) includes the published
GAFF2–AIM reference with original charges and Lennard-Jones parameters, alongside
current libraries and the HierGO/OPLS reference. New charge derivation protocols
focus on RESP; general assignment across all presets is ongoing.

General-purpose graphene oxide generator with explicit connectivity. This
directory contains code, examples, test suites, and parameter copies with full
provenance; it does not distribute Amber or Gaussian binaries. Legacy GOX/REGOX
corrected models and OLD parameters serve as reference points and are neither rebuilt
nor overwritten automatically. Alignment inputs are maintained with backups as
described in [ALIGNMENT.md](docs/ALIGNMENT.md).

The core workflow to **build and functionalize GO sheets** is:
`01_build_GO.py` → `02_open_epoxides.py` (optional) → `03_prepare_amber.py`.
`prepare_fragment.py` prepares grafting reagents once; the resulting recipe can be
reused across multiple loadings and random seeds.

## Installation & Environment Setup

Requires Python 3.11 or later. From the project directory:

```bash
uv sync
uv run python go.py --help
uv run pytest -q
```

Dependencies (`pyproject.toml` and `uv.lock`) include NumPy, ParmEd, and RDKit. For
environments without `uv`, a standard `requirements.txt` can be exported and installed
in a virtual environment. Amber/AmberTools (`tleap`, `sander`, `MMPBSA.py`, `espgen`,
`resp`) and Gaussian (`g16`) must be provided by the host environment or cluster modules.

This repository contains source code, examples, tests, reference assets, documentation,
and dependency specifications. It does not track `software/`, Vina/Open Babel/MGLTools
binaries, archives, or large `validation/` outputs. Recreate `.venv` on target machines:
macOS binary virtual environments cannot be transferred to Linux. Install system engines
separately as detailed in [Engine setup](docs/ENGINE_SETUP.md). To import modules directly
from source, use `PYTHONPATH=src uv run python ...`; root scripts resolve package paths
automatically.

## Scripts & Usage Workflow

| Script | Purpose |
|---|---|
| `prepare_fragment.py` | Optional: prepares a fragment from SDF/MOL/MOL2/SMILES and explicitly sets the reactive attachment site |
| `01_build_GO.py` | Generates basal sheet, vacancies/holes, functional groups, and edge H-terminations |
| `02_open_epoxides.py` | Optional: anti-opening of basal epoxides with the selected nucleophilic fragment |
| `03_prepare_amber.py` | Assigns GAFF2 templates and prepares LEaP inputs; `--run` executes external AmberTools |
| `04_prepare_docking_receptor.py` | Extracts native or equilibrated GO and constructs 3D bounding boxes without flattening |
| `05_dock_boxes.py` | Bounded docking using AutoDock Vina or GNINA (with Vina scoring and CNN disabled) |
| `06_assemble_docked_complex.py` | Assembles the complete receptor and pre-parameterized ligands in docked pose coordinates |
| `07_export_gromacs.py` | Exports a completed Amber topology into GROMACS `.top`, `.itp`, and `.gro` |
| `audit_charges.py` | Descriptive audit and comparison of CURRENT, OLD, and HierGO/OPLS charge assignments |
| `prepare_QM.py` | Generates Gaussian/MK ESP single-point inputs on explicitly capped finite fragments |
| `compare_ESP.py` | Evaluates electrostatic potential RMSD across candidate charge sets on identical grids |
| `import_Amber.py` | Alternative entry point to 01: imports existing topologies without guessing bonds |
| `check_model.py` | Inspects molecular graphs, GBSA prerequisites, and selection masks |
| `update_alignment.py` | Audits and migrates CPPTRAJ alignment inputs with automatic backups |

All commands are accessible via the `go.py` unified CLI: `fragment`, `build`,
`functionalize`, `amber`, `gromacs`, `audit`, `qm`, `esp-import`, `esp`, `import`,
`check`, `snapshot`, `dock`, `complex`, `alignment`. Run `go.py <command> --help`
to list specific options.

## Generating a Graphene Oxide Sheet

```bash
uv run python 01_build_GO.py --recipe examples/finite_GO.json \
  --output work/GO.json --pdb work/GO.pdb
uv run python 03_prepare_amber.py --model work/GO.json \
  --output-dir work/amber --radii mbondi2 --run
uv run python 07_export_gromacs.py --topology work/amber/go.prmtop \
  --coordinates work/amber/go.inpcrd --output-prefix work/gromacs/GO
```

Without `--run`, `03_prepare_amber.py` generates input files but does not execute
AmberTools. With `--run`, it requires `tleap` on `$PATH` and verifies topology creation.
Existing outputs are never silently overwritten: specify a new target directory.

JSON recipe fields:

| Field | Description |
|---|---|
| `nx`, `ny` | Number of rectangular unit cells along x and y (each initially containing 4 carbons) |
| `bond_length` | Initial C–C bond length in Å (default: `1.42`) |
| `pbc` | `[false, false, false]` for finite flakes; `[true, true, false]` for periodic graphs |
| `vacuum` | Box z-dimension in Å; does not add solvent |
| `vacancies` | Count or explicit list of carbon IDs to remove |
| `holes` | List of `{center: [x, y], radius: r}` in Å or `{size: n, seed_atom: id}` |
| `groups` | Counts for `OH`, `epoxide`, `carbonyl`, `COOH`, `COO`; supports explicit `count` and `id` |
| `faces` | Face assignment for basal groups: `both`, `up`, or `down` |
| `seed` | Reproducible random seed |

Unit cell dimensions are $3 \times \text{bond\_length} \times nx$ and $\sqrt{3} \times \text{bond\_length} \times ny$.
Pruning dangling carbons updates the final atom count. Hydroxyls and epoxides require
basal aromatic sites; carbonyls and carboxyls require edge or vacancy boundaries.
The generator reports actual group counts; requests exceeding site capacity or
unresolvable packing collisions fail gracefully with exit code 1 without truncating counts.

Experimental XPS data can inform group ratios, but automatic spectral interpretation
is out of scope. Counts represent chemical groups: an epoxide spans two carbons.

### Construction Details

The lattice originates from analytical hexagonal cells of four carbons (initial C–C = 1.42 Å).
Bonds are recorded during generation. Vacancies and holes remove atoms from the graph;
dangling carbons are pruned to keep the sheet fully connected. Functional groups are
assigned hierarchically: OH and epoxides on basal planes, carbonyl/COOH/COO along
outer edges and pore boundaries. Each epoxide spans one C–C pair; other groups occupy
single carbon sites without overlapping.

**Placement is randomized for standard groups unless explicit IDs are supplied.**
Each attempt generates an entire candidate layout, validating COOH, COO, carbonyls,
epoxides, and hydroxyls in sequence. If steric constraints are violated, up to 128
stochastic layout attempts are evaluated. The `seed` parameter ensures exact
reproducibility; `faces: both` randomly distributes basal groups between top and bottom
faces without forcing an exact 50/50 split. Unfunctionalized edge carbons are saturated
with hydrogen atoms to yield an initial unrelaxed starting geometry.

To anchor groups at specific sites, replace group counts with an object:

```json
"groups": {
  "OH": {"count": 1, "id": [58]},
  "epoxide": 2,
  "COOH": 1
}
```

Here, one hydroxyl is pinned to carbon 58, while remaining groups are placed randomly.
For epoxides, `id` specifies pairs of adjacent carbons (e.g., `{"count": 1, "id": [[58, 59]]}`).
IDs correspond to pre-decoration JSON carbon indices. Pinned sites are reserved prior
to stochastic sampling and are not shifted upon steric clashes. See
`examples/GO_mixed_sites.json` for a complete example.

## Covalent Functionalization (Amino Acids, PEG, Peptides, Cyclodextrins)

Functionalization can be **stochastic and reproducible**, mirroring standard group decoration.
The reactive attachment atom on the fragment is explicitly designated, while surface
attachment sites can be selected automatically.

```bash
uv run python prepare_fragment.py --smiles 'NCC(=O)O' --inspect --json
uv run python prepare_fragment.py --smiles 'NCC(=O)O' \
  --attachment-atom 1 --seed 17 --output work/open_glycine.json
uv run python 02_open_epoxides.py --model work/GO.json \
  --recipe work/open_glycine.json --output work/GO_glycine.json \
  --pdb work/GO_glycine.pdb
```

For peptides, amino-PEG chains, or amino-functionalized cyclodextrins, supply
`--input molecule.sdf` and inspect the atom mapping with `--inspect --json`.
Atom indices are **1-based** and preserved prior to hydrogen addition. Reactive atoms
are never inferred purely from molecule names: in lysine, for instance, $\alpha$-N and
$\epsilon$-N are distinct attachment points. MOL2 inputs require explicit hydrogens.

The implemented ring-opening mechanism cleaves one epoxide C–O bond, forms a new
C–nucleophile bond, and transfers a donor hydrogen to the remaining oxygen atom.
Nucleophilic attack occurs from the face opposite the cleaved C–O bond; the newly
formed hydroxyl remains on the opposite side of the graft. Atoms, formal charges,
provenance, and fragment stereochemistry are tracked throughout.

Example: generating GO and grafting two copies of an amino-PEG fragment across half of
the four available epoxides:

```bash
uv run python 01_build_GO.py --recipe examples/GO_for_random_grafting.json \
  --output work/GO_for_PEG.json --pdb work/GO_for_PEG.pdb
uv run python 02_open_epoxides.py --model work/GO_for_PEG.json \
  --recipe examples/open_epoxide_PEG2_amine.json \
  --count 2 --seed 17 \
  --output work/GO_PEG.json --pdb work/GO_PEG.pdb
```

Grafting options:

| Option | Description |
|---|---|
| `--count 2` | Number of grafts; available epoxide sites are sampled randomly by default |
| `--id 2 3` | Optional: targets specific epoxide rings identified in `metadata.groups` |
| `--seed 17` | Reproducible random seed for site selection and conformer orientation |

The grafting fraction is **computed from output data** (grafts achieved divided by
intact epoxides present prior to execution). In this phase, `--id` identifies
**epoxide group IDs**, distinct from construction recipe carbon IDs. Consecutive runs
with different fragments consume remaining intact epoxides, enabling heterogeneous
surface functionalization.

Fragment search samples torsional rotations around the newly formed bond while
preserving internal fragment geometry. Graft orientation faces away from the target
epoxide oxygen. If requested loadings cannot be accommodated without steric clashes,
the command issues a warning and exits with code 1 without saving a partial model.
With `--json`, structured diagnostic metadata is returned.

**Newly formed junctions require parameterization.** `03_prepare_amber.py` rejects
modified models until valid `--parameter-templates` are supplied, including atom types,
charges, and dedicated `frcmod` files. Epoxide partial charges are not reused on
opened products. See [FRAGMENTS.md](docs/FRAGMENTS.md) and [CHARGES.md](docs/CHARGES.md).

Pertinent literature is discussed in [FUNCTIONALIZATION_LITERATURE.md](docs/FUNCTIONALIZATION_LITERATURE.md).
The anti-addition mechanism represents the design choice implemented in this workbench.

## Amber Topologies, Periodicity, Docking & MM/GBSA

Supported Amber execution centers on **finite covalent flakes**, solvated within
periodic boundary solvent boxes. Covalent bonds spanning periodic boundaries in
infinite lattice graphs are preserved in JSON models, but export to Amber/GROMACS is
gated pending engine-level verification. PDB coordinates do not preserve image shifts.
`examples/periodic_GO.json` demonstrates graph connectivity, not production Amber MD.

Receptors for docking and MM/GBSA comprise the GO sheet **plus all covalent functionalizations**.
Retain the Amber topology as the definitive reference when preparing PDBQT files:
atom ordering, partial charges, and bonding in docking inputs must remain consistent
with MD topologies.

### Receptor Preparation and 3D Bounding Boxes

```bash
uv run python 04_prepare_docking_receptor.py \
  --topology compl_solv.parm7 --coordinates md_10ns.rst \
  --receptor-mask ':1-494' --box-size 20 20 20 --overlap 5 --padding 5 \
  --output-dir work/GO_snapshot
```

For REGOX models, use `--receptor-mask ':1-628'`. Freshly generated topologies may
differ in numbering: inspect masks produced by script 03. NetCDF trajectories can be
sampled with `--frame N` (1-based). The receptor mask must encompass the entire covalent
sheet and any attached grafts.

Script 04 stitches finite flake images across boundaries following explicit bonds while
preserving out-of-plane corrugations. It outputs `receptor.parm7`, `receptor.rst7`,
`receptor.pdb`, molecular graphs, and `boxes.json`. Grid boxes span the true 3D envelope
rather than assuming idealized planar sheets.

### AutoDock Vina & Empirical GNINA Docking

GNINA allows direct input of receptor PDB and ligand SDF files:

```bash
uv run python 05_dock_boxes.py --boxes work/GO_snapshot/boxes.json \
  --receptor work/GO_snapshot/receptor.pdb --ligand ligand.sdf \
  --engine gnina --gnina /path/to/gnina --cpus 8 \
  --runs 3 --exhaustiveness 16 --num-modes 10 \
  --output-dir work/dock_gnina --run
```

The execution worker sets `--scoring vina --cnn_scoring none --no_gpu --pose_sort_order Energy`.
For AutoDock Vina, use `--engine vina` with pre-prepared PDBQT files and specify `--vina /path/to/vina`.
Resource allocation (`--cpus`) and search depth (`--exhaustiveness`) are passed explicitly.
Without `--run`, the script writes execution plans; `--slurm --run` submits jobs to the scheduler.
Logs, random seeds, input hashes, and top poses are recorded. Details in [DOCKING_BOXES.md](docs/DOCKING_BOXES.md).

Standard Vina empirical scoring does not utilize supplied partial charges; Amber partial
charges remain governed by the topology. Refer to the [Official AutoDock Vina FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html).

### Complex Assembly for MD & MM/GBSA

```bash
uv run python 06_assemble_docked_complex.py \
  --receptor-topology work/amber/go.prmtop \
  --receptor-coordinates work/amber/go.inpcrd \
  --ligand-topology ligand.parm7 --ligand-coordinates docked_pose.rst7 \
  --output-dir work/complex_pose1
```

Coordinates must be complete and defined in the docked reference frame; missing hydrogens
are not reconstructed from PDBQT. Multiple ligands can be specified by repeating flags;
overlapping candidate poses require separate complex assemblies. Assembly preserves
intramolecular bonds and does not introduce covalent bonds between receptor and ligand.
Masks are extracted directly from the resulting topology. Details in [DOCKING.md](docs/DOCKING.md).

Script 03 supports `--radii bondi|mbondi|mbondi2|mbondi3`; select GB radii consistently
across receptor, ligand, and complex. Technical GBSA execution does not inherently validate
GO parameter accuracy. Covalent functionalization bonds reside within the receptor mask
and must not be severed during binding affinity evaluations.

`import_Amber.py` imports existing topologies with explicit coordinates, preserving unit
cell parameters (`null` if absent) without inferring missing bonds or formal charges.
It identifies basal `EPX` epoxides matching expected local graph patterns. Amber topologies
do not store bond orders: numerical import does not guess chemical hybridization.

## Charge Derivation: Audits & Validation

```bash
uv run python audit_charges.py --output-dir work/charge_audit
uv run python 01_build_GO.py --recipe examples/QM_sheet.json --output work/QM.json
uv run python prepare_QM.py --model work/QM.json --output-dir work/QM_reference
```

`prepare_QM.py` generates single-point HF/6-31G* inputs with Merz-Kollman (MK) ESP grids
and exploratory RESP fitting scripts. Gaussian (`g16`) and RESP calculations do not run
automatically.

Following QM completion and convergence verification, the generated `prepare_resp.sh`
script invokes AmberTools. Import ESP grids and evaluate candidate charges:

```bash
uv run python go.py esp-import --esp-file reference.esp --model work/QM.json \
  --total-charge 0 --output work/reference.json
uv run python compare_ESP.py --reference work/reference.json \
  --candidates candidates.json --output-dir work/ESP_comparison
```

`candidates.json` maps model identifiers to identical geometries and atom orderings,
contrasting candidate partial charge assignments. Total charge and protonation states
must match exactly. Common fragment topologies are required before benchmarking
CURRENT, OLD, and OPLS sets; isolated ESP scores alone do not determine optimal
adsorption force fields. See [CHARGES.md](docs/CHARGES.md).

## Local Execution & Slurm Integration

Graphene sheet generation via `01_build_GO.py` or `go.py build` runs locally:
it accepts custom recipes or presets and does not expose Slurm flags. `--dry-run`
prints execution commands without writing models or submission scripts.

Downstream pipeline stages share unified worker scripts for local and Slurm execution.
For example, preparing Amber inputs on an HPC cluster:

```bash
uv run python 03_prepare_amber.py --model work/GO.json \
  --output-dir work/amber_slurm --slurm --cpus 1 --memory 2G \
  --time 00:10:00 --partition ciamician_cpu
```

Common flags across downstream stages: `--slurm`, `--dry-run`, `--submission-dir`,
`--cpus`, `--memory`, `--time`, `--partition`, `--account`, `--setup`, `--python`, `--json`.
`--setup /path/to/site_setup.sh` initializes cluster modules; `--python` specifies
the worker Python interpreter. Slurm modules and partitions are never chosen automatically.
`--slurm --dry-run` generates inspectable job scripts without submitting; submission
receipts track executed commands, script hashes, and scheduler job IDs.

## Validation Status & Original Baselines

[VALIDATION.md](VALIDATION.md) delineates completed regression benchmarks from
ongoing scientific parameterization. [COMPARISON.md](docs/COMPARISON.md) compares
methodologies against HierGO and GOPY. Original reference files are indexed in
`resources/reference/SOURCES.json`; HierGO/OPLS tracks commit SHAs, file hashes,
and upstream licenses in `resources/reference/opls/PROVENANCE.json`.

## Citation & License

This project is licensed under the [Apache 2.0 License](LICENSE).

If you use this workbench or its parameterization methodologies in your research,
please cite this repository using the metadata provided in [CITATION.cff](CITATION.cff).
