# Docking boxes for graphene oxide snapshots

`go_builder.docking.dock_boxes` prepares reproducible AutoDock Vina or GNINA
runs against a selected receptor snapshot. The boxes JSON ties all box centers
to the exact `receptor.pdb` snapshot created by the snapshot workflow. This
supports both flat native structures and folded or undulating equilibrated
structures: centers and sizes come from the selected snapshot's analysis,
without assuming a global surface plane or fixed Z value.

The box input is JSON with `schema_version: 1`, a `receptor_pdb` path relative
to the JSON file, its `receptor_pdb_sha256`, and a nonempty `boxes` list. Each
box has a unique safe `id`, three-number `center`, three-positive-number
`size`, and a `config` path relative to the JSON file. A config contains exactly
the six Vina fields `center_x`, `center_y`, `center_z`, `size_x`, `size_y`, and
`size_z`. The function verifies that those values agree with the JSON and
writes its own configs in the docking output folder.

Before creating output files, the function checks the reference PDB hash and
matches every heavy atom in the supplied receptor to a heavy atom in that PDB
by element and coordinates within 0.02 Å. Atom order may differ and hydrogens
may be absent. Missing or changed heavy atoms are rejected. This is an identity
check only; it does not infer bonds, prepare receptor charges, or repair
coordinates.

Vina requires receptor and ligand PDBQT files. GNINA accepts the
snapshot-matching receptor PDB or PDBQT and either a one-record ligand SDF or
PDBQT. SDF ligand input is checked with RDKit for a valid molecular graph and
finite 3D coordinates. PDBQT input is checked for parseable atom records; this
worker does not infer or reconstruct its bond graph.

## Select the docking engine

The default remains Vina:

```python
manifest = dock_boxes(boxes, receptor_pdbqt, ligand_pdbqt, output_dir)
```

GNINA can use a snapshot PDB with an SDF ligand:

```python
manifest = dock_boxes(
    boxes,
    "snapshot/receptor.pdb",
    "ligand.sdf",
    "gnina_run",
    engine="gnina",
    gnina="gnina",
    cpus=8,
    seed=1001,
)
```

GNINA is invoked with Vina scoring, CNN scoring disabled, GPU disabled, and
energy pose sorting. It writes SDF poses, and the summary score is parsed from
the first output molecule's `minimizedAffinity` property. Vina writes PDBQT
poses and its `REMARK VINA RESULT` score. The manifest records which score
field was used. Both engines receive the requested CPU count and distinct
incrementing seeds, and runs execute sequentially. These options follow the
[GNINA command-line definitions](https://github.com/gnina/gnina/blob/master/gninasrc/main/main.cpp).

The supplied ligand must already have usable atom typing and hydrogens for the
selected engine. This function does not convert GAFF atom types to AutoDock
types, alter PDBQT charges, prepare a ligand, or replace an Amber topology or
its charges. Vina's standard scoring ignores user-supplied partial charges, so
a charge edit from one preparation workflow is not a general-purpose docking
correction. The engine's receptor graph is not asserted to match an Amber
topology; preserve the snapshot's Amber topology and coordinates for later MD
and analysis. See the [Vina FAQ on partial charges](https://autodock-vina.readthedocs.io/en/latest/faq.html).

## Prepare for review

```python
from go_builder.docking import dock_boxes

manifest = dock_boxes(
    "snapshot/boxes.json",
    "snapshot/receptor.pdbqt",
    "ligand.pdbqt",
    "docking_run",
    cpus=8,
    runs=3,
    seed=1001,
)
```

With `run=False` (the default), this writes `manifest.json`, `commands.json`,
`dock_boxes.sh`, and fresh configs under `boxes/`. The manifest says
`inputs_prepared` and `docking_completed: false`. Neither Vina nor GNINA needs
to be installed for this review step. The output directory must be new or
empty.

## Run docking

Pass `run=True` to execute. For each box, the function runs the requested
number of engine invocations sequentially. Every invocation receives the same
explicit CPU allocation and a distinct incrementing seed; boxes and runs do
not execute concurrently. It captures stdout and stderr for each run, checks
that a nonempty pose has the same heavy-atom count as the ligand input, parses
the engine-specific score, and writes `summary.csv` with the best-scoring pose
for each box plus a global best row. The manifest changes to
`docking_completed` only after all runs and summary generation succeed. An
execution or output error records `failed` and `docking_completed: false` in
the manifest.

Docking scores rank poses from the selected setup; they are not binding free
energies. The receptor and ligand inputs do not validate parameterization or
replace later coordinate reconstruction, Amber assembly, equilibration, or
free-energy analysis. This worker performs no MD, QM, GBSA, or implicit ligand
preparation.
