# Graphene and GO modeller comparison

This catalogue compares functionality. External tools were inspected through
manuals and pinned source, without running a common physical-accuracy or speed
benchmark. **ND** means not documented in the inspected material, not proof of
absence. **Research** identifies exercised prototypes outside the production CLI.

Sources: [GOPY](https://github.com/Iourarum/GOPY/tree/779c5ec45de09f81a00d8147b4889151e5151e12),
[HierGO](https://github.com/IFM-molecular-simulation-group/HierGO/tree/60382a1a1807ef0003ffcc3903392e913b1e5e9a),
[MakeGraphitics](https://github.com/velocirobbie/make-graphitics/tree/9d910ed525a8edf98f733551d616097288f9868f),
[GrapheneGUI](https://github.com/nicolas-loubet/GrapheneGUI/tree/3746ea6b7d2ca9330be46906244bfe51a33d2ba4).

The source inventory covers 23 files, 359 function definitions and 49 CLI
option definitions: all retrieved HierGO entry scripts, both GOPY variants,
and selected MakeGraphitics/GrapheneGUI implementation and manuals. These are
audit counts, not independent features. Legacy Python-2 source was inventoried
without executing or silently porting it.

The [downloadable source inventory](../research/benchmarks/2026-11-13/data/builder_source_inventory.json)
records per-file source URLs, commit IDs, hashes and detected entry points.

## Structure and chemistry

| Capability | GOPY | HierGO | MakeGraphitics | GrapheneGUI | Workbench (Blocco 1) |
|---|---|---|---|---|---|
| Pristine graphene | Yes | Yes | Yes | Yes | Yes |
| Shapes & physical sizing | x/y dimensions | Rectangular tiles | Rectangular / hexagonal flakes | Plate controls | Rectangular, hexagonal (`zigzag`/`armchair`), circular; `lx_nm`/`ly_nm`/`diameter_nm` or `nx`/`ny` |
| Periodic construction | ND | Yes | Yes | Box/export support; infinite-network behavior not audited | Explicit bond images (`[T,T,F]` / ribbons) |
| Infinite periodic MD topology | ND | OPLS tile workflow | Periodic LAMMPS workflow | ND | Planned in Blocco 5 (not yet enabled in Amber/GROMACS exporter) |
| Graphite / multiple layers | ND | In-plane tile stitching | Periodic graphite / structure combination | Multiple slabs / duplication | No dedicated stack constructor |
| Carbon nanotubes | ND | ND | ND in manual | Zigzag / armchair rolling | Not implemented |
| Vacancies / holes | Generation and cleanup | Percentage / size distributions | ND | Manual editing | HierGO-style distributions (`small`/`medium`/`large`/`mixed`/`random`), `percent_vacancies`/`count`, explicit circular/connected holes, 2D bare-C probe area |
| Stone–Wales defects | ND | Dedicated script | ND | ND | Not implemented |
| Mini grain boundaries | ND | Dedicated option | ND | ND | Not implemented |
| Rupture construction | ND | Dedicated script | ND | ND | Not implemented |
| N doping | Graphitic / pyridinic / pyrrolic | ND | ND | ND | Not implemented |
| Dangling-carbon cleanup | Hole cleanup | Trimming / removal scripts | No equivalent defect cleanup audited | Manual editing | Tip pruning + clash-screened H capping |
| Basal OH / epoxide | Yes | Yes | Yes | Yes | Yes (`HYD`, `EPX`) with `sp3` puckering & 12-azimuth screening |
| COOH / COO− & pH | Separate `GOPY_COO` variant | COOH yes; `COO−` not in audited charge table | Edge substitution | ND | Explicit `CBX`/`CBD` with formal charge, symmetric `COO−` (`1.26 Å`), Henderson–Hasselbalch `pH` / `deprotonation_fraction` |
| Other O chemistry | Main routine: OH/epoxide/COOH | Aged profiles include additional environments | Edge alcohol/carboxyl; basal OH/epoxide | OH/O documented | Edge epoxides (`EPH`, `EPT`) & carbonyls (`CAB`); lactones/phenols/quinones deferred |
| Reduced GO | Group removal routines | Oxygen / profile selection | Target oxidation ratio | Oxide removal | Direct low `O/C` (high `C/O`) recipes & `template-rgo`; no reduction trajectory |
| Pre-export geometric integrity gate | ND | Manual check workflow | ND | ND | Automatic $O(N)$ gate (`1-2`/`1-3`/`1-4+` clashes, motif-aware bond lengths, skew-cell minimum image, valence & charge check) |

## Patterns and construction controls

| Capability | GOPY | HierGO | MakeGraphitics | GrapheneGUI | Workbench (Blocco 1) |
|---|---|---|---|---|---|
| Random decoration | Yes | Yes | Reactivity-weighted growth | Random percentage | Seeded exact-count placement |
| Fresh / aged profiles | ND | Explicit profiles | Separate growth model | ND | No annealing profiles |
| Island / clustered oxidation | ND | `--bunching` (circular) | Island growth (not fixed disk) | Boolean spatial selection | Production CLI/config: `placement.mode = "island"` (`n_islands`, `centers`, `seed_atoms`, `radius`, `compactness`, `edge_mode`) |
| Kinetic / reactivity-driven islands | ND | Bunching is spatial | Random-forest reactivity / nucleation | ND | Not implemented |
| Deliberately dispersed pattern | ND | No decorator mode verified | ND | Expressions can impose regions | Production CLI/config: `placement.mode = "dispersed"` (max-min distance) |
| Individual group sites | ND | Region / edge selection | Reactivity controls | Interactive / Boolean selections | Exact IDs and epoxide bond pairs |
| Face control | Random up/down | Group-specific side rules | Site-dependent side rules | +Z / −Z / random | Both / up / down |
| Edge-first decoration | Edge groups treated separately | Explicit option | Edge substitution before basal oxidation | ND | Edge groups placed before basal; explicit `edge_mode` (`random` / `near_islands`) |
| Count guarantee | No comparable strict guarantee audited | Warnings for high loading / small regions | C/O target or site exhaustion | Percentage-based | Exact on success; explicit failure without silent truncation |
| Composition input | Group counts | Oxygen percentage / profiles | C/O target | Percentage / expressions | Explicit counts, `percent`/`fraction`, OR direct atomic `C/O` / `O/C` (scalar or interval) + `priority` + `pH` |
| Structural & QM diagnostics | ND | ND | ND | ND | Built-in `compute_chemistry_summary`: `sp2`/`sp3`, PBC-aware `sp2` domains, carbon chemical environments (`C-C sp2`, `C-O`, `C=O`, `O-C=O`), electron parity, sublattice counts |
| Reproducibility & config | Functions / GUI | CLI files/options; no seed flag audited | Scripts / global RNG | YAML / headless CLI | Direct CLI + `.toml` / `.json` configs (`--init-config`, `--write-config`), explicit `seed`, `--log`, provenance |

HierGO already provides circular bunching. Its decorator accepts a center and
diameter, filters eligible sites using periodic distances, and warns about
regions that are too small. The Workbench `island` mode is a controlled spatial
starting-structure option, not a kinetic oxidation model.
[Manual](https://github.com/IFM-molecular-simulation-group/HierGO/blob/60382a1a1807ef0003ffcc3903392e913b1e5e9a/README.md),
[source](https://github.com/IFM-molecular-simulation-group/HierGO/blob/60382a1a1807ef0003ffcc3903392e913b1e5e9a/Scripts/decorate-tile.py).

MakeGraphitics offers a different relevant strategy: local reactivity is
estimated with a random forest using published data, oxidation sites are
sampled by reactivity, and new islands can nucleate during growth. Reusing
this model would require a separate compatibility and transferability test;
it was not run here.
[Oxidiser](https://github.com/velocirobbie/make-graphitics/blob/9d910ed525a8edf98f733551d616097288f9868f/makegraphitics/reactors/oxidiser.py),
[paper](https://doi.org/10.1021/acs.jcim.9b00114).

## Export formats, topology and workflow

| Capability | GOPY | HierGO | MakeGraphitics | GrapheneGUI | Workbench (Blocco 1) |
|---|---|---|---|---|---|
| Coordinate / graph exports | PDB | PDB / GRO | LAMMPS data / PDB | PDB / MOL2 / GRO | JSON (authoritative graph), PDB (wwPDB `CONECT` + formal charges), MOL2 (contiguous residues, ParmEd-verified), XYZ, Extended XYZ (ASE-verified), SDF V3000 (Kekulé + `RAD=2`, RDKit `sanitize=True` verified), CIF (`P 1`, `>= 20 Å` non-periodic clearance) |
| PEG / covalent functionalization | Specific amino-PEG routine | ND | ND | ND | Atom-conserving anti-opening (`prepare_fragment.py` + `02_open_epoxides.py`); deferred until pure-GO FF is finalized |
| Connectivity | Coordinate-based identification | Coordinate/PBC neighbors | Connection graph / MD topology | Bond/topology export documented | Native explicit graph retained |
| Parameters | Amber tutorial / libraries | OPLS typing / topology workflow | OPLS / GraFF / ReaxFF documented | AMBER/GAFF-style export per manual | GAFF2 libraries (Blocco 2 `v3` + `L0_v2` in progress); multi-FF presets planned (Blocco 6) |
| Charge derivation / common ESP comparison | No common benchmark audited | Tabulated charges | Published parameter models | Charge export, not a new QM fit | Gaussian/RESP research; common-grid evaluation |
| Docking / complex assembly | ND | ND | ND | ND | Vina / bounded GNINA; topology-preserving assembly |
| Amber / GROMACS export | Amber tutorial / ParmEd | OPLS GROMACS workflow | ND in manual | MOL2 / GRO / TOP | Finite Amber (`tleap`) + GROMACS (via ParmEd); periodic MD export in Blocco 5 |
| LAMMPS export | ND | ND | Native data file | ND | Not implemented |
| Interface | Tkinter GUI / functions | CLI tools | Python library/scripts | GUI and headless CLI | Direct CLI (`01_build_GO.py` / `go.py build`) + `.toml`/`.json` configs; no modelling GUI |
| Evidence assessed here | Source/manual | Source/manual | Source/manual | Source/manual | Automated pytest suite (`368` tests), RDKit/ParmEd/ASE/Gemmi checks, and measured scaling below |

The GOPY PEG routine and Workbench anti-opening describe different chemical
operations. Neither a constructed geometry nor an exported file establishes
correct parameters for a new covalent junction.

## Measured Workbench scalability and geometric quality (Blocco 1)

All runs below used `build_sheet` + `check_geometry_integrity` + `compute_chemistry_summary`
on Apple Silicon (`CPython 3.13.15`, single thread, uninstrumented wall time; Python peak
memory measured in a separate `tracemalloc` pass, process RSS measured separately for `64×64`;
reproducible via `uv run python scripts/benchmark_build.py --blocco1`).
Finite cases use `O/C = 0.35`, `pH = 7.4`, `priority = ["OH", "epoxide", "COOH", "COO"]`;
periodic cases use `O/C = 0.30`, `priority = ["OH", "epoxide"]`.

| Benchmark case | Output atoms | Functional groups | Build (s) | Geometry gate (s) | Chemistry summary (s) | Python peak (MiB) | Realized O/C | Integrity gate |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| Finite `3×3` (`1.3×0.7 nm`) | 58 | 9 | 0.017 | 0.001 | 0.000 | 0.2 | 0.367 | PASS (`0` clashes) |
| Finite `6×6` (`2.6×1.5 nm`) | 231 | 40 | 0.059 | 0.011 | 0.001 | 0.9 | 0.348 | PASS (`0` clashes) |
| Hexagonal zigzag (`d = 4.0 nm`) | 638 | 119 | 0.112 | 0.020 | 0.002 | 3.2 | 0.350 | PASS (`0` clashes) |
| Finite `12×12` + 5% mixed holes + 2 islands | 899 | 161 | 0.246 | 0.027 | 0.002 | 4.1 | 0.348 | PASS (`0` clashes) |
| Finite `12×12` (`5.1×3.0 nm`) | 913 | 170 | 0.159 | 0.030 | 0.005 | 4.3 | 0.349 | PASS (`0` clashes) |
| Periodic `24×24` (`O/C = 0.30`) | 3,456 | 691 | 0.764 | 0.172 | 0.021 | 17.3 | 0.300 | PASS (`0` clashes) |
| Finite `24×24` (`10.2×5.9 nm`) | 3,666 | 729 | 0.671 | 0.116 | 0.025 | 17.5 | 0.350 | PASS (`0` clashes) |
| Finite `36×36` (`15.3×8.9 nm`) | 8,294 | 1,697 | 1.595 | 0.258 | 0.068 | 39.3 | 0.350 | PASS (`0` clashes) |
| Periodic `48×48` (`O/C = 0.30`) | 13,824 | 2,765 | 3.293 | 0.748 | 0.215 | 69.3 | 0.300 | PASS (`0` clashes) |
| Periodic `64×64` (`O/C = 0.30`) | 24,576 | 4,915 | 6.304 | 1.331 | 0.394 | ~203.4 RSS (process total) | 0.300 | PASS (`0` clashes) |

> [!NOTE]
> **Scope of this benchmark**: the table above measures Workbench alone (`O(N)` scaling across
> 420× in system size with zero clashes and zero recursion limits). A controlled head-to-head
> benchmark running external tools (`GOPY`, `HierGO`, `MakeGraphitics`, `GrapheneGUI`) on the
> same machine at matched sizes and oxidation levels (comparing speed, memory, initial geometry
> clashes/bond lengths, and topology/charge conservation) is scheduled for **Blocco 7**.

## Practical strengths and limitations

- **GOPY:** graphene/GO/rGO, N doping and a dedicated amino-PEG operation.
  Parametrization through its Amber tutorial is a separate step; raw motif charges
  do not conserve integer total charge on arbitrary oxidized flakes.
- **HierGO:** tile assembly, topological defects (Stone–Wales, grain boundaries,
  rupture), fresh/aged profiles and circular bunching; the documented workflow
  includes manual checks and omits `COO−` in the audited OPLS charge table.
- **MakeGraphitics:** graphite/LAMMPS systems and reactivity-driven oxidation.
  The inspected code uses Python 2 and reports substitutions for missing
  OPLS bonded environments.
- **GrapheneGUI:** interactive/Boolean selections, multiple slabs, CNT rolling
  and headless YAML configuration. File export alone does not validate parameters.
- **Workbench:** explicit connectivity, repeatable $O(N)$ construction up to
  `24,500+` atoms, clash-free constructed geometry (`constructed_unrelaxed`),
  motif-aware integrity gating, pH-aware `C/O` composition solving, multi-format
  QM/MD/docking exports (`PDB`/`MOL2`/`XYZ`/`EXTXYZ`/`SDF V3000`/`CIF`), and
  traceable numerical validation. It currently lacks N doping, Stone–Wales/grain
  boundaries, CNT rolling, a dedicated multi-layer stack constructor, kinetic
  oxidation, and periodic MD topology export (planned for Blocco 5).

The [13 November benchmark](../research/benchmarks/2026-11-13/README.md)
compares parameter combinations on common QM references. It does not rank
these external constructors by physical accuracy. Such a comparison needs
matched chemistry, composition, size, boundary conditions and ensembles.
