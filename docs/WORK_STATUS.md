> Historical baseline record. The current builder has 287 passing tests and local-only GO creation; current parameterization/preset status is documented in README, LAYERONE_PRESETS and research/parameterization.

> **Snapshot note:** this document preserves the work status recorded on 29 October 2026. The `validation/` folder, output archives, and binaries cited in the reports were excluded from the snapshot; the results summarized here therefore cannot be verified from the included files alone.

<a id="stato-di-lavoro--29-ottobre-2026"></a>
# Work status — 29 October 2026

The user's requested focus is GO creation and functionalization. The main
workflow is 01_build_GO.py → 02_open_epoxides.py (optional) →
03_prepare_amber.py, with prepare_fragment.py for reagents.

<a id="implementato-e-verificato"></a>
## Implemented and verified

- Explicit graph, finite/xy-periodic lattice, holes, vacancies, exact
  OH/epoxide/carbonyl/COOH/COO counts, H terminations, and reproducible seeds.
- Anti opening with a supported N/O/S fragment, transferred H, and preserved
  connectivity, formal charge, and identity; local normal and image handling.
- Simplified interface: --count with optional --id, --seed; random by
  default, with the fraction calculated from the output. Selection by fraction removed.
  The standard groups also accept {count,id}: carbon IDs reserved before
  random sites; lists of C-C pairs for epoxides.
- PlacementError for valid requests that cannot be fulfilled: WARNING, INCOMPLETE,
  known capacity, best count attempted and attempt count, suggestions; exit1,
  no partial model. Standard groups fall back to 128 arrangements; grafts try
  alternative sites/regioisomers/torsions. No physical maximum is inferred.
- Actual CLI runs: GO with standard groups, two PEG grafts on four epoxides, lysine
  through epsilon-N, and Gly-Gly through the N terminus. Data and logs in
  validation/final/generation_and_random_functionalization. New runs with
  count/id, standard sites, lysine warnings, capacity, and exhaustion of 128 attempts in
  validation/final/count_id_and_warnings; older records are historical.
- Independent review of the new selector: an unhashable ID input case
  fixed; no further defects in reproducibility, selectors, or accounting.
- 258 tests passed; Ruff clean (33 files). 14 entry points with help/Slurm; the Slurm worker
  run locally produces the same model as the serial command.
- Amber LEaP + Sander GB on the cluster: job 1374053 COMPLETED0:0, 179 atoms,
  232 bonds, q−1; graph/q/coordinates preserved. Numerical GROMACS round trip
  verified; no grompp or energy comparison between engines.
- Final CURRENT GAFF2/OLD/OPLS audit, provenance, and 51 original copies
  verified. Gaussian/RESP prepared, not run.

<a id="integrazioni-già-presenti"></a>
## Existing integrations

04 nonplanar snapshot and 3D box; 05 Vina/GNINA-Vina without CNN;
06 assembly of complete post-docking topologies; 07 GROMACS export.
These integrations follow the builder. Isolated external tests: GNINA 1.3.2
on actual GO + ethanol; MMPBSA 1 frame + idecomp1 + EO/DEO, job 1378877 COMPLETED0:0,
finite CSV values verified and downloaded. EO0/DEO1 indices require a frame map.
Assembly requires all H atoms; standard GNINA outputs may
require external reconstruction. No new development on these modules is
needed to use the generator.

<a id="manutenzione-allineamento-autorizzata-separatamente"></a>
## Separately authorized alignment maintenance

11533 HPC files (11520 campaign + 12 templates + 1 helper), 4 local templates updated
align→rms with backups. The 11520 previous inputs used align first without a mask.
CPPTRAJ 6.29.13 tested numerically: bare/masked align/rms all move
the entire frame, ligand error 1.57 × 10⁻⁶ Å. No reason to repeat GBSA solely because of
that command; no campaign trajectory or analysis recalculated.
HPC backups and helper downloaded to validation/final/alignment_migration.
README/SHA256SUMS from the old VS_POLL_WORKFLOW updated and verified (147 files);
the previous versions of the two documents are in the local backup. A valid autoimage
anchor is required before boxnobox; limitations for large/folded sheets are in ALIGNMENT.md.

<a id="limiti-scientifici-aperti"></a>
## Open scientific limitations

Parameterize new junctions; run a common QM/ESP benchmark before
choosing CURRENT/OLD/OPLS charges; relax/equilibrate the structures.
No preferred charge model selected. Winding covalent lattices
are preserved in the graph, but Amber export remains blocked; the tested workflow is finite.
Do not confuse format tests or a single GBSA snapshot with physical validation.

<a id="ambiente"></a>
## Environment

Isolated GO_BUILDER working folder; no Git at the time of the recorded status.
Python 3.14.10, uv lock; dependencies exported without editable installs to requirements.txt.
CLIs add src to the path; API access from this folder with PYTHONPATH=src.
Software copied byte-for-byte from VS_POLL_WORKFLOW: Vina/vina_split 1.2.7,
MGLTools 1.5.7p1 Linux x86_64 archive, wrappers, and licenses. Amber/Gaussian/GNINA
remain external. Do not extract MGLTools on a case-insensitive macOS filesystem.
User documentation: README.md, VALIDATION.md, docs/FRAGMENTS.md,
docs/COMPARISON.md, and docs/CHARGES.md.

<a id="confronto-richiesto-con-gopyhiergo"></a>
## Requested comparison with GOPY/HierGO

docs/COMPARISON.md updated with 7 functional components: 3 geometric building blocks
shared by both (lattice, defects/cleanup, decoration), periodicity related
to the HierGO approach, grafts/Amber-ParmEd related to the GOPY approach. Different algorithms and
parameters; no percentage of code attributed. The GOPY tutorial
uses AMBER99SB analogies, HierGO uses tabulated OPLS; this project uses 9 user GAFF2 libraries + frcmod.
Five HierGO originals archived for reference; only Oplsatom2 read by the audit,
no HierGO/GOPY program run by the generator. New junctions still require
parameterization, QM/RESP not run, no preferred charge model selected.
