# External engine setup

This source snapshot does not bundle molecular engine executables or software
archives. Install each engine separately under its own license and use a
version supported by your environment. The Python package dependencies do not
install these programs.

- **Vina:** install AutoDock Vina and pass its executable on `PATH` or with
  `--vina /path/to/vina` to `05_dock_boxes.py`. The docking code also expects
  `vina_split` when splitting output poses.
- **GNINA (optional):** install GNINA separately and pass its path with
  `--gnina /path/to/gnina`; the default workflow uses Vina scoring with CNN
  scoring disabled.
- **PDBQT preparation:** the Vina path requires prepared receptor and ligand
  PDBQT inputs. Open Babel and MGLTools/AutoDockTools are not included; install
  suitable tools separately if needed for preparation. The GNINA path can use
  receptor PDB and ligand SDF inputs as described in the main README.
- **Amber/AmberTools:** `tleap` is needed for topology preparation; Sander and
  `MMPBSA.py` are needed for the corresponding simulation/analysis actions.
  Gaussian (`g16`) plus AmberTools RESP utilities are needed only for the
  optional external QM/RESP workflow. These programs are not installed by
  `uv sync`.
- **HPC use:** Slurm, MPI, GPU drivers/runtime, and compatible compute-node
  installations are site-managed dependencies when those execution paths are
  used.

Confirm executable versions and local license/redistribution terms with their
respective projects or administrators before use. The snapshot records no tested
installation recipe for these external engines.
