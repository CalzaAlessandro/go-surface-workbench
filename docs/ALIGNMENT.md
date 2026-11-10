<a id="immagini-periodiche-e-allineamento-delle-traiettorie"></a>
# Periodic imaging and trajectory alignment

For finite GO in the classic workflow:

```text
parm compl_solv.parm7
trajin md_10ns.nc
strip :WAT
strip :Cl-
strip :Na+
autoimage anchor :1-494
box nobox
rms first :1-494
trajout md_10ns_imaged.nc
run
```

For REGOX, replace the mask with `:1-628`. These counts apply to the existing
corrected systems; new models require their own masks. The strip commands
shown do not renumber GO because water and ions follow it in the systems
examined. If the residue order changes, check the mask after stripping.

`autoimage :1-494` and `autoimage anchor :1-494` are equivalent syntax.
The anchor centers GO; other molecules are translated between periodic
images. The default behavior treats ligands as fixed molecules,
choosing the image closest to the anchor when bringing them closer. The command
requires the cell and therefore precedes `box nobox`. `rms first :1-494` uses GO
for the fit and applies the rotation/translation to the entire system; do not use
`nofit` for this purpose.
Sources: [autoimage](https://amberhub.chpc.utah.edu/autoimage/),
[RMSD](https://amberhub.chpc.utah.edu/rmsd/).

Autoimage does not flatten the sheet, nor does it by itself guarantee the reconstruction
of a covalent sheet already split across images. For very large or
folded GO, check bond continuity and ligand images: the criterion
relative to the anchor center does not always correspond to the minimum distance from
the surface. It does not automatically verify the contacts in every frame. Retain
the trajectory and topology with the cell for periodic analyses; the file without a box
cannot subsequently recover the periodic information that was removed.
The new extractor `04_prepare_docking_receptor.py` uses the covalent graph to
reassemble the receptor before preparing the boxes.

<a id="verifica-del-vecchio-align"></a>
## Verification of the previous align behavior

CPPTRAJ 6.29.13 on the cluster was tested with two frames related by a
known rotation and translation, using GOX + ligand, 1029 atoms. `align first`,
`align first :1-494` and `rms first :1-494` all moved the entire frame.
The residual error on the ligand is approximately 1.57 × 10⁻⁶ Å; the last two outputs are
numerically identical. The log for the command without a mask reports `*`,
1029 atoms. Evidence: `validation/final/alignment_behavior/verification.json`,
Slurm job 1378874 COMPLETED 0:0.

The [align page](https://amberhub.chpc.utah.edu/align/) describes a default
`move` behavior limited to the selection, which differs from the result for this
installed version. Incorrect GO–ligand geometry therefore cannot be attributed
solely to the use of the word align. In particular, the previous campaign inputs
used `align first` without a mask: GBSA does not need to be recalculated
solely because of that command. The GO-specific fit is now explicit for
future calculations; the resulting RMSD may differ from a previous fit over all
atoms without any change in internal distances.

<a id="manutenzione-applicata-il-29-ottobre-2026"></a>
## Maintenance applied on 29 October 2026

- Local: four active templates `align_{1,10,50,100}ns.in` in
  `VS_POLL_WORKFLOW/resources/LIB`, now using `rms first MASK`.
- HPC: 11520 campaign inputs (5760 GOX and 5760 REGOX), 12 runtime templates
  and one helper for free-ligand MD; 11533 files in total. The helper uses
  `rms first *`, because it contains no GO.
- All 11520 campaign inputs examined previously used the fit
  without a mask; no masked-only cases were found in this inventory.
- No campaign trajectories or analyses were recalculated.

Local backups: `validation/alignment_local_20261029`. Downloaded HPC backups and
manifests: `validation/final/alignment_migration`; the original archive is
`validation/final/alignment_backup_20261029.tar.gz`. The helper backup,
added after the archive, was downloaded separately into the migration folder.
Remote source: HPC validation area; the absolute path is omitted from the snapshot.

`update_alignment.py` is a maintenance utility: without `--apply` it produces
the inventory; with `--apply` it requires a new backup folder and preserves
the originals before making changes. It is not a reanalysis launcher.
