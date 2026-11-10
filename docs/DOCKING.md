<a id="complessi-dopo-docking-amber-e-mmgbsa"></a>
# Post-docking complexes: Amber and MM/GBSA

The receptor is the entire parameterized surface: GO, oxygen-containing groups,
and any linker or covalent functionalization. The adsorbed pollutant is a
separate component. This separation does not depend on names such as GOX or
on predefined masks.

<a id="input-e-comando"></a>
## Inputs and command

Complete Amber topologies for the receptor and ligand are required, together
with pose coordinates in the same Cartesian reference frame. All formats require every atom in the topology, including H atoms.
The command preserves Amber charges, types, masses, and parameters; any
charges in the docking file are not used.

| Coordinates | Matching to the Amber template |
|---|---|
| Amber restart (`.rst7`, `.inpcrd`, `.ncrst`, `.rst`, `.restrt`, `.crd`) | Complete topology order, explicitly provided by the user |
| PDB/ENT | Complete order; names, explicit elements, and residue boundaries checked |
| PDBQT | First MODEL; unique atom names, a single residue; Vina tree order may differ |
| SDF/MOL | First pose; explicit graph and elements matching the template, with all H atoms present |

```bash
uv run python 06_assemble_docked_complex.py \
  --receptor-topology GO/go.prmtop --receptor-coordinates GO/go.inpcrd \
  --ligand-topology contaminante.parm7 --ligand-coordinates posa.rst7 \
  --output-dir complesso_posa1
```

Alternatively, add `--slurm` and the partition/resource options: the scientific
logic is identical to local execution.

A PDBQT is accepted only if it contains **all** atoms in the template.
Names must be unique within a single molecule/residue; AutoDock types
are used to check the element. ROOT/BRANCH order does not determine Amber order.
Ring-closure pseudoatoms and hydration dummy atoms are not supported.
[Official Vina type definitions](https://github.com/ccsb-scripps/AutoDock-Vina/blob/develop/src/lib/atom_constants.h).

SDF/MOL files are read with RDKit without adding H atoms. The complete
connectivity, including the H atoms bonded to each atom, must match. Bond orders
are not inferred from Amber: the comparison uses adjacency and elements,
while the electronic/force-field model remains that of the Amber template.
This is not an independent check of aromaticity, electronic states, or E/Z stereochemistry.

Automatic mapping enumerates only heavy-atom matches
(max 4096); H atoms are assigned by parent atom. A symmetry that
exchanges atoms with different RESP charges, types, or masses is rejected. H atoms
with the same parent must have equivalent atomic parameters. To
resolve a known correspondence, add the property
`GO_BUILDER_ATOM_NAMES` to the **SDF**: a JSON array of all Amber names, in SDF
atom order. All template names must be unique. Example:

```text
> <GO_BUILDER_ATOM_NAMES>
["H2", "C1", "H1", "H3", "H4"]

```

Elements and the complete graph are checked even with an explicit map.
No SDF index is automatically assumed to be an Amber index.

For SDF files with potentially stereogenic tetrahedral C/Si centers, complete
reference coordinates from the template, obtained before docking, are required.
The pose geometry is compared with the reference geometry without
alignment: inversions and nearly planar centers are rejected. Provide
`--ligand-reference-coordinates contaminante.rst7` for each ligand.
In the Python API, the option is `ligand_reference_coordinates=[...]`, with one
file or `None` per ligand. The reference must be an Amber restart
or a PDB in the complete template order.

Standard GNINA outputs do not guarantee the required H atoms: the source
documentation for `stripH` states that nonpolar H atoms are always removed.
No flag is introduced here to retain them. A GNINA SDF lacking H atoms
must be completed **externally**, preserving the template's identity,
protonation, and stereochemical geometry, before this assembly step.
[Official GNINA source](https://github.com/gnina/gnina/blob/master/gninasrc/main/main.cpp).
The same constraint applies to PDBQT files that have lost H atoms.

This command does not perform alignment, assign charges,
reconstruct the ligand from distances, or generate H atoms.

## Output

| File | Contents |
|---|---|
| `compl.parm7`, `compl.rst7`, `compl.pdb` | Complete dry complex |
| `rec.parm7`, `rec.rst7`, `rec.pdb` | Receptor, including covalent modifications |
| `lig.parm7`, `lig.rst7`, `lig.pdb` | Ligand, or the set of simultaneous ligands |
| `lig_001.*`, etc. | Individual ligand components |
| `masks.json` | Masks for the receptor, ligand(s), and complex |
| `manifest.json` | Sources and hashes, charge, parameter-conservation checks |

The extensions `.parm7` and `.prmtop` denote the same Amber format. PDB files
include CONECT; the topology remains the reference for all parameters.

Assembly preserves bonds, angles, torsions, exclusions, 1–4 interactions,
LJ parameters, charges, and GB radii/screening factors; it does not add receptor–ligand bonds.
Atom names and geometric proximity do not determine connectivity.

<a id="pose-alternative-o-caricamento-simultaneo"></a>
## Alternative poses or simultaneous loading

Each alternative pose belongs in a separate complex. Repeating the ligand
options instead inserts several molecules **simultaneously** into the same
system. Coordinates must be compatible; severe overlaps
are rejected.

A `lig.parm7` containing multiple molecules describes a combined ligand group. A
single binding energy for that group is not the individual affinity of
each copy. Those analyses require consistent partitions and trajectories.

<a id="compatibilità-con-la-dinamica-e-il-gbsa"></a>
## Compatibility with dynamics and GBSA

The output files are dry. Solvation, counterions, and any relaxation
are subsequent steps in the MD workflow; this command does not
start them. The covalent solute must be finite, and the coordinates must
contain whole molecules.

The components must use a consistent protocol for GB radii, 1–4 factors, and
LJ types with the same names. Incompatible data are rejected: a
force field is not adapted on the fly. The numerical compatibility check does not validate
surface charges, pose correctness, or the thermodynamic
interpretation of GBSA.

The GO–PEG/amino-acid/peptide/cyclodextrin bond remains internal to the receptor.
Cutting it to treat the modification as a noncovalent ligand would change
the physical problem and is not the standard analysis of an adsorbed pollutant.
