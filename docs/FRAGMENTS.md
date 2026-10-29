# Molecular fragments and explicit attachment sites

`go_builder.fragments` converts an explicit molecular graph into an epoxide-opening
recipe. Amino acids, amine-terminated PEG linkers, peptides, cyclodextrins and other
molecules use the same interface. There is no molecule-name lookup or automatic
choice of an amino acid nitrogen, peptide terminus, linker end or hydroxyl group.
Choose the attachment atom explicitly from the atom mapping.

The current reaction uses a **neutral amine, alcohol or thiol with an attached
hydrogen**. That hydrogen transfers to the epoxide oxygen; the donor gains the new
C–N, C–O or C–S bond. Charged remote groups are preserved. Protonated ammonium,
carboxylate oxygen, tertiary amine and amide nitrogen are unsupported attachment
sites. No pH adjustment or deprotonation is performed. The operation constructs
an initial anti geometry and does not establish that a reaction is experimentally
accessible.

## Command-line use

```sh
# Inspect IDs before choosing alpha versus epsilon lysine N.
python go.py fragment --smiles 'N[C@@H](CCCCN)C(=O)O' --inspect --json

# Save a neutral epsilon-N attachment recipe.
python go.py fragment --smiles 'N[C@@H](CCCCN)C(=O)O' \
  --attachment-atom 7 --count 1 --seed 42 --output lysine_epsilon_recipe.json

# Existing 3D molecule: retain original coordinates.
python prepare_fragment.py --input peptide.sdf --attachment-atom 1 \
  --output peptide_recipe.json

python 02_open_epoxides.py --model GO.json --recipe peptide_recipe.json \
  --output grafted.json --pdb grafted.pdb
```

`--transfer-hydrogen <ID>` overrides the deterministic attached-H choice.
`--generate-3d` explicitly permits replacing file coordinates.
`--no-optimize` skips isolated-fragment MMFF optimization for generated geometry.
`--allow-undefined-stereo` explicitly permits unspecified input stereochemistry.
An input peptide must provide all required stereocenters unless that option is
deliberately used.

## Inspect before choosing a site

```python
from go_builder.fragments import inspect_fragment, prepare_fragment

# L-lysine, with an explicitly encoded stereocenter and neutral amino groups.
smiles = "N[C@@H](CCCCN)C(=O)O"
listing = inspect_fragment(smiles=smiles)
for atom in listing["atom_mapping"]:
    print(atom["id"], atom["name"], atom["element"], atom["formal_charge"],
          atom["neighbor_ids"], atom["hydrogen_ids"], atom["cip_code"])

# In this exact SMILES atom order: alpha N = 1; side-chain epsilon N = 7.
recipe = prepare_fragment(smiles=smiles, attachment_atom=7, count=1, seed=42)
```

Atom IDs are one-based **input atom order**, including original explicit
hydrogens. Hydrogens added from encoded valence are appended after all original
atoms. A different SMILES ordering or a file rewritten by another program can
have different IDs. The mapping includes `original_atom_id`, `added_hydrogen`,
neighbors, attached H IDs, formal charge, isotope and stereochemical labels.
`inspect_fragment` does not generate coordinates and can inspect 2D files.

If `transfer_hydrogen` is omitted, preparation selects the lowest-ID attached
hydrogen and records `lowest_id_attached_hydrogen` in the fragment metadata.
Supply `transfer_hydrogen=<ID>` to select another attached hydrogen explicitly.
Names such as `N1`, `C2` and `H33` are unique within the fragment; large IDs use
base-36 suffixes to retain four-character names where possible.

## Input files and geometry

```python
# Preserve every original coordinate of a trustworthy 3D input.
recipe = prepare_fragment("peptide.sdf", attachment_atom=1, seed=42)

# Explicit permission to replace a file's coordinates with generated 3D.
recipe = prepare_fragment("linker.mol", attachment_atom=1,
                          generate_3d=True, seed=42, optimize=False)
```

Accepted files are single-molecule SDF, MOL and MOL2 with explicit connectivity
and bond orders. Multi-record SDF/MOL2 and disconnected salts are rejected;
select and save one connected molecule first. PDB is deliberately unsupported:
the importer never guesses bonds or bond orders from distances.

SDF, MOL and SMILES can omit hydrogens when their valence is chemically encoded.
Adding those H atoms preserves the encoded formal charges and protonation.
For MOL2, provide complete explicit hydrogens and compatible atom types:
RDKit can otherwise interpret missing H as radicals or unexpected formal charges.
Such unsupported chemistry is rejected rather than silently corrected. MOL2
typing dialects unsupported by RDKit should be converted to a validated SDF/MOL.
MOL2 partial-charge columns are not adopted as force-field charges.

File inputs require declared 3D coordinates unless `generate_3d=True`. Adding
missing H to existing 3D inputs preserves all original coordinates, including
original H coordinates. An existing geometry is never minimized implicitly.
All graph coordinates are in Å.

SMILES always generate 3D using seeded RDKit ETKDGv3. Generated geometries receive
an optional isolated-fragment MMFF optimization (`optimize=True` by default),
when parameters exist. Metadata records convergence, iteration limits or missing
MMFF parameters. `optimize=False` retains the ETKDG result. This step does not
minimize a GO sheet, evaluate the attachment reaction, or produce an equilibrated
conformation. A large peptide or cyclodextrin has no special size cap, but
embedding can fail or be expensive; supplying a validated 3D structure avoids
embedding. Reproducibility applies to a fixed seed, input and RDKit version.

## Stereochemistry and chemical limits

Encoded stereocenters, isotopes, formal charges and bond stereochemistry are
retained in the preparation graph and provenance metadata. Undefined potential
stereocenters or double-bond stereochemistry are rejected by default. Set
`allow_undefined_stereo=True` explicitly to accept an arbitrary initial
stereoisomer; the unresolved input centers remain recorded in metadata. The
program does not imply that an unspecified input has a known stereoisomer.

For an amino-substituted cyclodextrin reagent, select its explicitly represented
amine nitrogen. Preparing an unmodified cyclodextrin OH fragment demonstrates
graph/geometry handling only; it does not establish OH-mediated grafting under
the conditions of an experiment using an amino-substituted cyclodextrin.

Aromatic bonds are exported as integral Kekulé orders, with the original
aromatic flags retained in bond metadata. This permits supported remote
heteroaromatic groups such as a histidine side chain without using approximate
aromatic valence sums as a chemistry model.

The reaction validator currently accepts complete H/C/N/O/S/halogen fragments
with its explicitly supported valences: neutral C, neutral or +1 N, neutral or
−1 O, and neutral divalent S. Phosphates, sulfoxides, sulfones, metals, radicals
and other unsupported valence/charge states require additional reaction and
parameterization support. General molecule input does not expand the validated
reaction chemistry. Even a supported large fragment may fail the sheet's
collision-free placement search for a particular starting conformation.

## Count, random sites and optional explicit sites

```bash
# Random reproducible selection; the same options work for any supported fragment.
python 02_open_epoxides.py --model GO.json --recipe peptide_recipe.json \
  --count 2 --seed 42 --output GO_peptide.json --pdb GO_peptide.pdb
# Explicit group IDs from GO.json metadata.groups, not atom IDs.
python 02_open_epoxides.py --model GO.json --recipe peptide_recipe.json \
  --count 2 --id 3 7 --seed 42 --output GO_selected.json
```

`count` controls loading, with random selection by default. `--id` is optional;
if supplied with count, its list length must agree. With only CLI `--id`, the
count is inferred from that list rather than inherited from the fragment
recipe. CLI count without IDs selects randomly, overriding any recipe sites.
With neither option the recipe is used; missing count defaults to one.
In API/JSON use `count` and optional `id`, plus `seed` and `fragment`.

The fraction is calculated as achieved grafts divided by intact epoxide
**rings at the start of this operation**, not carbon atoms or all oxygen
groups. It is output only. Two grafts on three rings report 2/3; two on four
report 1/2. No input fraction or rounding of a requested fraction is needed.

Random ordering and rigid attachment-torsion search are reproducible. Steric
rejection can bias accepted sites; this is a geometric starting model, not
equilibrium sampling or a randomization of every internal PEG torsion.
Unavailable or blocked requested loadings raise an error without saving a
partial model. `metadata.functionalization_history` records each operation,
seed, selected ring IDs, initial ring count and requested/achieved loading.
Sequential operations can attach different fragments to remaining rings.

Unmet counts raise `PlacementError`: CLI status `INCOMPLETE`, a human-readable
`[WARNING]`, requested/achieved counts, available rings, orientation attempts
and suggestions. The process exits 1 without a partial model. Available-ring
failures are detected before placement; steric failures try other sites and
rigid torsions. A different starting conformer may give a different result;
the recorded tried loading is not a proven physical maximum.

GOPY's `add_NH_PEG_NH2` chooses oxygen groups randomly according to requested
fractions, removes selected groups and builds a particular PEG chain. Here
random loading applies to a supplied molecular fragment, and the epoxide
oxygen is retained as OH during anti opening; OH/COOH groups are not silently
replaced. This preserves the atom balance of the chosen reaction.
Source: [GOPY.py, add_NH_PEG_NH2](https://github.com/Iourarum/GOPY/blob/master/GOPY.py).

## Apply or save the recipe through the API

```python
import json
from go_builder.functionalize import open_epoxides

with open("lysine_epsilon_recipe.json", "x") as handle:
    json.dump(recipe, handle, indent=2, allow_nan=False)

# model is an existing validated GO model with available epoxide groups.
opened = open_epoxides(model, recipe)
```

The result contains the exact atom/bond list, selected donor/H IDs, atom mapping,
input provenance, original and explicit-H isomeric SMILES and geometry method.
`charge` and `atom_type` remain unset. No QM program, Amber calculation or
parameter assignment is invoked; the new junction requires parameterization
before a production force-field simulation.
