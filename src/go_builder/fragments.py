"""Explicit molecular inputs for the existing neutral-donor epoxide reaction.

RDKit supplies file parsing, valence completion and optional initial geometry.
No atom selection, protonation, partial-charge or junction parameters are inferred.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .build import _integer
from .functionalize import _validated_fragment
from .model import canonical_bond


def _read(source, smiles, allow_undefined_stereo):
    if (source is None) == (smiles is None):
        raise ValueError("Provide exactly one source file or smiles string")
    if type(allow_undefined_stereo) is not bool:
        raise ValueError("allow_undefined_stereo must be boolean")
    if smiles is not None:
        if not isinstance(smiles, str) or not smiles.strip():
            raise ValueError("SMILES must be a nonempty string")
        parameters = Chem.SmilesParserParams()
        parameters.removeHs = False
        mol = Chem.MolFromSmiles(smiles, parameters)
        provenance = {"kind": "smiles", "value": smiles}
    else:
        path = Path(source)
        suffix = path.suffix.lower()
        if suffix not in {".mol", ".sdf", ".mol2"}:
            raise ValueError("Supported source formats are SDF, MOL and MOL2 with explicit bonds")
        if not path.is_file():
            raise ValueError(f"Fragment source does not exist: {path}")
        if suffix == ".sdf":
            records = list(Chem.SDMolSupplier(str(path), removeHs=False, strictParsing=True))
            if len(records) != 1:
                raise ValueError(
                    "SDF must contain exactly one molecule; split multiple records first"
                )
            mol = records[0]
        elif suffix == ".mol2":
            if path.read_text().count("@<TRIPOS>MOLECULE") != 1:
                raise ValueError("MOL2 must contain exactly one molecule")
            mol = Chem.MolFromMol2File(str(path), removeHs=False)
        else:
            mol = Chem.MolFromMolFile(str(path), removeHs=False, strictParsing=True)
        provenance = {"kind": suffix[1:], "path": str(path.resolve())}
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("Could not parse a nonempty chemically valid fragment")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Fragment must be one connected molecule; salts need explicit selection")
    if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        raise ValueError("Radical fragments are unsupported")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    undefined = [
        {
            "type": str(item.type),
            "centered_on": int(item.centeredOn) + 1,
            "index_kind": "bond" if "Bond" in str(item.type) else "atom",
        }
        for item in Chem.FindPotentialStereo(mol)
        if str(item.specified) == "Unspecified"
    ]
    if undefined and not allow_undefined_stereo:
        raise ValueError(
            "Undefined stereochemistry: provide explicit stereoisomers or set "
            "allow_undefined_stereo=True to accept an arbitrary initial stereoisomer"
        )
    provenance["rdkit_version"] = Chem.rdBase.rdkitVersion
    return mol, provenance, undefined


def _base36(value):
    text = ""
    while value:
        value, digit = divmod(value, 36)
        text = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[digit] + text
    return text or "0"


def _mapping(mol, original_count):
    names = set()
    mapping = []
    for atom in mol.GetAtoms():
        identifier = atom.GetIdx() + 1
        element = atom.GetSymbol()
        name = f"{element}{identifier}"
        if len(name) > 4:
            serial = identifier
            name = element + _base36(serial)
            while name in names:
                serial += 1
                name = element + _base36(serial)
        names.add(name)
        mapping.append(
            {
                "id": identifier,
                "original_atom_id": identifier if identifier <= original_count else None,
                "added_hydrogen": identifier > original_count,
                "element": element,
                "name": name,
                "formal_charge": atom.GetFormalCharge(),
                "isotope": atom.GetIsotope(),
                "chiral_tag": str(atom.GetChiralTag()),
                "cip_code": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
                "neighbor_ids": sorted(other.GetIdx() + 1 for other in atom.GetNeighbors()),
                "hydrogen_ids": sorted(
                    other.GetIdx() + 1 for other in atom.GetNeighbors() if other.GetAtomicNum() == 1
                ),
            }
        )
    return mapping


def inspect_fragment(source=None, *, smiles=None, allow_undefined_stereo=False):
    """Return stable one-based atom/site mapping without requiring or generating 3D.

    Original atoms retain input order, including any explicit H. Missing H atoms
    are appended from encoded valence; formal charges and protonation stay fixed.
    """
    mol, provenance, undefined = _read(source, smiles, allow_undefined_stereo)
    original_count = mol.GetNumAtoms()
    mol = Chem.AddHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return {
        "source": provenance,
        "atom_mapping": _mapping(mol, original_count),
        "original_atom_count": original_count,
        "atom_count": mol.GetNumAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "undefined_stereochemistry": undefined,
        "units": "angstrom",
    }


def prepare_fragment(
    source=None,
    *,
    smiles=None,
    attachment_atom: int,
    transfer_hydrogen: int | None = None,
    count=1,
    seed=0,
    generate_3d=False,
    allow_undefined_stereo=False,
    optimize=True,
):
    """Return an opening recipe with explicit coordinates, bonds and selected site.

    SMILES always require initial 3D generation. For files, existing 3D original
    coordinates remain exact unless generate_3d=True. Optional MMFF optimization
    concerns the isolated generated fragment only, never a GO reaction product.
    """
    _integer(attachment_atom, "attachment_atom", minimum=1)
    _integer(count, "count")
    _integer(seed, "seed")
    if transfer_hydrogen is not None:
        _integer(transfer_hydrogen, "transfer_hydrogen", minimum=1)
    if type(generate_3d) is not bool or type(optimize) is not bool:
        raise ValueError("generate_3d and optimize must be boolean")
    mol, provenance, undefined = _read(source, smiles, allow_undefined_stereo)
    original_count = mol.GetNumAtoms()
    original_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    generated = smiles is not None or generate_3d
    if not generated and (not mol.GetNumConformers() or not mol.GetConformer().Is3D()):
        raise ValueError("Input has no declared 3D geometry; set generate_3d=True explicitly")
    mol = Chem.AddHs(mol, addCoords=not generated)
    optimization = "not_requested" if not optimize else "not_applied_to_preserved_coordinates"
    if generated:
        mol.RemoveAllConformers()
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = seed % 2147483647
        if AllChem.EmbedMolecule(mol, parameters) != 0:
            raise ValueError("RDKit ETKDG failed to generate fragment 3D geometry")
        if optimize:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
                optimization = "converged" if status == 0 else "iteration_limit"
            else:
                optimization = "unavailable_parameters_etkdg_only"
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mapping = _mapping(mol, original_count)
    if attachment_atom > mol.GetNumAtoms():
        raise ValueError("attachment_atom ID does not exist; inspect the atom mapping")
    donor = mol.GetAtomWithIdx(attachment_atom - 1)
    hydrogen_ids = mapping[attachment_atom - 1]["hydrogen_ids"]
    if not hydrogen_ids:
        raise ValueError(
            "Selected donor has no explicit attached hydrogen after valence completion"
        )
    selected_automatically = transfer_hydrogen is None
    if selected_automatically:
        transfer_hydrogen = hydrogen_ids[0]
    if transfer_hydrogen not in hydrogen_ids:
        raise ValueError("transfer_hydrogen must identify H directly bonded to attachment_atom")
    # A peptide amide nitrogen is not a neutral amine reaction site.
    if donor.GetAtomicNum() == 7 and any(
        neighbor.GetAtomicNum() == 6
        and any(
            bond.GetBondTypeAsDouble() == 2
            and bond.GetOtherAtom(neighbor).GetAtomicNum() in (8, 16)
            for bond in neighbor.GetBonds()
        )
        for neighbor in donor.GetNeighbors()
    ):
        raise ValueError("Amide nitrogen is unsupported; select a neutral amine donor")
    aromatic_bonds = {bond.GetIdx(): bond.GetIsAromatic() for bond in mol.GetBonds()}
    # Integral Kekule orders allow chemically valid heteroaromatic remote groups
    # through the reaction validator while retaining aromatic provenance.
    export = Chem.Mol(mol)
    Chem.Kekulize(export, clearAromaticFlags=True)
    positions = mol.GetConformer().GetPositions()
    atoms = [
        {
            "id": row["id"],
            "element": row["element"],
            "name": row["name"],
            "resname": "NUC",
            "resid": 1,
            "position": positions[row["id"] - 1].tolist(),
            "formal_charge": row["formal_charge"],
            "charge": None,
            "atom_type": None,
            "metadata": {
                key: row[key]
                for key in (
                    "original_atom_id",
                    "added_hydrogen",
                    "isotope",
                    "chiral_tag",
                    "cip_code",
                )
            },
        }
        for row in mapping
    ]
    bonds = []
    for bond in export.GetBonds():
        entry = canonical_bond(
            bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1, bond.GetBondTypeAsDouble()
        )
        entry["metadata"] = {
            "aromatic_in_source": aromatic_bonds[bond.GetIdx()],
            "stereo": str(bond.GetStereo()),
            "stereo_atom_ids": [identifier + 1 for identifier in bond.GetStereoAtoms()],
        }
        bonds.append(entry)
    fragment = {
        "atoms": atoms,
        "bonds": bonds,
        "attachment_atom": attachment_atom,
        "transfer_hydrogen": transfer_hydrogen,
        "metadata": {
            "source": provenance,
            "atom_mapping": mapping,
            "original_atom_count": original_count,
            "original_isomeric_smiles": original_smiles,
            "explicit_h_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "undefined_stereochemistry": undefined,
            "transfer_hydrogen_selection": (
                "lowest_id_attached_hydrogen" if selected_automatically else "explicit_id"
            ),
            "units": "angstrom",
            "geometry": {
                "method": "ETKDGv3" if generated else "input_3d_with_missing_h_added",
                "original_coordinates_preserved": not generated,
                "mmff_optimization": optimization,
                "seed": seed,
                "status": "initial_fragment_geometry",
            },
            "parameterization_status": "unparameterized_junction",
        },
    }
    _validated_fragment(fragment)
    return {"count": count, "seed": seed, "fragment": fragment}
