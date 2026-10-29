"""Fragment preparation preserves explicit chemical identity and site choices."""

import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from go_builder import fragments
from go_builder.functionalize import _validated_fragment


def test_smiles_glycine_has_explicit_transfer_h_and_deterministic_geometry():
    recipe = fragments.prepare_fragment(smiles="NCC(=O)O", attachment_atom=1, seed=37)
    fragment = recipe["fragment"]
    assert len(fragment["atoms"]) == 10
    assert fragment["attachment_atom"] == 1
    assert fragment["transfer_hydrogen"] == 6
    assert recipe == fragments.prepare_fragment(smiles="NCC(=O)O", attachment_atom=1, seed=37)
    assert fragment["metadata"]["transfer_hydrogen_selection"] == "lowest_id_attached_hydrogen"
    assert all(a["charge"] is None and a["atom_type"] is None for a in fragment["atoms"])
    assert len({a["name"] for a in fragment["atoms"]}) == 10
    json.dumps(recipe, allow_nan=False)
    _validated_fragment(fragment)


def test_site_mapping_distinguishes_alpha_and_lysine_epsilon_nitrogen():
    smiles = "N[C@@H](CCCCN)C(=O)O"
    listing = fragments.inspect_fragment(smiles=smiles)
    assert [(a["id"], a["element"]) for a in listing["atom_mapping"] if a["element"] == "N"] == [
        (1, "N"),
        (7, "N"),
    ]
    alpha = fragments.prepare_fragment(smiles=smiles, attachment_atom=1)
    epsilon = fragments.prepare_fragment(smiles=smiles, attachment_atom=7)
    assert alpha["fragment"]["attachment_atom"] == 1
    assert epsilon["fragment"]["attachment_atom"] == 7
    assert epsilon["fragment"]["metadata"]["atom_mapping"][1]["cip_code"] == "S"


def test_3d_input_original_atoms_preserved_when_adding_hydrogens(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("NCCO"))
    AllChem.EmbedMolecule(mol, randomSeed=12)
    mol = Chem.RemoveHs(mol)
    path = tmp_path / "fragment.mol"
    Chem.MolToMolFile(mol, str(path))
    expected = Chem.MolFromMolFile(str(path), removeHs=False).GetConformer().GetPositions()
    recipe = fragments.prepare_fragment(path, attachment_atom=1)
    observed = np.array([a["position"] for a in recipe["fragment"]["atoms"][:4]])
    np.testing.assert_array_equal(observed, expected)
    assert recipe["fragment"]["metadata"]["geometry"]["original_coordinates_preserved"]


def test_two_dimensional_file_requires_explicit_generation(tmp_path):
    mol = Chem.MolFromSmiles("NCCO")
    AllChem.Compute2DCoords(mol)
    path = tmp_path / "fragment.sdf"
    Chem.MolToMolFile(mol, str(path))
    with pytest.raises(ValueError, match="3D|3d"):
        fragments.prepare_fragment(path, attachment_atom=1)
    _validated_fragment(
        fragments.prepare_fragment(path, attachment_atom=1, generate_3d=True)["fragment"]
    )


def test_multi_record_sdf_is_rejected(tmp_path):
    path = tmp_path / "ambiguous.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(Chem.MolFromSmiles("CN"))
        writer.write(Chem.MolFromSmiles("CO"))
    with pytest.raises(ValueError, match="single|one|multiple"):
        fragments.prepare_fragment(path, attachment_atom=2, generate_3d=True)


def test_undefined_stereocenters_require_opt_in():
    with pytest.raises(ValueError, match="stereo"):
        fragments.prepare_fragment(smiles="NC(C)C(=O)O", attachment_atom=1)
    recipe = fragments.prepare_fragment(
        smiles="NC(C)C(=O)O", attachment_atom=1, allow_undefined_stereo=True
    )
    assert recipe["fragment"]["metadata"]["undefined_stereochemistry"]


def test_remote_charge_and_encoded_protonation_are_preserved():
    recipe = fragments.prepare_fragment(smiles="NCC(=O)[O-]", attachment_atom=1)
    atoms = recipe["fragment"]["atoms"]
    assert sum(a["formal_charge"] for a in atoms) == -1
    assert atoms[4]["formal_charge"] == -1
    assert len(atoms) == 9
    with pytest.raises(ValueError, match="neutral|donor|Supported"):
        fragments.prepare_fragment(smiles="[NH3+]CC(=O)[O-]", attachment_atom=1)


@pytest.mark.parametrize(
    "smiles,donor",
    [
        ("NCCOCCOCCO", 1),
        ("NCC(=O)N[C@@H](C)C(=O)O", 1),
        ("OCC(O)CO", 1),
        ("SCCO", 1),
        ("NCCc1c[nH]cn1", 1),
    ],
)
def test_generic_peg_peptide_polyol_thiol_and_heteroaromatic_graphs(smiles, donor):
    _validated_fragment(
        fragments.prepare_fragment(smiles=smiles, attachment_atom=donor)["fragment"]
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"attachment_atom": 0},
        {"attachment_atom": True},
        {"attachment_atom": 99},
        {"attachment_atom": 2},
        {"attachment_atom": 1, "transfer_hydrogen": 2},
        {"attachment_atom": 1, "count": -1},
        {"attachment_atom": 1, "seed": True},
    ],
)
def test_invalid_site_hydrogen_and_integer_options_reject(kwargs):
    with pytest.raises(ValueError):
        fragments.prepare_fragment(smiles="NCCO", **kwargs)


def test_file_formats_do_not_infer_bonds_from_pdb(tmp_path):
    path = tmp_path / "molecule.pdb"
    path.write_text("END\n")
    with pytest.raises(ValueError, match="SDF|MOL|format"):
        fragments.prepare_fragment(path, attachment_atom=1)


def test_explicit_hydrogen_ids_remain_in_original_input_order():
    recipe = fragments.prepare_fragment(
        smiles="[H]N([H])CCO", attachment_atom=2, transfer_hydrogen=3
    )
    fragment = recipe["fragment"]
    assert [a["element"] for a in fragment["atoms"][:4]] == ["H", "N", "H", "C"]
    assert fragment["transfer_hydrogen"] == 3
    assert fragment["metadata"]["atom_mapping"][2]["original_atom_id"] == 3


def test_missing_or_conflicting_sources_reject():
    with pytest.raises(ValueError):
        fragments.prepare_fragment(attachment_atom=1)
    with pytest.raises(ValueError):
        fragments.prepare_fragment("x.sdf", smiles="CN", attachment_atom=2)


def test_mol2_explicit_connectivity_and_original_atom_order(tmp_path):
    path = tmp_path / "ethylamine.mol2"
    path.write_text("""@<TRIPOS>MOLECULE
ethylamine
10 9 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 CA 0.0000 0.0000 0.0000 C.3 1 NUC 0.0
2 CB 1.5000 0.0000 0.0000 C.3 1 NUC 0.0
3 NZ 2.0000 1.3000 0.2000 N.3 1 NUC 0.0
4 H1 -1.0000 0.0000 0.0000 H 1 NUC 0.0
5 H2 0.0000 1.0000 0.0000 H 1 NUC 0.0
6 H3 0.0000 0.0000 1.0000 H 1 NUC 0.0
7 H4 1.5000 0.0000 -1.0000 H 1 NUC 0.0
8 H5 1.5000 0.0000 1.0000 H 1 NUC 0.0
9 H6 2.0000 2.2000 0.2000 H 1 NUC 0.0
10 H7 2.8000 1.3000 0.2000 H 1 NUC 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 1 4 1
4 1 5 1
5 1 6 1
6 2 7 1
7 2 8 1
8 3 9 1
9 3 10 1
""")
    recipe = fragments.prepare_fragment(path, attachment_atom=3)
    fragment = recipe["fragment"]
    assert [atom["element"] for atom in fragment["atoms"][:3]] == ["C", "C", "N"]
    assert fragment["atoms"][2]["position"] == [2.0, 1.3, 0.2]
    assert len(fragment["atoms"]) == 10


def test_peptide_amide_is_not_selected_as_an_amine():
    with pytest.raises(ValueError, match="Amide|amine"):
        fragments.prepare_fragment(smiles="NCC(=O)NCC(=O)O", attachment_atom=5)


def test_cyclodextrin_sized_macrocycle_uses_same_explicit_oh_site():
    # Geometry/graph handling only; this is not evidence for an OH grafting reaction.
    smiles = (
        "OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)"
        "O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)"
        "[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)"
        "[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)"
        "O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]2[C@H](O)[C@@H]1O"
    )
    recipe = fragments.prepare_fragment(smiles=smiles, attachment_atom=1, optimize=False, seed=9)
    fragment = recipe["fragment"]
    assert len(fragment["atoms"]) == 147
    assert sum(atom["element"] == "C" for atom in fragment["atoms"]) == 42
    assert sum(atom["element"] == "O" for atom in fragment["atoms"]) == 35
    _validated_fragment(fragment)


def test_stereochemistry_is_consistent_with_generated_coordinates():
    recipe = fragments.prepare_fragment(smiles="N[C@@H](C)C(=O)O", attachment_atom=1)
    fragment = recipe["fragment"]
    molecule = Chem.AddHs(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    conformer = Chem.Conformer(len(fragment["atoms"]))
    for atom in fragment["atoms"]:
        conformer.SetAtomPosition(atom["id"] - 1, atom["position"])
    molecule.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    assert molecule.GetAtomWithIdx(1).GetProp("_CIPCode") == "S"
