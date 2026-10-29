"""Fragment chemistry tests: stoichiometry, stable maps, functionality and geometry."""
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from go_param import fragments

EXPECTED = {
    'pyrene_pristine': ('C16H10', 0),
    'pyrene_epoxide': ('C16H10O', 0),
    'pyrene_diol_cis': ('C16H12O2', 0),
    'pyrene_diol_trans': ('C16H12O2', 0),
    'pyrene_dione': ('C16H8O2', 0),
    'pyrene_COOH': ('C17H10O2', 0),
    'pyrene_COO': ('C17H9O2-', -1),
    'pyrene_diol_epoxide': ('C16H12O3', 0),
}


@pytest.fixture(scope='module')
def molecules():
    return {name: fragments.build_fragment(name) for name in EXPECTED}


@pytest.mark.parametrize('name', EXPECTED)
def test_stoichiometry_closed_shell_and_explicit_finite_geometry(molecules, name):
    mol = molecules[name]
    formula, charge = EXPECTED[name]
    assert rdMolDescriptors.CalcMolFormula(mol) == formula
    assert Chem.GetFormalCharge(mol) == charge
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 0
    assert (sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge) % 2 == 0
    assert all(a.GetNumImplicitHs() == 0 for a in mol.GetAtoms())
    assert mol.GetNumConformers() == 1
    xyz = mol.GetConformer().GetPositions()
    assert np.isfinite(xyz).all()
    for bond in mol.GetBonds():
        length = np.linalg.norm(xyz[bond.GetBeginAtomIdx()] - xyz[bond.GetEndAtomIdx()])
        assert 0.85 < length < 1.9


@pytest.mark.parametrize('name', ['pyrene_epoxide', 'pyrene_diol_cis', 'pyrene_diol_trans', 'pyrene_diol_epoxide'])
def test_oxygenated_basal_carbons_are_three_carbon_junctions(molecules, name):
    mol = molecules[name]
    functional = [a for a in mol.GetAtoms() if a.GetProp('role') in ('epoxide_C', 'hydroxyl_C')]
    assert len(functional) >= 2
    for carbon in functional:
        assert sum(n.GetSymbol() == 'C' for n in carbon.GetNeighbors()) == 3
        assert sum(n.GetSymbol() == 'H' for n in carbon.GetNeighbors()) == 0
        assert carbon.GetHybridization() == Chem.HybridizationType.SP3


@pytest.mark.parametrize('name,sign', [('pyrene_diol_cis', 1), ('pyrene_diol_trans', -1)])
def test_diol_faces_and_stereocenters_are_explicit(molecules, name, sign):
    mol = molecules[name]
    oxygens = [a.GetIdx() for a in mol.GetAtoms() if a.GetProp('role') == 'hydroxyl_O']
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C' and a.GetProp('role') != 'hydroxyl_C']
    xyz = mol.GetConformer().GetPositions()
    center = xyz[carbons].mean(axis=0)
    normal = np.linalg.svd(xyz[carbons] - center)[2][-1]
    heights = (xyz[oxygens] - center) @ normal
    assert abs(heights).min() > 0.3
    assert np.sign(np.prod(heights)) == sign
    for atom in mol.GetAtoms():
        if atom.GetProp('role') == 'hydroxyl_C':
            assert atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED


def test_carboxylate_oxygen_equivalence_is_declared_without_fractional_qm_bonds(molecules):
    record = fragments.model_record('pyrene_COO', molecules['pyrene_COO'])
    oxygens = [a for a in record['atoms'] if a['role'] == 'COO_O']
    assert len(oxygens) == 2
    assert len({a['symmetry_class'] for a in oxygens}) == 1
    oxygen_ids = {a['id'] for a in oxygens}
    orders = sorted(b['order'] for b in record['bonds'] if b['a'] in oxygen_ids or b['b'] in oxygen_ids)
    assert orders == [1.0, 2.0]
    assert record['metadata']['carboxylate_resonance_equivalence'] == sorted(oxygen_ids)


def test_export_preserves_atom_ids_roles_and_sdf_graph(molecules, tmp_path):
    catalog = fragments.export_models(tmp_path, names=['pyrene_pristine', 'pyrene_epoxide'])
    assert len(catalog['models']) == 2
    for entry in catalog['models']:
        folder = tmp_path / entry['model_id']
        record = json.loads((folder / 'model.json').read_text())
        sdf = Chem.SDMolSupplier(str(folder / 'initial.sdf'), removeHs=False)[0]
        assert sdf is not None
        assert len(record['atoms']) == sdf.GetNumAtoms()
        assert [a['id'] for a in record['atoms']] == list(range(1, sdf.GetNumAtoms() + 1))
        assert all(a['role'] and a['roles'] and a['name'] for a in record['atoms'])
        assert [a.GetAtomMapNum() for a in sdf.GetAtoms()] == list(range(1, sdf.GetNumAtoms() + 1))
        assert int((folder / 'initial.xyz').read_text().splitlines()[0]) == sdf.GetNumAtoms()
        assert record['multiplicity'] == 1
        assert record['metadata']['geometry_status'] == 'initial_force_field_geometry_not_QM_optimized'


def test_unknown_model_is_rejected():
    with pytest.raises(ValueError, match='Unknown model'):
        fragments.build_fragment('singleOH_radical')
