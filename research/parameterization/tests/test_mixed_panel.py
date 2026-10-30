"""Mixed panels must retain explicit chemistry and report measured separations."""

import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from go_param import mixed_panel
from go_param.modular import partition_model


@pytest.mark.parametrize('shell,formula,carbons', [(2, 'C24H12', 24), (3, 'C54H18', 54)])
def test_analytic_parents_are_closed_shell_compact_pahs(shell, formula, carbons):
    mol = mixed_panel.build_parent(shell)
    assert rdMolDescriptors.CalcMolFormula(mol) == formula
    assert sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms()) == carbons
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    assert all(a.GetAtomMapNum() == a.GetIdx() + 1 for a in mol.GetAtoms())
    assert len(Chem.GetMolFrags(mol)) == 1


@pytest.mark.parametrize('category,oxygen,hydrogen', [('OH_OH', 2, 14),
                                                    ('EPX_EPX', 2, 12),
                                                    ('OH_EPX', 3, 14)])
def test_candidates_preserve_counts_and_perfect_pi_matching(category, oxygen, hydrogen):
    candidates = mixed_panel.enumerate_candidates(2, category)
    assert candidates
    for candidate in candidates:
        sites = [i for motif in candidate['motifs'] for i in motif['carbon_ids']]
        assert len(sites) == len(set(sites))
        assert candidate['sublattice_imbalance'] == 0
        assert candidate['pi_matching_complete']
    candidate = candidates[0]
    mol, record = mixed_panel.build_configuration(2, category, candidate, 'cis')
    assert record['formula'] == f'C24H{hydrogen}O{oxygen}'
    assert record['formal_charge'] == 0
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    assert record['metadata']['multiplicity_basis'].endswith('stability unverified')


@pytest.mark.parametrize('face', ['cis', 'trans'])
def test_reported_faces_and_distances_match_actual_graph_and_coordinates(face):
    candidate = mixed_panel.enumerate_candidates(2, 'OH_EPX')[0]
    mol, record = mixed_panel.build_configuration(2, 'OH_EPX', candidate, face)
    metadata = record['metadata']
    assert metadata['relative_face'] == face
    assert metadata['observed_relative_face'] == face
    motifs = metadata['functional_groups']
    a, b = metadata['measured_pair_motif_indices']
    xyz = mol.GetConformer().GetPositions()
    center_a = xyz[[i - 1 for i in motifs[a]['carbon_ids']]].mean(axis=0)
    center_b = xyz[[i - 1 for i in motifs[b]['carbon_ids']]].mean(axis=0)
    assert metadata['pair_centroid_distance_angstrom'] == pytest.approx(
        float(np.linalg.norm(center_a - center_b)))
    assert metadata['pair_graph_distance'] == candidate['graph_distance']
    assert len(motifs) == 3  # two OH are essential to the balanced pi graph
    assert metadata['geometry_status'] == 'initial_constrained_force_field_geometry_not_QM_optimized'
    assert all(atom.get('charge') is None for atom in record['atoms'])


def test_panel_audits_full_factorial_and_contains_both_parent_sizes(tmp_path):
    catalog = mixed_panel.export_panel(tmp_path)
    models = catalog['models']
    assert len(models) == 34
    assert {entry['parent_carbon_count'] for entry in models} == {24, 54}
    assert {entry['relative_face'] for entry in models} == {'cis', 'trans'}
    assert {entry['distance_bin'] for entry in models} == {'near', 'intermediate', 'far'}
    for category in mixed_panel.CATEGORIES:
        for distance_bin in mixed_panel.BINS:
            selected = [m for m in models if m['category'] == category and m['distance_bin'] == distance_bin]
            assert len({m['graph_distance'] for m in selected}) == 1
    assert catalog['omitted_cells']
    assert all('unattainable' in cell['reason'] for cell in catalog['omitted_cells'])
    assert {entry['control_class'] for entry in catalog['edge_controls']} == {'EPH', 'EPT', 'CAB', 'COOH', 'COO'}
    assert catalog['status'] == 'initial_models_only_no_QM_results_or_fitted_charges'
    loaded = json.loads((tmp_path / 'catalog.json').read_text())
    assert loaded == catalog
    for item in catalog['models']:
        record = json.loads((tmp_path / item['path']).read_text())
        mol = Chem.SDMolSupplier(str(tmp_path / item['model_id'] / 'initial.sdf'), removeHs=False)[0]
        assert mol.GetNumAtoms() == len(record['atoms'])
        assert [a.GetAtomMapNum() for a in mol.GetAtoms()] == [a['id'] for a in record['atoms']]
        assert partition_model(record)['full_coverage']
    for item in catalog['edge_controls']:
        record = json.loads((tmp_path / item['path']).read_text())
        partition = partition_model(record)
        assert partition['full_coverage']
        assert item['control_class'] in {group['block_type'] for group in partition['groups']}
        assert sum(group['target_charge_e'] for group in partition['groups']) == record['formal_charge']
        assert all(atom.get('charge') is None for atom in record['atoms'])


@pytest.mark.parametrize('subtype,counts', [('EPH', [0, 1]), ('EPT', [1, 1])])
def test_terminal_epoxide_control_subtype_is_measured_from_hydrogen_neighbors(subtype, counts):
    mol, record = mixed_panel.build_edge_control(2, subtype)
    group = record['metadata']['functional_groups'][0]
    observed = sorted(sum(n.GetAtomicNum() == 1 for n in mol.GetAtomWithIdx(i - 1).GetNeighbors())
                      for i in group['carbon_ids'])
    assert observed == counts
    assert record['metadata']['epoxide_hydrogen_subtype'] == subtype
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    partition = partition_model(record)
    groups = [group for group in partition['groups'] if group['block_type'] == subtype]
    assert len(groups) == 1
    expected_roles = (['EPH_C', 'EPH_CH', 'EPH_H', 'EPH_O'] if subtype == 'EPH'
                      else ['EPT_C', 'EPT_C', 'EPT_H', 'EPT_H', 'EPT_O'])
    assert sorted(groups[0]['roles']) == sorted(expected_roles)
    assert groups[0]['target_charge_e'] == 0


def test_invalid_or_overlapping_sites_are_rejected():
    candidate = mixed_panel.enumerate_candidates(2, 'OH_OH')[0]
    candidate['motifs'][1]['carbon_ids'] = candidate['motifs'][0]['carbon_ids'][:]
    with pytest.raises(ValueError, match='overlap|distinct'):
        mixed_panel.build_configuration(2, 'OH_OH', candidate, 'cis')
