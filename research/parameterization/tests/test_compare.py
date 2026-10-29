"""Charge transfer guards: prevent invented coverage and mismatched ESP ranking."""
import json
from pathlib import Path
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((ROOT / 'models' / name / 'model.json').read_text())


def test_current_mapping_preserves_ids_and_functional_precedence():
    from go_param.compare import assign_candidate
    result = assign_candidate(ROOT / 'models/pyrene_COOH/model.json', 'current_gaff2')
    charges = {a['atom_id']: a['charge_e'] for a in result['atoms']}
    assert charges[1] == -0.115
    assert charges[17] == 0.735
    assert charges[29] == 0.46
    assert result['coverage']['complete']
    assert result['charge_error_e'] == pytest.approx(0, abs=1e-12)
    assert result['primary_esp_eligible']


@pytest.mark.parametrize('candidate', ['old_whole_sheet', 'pinto_gaff2_aim'])
def test_context_specific_sheets_cannot_be_averaged_onto_fragment(candidate):
    from go_param.compare import assign_candidate
    result = assign_candidate(load('pyrene_diol_cis'), candidate)
    assert not result['coverage']['complete']
    assert all(a['charge_e'] is None for a in result['atoms'])
    assert result['total_charge_e'] is None
    assert not result['primary_esp_eligible']


def test_gopy_raw_nonzero_motif_is_retained_but_not_ranked():
    from go_param.compare import assign_candidate
    result = assign_candidate(load('pyrene_epoxide'), 'gopy_raw')
    assert result['charge_error_e'] == pytest.approx(-0.3691)
    assert result['coverage']['complete']
    assert result['model_kind'] == 'hybrid_control'
    assert not result['primary_esp_eligible']


def test_unsupported_carboxylate_is_not_scaffold_fallback():
    from go_param.compare import assign_candidate
    result = assign_candidate(load('pyrene_COO'), 'hiergo_opls')
    assert result['coverage']['unsupported_roles'] == ['COO_C', 'COO_O', 'COO_anchor']
    assert result['total_charge_e'] is None
    assert result['charge_error_e'] is None


def test_stauffer_pure_caps_missing_hybrid_explicit_and_finite():
    from go_param.compare import assign_candidate
    pure = assign_candidate(load('pyrene_diol_cis'), 'stauffer_esp')
    assert not pure['coverage']['complete']
    hybrid = assign_candidate(load('pyrene_diol_cis'), 'stauffer_sesp_hybrid')
    assert hybrid['coverage']['complete']
    assert hybrid['model_kind'] == 'hybrid_control'
    assert hybrid['total_charge_e'] == pytest.approx(0, abs=1e-12)
    assert next(a for a in hybrid['atoms'] if a['role'] == 'hydroxyl_O')['charge_e'] == -0.57


def reference(model):
    return {'atom_ids': [a['id'] for a in model['atoms']],
            'atom_positions_angstrom': [a['position'] for a in model['atoms']],
            'formal_charge': model['formal_charge'], 'points_angstrom': [[20., 0., 0.], [0., 20., 0.]],
            'potential_au': [0., 0.], 'dipole_origin_angstrom': [0., 0., 0.]}


def test_esp_exact_zero_charges_and_chunking_and_mismatch_guard():
    from go_param.compare import compare_esp
    model = load('pyrene_pristine')
    ref = reference(model)
    candidate = {'candidate': 'synthetic_zero', 'atoms': [{'atom_id': a['id'], 'charge_e': 0.} for a in model['atoms']],
                 'primary_esp_eligible': True, 'coverage': {'complete': True}}
    output = compare_esp(model, ref, [candidate], chunk_size=1)
    assert output['results'][0]['rmse_au'] == 0
    assert output['results'][0]['dipole_e_angstrom'] == [0., 0., 0.]
    ref['atom_ids'] = list(reversed(ref['atom_ids']))
    with pytest.raises(ValueError, match='identity'):
        compare_esp(model, ref, [candidate])


def test_esp_charge_mismatch_only_explicit_diagnostic():
    from go_param.compare import compare_esp
    model = load('pyrene_epoxide')
    output = compare_esp(model, reference(model), ['gopy_raw'])
    assert output['results'][0]['status'] == 'skipped'
    diagnostic = compare_esp(model, reference(model), ['gopy_raw'], allow_charge_mismatch_diagnostic=True)
    assert diagnostic['results'][0]['status'] == 'diagnostic_charge_mismatch'
    assert np.isfinite(diagnostic['results'][0]['rmse_au'])


def test_reject_unknown_roles_and_nonfinite_geometry():
    from go_param.compare import assign_candidate, compare_esp
    model = load('pyrene_pristine')
    model['atoms'][0]['role'] = 'unpublished_C'
    model['atoms'][0]['roles'].append('unpublished_C')
    assert not assign_candidate(model, 'current_gaff2')['coverage']['complete']
    model['atoms'][0]['position'][0] = float('nan')
    with pytest.raises(ValueError, match='finite'):
        compare_esp(model, reference(model), ['current_gaff2'])


def test_coulomb_units_nonzero_charge_and_grid_intersection():
    from go_param.compare import compare_esp
    model = {'model_id': 'synthetic', 'formal_charge': 1,
             'atoms': [{'id': 8, 'element': 'H', 'position': [0., 0., 0.]}]}
    ref = {'atom_ids': [8], 'atom_positions_angstrom': [[0., 0., 0.]], 'formal_charge': 1,
           'points_angstrom': [[0.529177210903, 0., 0.]], 'potential_au': [1.],
           'dipole_origin_angstrom': [1., 0., 0.]}
    candidate = {'candidate': 'unit_charge', 'atoms': [{'atom_id': 8, 'charge_e': 1.}],
                 'coverage': {'complete': True}}
    result = compare_esp(model, ref, [candidate], chunk_size=1)['results'][0]
    assert result['rmse_au'] == pytest.approx(0.)
    assert result['dipole_e_angstrom'] == [-1., 0., 0.]
    ref['points_angstrom'] = [[0., 0., 0.]]
    with pytest.raises(ValueError, match='intersects'):
        compare_esp(model, ref, [candidate])


def test_export_includes_literal_role_charges_with_source_urls(tmp_path):
    from go_param.compare import export_reference_summary
    (tmp_path/'models/pyrene_epoxide').mkdir(parents=True)
    (tmp_path/'models/pyrene_epoxide/model.json').write_text(json.dumps(load('pyrene_epoxide')))
    (tmp_path/'references').mkdir()
    (tmp_path/'references/reference_catalog.json').write_text((ROOT/'references/reference_catalog.json').read_text())
    result = export_reference_summary(tmp_path)
    row = next(r for r in result['role_charges'] if r['candidate'] == 'hiergo_opls' and r['role'] == 'hydroxyl_O')
    assert row['charge_e'] == -0.683
    assert row['source_url'] == 'https://github.com/IFM-molecular-simulation-group/HierGO'
    assert next(r for r in result['role_charges'] if r['candidate'] == 'hiergo_opls' and r['role'] == 'ether_O')['charge_e'] == -0.17
    assert (tmp_path/'reports/reference_charge_summary.roles.csv').exists()


@pytest.mark.parametrize('tolerance', [0, -1, float('nan'), float('inf'), True])
def test_explicit_charge_tolerance_must_be_finite_positive(tolerance):
    from go_param.compare import compare_esp
    model = load('pyrene_pristine')
    with pytest.raises(ValueError, match='tolerance'):
        compare_esp(model, reference(model), [], charge_tolerance_e=tolerance)


def test_recorded_resp_rounding_tolerance_never_renormalizes():
    from go_param.compare import compare_esp
    model = load('pyrene_pristine')
    q = [0.]*len(model['atoms']); q[0] = 3e-6
    candidate = {'candidate': 'rounding_fixture', 'atoms': [
        {'atom_id': a['id'], 'charge_e': v} for a, v in zip(model['atoms'], q)], 'coverage': {'complete': True}}
    strict = compare_esp(model, reference(model), [candidate])
    assert strict['results'][0]['status'] == 'skipped'
    scored = compare_esp(model, reference(model), [candidate], charge_tolerance_e=5e-6)
    assert scored['charge_tolerance_e'] == 5e-6
    assert scored['results'][0]['status'] == 'matched_charge_control'
    assert scored['results'][0]['charge_error_e'] == 3e-6
    assert candidate['atoms'][0]['charge_e'] == 3e-6
    assert scored['results'][0]['rmse_au'] > 0
