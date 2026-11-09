"""Chemical-context invariants and conservation; no synthetic fit claims."""
import copy
import json
import random
from pathlib import Path

import numpy as np
import pytest

from go_param.modular import BLOCK_REGISTRY, make_resp_input, validate_fitted_charges

ROOT = Path(__file__).resolve().parents[1]


def model(name):
    return json.loads((ROOT / 'models' / name / 'model.json').read_text())


def roles(result):
    return {r['atom_id']: r['fit_role'] for r in result['atoms']}


def test_numbering_model_name_and_geometry_do_not_change_classes():
    from go_param import environments
    m = model('pyrene_diol_epoxide')
    original = roles(environments.classify_model(m))
    changed = copy.deepcopy(m)
    order = list(range(1, len(m['atoms']) + 1))
    random.Random(17).shuffle(order)
    mapping = {old: new for new, old in enumerate(order, 1)}
    changed['model_id'] = 'unrelated_name'
    changed['atoms'] = [copy.deepcopy(m['atoms'][old - 1]) for old in order]
    for a in changed['atoms']:
        a['id'] = mapping[a['id']]
        x, y, z = a['position']
        a['position'] = [-y + 19, x - 7, z + 4]
    changed['bonds'] = [{'a': mapping[b['a']], 'b': mapping[b['b']],
                         'order': b['order']} for b in reversed(m['bonds'])]
    result = roles(environments.classify_model(changed))
    assert all(result[mapping[old]] == value for old, value in original.items())


def test_coo_resonance_drawing_preserves_classes_and_oxygen_equivalence():
    from go_param import environments
    m = model('pyrene_COO')
    before = roles(environments.classify_model(m))
    oxy = [a for a in m['atoms'] if a['role'] == 'COO_O']
    assert before[oxy[0]['id']] == before[oxy[1]['id']]
    oxy[0]['formal_charge'], oxy[1]['formal_charge'] = oxy[1]['formal_charge'], oxy[0]['formal_charge']
    for b in m['bonds']:
        if b['a'] in {a['id'] for a in oxy} or b['b'] in {a['id'] for a in oxy}:
            b['order'] = 3 - b['order']
    assert roles(environments.classify_model(m)) == before


def test_local_chemical_change_creates_a_distinct_hydroxyl_environment():
    from go_param import environments
    cases = json.loads((ROOT / 'calculations/nagl_benchmark_20261105_attempt_001/inputs.json').read_text())['cases']
    ms = {c['model']['model_id']: c['model'] for c in cases}
    def oh_labels(mid):
        return {r['fit_role'] for r in environments.classify_model(ms[mid])['atoms']
                if r['base_role'] == 'hydroxyl_O'}
    assert oh_labels('C24_OH_OH_near_trans') != oh_labels('C24_OH_OH_far_trans')


def test_cis_trans_are_not_spurious_graph_classes():
    from go_param import environments
    a = environments.classify_model(model('pyrene_diol_cis'))
    b = environments.classify_model(model('pyrene_diol_trans'))
    assert sorted(roles(a).values()) == sorted(roles(b).values())


def test_shared_block_solutions_pass_native_card_and_every_total_check():
    from go_param import environments
    ms = [model(n) for n in ['pyrene_pristine', 'pyrene_diol_epoxide', 'pyrene_COO']]
    preserved = copy.deepcopy(ms)
    p = environments.build_problem(ms)
    assert ms == preserved
    assert np.linalg.matrix_rank([r['coefficients'] for r in p['constraints']]) == len(p['constraints'])
    base = {r: q for block in BLOCK_REGISTRY.values()
            for r, q in block['current_role_charges_e'].items()}
    q = {'scaffold_C': 0., **{r['fit_role']: base[r['base_role']]
                             for c in p['classifications'].values() for r in c['atoms']}}
    records = [{'model_id': m['model_id'], 'grid_points': 100} for m in ms for _ in range(3)]
    _, plan = make_resp_input(p, records)
    values = [q[a['role']] for m in p['models'] for _ in range(3) for a in m['atoms']]
    fitted = validate_fitted_charges(p, plan, values)
    assert all(abs(r['charge_error_e']) < 1e-9 for r in fitted['evaluation']['groups'])
    for m in ms:
        pred = environments.predict(m, q)
        assert abs(sum(pred) - m['formal_charge']) < 1e-9


def test_unseen_context_is_rejected_without_fallback():
    from go_param import environments
    m = model('pyrene_pristine')
    p = environments.build_problem([m])
    target = model('pyrene_diol_epoxide')
    assert not environments.coverage(target, p)['fully_supported']
    q = {r: 0. for r in p['roles']}; q['scaffold_C'] = 0.
    with pytest.raises(ValueError, match='unsupported'):
        environments.predict(target, q)


def test_nonfinite_wrong_total_and_unrecognized_chemistry_fail():
    from go_param import environments
    m = model('pyrene_COO')
    c = environments.classify_model(m)
    q = {r['fit_role']: 0. for r in c['atoms']}
    with pytest.raises(ValueError, match='charge|total|block'):
        environments.predict(m, q)
    q[next(iter(q))] = float('nan')
    with pytest.raises(ValueError, match='finite'):
        environments.predict(m, q)
    m['atoms'][0]['role'] = 'unknown_C'
    with pytest.raises(ValueError):
        environments.classify_model(m)
