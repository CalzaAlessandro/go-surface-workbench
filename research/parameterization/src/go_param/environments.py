"""First-neighbor RESP research classes, not the published GAFF2-AIM preset.

Chemical identity comes from the validated GO role/graph schema. Coordinates,
model names, atom numbers and arbitrary Kekule/resonance bond orders never enter
a class key. This controlled pilot retains L0's block sums and zero background.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from collections import Counter

import numpy as np

from . import modular

POLICY = {
    'name': 'GO_first_neighbor_group_context_v1',
    'native_GAFF2_AIM': False,
    'background_carbon': 'fixed_zero_in_this_controlled_pilot',
    'edge_caps': 'shared_GRT_finite_H_caps',
    'unseen_context': 'unsupported_no_fallback',
    'relative_face': 'not_encoded_graph_only',
    'bond_orders': 'validated_upstream_not_used_as_resonance_sensitive_features',
}


def classify_model(value):
    """Return a deterministic role-to-environment assignment without mutation."""
    model = modular._read(value)
    partition = modular.partition_model(model)
    atoms = {a['id']: a for a in model['atoms']}
    neighbors = {i: [] for i in atoms}
    for b in model['bonds']:
        neighbors[b['a']].append(b['b'])
        neighbors[b['b']].append(b['a'])

    def carbon_environment(i):
        return [atoms[i]['role'], sorted(atoms[j]['role'] for j in neighbors[i])]

    catalog = {}
    rows = {}
    for group in partition['groups']:
        context = sorted([carbon_environment(i) for i in group['atom_ids']
                          if atoms[i]['element'] == 'C'], key=lambda x: json.dumps(x))
        for i in group['atom_ids']:
            a = atoms[i]
            descriptor = {'policy': POLICY['name'], 'base_role': a['role'],
                          'element': a['element'], 'block': group['block_type']}
            if group['block_type'] != 'GRT':
                descriptor['group_carbon_first_neighbors'] = context
                if a['element'] == 'C':
                    descriptor['carbon_site_first_neighbors'] = carbon_environment(i)
            canonical = json.dumps(descriptor, sort_keys=True, separators=(',', ':'))
            label = a['role'] if group['block_type'] == 'GRT' else (
                a['role'] + '@' + hashlib.sha256(canonical.encode()).hexdigest()[:16])
            if label in catalog and catalog[label] != descriptor:
                raise ValueError('Environment class collision')
            catalog[label] = descriptor
            rows[i] = {'atom_id': i, 'base_role': a['role'], 'fit_role': label,
                       'block_type': group['block_type']}
    for i in partition['scaffold_atom_ids']:
        rows[i] = {'atom_id': i, 'base_role': 'scaffold_C',
                   'fit_role': 'scaffold_C', 'block_type': 'background'}
    assert set(rows) == set(atoms)
    return {'model_id': model['model_id'], 'policy': copy.deepcopy(POLICY),
            'atoms': [rows[i] for i in sorted(rows)], 'classes': catalog,
            'partition': partition}


def build_problem(values):
    """Build a class-shared constrained problem for the native RESP card writer.

    Every distinct block equation is retained for output validation; only a
    consistent independent subset is emitted to RESP to avoid singular cards.
    """
    originals = [modular._read(v) for v in values]
    # Reuse the original chemistry, exclusion and no-stage-2 gates.
    modular.build_problem(originals)
    classifications = {m['model_id']: classify_model(m) for m in originals}
    models = copy.deepcopy(originals)
    partitions = {}
    elements = {}
    support = {}
    for m in models:
        classification = classifications[m['model_id']]
        assignment = {r['atom_id']: r['fit_role'] for r in classification['atoms']}
        for a in m['atoms']:
            role = assignment[a['id']]
            if role in elements and elements[role] != a['element']:
                raise ValueError('Environment class mixes elements')
            elements[role] = a['element']
            a['base_role'] = a['role']
            a['role'] = role
            support.setdefault(role, set()).add(m['model_id'])
        partition = copy.deepcopy(classification['partition'])
        for g in partition['groups']:
            g['base_roles'] = g['roles']
            g['roles'] = [assignment[i] for i in g['atom_ids']]
        partitions[m['model_id']] = partition
    roles = sorted(r for r in elements if r != 'scaffold_C')
    constraints = []
    all_constraints = []
    seen = {}
    a_kept = np.empty((0, len(roles)))
    b_kept = np.empty(0)
    for m in models:
        for g in partitions[m['model_id']]['groups']:
            counts = Counter(g['roles'])
            coefficients = [counts[r] for r in roles]
            key = tuple(coefficients)
            target = g['target_charge_e']
            if key in seen:
                if seen[key] != target:
                    raise ValueError('Infeasible shared class constraints')
                continue
            seen[key] = target
            row = {'block_type': g['block_type'], 'coefficients': coefficients,
                   'target_charge_e': target, 'representative_model_id': m['model_id'],
                   'representative_atom_ids': g['atom_ids']}
            all_constraints.append(row)
            test_a = np.vstack((a_kept, coefficients))
            test_b = np.append(b_kept, target)
            rank = int(np.linalg.matrix_rank(test_a))
            if np.linalg.matrix_rank(np.column_stack((test_a, test_b))) != rank:
                raise ValueError('Infeasible environment block constraints')
            if rank > len(constraints):
                constraints.append(row)
                a_kept, b_kept = test_a, test_b
    observed = {g['block_type'] for p in partitions.values() for g in p['groups']}
    return {
        'status': 'L1_RESP_RESEARCH_NOT_PRODUCTION', 'scope': 'first_neighbor_classes',
        'policy': copy.deepcopy(POLICY), 'models': models,
        'model_order': [m['model_id'] for m in models], 'partitions': partitions,
        'classifications': classifications, 'roles': roles, 'role_elements': elements,
        'constraints': constraints, 'all_distinct_constraints': all_constraints,
        'constraint_rank': len(constraints),
        'free_dimensions_after_constraints': len(roles) - len(constraints),
        'fixed_roles': {'scaffold_C': 0.0},
        'class_support_models': {k: sorted(v) for k, v in support.items()},
        'unfitted_terminal_blocks': sorted({'EPH', 'EPT'} - observed),
        'unobserved_library_units': [k for k, v in modular.BLOCK_REGISTRY.items()
                                   if k != 'GGG' and v['block_type'] not in observed],
    }


def coverage(model, problem):
    """Report strict class support; atoms and models are different sample units."""
    c = classify_model(model)
    known = set(problem['roles']) | {'scaffold_C'}
    unknown = sorted({r['fit_role'] for r in c['atoms']} - known)
    return {'model_id': c['model_id'], 'fully_supported': not unknown,
            'unsupported_classes': unknown, 'atoms': len(c['atoms']),
            'supported_atoms': sum(r['fit_role'] in known for r in c['atoms'])}


def predict(model, charges):
    """Map a fitted class library with no guessed classes or normalization."""
    c = classify_model(model)
    missing = {r['fit_role'] for r in c['atoms']} - set(charges)
    if missing:
        raise ValueError('unsupported environment classes: ' + ', '.join(sorted(missing)))
    q = [float(charges[r['fit_role']]) for r in c['atoms']]
    if not all(math.isfinite(x) for x in q):
        raise ValueError('Nonfinite environment charges')
    if any(q[r['atom_id'] - 1] != 0 for r in c['atoms'] if r['fit_role'] == 'scaffold_C'):
        raise ValueError('Frozen background charge is not zero')
    for g in c['partition']['groups']:
        if abs(sum(q[i - 1] for i in g['atom_ids']) - g['target_charge_e']) > len(g['atom_ids']) * 5.1e-7 + 2e-6:
            raise ValueError('Environment block charge constraint failed')
    formal = modular._read(model)['formal_charge']
    if abs(sum(q) - formal) > len(q) * 5.1e-7 + 2e-6:
        raise ValueError('Environment total charge constraint failed')
    return q
