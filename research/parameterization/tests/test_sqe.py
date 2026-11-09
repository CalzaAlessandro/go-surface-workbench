"""Numerical gates for a RESP-reference SQE research variant, not FF validation."""
from copy import deepcopy
import json
import math
from pathlib import Path

import numpy as np
import pytest


def test_two_atom_solution_matches_analytic_equation_with_initial_charges():
    from go_param.sqe import solve_sqe

    distance=1.5
    interaction=math.erf(distance/math.sqrt(2*.7**2+2*.5**2))/distance
    transfer=(-.5-.8*interaction)/(4+3-2*interaction+2)
    q=solve_sqe([[0,0,0],[distance,0,0]],[[0,1]],chi=[.3,-.2],
                hardness=[4,3],kappa=[2],widths=[.7,.5],reference=[-.4,.4])
    np.testing.assert_allclose(q,[-.4+transfer,.4-transfer],atol=1e-13)


def test_permutation_bond_direction_rotation_and_gauge_invariance():
    from go_param.sqe import solve_sqe

    xyz=np.array([[0,0,0],[1.4,0,0],[2.1,1.2,0],[5,0,0.]])
    bonds=np.array([[0,1],[1,2]])
    chi=np.array([0.,1.,2.,9.]);eta=np.array([4.,5.,6.,7.]);width=np.full(4,.7)
    reference=np.array([-.8,-.1,-.1,1.])
    def calc(x,b,c,h,w,r):
        return solve_sqe(x,b,chi=c,hardness=h,kappa=[3,3],widths=w,reference=r)
    q=calc(xyz,bonds,chi,eta,width,reference)
    assert q[:3].sum()==pytest.approx(-1,abs=1e-13)
    assert q[3]==1
    np.testing.assert_allclose(calc(xyz,bonds[:,::-1],chi+31,eta,width,reference),q,atol=1e-12)
    order=np.array([3,1,0,2]);reverse=np.argsort(order)
    permuted=calc(xyz[order],reverse[bonds],chi[order],eta[order],width[order],reference[order])
    np.testing.assert_allclose(permuted,q[order],atol=1e-12)
    rot=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    np.testing.assert_allclose(calc(xyz@rot+[3,-2,4],bonds,chi,eta,width,reference),q,atol=1e-12)


def test_invalid_geometry_and_unstable_system_are_rejected():
    from go_param.sqe import solve_sqe

    with pytest.raises(ValueError,match='coincident'):
        solve_sqe([[0,0,0],[0,0,0]],[[0,1]],chi=[0,1],hardness=[4,4],
                  kappa=[2],widths=[.7,.7],reference=[0,0])
    with pytest.raises(ValueError,match='positive definite'):
        solve_sqe([[0,0,0],[1,0,0]],[[0,1]],chi=[0,1],hardness=[.01,.01],
                  kappa=[-9],widths=[.7,.7],reference=[0,0])


def test_role_response_is_resonance_invariant_and_unsupported_roles_fail():
    from go_param.sqe import prepare_response

    path=Path(__file__).resolve().parents[1]/'calculations/espaloma_benchmark_20261103_attempt_001/inputs.json'
    cases=json.loads(path.read_text())['cases']
    case=next(c for c in cases if c['model']['model_id']=='pyrene_COO')
    model=case['model'];baseline=case['references']['current_library']
    a=prepare_response(model,baseline)
    changed=deepcopy(model)
    oxygens=[v for v in changed['atoms'] if v['role']=='COO_O']
    assert len(oxygens)==2
    oxygens[0]['formal_charge'],oxygens[1]['formal_charge']=oxygens[1]['formal_charge'],oxygens[0]['formal_charge']
    ids={v['id'] for v in oxygens}
    for bond in changed['bonds']:
        if bond['a'] in ids or bond['b'] in ids:bond['order']=3-bond['order']
    b=prepare_response(changed,baseline)
    assert a['roles']==b['roles']
    np.testing.assert_allclose(a['offset'],b['offset'],atol=1e-13)
    np.testing.assert_allclose(a['design'],b['design'],atol=1e-13)
    changed['atoms'][0]['role']='unknown_chemical_role'
    with pytest.raises(ValueError,match='unsupported'):
        prepare_response(changed,baseline)


def test_response_linear_fit_recovers_synthetic_observable_charges():
    from go_param.sqe import prepare_response, fit_response, predict_response

    path=Path(__file__).resolve().parents[1]/'calculations/espaloma_benchmark_20261103_attempt_001/inputs.json'
    cases=json.loads(path.read_text())['cases'][:7]
    examples=[prepare_response(c['model'],c['references']['current_library']) for c in cases]
    roles=examples[0]['roles']
    theta=np.array([.15*math.sin(i) for i in range(len(roles))])
    targets=[e['offset']+e['design']@theta for e in examples]
    fitted=fit_response(examples,targets,ridge=1e-12)
    for e,target in zip(examples,targets):
        np.testing.assert_allclose(predict_response(e,fitted),target,atol=1e-7)
    assert fitted['gauge_reference_role']=='scaffold_C'


def test_native_chargefw_initial_charge_equations_agree_with_recorded_controls():
    from go_param.sqe import solve_sqe

    fixture=Path(__file__).parent/'data/sqe_native_controls.json'
    controls=json.loads(fixture.read_text())
    for c in controls:
        actual=solve_sqe(c['positions'],c['bonds'],chi=c['chi'],hardness=c['hardness'],
                         kappa=c['kappa'],widths=c['widths'],reference=c['reference'])
        np.testing.assert_allclose(actual,c['expected_native_charges'],atol=2e-10)
