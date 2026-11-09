"""Research SQE variant with fixed RESP reference charges and invariant roles.

The core follows the SQE initial-charge equations (Schindler et al., 2021).
Role coefficients fitted here are not the published NPA parameterization.
No orbital/density calculation or dynamic charge update is performed.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.special import erf

ROLES = tuple(sorted({
    'scaffold_C','edge_C','edge_H','hydroxyl_C','hydroxyl_O','hydroxyl_H',
    'epoxide_C','epoxide_O','carbonyl_C','carbonyl_O','COOH_anchor','COOH_C',
    'COOH_O','COOH_OH','COOH_H','COO_anchor','COO_C','COO_O',
    'EPH_C','EPH_CH','EPH_O','EPH_H','EPT_C','EPT_O','EPT_H'}))
FREE_ROLES=tuple(r for r in ROLES if r!='scaffold_C')
DEFAULTS={'hardness':4.0,'bond_hardness':10.0,'width_scale':1.0}
WIDTHS={'C':.7,'O':.6,'H':.4}


def _operators(positions,bonds,hardness,kappa,widths):
    xyz=np.asarray(positions,dtype=float)
    if xyz.ndim!=2 or xyz.shape[1]!=3 or len(xyz)==0 or not np.isfinite(xyz).all():
        raise ValueError('positions must be a finite nonempty N by 3 array')
    n=len(xyz);edges=np.asarray(bonds,dtype=int).reshape((-1,2))
    eta=np.asarray(hardness,dtype=float);width=np.asarray(widths,dtype=float)
    kap=np.asarray(kappa,dtype=float)
    if eta.shape!=(n,) or width.shape!=(n,) or kap.shape!=(len(edges),):
        raise ValueError('inconsistent SQE parameter dimensions')
    if not all(np.isfinite(v).all() for v in [eta,width,kap]) or np.any(width<0):
        raise ValueError('invalid SQE parameter values')
    if len(edges) and (edges.min()<0 or edges.max()>=n or np.any(edges[:,0]==edges[:,1])):
        raise ValueError('invalid bond endpoints')
    if len({tuple(sorted(e)) for e in edges})!=len(edges):
        raise ValueError('duplicate bond endpoints')
    distance=np.linalg.norm(xyz[:,None,:]-xyz[None,:,:],axis=2)
    offdiag=~np.eye(n,dtype=bool)
    if np.any(distance[offdiag]<1e-10):raise ValueError('coincident atom positions')
    denominator=np.sqrt(2*width[:,None]**2+2*width[None,:]**2)
    argument=np.divide(distance,denominator,out=np.full((n,n),np.inf),where=denominator>0)
    interaction=np.divide(erf(argument),distance,out=np.zeros((n,n)),where=offdiag)
    matrix=interaction+np.diag(eta)
    incidence=np.zeros((n,len(edges)))
    for i,(a,b) in enumerate(edges):incidence[a,i]=1;incidence[b,i]=-1
    split=incidence.T@matrix@incidence+np.diag(kap)
    if not len(edges):return interaction,np.zeros((n,n))
    try:factor=cho_factor(split,lower=True,check_finite=True)
    except np.linalg.LinAlgError as exc:
        raise ValueError('SQE split-charge system is not positive definite') from exc
    response=incidence@cho_solve(factor,incidence.T)
    return interaction,response


def solve_sqe(positions,bonds,*,chi,hardness,kappa,widths,reference):
    """Return q0 - R(chi + K_off q0), conserving every component's q0 sum."""
    interaction,response=_operators(positions,bonds,hardness,kappa,widths)
    q0=np.asarray(reference,dtype=float);electro=np.asarray(chi,dtype=float)
    n=len(response)
    if q0.shape!=(n,) or electro.shape!=(n,) or not np.isfinite(q0).all() or not np.isfinite(electro).all():
        raise ValueError('invalid initial charges or electronegativities')
    charges=q0-response@(electro+interaction@q0)
    if not np.isfinite(charges).all():raise ValueError('nonfinite SQE charges')
    return charges


def prepare_response(model,baseline,*,hardness=4.,bond_hardness=10.,width_scale=1.):
    """Build an affine response independent of atom IDs and resonance labels.

    Bond connectivity and geometry matter. Numeric bond orders and localized
    COO formal-charge labels do not define coefficients in this research variant.
    """
    atoms=model['atoms'];n=len(atoms);ids={a['id']:i for i,a in enumerate(atoms)}
    if len(ids)!=n:raise ValueError('duplicate atom IDs')
    atom_roles=[a['role'] for a in atoms]
    unknown=set(atom_roles)-set(ROLES)
    if unknown:raise ValueError('unsupported chemical roles: '+','.join(sorted(unknown)))
    q0=np.asarray(baseline,dtype=float)
    tolerance=max(1e-6,n*5.1e-7)
    if q0.shape!=(n,) or abs(q0.sum()-model['formal_charge'])>tolerance:
        raise ValueError('reference library does not preserve formal charge')
    if hardness<=0 or bond_hardness<=0 or width_scale<=0:
        raise ValueError('fitted-model response parameters must be positive')
    widths=np.array([WIDTHS[a['element']]*width_scale for a in atoms])
    # The Gaussian Coulomb kernel including self terms is PSD. This lower
    # bound on eta makes its diagonal replacement positive for any sheet size.
    if hardness<=float(np.max(1/(np.sqrt(np.pi)*widths))):
        raise ValueError('atomic hardness is too small for a stable global kernel')
    edges=[[ids[b['a']],ids[b['b']]] for b in model['bonds']]
    interaction,response=_operators([a['position'] for a in atoms],edges,np.full(n,hardness),
                                    np.full(len(edges),bond_hardness),widths)
    features=np.array([[float(role==candidate) for candidate in FREE_ROLES] for role in atom_roles])
    return {'model_id':model['model_id'],'roles':list(FREE_ROLES),'atom_roles':atom_roles,
            'offset':q0-response@interaction@q0,'design':-response@features,
            'baseline':q0,'formal_charge':model['formal_charge'],
            'response_parameters':{'hardness':hardness,'bond_hardness':bond_hardness,'width_scale':width_scale}}


def fit_response(examples,targets,*,ridge=1e-4):
    """Fit role chi with equal molecule weight and equal present-role weight."""
    if not examples or len(examples)!=len(targets) or ridge<0:
        raise ValueError('invalid SQE fit inputs')
    roles=examples[0]['roles'];blocks=[];values=[];observed=set()
    for example,target in zip(examples,targets):
        if example['roles']!=roles:raise ValueError('inconsistent role basis')
        target=np.asarray(target,dtype=float)
        if target.shape!=example['offset'].shape or not np.isfinite(target).all():
            raise ValueError('invalid reference target')
        tolerance=max(1e-6,len(target)*5.1e-7)
        if abs(target.sum()-example['formal_charge'])>tolerance:
            raise ValueError('target does not preserve formal charge')
        counts=Counter(example['atom_roles']);observed.update(counts)
        weight=np.sqrt([1/(len(examples)*len(counts)*counts[r]) for r in example['atom_roles']])
        blocks.append(example['design']*weight[:,None]);values.append((target-example['offset'])*weight)
    design=np.vstack(blocks);residual=np.concatenate(values)
    # Ridge fixes poorly identified relative coefficients without inventing
    # coefficients for chemistry absent from the training molecules.
    active=np.array([r in observed for r in roles])
    a=design[:,active]
    augmented=np.vstack([a,np.sqrt(ridge)*np.eye(a.shape[1])])
    y=np.concatenate([residual,np.zeros(a.shape[1])])
    fitted,_,rank,singular=np.linalg.lstsq(augmented,y,rcond=None)
    theta=np.zeros(len(roles));theta[active]=fitted
    return {'method':'SQE_fixed_RESP_reference_role_fit','roles':roles,'theta':theta.tolist(),
            'observed_roles':sorted(observed),'gauge_reference_role':'scaffold_C',
            'ridge':ridge,'training_models':[e['model_id'] for e in examples],
            'augmented_rank':int(rank),'singular_values':singular.tolist(),
            'response_parameters':examples[0]['response_parameters'],
            'production_parameter_validation':False}


def predict_response(example,fit):
    if example['roles']!=fit['roles'] or example['response_parameters']!=fit['response_parameters']:
        raise ValueError('response settings differ from fitted model')
    unsupported=set(example['atom_roles'])-set(fit['observed_roles'])
    if unsupported:raise ValueError('unsupported roles absent from training: '+','.join(sorted(unsupported)))
    q=example['offset']+example['design']@np.asarray(fit['theta'])
    if abs(q.sum()-example['formal_charge'])>max(1e-6,len(q)*5.1e-7):
        raise ValueError('predicted charge differs from formal total')
    return q
