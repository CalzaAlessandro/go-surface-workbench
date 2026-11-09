"""Count-, distance- and protonation-controlled finite GO charge references.

Original graph construction reuses the project's honeycomb and exact residual
pi matching. New geometries are constrained MMFF starting points, never QM data.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations
import random

import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,rdMolDescriptors

from .mixed_panel import _lattice,_matching

TYPES={'OH','EPX','EPH','EPT','CAB','COOH','COO'}


def _placements(shell,counts,distance_bin,seed):
    points,edges,neighbors,_,distances=_lattice(shell)
    basal=[i for i,ns in neighbors.items() if len(ns)==3]
    boundary=[i for i,ns in neighbors.items() if len(ns)==2]
    if counts.get('OH',0)>len(basal) or sum(counts.get(k,0) for k in ['COOH','COO','CAB'])>len(boundary):
        raise ValueError('requested group count exceeds available sites')
    rng=random.Random(seed)
    for attempt in range(2000):
        used=set();motifs=[];acids=[];failed=False
        for kind in ['EPX','EPH','EPT']:
            for _ in range(counts.get(kind,0)):
                candidates=[(a,b) for a,b in edges if a not in used and b not in used and
                    sum(len(neighbors[i])==2 for i in [a,b])=={'EPX':0,'EPH':1,'EPT':2}[kind]]
                if not candidates:failed=True;break
                pair=rng.choice(candidates);used.update(pair)
                motifs.append({'type':kind,'carbon_ids':[i+1 for i in pair]})
            if failed:break
        if failed:continue
        for kind,allowed in [('OH',basal),('CAB',boundary)]:
            available=[i for i in allowed if i not in used]
            if len(available)<counts.get(kind,0):failed=True;break
            for i in rng.sample(available,counts.get(kind,0)):
                used.add(i);motifs.append({'type':kind,'carbon_ids':[i+1]})
        if failed or _matching(shell,used) is None:continue
        available=[i for i in boundary if i not in used]
        acid_types=[kind for kind in ['COOH','COO'] for _ in range(counts.get(kind,0))]
        if len(available)<len(acid_types):continue
        for kind,i in zip(acid_types,rng.sample(available,len(acid_types))):
            acids.append({'type':kind,'carbon_ids':[i+1]})
        selected=(acids[:2] if len(acids)>1 else [acids[0],motifs[0]] if acids and motifs else motifs[:2])
        distance=None if len(selected)<2 else min(distances[a-1][b-1] for a in selected[0]['carbon_ids'] for b in selected[1]['carbon_ids'])
        if distance_bin=='near' and (distance is None or distance>2):continue
        if distance_bin=='far' and (distance is None or distance<4):continue
        return motifs,acids,attempt,distance
    raise ValueError('unable to place exact counts with requested distance and closed-shell pi matching')


def build_expanded_record(shell,counts,*,distance_bin=None,seed=20261104):
    if not counts or set(counts)-TYPES or any(type(v) is not int or v<1 for v in counts.values()):
        raise ValueError('positive integer counts of supported motifs required')
    if distance_bin not in (None,'near','far'):raise ValueError('unknown distance bin')
    points,edges,neighbors,_,distances=_lattice(shell);nc=len(points)
    motifs,acids,attempt,distance=_placements(shell,counts,distance_bin,seed)
    removed={i-1 for m in motifs for i in m['carbon_ids']}
    matching=_matching(shell,removed);rw=Chem.RWMol()
    for i in range(nc):
        atom=Chem.Atom('C');atom.SetProp('role','edge_C' if len(neighbors[i])==2 else 'scaffold_C');rw.AddAtom(atom)
    for a,b in edges:rw.AddBond(a,b,Chem.BondType.DOUBLE if (a,b) in matching else Chem.BondType.SINGLE)
    positions=[np.array([*p,0.]) for p in points]
    groups=[]
    def add(element,role,position,charge=0):
        atom=Chem.Atom(element);atom.SetProp('role',role);atom.SetFormalCharge(charge)
        i=rw.AddAtom(atom);positions.append(np.asarray(position,dtype=float));return i
    def outward(i):
        v=np.array(points[i])-np.mean(points,axis=0);v/=np.linalg.norm(v)
        return np.array([*v,0.]),np.array([-v[1],v[0],0.])
    for k,m in enumerate(motifs):
        kind=m['type'];ids=[i-1 for i in m['carbon_ids']];face=1 if k%2==0 else -1
        if kind=='CAB':
            normal,_=outward(ids[0]);pos=positions[ids[0]]+1.23*normal
            oxygen=add('O','carbonyl_O',pos);rw.GetAtomWithIdx(ids[0]).SetProp('role','carbonyl_C')
        else:
            for i in ids:positions[i][2]=.2*face
            pos=np.mean([positions[i] for i in ids],axis=0);pos[2]=(1.63 if kind=='OH' else 1.43)*face
            oxygen=add('O',{'OH':'hydroxyl_O','EPX':'epoxide_O','EPH':'EPH_O','EPT':'EPT_O'}[kind],pos)
            for i in ids:
                role=('hydroxyl_C' if kind=='OH' else 'epoxide_C' if kind=='EPX' else
                      'EPT_C' if kind=='EPT' else 'EPH_CH' if len(neighbors[i])==2 else 'EPH_C')
                rw.GetAtomWithIdx(i).SetProp('role',role)
        for i in ids:rw.AddBond(i,oxygen,Chem.BondType.DOUBLE if kind=='CAB' else Chem.BondType.SINGLE)
        groups.append({**m,'oxygen_ids':[oxygen+1],'requested_face':face if kind!='CAB' else None})
    for m in acids:
        kind=m['type'];i=m['carbon_ids'][0]-1;normal,tangent=outward(i)
        rw.GetAtomWithIdx(i).SetProp('role',kind+'_anchor')
        carbon=add('C',kind+'_C',positions[i]+1.48*normal);rw.AddBond(i,carbon,Chem.BondType.SINGLE)
        o1=add('O',kind+'_O',positions[carbon]+.6*normal+1.065*tangent)
        o2=add('O','COOH_OH' if kind=='COOH' else 'COO_O',positions[carbon]+.65*normal-1.2*tangent,-1 if kind=='COO' else 0)
        rw.AddBond(carbon,o1,Chem.BondType.DOUBLE);rw.AddBond(carbon,o2,Chem.BondType.SINGLE)
        groups.append({**m,'carboxyl_carbon_id':carbon+1,'oxygen_ids':[o1+1,o2+1]})
    mol=rw.GetMol();Chem.SanitizeMol(mol)
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,p in enumerate(positions):conf.SetAtomPosition(i,p.tolist())
    mol.AddConformer(conf);mol=Chem.AddHs(mol,addCoords=True)
    for a in mol.GetAtoms():
        if a.GetAtomicNum()==1:
            parent=a.GetNeighbors()[0].GetProp('role')
            a.SetProp('role',{'hydroxyl_O':'hydroxyl_H','COOH_OH':'COOH_H','EPH_CH':'EPH_H','EPT_C':'EPT_H'}.get(parent,'edge_H'))
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):raise ValueError('radical starting graph rejected')
    if not AllChem.MMFFHasAllMoleculeParams(mol):raise ValueError('initial geometry has unsupported MMFF chemistry')
    ff=AllChem.MMFFGetMoleculeForceField(mol,AllChem.MMFFGetMoleculeProperties(mol))
    for i in range(nc):ff.MMFFAddPositionConstraint(i,.35,100.)
    for g in groups:
        if g['type'] in ['OH','EPX','EPH','EPT']:
            for i in g['oxygen_ids']:ff.MMFFAddPositionConstraint(i-1,.4,100.)
    ff.Initialize();status=ff.Minimize(maxIts=4000)
    Chem.AssignAtomChiralTagsFromStructure(mol,replaceExistingTags=True);Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    xyz=mol.GetConformer().GetPositions()
    if not np.isfinite(xyz).all():raise ValueError('nonfinite initial geometry')
    for b in mol.GetBonds():
        length=np.linalg.norm(xyz[b.GetBeginAtomIdx()]-xyz[b.GetEndAtomIdx()])
        if not .75<length<2.0:raise ValueError('implausible initial bond length')
    classes=list(Chem.CanonicalRankAtoms(mol,breakTies=False));next_class=max(classes)+1
    for group in groups:
        if group['type']=='COO':
            for i in group['oxygen_ids']:classes[i-1]=next_class
            next_class+=1
    atoms=[]
    for a in mol.GetAtoms():
        i=a.GetIdx();role=a.GetProp('role');roles={role}
        if i<nc:
            roles.add('scaffold_C')
            if len(neighbors[i])==2:roles.add('edge_C')
        atoms.append({'id':i+1,'name':f'{a.GetSymbol()}{i+1}','element':a.GetSymbol(),
            'position':xyz[i].tolist(),'role':role,'roles':sorted(roles),'formal_charge':a.GetFormalCharge(),
            'charge':None,'symmetry_class':int(classes[i]),'chiral_tag':str(a.GetChiralTag()),
            'cip_code':a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None})
    pair_distances=[{'group_indices':[i,j],'distance':min(distances[a-1][b-1] for a in g['carbon_ids'] for b in h['carbon_ids'])}
                    for (i,g),(j,h) in combinations(enumerate(groups),2)]
    q=Chem.GetFormalCharge(mol);electron_count=sum(a.GetAtomicNum() for a in mol.GetAtoms())-q
    if q!=-counts.get('COO',0) or electron_count%2:raise ValueError('formal-charge/electron parity mismatch')
    model_id=f'C{nc}_'+ '_'.join(f'{k}{v}' for k,v in sorted(counts.items()))+'_'+(distance_bin or 'loading')+f'_s{seed}'
    record={'schema_version':1,'model_id':model_id,'formula':rdMolDescriptors.CalcMolFormula(mol),
        'formal_charge':q,'multiplicity':1,'electron_count':electron_count,'atoms':atoms,
        'bonds':[{'a':b.GetBeginAtomIdx()+1,'b':b.GetEndAtomIdx()+1,'order':b.GetBondTypeAsDouble()} for b in mol.GetBonds()],
        'metadata':{'scope':'Expanded finite-fragment charge reference; starting structure only',
            'parent_carbon_count':nc,'scaffold_carbon_ids':list(range(1,nc+1)),
            'functional_groups':groups,'actual_group_counts':dict(Counter(g['type'] for g in groups)),
            'requested_group_counts':dict(counts),'pair_graph_distances':pair_distances,
            'requested_distance_bin':distance_bin,'selected_pair_distance':distance,
            'graph_distance_definition':'minimum C-C scaffold bond distance between group anchors',
            'pi_matching_complete':True,'placement_attempt':attempt,'seed':seed,'rdkit_version':rdBase.rdkitVersion,
            'geometry_status':'initial_constrained_force_field_geometry_not_QM_optimized',
            'initial_geometry_method':'analytic honeycomb plus constrained MMFF94',
            'initial_MMFF_converged':status==0,'coordinate_units':'angstrom',
            'canonical_isomeric_SMILES':Chem.MolToSmiles(mol,isomericSmiles=True),
            'wavefunction_stability_validated':False,'fitted_partial_charges_available':False}}
    return mol,record
