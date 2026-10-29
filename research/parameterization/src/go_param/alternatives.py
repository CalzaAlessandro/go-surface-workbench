"""System AmberTools AM1-BCC worker; standard library only, no charge substitute."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

ELEMENTS = {'H': 1, 'C': 6, 'O': 8}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def atom_name(atom):
    """Element-compatible names fit Amber's four-character naming convention."""
    if atom['element'] not in ELEMENTS or type(atom['id']) is not int or not 1 <= atom['id'] <= 999:
        raise ValueError('Atom name requires a C/H/O element and integer ID from 1 to 999')
    return f'{atom["element"]}{atom["id"]:03d}'


def validate_model(model):
    atoms = model['atoms']
    if not atoms or [a['id'] for a in atoms] != list(range(1, len(atoms)+1)) or any(type(a['id']) is not int for a in atoms):
        raise ValueError('Atom IDs must be contiguous ordered integers from 1')
    if len(atoms) > 999:
        raise ValueError('Atom IDs exceed supported 999-name/V2000 limit')
    if type(model['formal_charge']) is not int or model.get('multiplicity') != 1:
        raise ValueError('AM1-BCC campaign requires integer formal charge and singlet multiplicity')
    if any(a['element'] not in ELEMENTS or len(a['position']) != 3 or not all(math.isfinite(x) for x in a['position']) for a in atoms):
        raise ValueError('Finite C/H/O model coordinates are required')
    if sum(a['formal_charge'] for a in atoms) != model['formal_charge']:
        raise ValueError('Atomic and whole-model formal charge disagree')
    if (sum(ELEMENTS[a['element']] for a in atoms)-model['formal_charge']) % 2:
        raise ValueError('Closed-shell electron parity mismatch')
    seen = set()
    for b in model['bonds']:
        pair = tuple(sorted((b['a'], b['b'])))
        if pair in seen or not 1 <= pair[0] < pair[1] <= len(atoms) or b['order'] not in (1, 1.5, 2, 3):
            raise ValueError('Invalid model bond graph')
        seen.add(pair)


def sdf_text(model):
    """Reviewable V2000 graph with atom-map IDs and explicit formal charges."""
    n, nb = len(model['atoms']), len(model['bonds'])
    if max(n, nb) > 999:
        raise ValueError('V2000 counts exceed supported limit')
    lines = [model['model_id'], 'go_param            3D', '', f'{n:3d}{nb:3d}  0  0  0  0            999 V2000']
    for atom in model['atoms']:
        x, y, z = atom['position']
        lines.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {atom["element"]:<3s}' + ' 0  0  0  0  0  0  0  0  0' + f'{atom["id"]:3d}' + '  0  0')
    for b in model['bonds']:
        lines.append(f'{b["a"]:3d}{b["b"]:3d}{4 if b["order"] == 1.5 else int(b["order"]):3d}  0  0  0  0')
    charged = [(a['id'], a['formal_charge']) for a in model['atoms'] if a['formal_charge']]
    for start in range(0, len(charged), 8):
        block = charged[start:start+8]
        lines.append(f'M  CHG{len(block):3d}' + ''.join(f'{i:4d}{q:4d}' for i, q in block))
    return '\n'.join(lines+['M  END', '$$$$', ''])


def input_mol2(model):
    """Names survive MOL2 ingestion; original graph/coordinates are preserved."""
    lines = ['@<TRIPOS>MOLECULE', 'GOF', f'{len(model["atoms"])} {len(model["bonds"])} 1 0 0', 'SMALL', 'USER_CHARGES', '@<TRIPOS>ATOM']
    for a in model['atoms']:
        bonds = [b['order'] for b in model['bonds'] if a['id'] in (b['a'], b['b'])]
        typ = a['element']
        if typ == 'C': typ += '.ar' if 1.5 in bonds else '.2' if 2 in bonds else '.1' if 3 in bonds else '.3'
        elif typ == 'O': typ += '.2' if 2 in bonds else '.3'
        x, y, z = a['position']
        lines.append(f'{a["id"]} {atom_name(a)} {x:.12f} {y:.12f} {z:.12f} {typ} 1 GOF {a["formal_charge"]:.8f}')
    lines.append('@<TRIPOS>BOND')
    lines += [f'{i} {b["a"]} {b["b"]} {"ar" if b["order"] == 1.5 else int(b["order"])}' for i, b in enumerate(model['bonds'], 1)]
    return '\n'.join(lines+['@<TRIPOS>SUBSTRUCTURE', '1 GOF 1', ''])


def read_mol2(path):
    section = None; atoms = []; bonds = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            section = line[9:]; continue
        if not line.strip() or line.lstrip().startswith('#'): continue
        f = line.split()
        if section == 'ATOM':
            if len(f) < 9: raise ValueError('MOL2 atom lacks charges or fields')
            atoms.append({'id': int(f[0]), 'name': f[1], 'position': list(map(float, f[2:5])), 'type': f[5], 'charge_e': float(f[8])})
        elif section == 'BOND':
            if len(f) < 4: raise ValueError('Malformed MOL2 bond')
            bonds.append((tuple(sorted((int(f[1]), int(f[2])))), f[3].lower()))
    return atoms, bonds


def validate_geometry(model, positions):
    radii = {'H': .31, 'C': .76, 'O': .66}
    declared = {tuple(sorted((b['a'], b['b']))) for b in model['bonds']}
    for i, atom in enumerate(model['atoms']):
        for j in range(i+1, len(positions)):
            r = radii[atom['element']]+radii[model['atoms'][j]['element']]
            d = math.dist(positions[i], positions[j])
            if (i+1,j+1) in declared:
                if not .55*r < d < 1.30*r:
                    raise ValueError(f'AM1 geometry broken/short bond {i+1}-{j+1}')
            elif d < 1.15*r:
                raise ValueError(f'AM1 geometry new bond/clash {i+1}-{j+1}')


def graph_digest(model):
    graph={'atoms':[(a['id'],a['element'],a.get('formal_charge',0)) for a in model['atoms']],
           'bonds':sorted((min(b['a'],b['b']),max(b['a'],b['b']),b['order']) for b in model['bonds'])}
    return hashlib.sha256(json.dumps(graph,sort_keys=True).encode()).hexdigest()


def prepare_resonance_model(model_path, output_path):
    """LOCAL RDKit preparation only; worker subsequently needs only this metadata."""
    from rdkit import Chem, rdBase
    model=json.loads(Path(model_path).read_text()); validate_model(model)
    mol=Chem.MolFromMolBlock(sdf_text(model),removeHs=False)
    if mol is None: raise ValueError('RDKit cannot validate source graph')
    model.setdefault('metadata',{})['am1_bcc_resonance']={
        'graph_sha256':graph_digest(model),'rdkit_version':rdBase.rdkitVersion,
        'aromatic_bond_pairs':sorted([sorted([b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1]) for b in mol.GetBonds() if b.GetIsAromatic()])}
    write_json(output_path,model)
    return model


def validate_bond_representation(model, actual, output_atoms):
    expected={tuple(sorted((b['a'],b['b']))):b['order'] for b in model['bonds']}
    if set(actual)!=set(expected):raise ValueError('AM1-BCC bond connectivity changed')
    metadata=model.get('metadata',{}).get('am1_bcc_resonance',{})
    if metadata and metadata.get('graph_sha256')!=graph_digest(model):raise ValueError('Resonance metadata graph digest mismatch')
    aromatic={tuple(p) for p in metadata.get('aromatic_bond_pairs',[])}
    types={a['id']:a['type'] for a in output_atoms}
    coo=set()
    neighbors={a['id']:[] for a in model['atoms']}
    for a,b in expected:neighbors[a].append(b);neighbors[b].append(a)
    for carbon in model['atoms']:
        if carbon['element']!='C':continue
        oxy=[i for i in neighbors[carbon['id']] if model['atoms'][i-1]['element']=='O' and len(neighbors[i])==1]
        if (len(oxy)==2 and len(neighbors[carbon['id']])==3
                and sum(model['atoms'][i-1].get('formal_charge',0) for i in oxy)==-1
                and sorted(expected[tuple(sorted((carbon['id'],i)))] for i in oxy)==[1,2]
                and all(types.get(i) in ('o','o2') for i in oxy)
                and model['atoms'][oxy[0]-1].get('symmetry_class') is not None
                and model['atoms'][oxy[0]-1].get('symmetry_class')==model['atoms'][oxy[1]-1].get('symmetry_class')):
            coo.update(tuple(sorted((carbon['id'],i))) for i in oxy)
    changes=[]; cyclic={}
    for pair,old in expected.items():
        value=actual[pair]
        if value==str(int(old)) and old!=1.5:continue
        if old==1 and value=='am':continue
        if old==1.5 and value in ('ar','1','2'):continue
        if pair in aromatic and value in ('ar','1','2'):
            changes.append({'atom_ids':list(pair),'kind':'RDKit_aromatic_representation','input_order':old,'output_order':value});continue
        if pair in coo and value in ('1','2') and all(actual[p] in ('1','2') for p in coo):
            changes.append({'atom_ids':list(pair),'kind':'terminal_carboxylate_delocalized_MOL2_representation','input_order':old,'output_order':value});continue
        if (old in (1,2) and value in ('1','2') and int(value)==3-old
                and all(model['atoms'][i-1]['element']=='C' for i in pair)):
            cyclic[pair]=(old,int(value));continue
        raise ValueError(f'AM1-BCC fixed bond order changed: {pair}: {old}->{value}')
    # A numeric Kekule change must be a closed even alternating C-C cycle and
    # preserve each atom's integer bond-valence. Isolated sigma/FG changes fail.
    adjacency={}
    for pair,(old,new) in cyclic.items():
        for i in pair:adjacency.setdefault(i,[]).append((pair,old,new))
    if any(len(edges)!=2 or sum(e[2]-e[1] for e in edges)!=0 for edges in adjacency.values()):
        raise ValueError('Changed numeric bond orders do not form valence-preserving alternating cycles')
    cycles=[]; unseen=set(adjacency)
    while unseen:
        stack=[min(unseen)]; component=set()
        while stack:
            i=stack.pop()
            if i in component:continue
            component.add(i)
            stack.extend(j for edge,_,_ in adjacency[i] for j in edge if j not in component)
        unseen-=component
        if len(component)<6 or len(component)%2:raise ValueError('Invalid conjugated carbon cycle')
        cycles.append(sorted(component))
    return {'representation_changes':changes,'alternating_cycle_atom_ids':cycles,
            'validated_aromatic_metadata':bool(metadata),'fixed_sigma_and_functional_group_bonds_preserved':True}


def constrain_rounding(raw, formal_charge, *, source_decimals):
    if source_decimals!=3 or not raw or not all(math.isfinite(q) for q in raw):
        raise ValueError('Constrained-rounding policy requires proven finite three-decimal SQM populations')
    residual=sum(raw)-formal_charge; bound=len(raw)*5.1e-4
    if abs(residual)>bound:raise ValueError('Raw AM1-BCC charge residual exceeds propagated SQM rounding bound')
    delta=-residual/len(raw); corrected=[q+delta for q in raw]
    return corrected,{'policy':'constrained-rounding','raw_charge_sum_e':sum(raw),'raw_charge_error_e':residual,
        'corrected_charge_sum_e':sum(corrected),'per_atom_correction_e':[delta]*len(raw),
        'max_abs_correction_e':abs(delta),'source_mulliken_decimals':3,'propagated_rounding_bound_e':bound,
        'method':'uniform minimum-L2 correction subject to the integer total; explicit derived candidate, not unmodified standard AM1-BCC'}


def _bcc_parameter_file(programs):
    import os
    roots=[Path(os.environ[k]) for k in ('AMBERHOME','AMBERCLASSICHOME') if k in os.environ]
    roots.append(Path(programs['antechamber']['path']).parent.parent)
    for root in roots:
        path=root/'dat/antechamber/BCCPARM.DAT'
        if path.is_file():return path
    raise FileNotFoundError('Cannot locate active system BCCPARM.DAT for numerical BCC verification')


def _read_ac(path, model):
    atoms=[];bonds=[]
    for line in Path(path).read_text().splitlines():
        f=line.split()
        if line.startswith('ATOM'):
            atoms.append({'id':int(f[1]),'name':line[13:17].strip(),'q':float(f[-2]),'type':f[-1]})
        elif line.startswith('BOND'):bonds.append((int(f[2]),int(f[3]),int(f[4])))
    if any(not math.isfinite(a['q']) for a in atoms):raise ValueError('Nonfinite native AC charges')
    if ([a['id'] for a in atoms]!=[a['id'] for a in model['atoms']]
            or [a['name'] for a in atoms]!=[atom_name(a) for a in model['atoms']]):
        raise ValueError('Native AM1-BCC AC atom identity mismatch')
    return atoms,bonds


def validate_bcc_evidence(folder, model, output_charges):
    text=(folder/'sqm.out').read_text()
    headers=list(re.finditer(r'^\s*Atom\s+Element\s+Mulliken Charge\s*$',text,re.M))
    if not headers:raise ValueError('SQM Mulliken population table missing')
    rows=[]
    for line in text[headers[-1].end():].splitlines():
        f=line.split()
        if not f:continue
        if len(f)!=3 or not f[0].isdigit():break
        rows.append(f)
    if ([int(r[0]) for r in rows]!=[a['id'] for a in model['atoms']]
            or [r[1] for r in rows]!=[a['element'] for a in model['atoms']]):
        raise ValueError('SQM Mulliken atom mapping mismatch')
    qmul=[float(r[2]) for r in rows]
    decimals=[len(r[2].split('.')[1]) if '.' in r[2] else 0 for r in rows]
    if len(set(decimals))!=1:raise ValueError('Inconsistent SQM Mulliken precision')
    pre,_=_read_ac(folder/'ANTECHAMBER_AM1BCC_PRE.AC',model)
    post,bonds=_read_ac(folder/'ANTECHAMBER_AM1BCC.AC',model)
    if len(output_charges)!=len(post) or any(abs(a['q']-q)>1.1e-6 for a,q in zip(post,output_charges)):
        raise ValueError('MOL2 charges differ from native post-BCC AC output')
    if abs(sum(a['q'] for a in pre)-sum(qmul))>len(pre)*5.1e-7+2e-6:
        raise ValueError('Native pre-BCC charge sum differs from rounded SQM table')
    parameters={};corrections={}
    for line in (folder/'BCCPARM.DAT').read_text().splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='CORR':corrections[f[1]]=f[2]
        elif len(f)==5:parameters[(f[1],f[2],int(f[3]))]=float(f[4])
    if not all(math.isfinite(v) for v in parameters.values()):raise ValueError('Nonfinite BCC parameters')
    predicted=[0.]*len(post);covered=0;nonzero=0
    for a,b,kind in bonds:
        ta,tb=post[a-1]['type'],post[b-1]['type'];value=None
        for x,y in [(ta,tb),(corrections.get(ta,ta),tb),(ta,corrections.get(tb,tb)),(corrections.get(ta,ta),corrections.get(tb,tb))]:
            if (x,y,kind) in parameters:value=parameters[x,y,kind];break
            if (y,x,kind) in parameters:value=-parameters[y,x,kind];break
        if value is None:raise ValueError(f'BCC parameter missing for native types {ta}-{tb}-{kind}')
        predicted[a-1]+=value;predicted[b-1]-=value;covered+=1;nonzero+=bool(value)
    if any(abs((b['q']-a['q'])-d)>2.1e-6 for a,b,d in zip(pre,post,predicted)):
        raise ValueError('Native BCC correction does not match active parameter table')
    expected_edges={tuple(sorted((b['a'],b['b']))) for b in model['bonds']}
    if {tuple(sorted((a,b))) for a,b,_ in bonds}!=expected_edges or len(bonds)!=len(expected_edges):
        raise ValueError('Native BCC graph connectivity changed')
    return {'native_BCC_verified':True,'parameterized_bonds':covered,'nonzero_BCC_bonds':nonzero,
            'max_abs_BCC_delta_e':max(abs(d) for d in predicted),'source_mulliken_decimals':decimals[0],
            'raw_mulliken_sum_e':sum(qmul),'pre_BCC_sum_e':sum(a['q'] for a in pre),
            'post_BCC_sum_e':sum(a['q'] for a in post),'parameter_sha256':digest(folder/'BCCPARM.DAT')}


def validate_mol2(path, model, *, check_total=True):
    atoms, bonds = read_mol2(path)
    if [a['id'] for a in atoms] != [a['id'] for a in model['atoms']]:
        raise ValueError('AM1-BCC atom count/order/IDs differ')
    if len({a['name'] for a in atoms}) != len(atoms):
        raise ValueError('AM1-BCC output atom names must be unique; duplicate names found')
    for a, expected in zip(atoms, model['atoms']):
        if a['name'] != atom_name(expected):
            raise ValueError('AM1-BCC atom name/identity mapping changed')
        # GAFF2 campaign atom types begin with their element letter.
        if a['type'][0].upper() != expected['element']:
            raise ValueError('AM1-BCC atom element mapping changed')
        if not all(math.isfinite(v) for v in [*a['position'], a['charge_e']]):
            raise ValueError('AM1-BCC nonfinite charges/coordinates')
    validate_geometry(model, [a['position'] for a in atoms])
    expected = {tuple(sorted((b['a'], b['b']))): b['order'] for b in model['bonds']}
    actual = dict(bonds)
    if len(actual) != len(bonds) or set(actual) != set(expected):
        raise ValueError('AM1-BCC bond connectivity changed')
    representations = validate_bond_representation(model, actual, atoms)
    total = sum(a['charge_e'] for a in atoms)
    tolerance = len(atoms)*5.1e-7+2e-6  # Amber wmol2 writes charges with six decimals.
    if check_total and abs(total-model['formal_charge']) > tolerance:
        raise ValueError('AM1-BCC whole-model charge sum differs from integer target')
    return atoms, total, tolerance, representations


def validate_sqm(folder, model):
    text = (folder/'sqm.out').read_text()
    inp = (folder/'sqm.in').read_text()
    if not re.search(r"qm_theory\s*=\s*['\"]AM1['\"]", inp, re.I):
        raise ValueError('SQM input is not AM1')
    if re.search(r'failed to converge|SCF.*(?:not converged|failure)|SANDER BOMB|fatal|NaN|unable to converge|maximum.*(?:steps|cycles).*reached', text, re.I):
        raise ValueError('SQM failure/nonconvergence reported')
    if not re.search(r'Calculation Completed', text) or not re.search(r'geometry converged', text, re.I):
        raise ValueError('SQM AM1 optimization/completion evidence missing')
    final = text.rsplit('Final Structure', 1)
    if len(final) != 2: raise ValueError('SQM final AM1 geometry missing')
    rows = re.findall(r'^\s*QMMM:\s+(\d+)\s+(\d+)\s+([A-Za-z]+)\s+([-+\d.EeDd]+)\s+([-+\d.EeDd]+)\s+([-+\d.EeDd]+)', final[1], re.M)
    if ([int(r[0]) for r in rows] != list(range(1,len(model['atoms'])+1))
            or [r[2] for r in rows] != [a['element'] for a in model['atoms']]):
        raise ValueError('SQM final geometry atom mapping/order mismatch')
    positions = [[float(v.replace('D','E').replace('d','e')) for v in r[3:]] for r in rows]
    if not all(math.isfinite(v) for p in positions for v in p): raise ValueError('Nonfinite SQM final geometry')
    validate_geometry(model, positions)
    displacement = [math.dist(p,a['position']) for p,a in zip(positions,model['atoms'])]
    return {'sqm_optimization_converged': True, 'method': 'standard antechamber AM1-BCC defaults; sqm.in retained verbatim',
            'am1_final_positions_angstrom': positions, 'max_input_to_am1_displacement_angstrom': max(displacement),
            'am1_geometry_changed': max(displacement)>1e-4,
            'scoring_geometry': 'assign these charges by stable atom ID onto the common HF reference geometry; do not substitute AM1 coordinates'}


def _program(name):
    path = shutil.which(str(name))
    if path is None: raise FileNotFoundError(f'System AmberTools executable not found: {name}')
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': digest(path)}


def _attempt(root):
    root = Path(root).resolve(); root.mkdir(parents=True, exist_ok=True)
    i = 1
    while True:
        path = root/f'attempt_{i:03d}'
        try: path.mkdir()
        except FileExistsError: i += 1
        else: return path


def run_am1_bcc(model_path, output_dir, *, antechamber='antechamber', charge_policy='strict'):
    out = _attempt(output_dir)
    try:
        if charge_policy not in ('strict','constrained-rounding'):raise ValueError('Unknown charge policy')
        model_path = Path(model_path).resolve()
        (out/'model.json').write_bytes(model_path.read_bytes())
        model = json.loads((out/'model.json').read_text()); validate_model(model)
        (out/'input.sdf').write_text(sdf_text(model))
        (out/'input.mol2').write_text(input_mol2(model))
        write_json(out/'atom_mapping.json', [{'atom_id': a['id'], 'name': atom_name(a), 'element': a['element']} for a in model['atoms']])
        programs = {'antechamber': _program(antechamber), 'sqm': _program('sqm'), 'am1bcc': _program('am1bcc')}
        parameter_source=_bcc_parameter_file(programs)
        shutil.copy2(parameter_source,out/'BCCPARM.DAT')
        args = [programs['antechamber']['path'], '-i', 'input.mol2', '-fi', 'mol2', '-o', 'am1_bcc.mol2', '-fo', 'mol2', '-c', 'bcc', '-at', 'gaff2', '-nc', str(model['formal_charge']), '-rn', 'GOF', '-s', '2']
        started = datetime.now(timezone.utc).isoformat()
        write_json(out/'command.json', {'argv': args, 'started_utc': started,'bcc_parameter_source_path':str(parameter_source.resolve()),
                   'input_format_reason': 'SDF is archived; named MOL2 working input preserves element-prefixed four-character atom identity names that the SDF reader discards.'})
        with (out/'antechamber.log').open('x') as handle:
            process = subprocess.run(args, cwd=out, stdout=handle, stderr=subprocess.STDOUT, check=False)
        write_json(out/'execution.json', {'returncode': process.returncode, 'finished_utc': datetime.now(timezone.utc).isoformat()})
        if process.returncode: raise RuntimeError(f'antechamber failed ({process.returncode})')
        for name in ('sqm.in','sqm.out','am1_bcc.mol2'):
            if not (out/name).is_file() or not (out/name).stat().st_size: raise ValueError(f'Required AM1-BCC output missing: {name}')
        geometry = validate_sqm(out,model)
        atoms, total, tolerance, representation = validate_mol2(out/'am1_bcc.mol2',model,check_total=False)
        raw=[a['charge_e'] for a in atoms]
        evidence=validate_bcc_evidence(out,model,raw)
        write_json(out/'raw_am1_bcc.json',{'status':'raw_native_AM1_BCC_diagnostic','candidate':'raw_native_am1_bcc','model_id':model['model_id'],
            'coverage':{'complete':True},'final_parameter_validation':False,'atoms':[
            {'atom_id':a['id'],'name':a['name'],'element':m['element'],'role':m.get('role'),'charge_e':a['charge_e']} for a,m in zip(atoms,model['atoms'])],
            'charge_sum_e':total,'formal_charge':model['formal_charge'],'bcc_evidence':evidence})
        correction={'policy':'strict','raw_charge_sum_e':total,'corrected_charge_sum_e':total}
        if charge_policy=='constrained-rounding':
            fixed,correction=constrain_rounding(raw,model['formal_charge'],source_decimals=evidence['source_mulliken_decimals'])
            for a,q in zip(atoms,fixed):a['charge_e']=q
            total=sum(fixed)
        if abs(total-model['formal_charge'])>tolerance:
            raise ValueError('AM1-BCC whole-model charge sum differs from integer target; raw diagnostic retained')
        meta = model.get('metadata', {})
        scope = meta.get('scope', meta.get('scientific_status','input_reference_scope_unspecified'))
        records = [dict(atom_id=m['id'], name=a['name'], element=m['element'], role=m.get('role'),
                        charge_e=a['charge_e'], amber_type=a['type']) for a,m in zip(atoms,model['atoms'])]
        candidate = {'candidate': ('new_am1_bcc_charge_constrained_vsraw' if charge_policy=='constrained-rounding' else 'new_am1_bcc_gaff2'), 'model_id': model['model_id'], 'model_kind': ('derived_charge_constrained_AM1_BCC_control' if charge_policy=='constrained-rounding' else 'computed_AM1_BCC_charge_control'),
                     'source_reference_scope': scope, 'atoms': records, 'coverage': {'complete': True, 'covered_atoms': len(atoms), 'total_atoms': len(atoms), 'unsupported_atom_ids': []},
                     'formal_charge_e': model['formal_charge'], 'total_charge_e': total, 'charge_error_e': total-model['formal_charge'], 'charge_tolerance_e': tolerance,
                     'final_parameter_validation': False, 'geometry': geometry,'charge_policy':correction,'bcc_evidence':evidence}
        write_json(out/'candidate.json',candidate)
        write_json(out/'charges.json',{'model_id':model['model_id'],'formal_charge':model['formal_charge'],'atoms':records})
        with (out/'charges.csv').open('x', newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=['atom_id','name','element','role','charge_e','amber_type']); writer.writeheader(); writer.writerows(records)
        result = {'status':('DERIVED_AM1_BCC_CHARGE_CONSTRAINED_CONTROL_COMPLETE' if charge_policy=='constrained-rounding' else 'AM1_BCC_CHARGE_CONTROL_COMPLETE'),'candidate':candidate['candidate'],'model_id':model['model_id'],'attempt_dir':str(out),
                  'candidate_json':str(out/'candidate.json'),'charges_json':str(out/'charges.json'),'charges_csv':str(out/'charges.csv'),
                  'sdf_coordinate_precision_angstrom':0.0001,'working_mol2_coordinate_decimals':12,
                  'source_reference_scope':scope,'input_geometry_status':meta.get('geometry_status','unspecified'),
                  'final_parameter_validation':False,'geometry':geometry,'aromatic_bond_representation_changes':representation,
                  'charge_sum_e':total,'charge_tolerance_e':tolerance,'charge_policy':correction,'bcc_evidence':evidence,'provenance':{'model_sha256':digest(out/'model.json'),'code_sha256':digest(__file__),'executables':programs}}
        result['artifacts_sha256']={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()}
        write_json(out/'AM1_BCC_COMPLETE.json',result)
        return result
    except Exception as exc:
        write_json(out/'AM1_BCC_FAILED.json',{'exception':type(exc).__name__,'message':str(exc)})
        raise


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',required=True); p.add_argument('--output-dir',required=True)
    p.add_argument('--antechamber',default='antechamber'); p.add_argument('--slurm',action='store_true')
    p.add_argument('--partition',default='ciamician_cpu'); p.add_argument('--dependency')
    p.add_argument('--charge-policy',choices=('strict','constrained-rounding'),default='strict')
    args=p.parse_args(argv)
    if not args.slurm:
        result=run_am1_bcc(args.model,args.output_dir,antechamber=args.antechamber,charge_policy=args.charge_policy)
        print(json.dumps(result,indent=2)); return
    output=Path(args.output_dir).resolve(); output.parent.mkdir(parents=True,exist_ok=True)
    worker=output.with_name(output.name+'_am1bcc_worker.slurm')
    root=Path(__file__).resolve().parents[2]
    command=['python','-m','go_param.alternatives','--model',str(Path(args.model).resolve()),'--output-dir',str(output),'--antechamber',args.antechamber,'--charge-policy',args.charge_policy]
    with worker.open('x') as h:
        h.write('#!/bin/bash -l\nset -eo pipefail\nmodule load amber/25\nset -u\nexport OMP_NUM_THREADS=1\nexport PYTHONPATH='+shlex.quote(str(root/'src'))+'\n'+shlex.join(command)+'\n')
    cmd=['sbatch','--parsable','--partition',args.partition,'--nodes=1','--ntasks=1','--cpus-per-task=1','--mem=4G','--time=02:00:00','--job-name=GO-AM1BCC','--output',str(worker)+'.out','--error',str(worker)+'.err']
    if args.dependency:cmd+=['--dependency',args.dependency]
    receipt=subprocess.run(cmd+[str(worker)],capture_output=True,text=True,check=True)
    write_json(worker.with_suffix('.receipt.json'),{'status':'SUBMITTED','job_id':receipt.stdout.strip(),'command':cmd+[str(worker)]})
    print(receipt.stdout.strip())


if __name__=='__main__':main()
