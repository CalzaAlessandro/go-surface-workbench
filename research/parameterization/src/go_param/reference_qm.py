"""Explicit electronic-reference controls; no MD charge assignment or FF changes."""
from __future__ import annotations

import math
import re

Z={'H':1,'C':6,'O':8}
AU_DIPOLE_DEBYE=2.541746473
FIELD_UNIT_AU=1e-4
METHODS={
    'hf_631gstar':{'method':'HF','basis':'6-31G*','kind':'hf'},
    'wb97xd_pople':{'method':'wB97XD','basis':'6-311++G(d,p)','kind':'dft'},
    'wb97xd_avtz':{'method':'wB97XD','basis':'aug-cc-pVTZ','kind':'dft'},
    'pbe0_avtz':{'method':'PBE1PBE','basis':'aug-cc-pVTZ','kind':'dft'},
    'pw6b95_avtz':{'method':'PW6B95','basis':'aug-cc-pVTZ','kind':'dft'},
    'pw6b95_avdz':{'method':'PW6B95','basis':'aug-cc-pVDZ','kind':'dft'},
    'hf_avdz':{'method':'HF','basis':'aug-cc-pVDZ','kind':'hf'},
    'hf_avtz':{'method':'HF','basis':'aug-cc-pVTZ','kind':'hf'},
    'mp2_avdz':{'method':'MP2','basis':'aug-cc-pVDZ','kind':'mp2'},
    'mp2_avtz':{'method':'MP2','basis':'aug-cc-pVTZ','kind':'mp2'},
    'ccsd_avdz':{'method':'CCSD','basis':'aug-cc-pVDZ','kind':'ccsd'},
    'ccsdt_avdz':{'method':'CCSD(T)','basis':'aug-cc-pVDZ','kind':'triples_energy'},
}


def make_sp_input(model,method_key,*,cpus,mem_gb,field=None,include_stability=True):
    if method_key not in METHODS:raise ValueError('unsupported reference method')
    spec=METHODS[method_key];atoms=model['atoms'];q=model['formal_charge'];m=model['multiplicity']
    if not atoms or type(q) is not int or type(m) is not int or m<1:raise ValueError('explicit atoms, Q and multiplicity required')
    electrons=sum(Z[a['element']] for a in atoms)-q
    if electrons<1 or electrons<m-1 or (electrons-m+1)%2:raise ValueError('electron/multiplicity parity mismatch')
    if type(cpus) is not int or cpus<1 or type(mem_gb) is not int or mem_gb<1:raise ValueError('invalid resources')
    xyz=[]
    for atom in atoms:
        p=atom['position']
        if len(p)!=3 or not all(math.isfinite(v) for v in p):raise ValueError('invalid fixed coordinates')
        xyz.append(atom['element']+' '+' '.join(f'{v:.12f}' for v in p))
    energy_only=spec['kind']=='triples_energy' or field is not None
    extra=' Integral=SuperFineGrid' if spec['kind']=='dft' else ''
    if not energy_only:extra+=' Density=Current Pop=MK IOp(6/33=2,6/41=10,6/42=17)'
    if field is not None:
        axis,n=field
        if axis not in 'XYZ' or len(axis)!=1 or type(n) is not int or n==0 or abs(n)>10:raise ValueError('invalid field specification')
        extra+=f' Field={axis}{n:+d}'
    route=f'{spec["method"]}/{spec["basis"]}'
    header=[f'%NProcShared={cpus}',f'%Mem={mem_gb}GB','%Chk=reference.chk']
    lines=header+[f'#P {route} NoSymm SCF=(VeryTight,XQC,MaxCycle=512){extra}',
                  '',f'{model["model_id"]}: fixed-geometry reference {method_key}', '',f'{q} {m}',*xyz,'']
    if include_stability and field is None:
        reference=spec['method'] if spec['kind']=='dft' else 'HF'
        integration=' Integral=SuperFineGrid' if spec['kind']=='dft' else ''
        lines+=['--Link1--',*header,f'#P {reference}/{spec["basis"]} NoSymm Geom=AllCheck Guess=Read Stable SCF=(VeryTight,XQC,MaxCycle=512){integration}','','']
    else:lines+=['']
    return '\n'.join(lines)


def finite_field_dipole(energies,*,operator_sign):
    """Dipole from h=1e-4 and 2h energies with a runtime-verified field sign."""
    if operator_sign not in (-1,1):raise ValueError('field operator sign must be explicitly verified as -1 or +1')
    dipoles=[]
    for n in (1,2):
        mu=[]
        for axis in 'XYZ':
            try:plus=float(energies[axis][n]);minus=float(energies[axis][-n])
            except (KeyError,TypeError,ValueError) as exc:raise ValueError('missing finite-field energies') from exc
            if not math.isfinite(plus) or not math.isfinite(minus):raise ValueError('nonfinite field energy')
            mu.append(operator_sign*(plus-minus)/(2*n*FIELD_UNIT_AU))
        dipoles.append(mu)
    extrapolated=[(4*a-b)/3 for a,b in zip(*dipoles)]
    return {'dipole_au':extrapolated,'dipole_debye':[v*AU_DIPOLE_DEBYE for v in extrapolated],
            'dipole_h_debye':[v*AU_DIPOLE_DEBYE for v in dipoles[0]],
            'dipole_2h_debye':[v*AU_DIPOLE_DEBYE for v in dipoles[1]],
            'max_step_difference_debye':max(abs(a-b)*AU_DIPOLE_DEBYE for a,b in zip(*dipoles)),
            'field_step_au':FIELD_UNIT_AU,'field_operator_sign':operator_sign,
            'method':'central differences, two steps, Richardson extrapolation; fixed coordinates and origin'}


def property_job_text(text):
    """Only the first terminated property job; a later HF link is not CC data."""
    end=re.search(r'Normal termination of Gaussian[^\n]*',text)
    return text[:end.end()]+'\n' if end else text


def inspect_sp_output(text,model,method_key,*,field=None,require_stability=True):
    """Gate fixed-geometry properties using actual G16 density evidence.

    Acceptance is for this comparison, never a ground-state or MD validation.
    Field jobs require a separately accepted, zero-field determinant check.
    """
    from .qm import parse_gaussian_log
    if method_key not in METHODS:raise ValueError('unsupported reference method')
    spec=METHODS[method_key];first=property_job_text(text)
    evidence=parse_gaussian_log(first);whole=parse_gaussian_log(text)
    issues=[]
    number=lambda s:float(s.replace('D','E').replace('d','e'))
    if evidence['normal_terminations']!=1 or whole['error_lines'] or not re.search(r'Normal termination of Gaussian[^\n]*\s*\Z',text):
        issues.append('incomplete job or Gaussian error')
    route=''.join(evidence['routes'][0].lower().split()) if evidence['routes'] else ''
    if (spec['method']+'/'+spec['basis']).lower() not in route:
        issues.append('requested method/basis not evidenced')
    if 'nosymm' not in route or re.search(r'\b(?:opt|freq)\b',route):
        issues.append('fixed geometry/orientation protocol not evidenced')
    qm=evidence['charge_multiplicity_records']
    if not qm or any(v!=(model['formal_charge'],model['multiplicity']) for v in qm):
        issues.append('charge or multiplicity mismatch')
    geom=evidence['final_orientation']
    if not geom or geom['parse_error'] or len(geom['atoms'])!=len(model['atoms']):
        issues.append('atomic identity/geometry missing')
    else:
        for row,a in zip(geom['atoms'],model['atoms']):
            if row['atomic_number']!=Z[a['element']] or math.dist(row['position'],a['position'])>2e-6:
                issues.append('fixed atomic identity or coordinates mismatch');break
    if require_stability and (whole['normal_terminations']!=2 or not whole['wavefunction_stable'] or not any('stable' in r.lower() for r in whole['routes'][1:])):
        issues.append('separate determinant stability evidence missing or unstable')
    if field is not None:
        axis,n=field
        if f'field={axis}{n:+d}'.lower() not in route:issues.append('field specification mismatch')
    patterns={'dft':r'SCF Done:\s+E\([^)]+\)\s*=\s*(\S+)',
              'hf':r'SCF Done:\s+E\([^)]+\)\s*=\s*(\S+)',
              'mp2':r'EUMP2\s*=\s*(\S+)',
              'ccsd':r'E\(CORR\)\s*=\s*(\S+)',
              'triples_energy':r'^\s*CCSD\(T\)\s*=\s*(\S+)'}
    matches=re.findall(patterns[spec['kind']],first,re.M|re.I)
    try:energy=number(matches[-1]) if matches else None
    except ValueError:energy=None
    if energy is None or not math.isfinite(energy):issues.append('finite method-specific energy missing')
    if spec['kind'] in ('ccsd','triples_energy') and 'Wavefunction amplitudes converged.' not in first:
        issues.append('correlated amplitude convergence missing')
    energy_only=field is not None or spec['kind']=='triples_energy'
    density=None;dipole=None
    if not energy_only:
        expected={'dft':'SCF','hf':'SCF','mp2':'MP2','ccsd':'CC'}[spec['kind']]
        populations=re.findall(r'Population analysis using the (\w+) Density',first,re.I)
        properties=re.findall(r'Electrostatic Properties Using The (\w+) Density',first,re.I)
        density=populations[-1].upper() if populations else None
        if 'density=current' not in route or density!=expected or not properties or properties[-1].upper()!=expected:
            issues.append('selected property density missing or inconsistent with requested method')
        dipoles=re.findall(r'Dipole moment \(field-independent basis, Debye\):\s*X=\s*(\S+)\s+Y=\s*(\S+)\s+Z=\s*(\S+)',first)
        try:dipole=[number(v) for v in dipoles[-1]] if dipoles else None
        except ValueError:dipole=None
        if dipole is None or not all(math.isfinite(v) for v in dipole):issues.append('finite property dipole missing')
    diagnostic=re.findall(r'T1 Diagnostic\s*=\s*(\S+)',first,re.I)
    return {'accepted':not issues,'issues':issues,'method_key':method_key,
            'energy_hartree':energy,'density_label':density,'dipole_debye':dipole,
            'origin_angstrom':[0.,0.,0.],'ionic_dipole_origin_dependent':model['formal_charge']!=0,
            'determinant_stability_checked':require_stability and whole['wavefunction_stable'],
            'T1_diagnostic':number(diagnostic[-1]) if diagnostic else None,
            'field':list(field) if field else None,
            'interpretation':'comparison gate only; no proof of global ground state, multireference safety or MD transferability'}
