"""System-ORCA whole-sheet DFT controls; never label populations as RESP."""
from __future__ import annotations

import math
import re

Z={'H':1,'C':6,'O':8}
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


def make_orca_input(model,*,multiplicity,processes,maxcore_mb):
    atoms=model['atoms'];charge=model['formal_charge']
    if not atoms or type(charge) is not int or type(multiplicity) is not int:
        raise ValueError('explicit atoms, integer charge and multiplicity required')
    electrons=sum(Z[a['element']] for a in atoms)-charge
    if multiplicity<1 or electrons<multiplicity-1 or (electrons-multiplicity+1)%2:
        raise ValueError('electron count and multiplicity have incompatible parity')
    if processes<1 or maxcore_mb<1:raise ValueError('invalid process/memory request')
    xyz=[]
    for atom in atoms:
        p=atom['position']
        if len(p)!=3 or not all(math.isfinite(v) for v in p):raise ValueError('invalid coordinates')
        xyz.append(atom['element']+' '+' '.join(f'{v:.10f}' for v in p))
    guess=' GuessMix' if multiplicity==1 else ''
    return (f'! UKS PBE ma-def2-SVP RI def2/J TightSCF SlowConv NoAutoStart{guess} CPCM(Water) Hirshfeld\n'
            f'%pal nprocs {processes} end\n%maxcore {maxcore_mb}\n'
            '%scf\n MaxIter 500\n AutoTRAH true\nend\n'
            f'* xyz {charge} {multiplicity}\n'+'\n'.join(xyz)+'\n*\n')


def inspect_orca(text,elements,*,formal_charge,multiplicity):
    """Gate the single point, preserving spin/population warnings separately."""
    energy=re.findall(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+')',text)
    s2=re.findall(r'Expectation value of <S\*\*2>\s*:\s*('+NUMBER+')',text)
    nel=re.findall(r'Number of electrons\s*=\s*(\d+)',text,re.I)
    mult=re.findall(r'Multiplicity\s*=\s*(\d+)',text,re.I)
    expected_electrons=sum(Z[e] for e in elements)-formal_charge
    rows=[]
    if 'HIRSHFELD ANALYSIS' in text:
        section=text.rsplit('HIRSHFELD ANALYSIS',1)[1].split('TOTAL',1)[0]
        for i,s,q,spin in re.findall(r'^\s*(\d+)\s+([A-Z][a-z]?)\s+('+NUMBER+r')\s+('+NUMBER+r')\s*$',section,re.M):
            rows.append({'atom_index':int(i),'element':s,'charge_e':float(q),'spin_population':float(spin)})
    qsum=sum(r['charge_e'] for r in rows)
    tolerance=max(1e-3,len(elements)*2e-5)
    checks={'normal_termination':'ORCA TERMINATED NORMALLY' in text,
            'SCF_converged':'SCF CONVERGED AFTER' in text and 'SCF NOT CONVERGED' not in text,
            'finite_energy':bool(energy) and math.isfinite(float(energy[-1])),
            'electron_count_matches':bool(nel) and int(nel[-1])==expected_electrons,
            'trial_multiplicity_matches':bool(mult) and int(mult[-1])==multiplicity,
            'population_atom_mapping_matches':len(rows)==len(elements)
                and [r['atom_index'] for r in rows]==list(range(len(elements)))
                and [r['element'] for r in rows]==elements,
            'population_total_within_integration_tolerance':bool(rows) and abs(qsum-formal_charge)<=tolerance,
            'finite_populations':all(math.isfinite(r['charge_e']) and math.isfinite(r['spin_population']) for r in rows)}
    spin=(multiplicity-1)/2;observed_s2=float(s2[-1]) if s2 else None
    return {'accepted_single_point_diagnostic':all(checks.values()),'checks':checks,
            'energy_hartree':float(energy[-1]) if energy else None,
            'expected_electrons':expected_electrons,'trial_multiplicity':multiplicity,
            'S_squared':observed_s2,'ideal_pure_spin_S_squared':spin*(spin+1),
            'spin_contamination_flag':observed_s2 is None or abs(observed_s2-spin*(spin+1))>.2,
            'hirshfeld_charge_sum_e':qsum,'integration_tolerance_e':tolerance,
            'hirshfeld_atoms':rows,'charges_are_RESP':False,
            'wavefunction_stability_validated':False,'ground_spin_state_established':False,
            'geometry_optimized_at_QM_level':False,'production_parameter_validation':False}
