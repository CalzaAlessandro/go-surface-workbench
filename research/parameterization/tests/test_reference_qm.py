"""Meaningful gates for reference calculations, separate from charge fitting."""
import pytest
from pathlib import Path


def water():
    return {'model_id':'runtime_water','formal_charge':0,'multiplicity':1,
            'atoms':[{'id':1,'element':'O','position':[0.,0.,0.]},
                     {'id':2,'element':'H','position':[.758,0.,.504]},
                     {'id':3,'element':'H','position':[-.758,0.,.504]}],
            'bonds':[{'a':1,'b':2,'order':1},{'a':1,'b':3,'order':1}]}


def test_fixed_geometry_explicit_density_and_reference_stability():
    from go_param.reference_qm import make_sp_input
    s=make_sp_input(water(),'ccsd_avdz',cpus=4,mem_gb=6)
    assert 'CCSD/aug-cc-pVDZ' in s and 'Density=Current' in s and 'Pop=MK' in s
    assert 'HF/aug-cc-pVDZ' in s and 'Stable' in s
    assert 'Opt' not in s and 'Freq' not in s and '0.758000000000' in s
    assert 'Integral=SuperFineGrid' in make_sp_input(water(),'pw6b95_avtz',cpus=4,mem_gb=6)
    bad=water();bad['formal_charge']=-1
    with pytest.raises(ValueError,match='parity'):make_sp_input(bad,'pbe0_avtz',cpus=4,mem_gb=6)
    with pytest.raises(ValueError,match='method'):make_sp_input(water(),'invented',cpus=4,mem_gb=6)


def test_field_units_and_no_invented_triples_density():
    from go_param.reference_qm import make_sp_input,finite_field_dipole
    s=make_sp_input(water(),'ccsdt_avdz',cpus=4,mem_gb=6,field=('Z',-2))
    assert 'CCSD(T)/aug-cc-pVDZ' in s and 'Field=Z-2' in s
    assert 'Pop=MK' not in s and 'Density=' not in s and 'Opt' not in s
    mu=[.2,-.5,.9];h=1e-4
    values={axis:{1:-10-m*h, -1:-10+m*h, 2:-10-m*2*h, -2:-10+m*2*h} for axis,m in zip('XYZ',mu)}
    result=finite_field_dipole(values,operator_sign=-1)
    assert result['dipole_au']==pytest.approx(mu,abs=1e-10)
    assert result['max_step_difference_debye']<1e-9
    with pytest.raises(ValueError,match='sign'):finite_field_dipole(values,operator_sign=0)
    with pytest.raises(ValueError,match='field'):make_sp_input(water(),'ccsdt_avdz',cpus=4,mem_gb=6,field=('Z',0))


def test_actual_gaussian_density_gates_and_first_link_dipole():
    from go_param.reference_qm import inspect_sp_output
    p=Path(__file__).parent/'fixtures/reference_qm'
    for case,key,mu,density in [('pw6b95_avtz','pw6b95_avtz',None,'SCF'),
                               ('mp2_avdz','mp2_avdz',1.7207,'MP2'),
                               ('ccsd_avdz','ccsd_avdz',1.7242,'CC')]:
        t=(p/(case+'.log')).read_text()
        r=inspect_sp_output(t,water(),key)
        assert r['accepted'],r
        assert r['density_label']==density
        if mu:assert r['dipole_debye'][2]==pytest.approx(mu,abs=1e-6)
        assert not inspect_sp_output(t.replace('is stable under the perturbations considered.','has an RHF -> UHF instability.'),water(),key)['accepted']
        wrong=water();wrong['atoms'][1]['position'][0]+=.02
        assert not inspect_sp_output(t,wrong,key)['accepted']
        assert not inspect_sp_output(t.split('Normal termination')[0],water(),key)['accepted']
    negative=inspect_sp_output((p/'ccsd_deliberate_scf_density.log').read_text(),water(),'ccsd_avdz')
    assert not negative['accepted'] and any('density' in i for i in negative['issues'])
    energy=inspect_sp_output((p/'ccsdt_avdz_Z+1.log').read_text(),water(),'ccsdt_avdz',field=('Z',1),require_stability=False)
    assert energy['accepted'] and energy['dipole_debye'] is None
    assert energy['energy_hartree']==pytest.approx(-76.266203501,abs=1e-9)
