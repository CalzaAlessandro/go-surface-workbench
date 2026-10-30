"""Synthetic evidence tests; no Gaussian calculation is simulated as science."""
import copy
import json

import pytest

from test_qm import model, synthetic_log


def fixed_log(record):
    numbers = {'C': 6, 'O': 8, 'H': 1}
    lines = ['SYNTHETIC HF FIXED-GEOMETRY PARSER TEST',
             ' #P HF/6-31G* Stable SCF=(Tight,XQC) Symmetry=None',
             ' ---------------------------------------------------------------------',
             f" Charge = {record['formal_charge']} Multiplicity = {record['multiplicity']}",
             ' SCF Done:  E(RHF) = -765.123456789 A.U. after 12 cycles',
             ' Input orientation:',
             ' ---------------------------------------------------------------------',
             ' Center Atomic Atomic Coordinates (Angstroms)',
             ' Number Number Type X Y Z',
             ' ---------------------------------------------------------------------']
    for a in record['atoms']:
        lines.append(f" {a['id']} {numbers[a['element']]} 0 "+' '.join(f'{x:.10f}' for x in a['position']))
    lines += [' ---------------------------------------------------------------------',
              ' The wavefunction is stable under the perturbations considered.',
              ' Normal termination of Gaussian 16']
    return '\n'.join(lines)+'\n'


def test_probe_preserves_dft_geometry_and_has_separate_hf_target(tmp_path):
    from go_param import hf_fixed
    record=model(); original=copy.deepcopy(record)
    log=synthetic_log(record).replace('HF/6-31G*','wB97XD/6-311++G(d,p)').replace('E(RHF)','E(RwB97XD)')
    result=hf_fixed.prepare_probe(record,log,tmp_path/'probe',cpus=4,mem_gb=6)
    assert record==original
    prepared=json.loads((tmp_path/'probe/geometry_model.json').read_text())
    assert [a['id'] for a in prepared['atoms']]==[a['id'] for a in record['atoms']]
    assert prepared['metadata']['qm_method']=='wB97XD/6-311++G(d,p)'
    assert result['geometry_method']=='wB97XD/6-311++G(d,p)'
    assert result['probe_method']=='HF/6-31G*'
    assert result['hf_reference_accepted'] is False
    assert result['final_parameter_validation'] is False


def test_probe_rejects_unaccepted_dft_parent_before_creating_output(tmp_path):
    from go_param import hf_fixed
    with pytest.raises(ValueError,match='DFT parent'):
        hf_fixed.prepare_probe(model(),synthetic_log(model()),tmp_path/'rejected')
    assert not (tmp_path/'rejected').exists()


def test_valid_fixed_hf_reference_passes_without_claiming_hf_optimization():
    from go_param import hf_fixed
    report=hf_fixed.validate_probe(model(),fixed_log(model()))
    assert report['accepted'] and report['issues']==[]
    assert report['hf_wavefunction_stable'] is True
    assert report['hf_geometry_optimization_validated'] is False
    assert report['final_parameter_validation'] is False


@pytest.mark.parametrize('old,new,issue',[
    ('The wavefunction is stable under the perturbations considered.','The wavefunction has an RHF -> UHF instability.','hf_stability_not_confirmed'),
    ('E(RHF)','E(RwB97XD)','hf_energy_method_mismatch'),
    ('-765.123456789','NaN','nonfinite_hf_energy'),
    ('HF/6-31G*','wB97XD/6-311++G(d,p)','hf_route_mismatch'),
    (' Stable ',' Stable=Opt ','hf_route_mismatch'),
    ('Normal termination of Gaussian 16','Error termination','hf_incomplete_or_error'),
    ('Multiplicity = 1','Multiplicity = 3','charge_or_multiplicity_mismatch'),
])
def test_fixed_hf_reference_rejects_unsafe_evidence(old,new,issue):
    from go_param import hf_fixed
    record=model()
    report=hf_fixed.validate_probe(record,fixed_log(record).replace(old,new))
    assert not report['accepted'] and issue in report['issues']


def test_probe_rejects_movement_and_element_mapping_changes():
    from go_param import hf_fixed
    record=model(); moved=copy.deepcopy(record)
    moved['atoms'][0]['position'][0]+=0.1
    report=hf_fixed.validate_probe(record,fixed_log(moved))
    assert not report['accepted'] and 'fixed_geometry_changed' in report['issues']
    altered=copy.deepcopy(record); altered['atoms'][0]['element']='O'
    report=hf_fixed.validate_probe(record,fixed_log(altered))
    assert not report['accepted'] and 'atom_identity_mismatch' in report['issues']
