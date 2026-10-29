"""Synthetic Gaussian text fixtures exercise gates; these are not QM results."""
import json
import math
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]


def model():
    return json.loads((ROOT/'models/pyrene_epoxide/model.json').read_text())


def synthetic_log(record, *, frequencies=None, stable=True):
    frequencies = frequencies if frequencies is not None else [100.] * (3*len(record['atoms'])-6)
    elements = {'H': 1, 'C': 6, 'O': 8}
    lines = ['SYNTHETIC PARSER FIXTURE - NOT A QM RESULT',
             ' #P HF/6-31G* Opt=(Tight,CalcFC,MaxCycles=200) Freq Pop=MK SCF=(Tight,XQC) NoSymm',
             ' ---------------------------------------------------------------------',
             f" Charge = {record['formal_charge']} Multiplicity = {record['multiplicity']}",
             ' SCF Done:  E(RHF) =  -765.123456789 A.U. after 12 cycles',
             ' Optimization completed.', ' Harmonic frequencies (cm**-1)']
    for i in range(0, len(frequencies), 3):
        lines.append(' Frequencies -- ' + ' '.join(str(x) for x in frequencies[i:i+3]))
    lines += [' Normal termination of Gaussian 16', ' Input orientation:',
              ' ---------------------------------------------------------------------',
              ' Center     Atomic      Atomic             Coordinates (Angstroms)',
              ' Number     Number       Type             X           Y           Z',
              ' ---------------------------------------------------------------------']
    for i, a in enumerate(record['atoms'], 1):
        x,y,z = a['position']
        lines.append(f' {i} {elements[a["element"]]} 0 {x:.10f} {y:.10f} {z:.10f}')
    lines += [' ---------------------------------------------------------------------',
              ' #P HF/6-31G* Geom=AllCheck Guess=Read Stable NoSymm SCF=(Tight,XQC)',
              ' ---------------------------------------------------------------------',
              f" Charge = {record['formal_charge']} Multiplicity = {record['multiplicity']}",
              ' SCF Done:  E(RHF) =  -765.123456790 A.U. after 1 cycles']
    if stable:
        lines.append(' The wavefunction is stable under the perturbations considered.')
    lines.append(' Normal termination of Gaussian 16')
    return '\n'.join(lines)+'\n'


def test_gaussian_input_contains_explicit_stability_check_separate_from_opt():
    from go_param.qm import gaussian_input
    text = gaussian_input(model(), cpus=4, mem_gb=3)
    assert '%NProcShared=4' in text
    assert '%Mem=3GB' in text
    assert 'Opt=(Tight,CalcFC,MaxCycles=200) Freq Pop=MK' in text
    assert 'IOp(6/33=2,6/41=10,6/42=17)' in text
    assert '--Link1--' in text
    assert 'Geom=AllCheck Guess=Read Stable Symmetry=None' in text
    assert text.count('Symmetry=None') == 2 and 'NoSymm' not in text
    assert 'Stable=Opt' not in text


def test_complete_synthetic_evidence_passes_and_preserves_atom_identity():
    from go_param.qm import validate_gaussian_log
    record = model()
    result = validate_gaussian_log(record, synthetic_log(record))
    assert result['accepted']
    assert result['issues'] == []
    assert len(result['frequencies_cm1']) == 75
    assert result['final_hf_energy_hartree'] == -765.123456790
    assert [a['id'] for a in result['optimized_atoms']] == [a['id'] for a in record['atoms']]
    assert result['topology']['passed']


@pytest.mark.parametrize('old,new,issue', [
    (' Optimization completed.', '', 'optimization_not_completed'),
    (' The wavefunction is stable under the perturbations considered.', '', 'stability_not_confirmed'),
    (' Normal termination of Gaussian 16', ' Error termination', 'normal_terminations_missing'),
    (' Frequencies -- 100.0 100.0 100.0', ' Frequencies -- -50.0 100.0 100.0', 'significant_imaginary_frequencies'),
    (' Charge = 0 Multiplicity = 1', ' Charge = 1 Multiplicity = 1', 'charge_or_multiplicity_mismatch'),
])
def test_missing_or_negative_evidence_is_rejected(old, new, issue):
    from go_param.qm import validate_gaussian_log
    record = model()
    result = validate_gaussian_log(record, synthetic_log(record).replace(old, new))
    assert not result['accepted']
    assert issue in result['issues']


def test_small_negative_mode_requires_manual_review_and_nonfinite_is_rejected():
    from go_param.qm import validate_gaussian_log
    record = model()
    f = [100.] * 75
    f[0] = -10.
    result = validate_gaussian_log(record, synthetic_log(record, frequencies=f))
    assert not result['accepted']
    assert result['manual_review_required']
    f[0] = float('nan')
    result = validate_gaussian_log(record, synthetic_log(record, frequencies=f))
    assert not result['accepted']
    assert 'nonfinite_frequencies' in result['issues']


def test_identity_and_topology_changes_fail_and_absent_evidence_never_passes():
    from go_param.qm import validate_gaussian_log
    record = model()
    bad = json.loads(json.dumps(record))
    bad['atoms'][0]['position'] = [100., 0., 0.]
    result = validate_gaussian_log(record, synthetic_log(bad))
    assert not result['accepted']
    assert 'topology_distance_check_failed' in result['issues']
    assert result['topology']['lost_bonds']
    bad = json.loads(json.dumps(record))
    bad['atoms'][0]['element'] = 'O'
    assert 'atom_identity_mismatch' in validate_gaussian_log(record, synthetic_log(bad))['issues']
    assert not validate_gaussian_log(record, ' Normal termination of Gaussian 16\n')['accepted']


def test_prepare_is_gated_and_three_rotations_preserve_distances(tmp_path):
    from go_param.qm import prepare_orientations
    record = model()
    with pytest.raises(ValueError, match='not accepted'):
        prepare_orientations(record, 'incomplete log', tmp_path/'blocked')
    assert not (tmp_path/'blocked').exists()
    manifest = prepare_orientations(record, synthetic_log(record), tmp_path/'prepared')
    assert len(manifest['orientations']) == 3
    original_distance = math.dist(record['atoms'][0]['position'], record['atoms'][1]['position'])
    for entry in manifest['orientations']:
        folder = tmp_path/'prepared'/entry['directory']
        rotated = json.loads((folder/'model.json').read_text())
        assert math.dist(rotated['atoms'][0]['position'], rotated['atoms'][1]['position']) == pytest.approx(original_distance, abs=1e-8)
        assert len(entry['input_sha256']) == 64
        text = (folder/'mk.gjf').read_text()
        assert 'Pop=MK' in text and 'Symmetry=None' in text
        assert 'NoSymm' not in text
        assert 'Opt=' not in text and 'Freq' not in text
    assert (tmp_path/'prepared/model.json').exists()
    with pytest.raises(FileExistsError):
        prepare_orientations(record, synthetic_log(record), tmp_path/'prepared')


def test_controller_inspection_and_absolute_prepared_paths(tmp_path):
    from go_param.qm import inspect_reference, prepare_orientations
    import hashlib
    record = model()
    source = tmp_path/'synthetic.log'
    source.write_bytes(synthetic_log(record).replace('\n', '\r\n').encode())
    report = inspect_reference(record, source)
    assert report['accepted']
    assert report['optimized_model']['metadata']['source_log_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    prepared = prepare_orientations(record, source, tmp_path/'prepared')
    assert Path(prepared['optimized_model_path']).is_file()
    for entry in prepared['orientations']:
        assert Path(entry['input_path']).is_file()
        assert Path(entry['model_path']).is_file()


def test_later_instability_cannot_reuse_an_earlier_stable_phrase():
    from go_param.qm import validate_gaussian_log
    record = model()
    text = synthetic_log(record).replace(' The wavefunction is stable under the perturbations considered.',
        ' The wavefunction is stable under the perturbations considered.\n The wavefunction has an internal instability.')
    result = validate_gaussian_log(record, text)
    assert not result['accepted']
    assert 'stability_not_confirmed' in result['issues']


def test_explicit_rhf_diagnostic_retains_rejection_and_never_bypasses_other_issues(tmp_path):
    from go_param.qm import prepare_orientations
    record = model()
    log = synthetic_log(record).replace('The wavefunction is stable under the perturbations considered.',
        'The wavefunction has an RHF -> UHF instability.')
    with pytest.raises(ValueError, match='not accepted'):
        prepare_orientations(record, log, tmp_path/'strict')
    result = prepare_orientations(record, log, tmp_path/'diagnostic', allow_rhf_instability_diagnostic=True)
    assert result['scope'] == 'rhf_instability_diagnostic_control'
    validation = json.loads((tmp_path/'diagnostic/validation.json').read_text())
    assert validation['accepted'] is False
    assert validation['issues'] == ['stability_not_confirmed']
    prepared = json.loads((tmp_path/'diagnostic/model.json').read_text())
    assert prepared['metadata']['wavefunction_stability_validated'] is False
    with pytest.raises(ValueError, match='not accepted'):
        prepare_orientations(record, log.replace('Optimization completed.', ''), tmp_path/'bad',
            allow_rhf_instability_diagnostic=True)


def test_route_provenance_cannot_combine_dft_optimization_with_hf_stability():
    from go_param.qm import inspect_reference
    record = model()
    text = synthetic_log(record).replace('HF/6-31G* Opt=', 'B3LYP/6-31G* Opt=')
    result = inspect_reference(record, text)
    assert not result['accepted']
    assert 'route_provenance_mismatch' in result['issues']


@pytest.mark.parametrize('keyword', ['NoSymm', 'Symmetry=None', 'Symmetry = None'])
def test_reference_gate_accepts_both_explicit_symmetry_disabling_forms(keyword):
    from go_param.qm import inspect_reference
    record = model()
    assert inspect_reference(record, synthetic_log(record).replace('NoSymm', keyword))['accepted']


def test_reference_gate_rejects_symmetry_enabled_routes():
    from go_param.qm import inspect_reference
    record = model()
    report = inspect_reference(record, synthetic_log(record).replace('NoSymm', 'Symmetry=On'))
    assert not report['accepted']
    assert 'route_provenance_mismatch' in report['issues']


# Exact route excerpts from the real pilot log, without geometry, energies or
# scientific results. The full log is not embedded in the test suite.
PILOT_ROUTE_EXCERPT = '''
 ----------------------------------------------------------------------
 #P HF/6-31G* Opt=(Tight,CalcFC,MaxCycles=200) Freq Pop=MK SCF=(Tight,X
 QC) NoSymm IOp(6/33=2,6/41=10,6/42=17)
 ----------------------------------------------------------------------
 #P Geom=AllCheck Guess=TCheck SCRF=Check Test GenChk RHF/6-31G(d) Freq
 ----------------------------------------------------------------------
 #P HF/6-31G* Geom=AllCheck Guess=Read Stable NoSymm
 ----------------------------------------------------------------------
'''


def pilot_route_parser_fixture(record, internal_change=None):
    from go_param.qm import parse_gaussian_log
    routes = parse_gaussian_log(PILOT_ROUTE_EXCERPT)['routes']
    assert len(routes) == 3
    if internal_change:
        routes[1] = internal_change(routes[1])
    text = synthetic_log(record, stable=False)
    original = parse_gaussian_log(text)['routes']
    text = text.replace(original[0], routes[0]).replace(original[1], routes[2])
    text = text.replace(' Harmonic frequencies', ' '+routes[1]+'\n ----------------------------------------------------------------------\n Harmonic frequencies')
    return text.replace(' Normal termination of Gaussian 16\n',
        ' The wavefunction has an RHF -> UHF instability.\n Normal termination of Gaussian 16\n')


def test_real_pilot_internal_frequency_route_retains_only_stability_rejection(tmp_path):
    from go_param.qm import inspect_reference, prepare_orientations
    record = model()
    text = pilot_route_parser_fixture(record)
    report = inspect_reference(record, text)
    assert report['accepted'] is False
    assert report['issues'] == ['stability_not_confirmed']
    manifest = prepare_orientations(record, text, tmp_path/'diagnostic', allow_rhf_instability_diagnostic=True)
    assert manifest['reference_accepted'] is False


@pytest.mark.parametrize('alter', [
    lambda r: r.replace('RHF/6-31G(d)', 'B3LYP/6-31G(d)'),
    lambda r: r.replace('GenChk', ''),
    lambda r: r.replace('Guess=TCheck', 'Guess=Read'),
])
def test_internal_frequency_exception_cannot_bypass_method_or_checkpoint_provenance(alter):
    from go_param.qm import inspect_reference
    record = model()
    report = inspect_reference(record, pilot_route_parser_fixture(record, alter))
    assert report['accepted'] is False
    assert 'route_provenance_mismatch' in report['issues']


def test_internal_frequency_exception_requires_prior_explicit_opt_freq_symmetry():
    from go_param.qm import inspect_reference
    record = model()
    text = pilot_route_parser_fixture(record).replace('NoSymm IOp', 'Symmetry=On IOp')
    assert 'route_provenance_mismatch' in inspect_reference(record, text)['issues']
