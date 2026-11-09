"""Explicit QM method controls; synthetic fixtures are not scientific results."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from go_param import pipeline, qm, resp
from test_qm import model, synthetic_log
from test_resp import formate, write_case, external_fixtures
from test_pipeline import source

DFT = 'wB97XD/6-311++G(d,p)'
HF = 'HF/6-31G*'


def dft_log(record):
    return synthetic_log(record).replace(HF, DFT).replace('E(RHF)', 'E(RwB97XD)')


def test_explicit_dft_inputs_and_three_orientations_retain_method(tmp_path):
    record = model()
    text = qm.gaussian_input(record, method=DFT)
    assert text.count('#P '+DFT+' ') == 2
    assert text.count('Integral=UltraFine') == 2
    report = qm.inspect_reference(record, dft_log(record), method=DFT)
    assert report['accepted'] and report['issues'] == []
    assert report['qm_method'] == DFT
    assert report['final_scf_energy_hartree'] == -765.123456790
    assert report['final_scf_energy_method'] == 'RwB97XD'
    assert report['final_hf_energy_hartree'] is None
    metadata = report['optimized_model']['metadata']
    assert metadata['qm_method'] == DFT
    assert metadata['geometry_status'] == 'wB97XD_6-311++G(d,p)_optimized_frequency_and_stability_validated'
    manifest = qm.prepare_orientations(record, dft_log(record), tmp_path/'prepared', method=DFT)
    assert manifest['qm_method'] == DFT
    assert manifest['scope'] == 'dft_resp_control'
    assert manifest['final_parameter_validation'] is False
    assert len(manifest['orientations']) == 3
    for entry in manifest['orientations']:
        assert '#P '+DFT+' Pop=MK ' in Path(entry['input_path']).read_text()
        assert 'Integral=UltraFine' in Path(entry['input_path']).read_text()
        assert json.loads(Path(entry['model_path']).read_text())['metadata']['qm_method'] == DFT


def test_hf_default_is_preserved_and_dft_requires_explicit_selection():
    record = model()
    report = qm.inspect_reference(record, synthetic_log(record))
    assert report['accepted']
    assert report['optimized_model']['metadata']['geometry_status'] == 'HF_6-31Gstar_optimized_frequency_and_stability_validated'
    assert not qm.inspect_reference(record, dft_log(record))['accepted']
    assert not qm.inspect_reference(record, synthetic_log(record), method=DFT)['accepted']


@pytest.mark.parametrize('alter', [
    lambda text: text.replace(DFT+' Geom=', HF+' Geom='),
    lambda text: text.replace(DFT+' Opt=', 'B3LYP/6-311++G(d,p) Opt='),
    lambda text: text.replace(DFT, DFT+'Extra'),
    lambda text: text.replace('E(RwB97XD)', 'E(RB3LYP)'),
    lambda text: text.replace('E(RwB97XD)', 'E(RHF)', 1),
])
def test_dft_reference_rejects_mixed_or_wrong_route_and_energy_provenance(alter):
    record = model()
    report = qm.inspect_reference(record, alter(dft_log(record)), method=DFT)
    assert not report['accepted']
    assert any('provenance' in issue or 'energy' in issue for issue in report['issues'])


def test_dft_internal_frequency_route_accepts_only_matching_parent():
    record = model()
    text = dft_log(record).replace(' Harmonic frequencies',
        ' #P Geom=AllCheck Guess=TCheck SCRF=Check GenChk RwB97XD/6-311++G(d,p) Freq\n'
        ' ---------------------------------------------------------------------\n Harmonic frequencies')
    assert qm.inspect_reference(record, text, method=DFT)['accepted']
    bad = text.replace('GenChk RwB97XD', 'GenChk RB3LYP')
    assert 'route_provenance_mismatch' in qm.inspect_reference(record, bad, method=DFT)['issues']
    for old, new in [('GenChk', ''), ('Guess=TCheck', 'Guess=Read'), ('NoSymm', 'Symmetry=On')]:
        assert not qm.inspect_reference(record, text.replace(old, new), method=DFT)['accepted']


@pytest.mark.parametrize('alter', [
    lambda text: text.replace('Optimization completed.', ''),
    lambda text: text.replace('Frequencies -- 100.0', 'Frequencies -- -30.0', 1),
    lambda text: text.replace('stable under the perturbations considered.', 'unstable'),
    lambda text: text.replace('Charge = 0', 'Charge = 1'),
])
def test_dft_method_preserves_strict_scientific_gates(alter):
    record = model()
    assert not qm.inspect_reference(record, alter(dft_log(record)), method=DFT)['accepted']


def test_rhf_diagnostic_flag_is_hf_only_even_for_stable_dft(tmp_path):
    record = model()
    with pytest.raises(ValueError, match='HF.only|HF/6-31G'):
        qm.prepare_orientations(record, dft_log(record), tmp_path/'blocked', method=DFT,
                                allow_rhf_instability_diagnostic=True)
    assert not (tmp_path/'blocked').exists()
    source = tmp_path/'model.json'; source.write_text(json.dumps(record))
    log = tmp_path/'reference.log'; log.write_text(dft_log(record))
    with pytest.raises(ValueError, match='HF.only|HF/6-31G'):
        pipeline.run(source, log, tmp_path/'pipeline_blocked', method=DFT, rhf_diagnostic=True)


@pytest.mark.parametrize('method', ['B3LYP/6-31G*', 'wB97XD/6-31G*', 'HF/6-311++G(d,p)', ''])
def test_only_two_explicit_methods_are_supported(method, tmp_path):
    with pytest.raises(ValueError, match='method|Method'):
        qm.gaussian_input(model(), method=method)
    with pytest.raises(ValueError, match='method|Method'):
        resp.fit_atomic('absent.json', [], tmp_path/'fits', method=method)


def test_dft_atomic_fit_keeps_resp_prior_and_explicit_provenance(tmp_path, formate):
    path, log, _ = write_case(tmp_path, formate)
    log.write_text(log.read_text().replace(HF, DFT).replace('E(RHF)', 'E(RwB97XD)'))
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    result = resp.fit_atomic(path, [log], tmp_path/'fits', espgen=eg, resp=rp, method=DFT)
    assert result['protocol']['method'] == DFT+' MK ESP; Amber hyperbolic RESP'
    assert result['protocol']['qm_method'] == DFT
    assert result['protocol']['scope'].startswith('DFT-RESP control;')
    assert result['protocol']['stage1_qwt'] == .0005
    assert result['protocol']['stage2_qwt'] == .001
    assert result['protocol']['ihfree'] == 1
    assert result['orientation_weights'] == [1.]
    assert result['final_parameter_validation'] is False
    assert result['provenance']['qm_method'] == DFT
    assert result['provenance']['gaussian_logs'][0]['qm_method'] == DFT
    assert result['metrics'][0]['esp_rrms'] < 1e-6


@pytest.mark.parametrize('other', [HF, 'B3LYP/6-311++G(d,p)', DFT+'Extra'])
def test_dft_mk_method_gate_rejects_other_methods_and_mixed_jobs(tmp_path, formate, other):
    _, log, _ = write_case(tmp_path, formate)
    text = log.read_text().replace(HF, other)
    with pytest.raises(ValueError, match='protocol|method'):
        resp._gaussian_geometry(text, formate, method=DFT)
    # An earlier foreign job must not be hidden behind the final selected route.
    text = ' #p '+other+' NoSymm Pop=MK\n ----------------\n'+log.read_text().replace(HF, DFT)
    with pytest.raises(ValueError, match='protocol|method'):
        resp._gaussian_geometry(text, formate, method=DFT)


def test_atomic_fit_rejects_method_mismatched_optimized_model(tmp_path, formate):
    formate['metadata'] = {'qm_method': HF}
    path, log, _ = write_case(tmp_path, formate)
    log.write_text(log.read_text().replace(HF, DFT).replace('E(RHF)', 'E(RwB97XD)'))
    with pytest.raises(ValueError, match='method'):
        resp.fit_atomic(path, [log], tmp_path/'fits', method=DFT)


def test_pipeline_propagates_dft_method_to_gate_orientations_and_fit(tmp_path, monkeypatch, formate):
    path, log, _ = write_case(tmp_path, formate)
    source = tmp_path/'parent.log'; source.write_text('SYNTHETIC CONTROLLER FIXTURE')
    seen = []
    def gate(*args, **kwargs):
        seen.append(('gate', kwargs.get('method', args[2] if len(args) > 2 else None)))
        return {'accepted': True, 'issues': []}
    def prepare(*args, **kwargs):
        seen.append(('prepare', kwargs['method']))
        folder = Path(args[2]); folder.mkdir()
        inp = folder/'mk.gjf'; inp.write_text('CONTROLLER TRANSPORT FIXTURE')
        return {'optimized_model_path': str(path), 'orientations': [{'input_path': str(inp), 'input_sha256': 'fixture'}]}
    def fit(*args, **kwargs):
        seen.append(('fit', kwargs['method']))
        raise RuntimeError('stopped after method transport; no fit executed')
    monkeypatch.setattr(pipeline, 'inspect_reference', gate)
    monkeypatch.setattr(pipeline, 'prepare_orientations', prepare)
    monkeypatch.setattr(pipeline, 'fit_atomic', fit)
    monkeypatch.setattr(pipeline.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=0))
    with pytest.raises(RuntimeError, match='method transport'):
        pipeline.run(path, source, tmp_path/'output', method=DFT)
    assert seen == [('gate', DFT), ('prepare', DFT), ('fit', DFT)]
    scope = json.loads((tmp_path/'output/execution_scope.json').read_text())
    assert scope['qm_method'] == DFT and scope['scope'] == 'dft_resp_control'
    assert scope['final_parameter_validation'] is False


def test_qm_cli_accepts_explicit_dft_method(tmp_path, capsys):
    record = model()
    path = tmp_path/'model.json'; path.write_text(json.dumps(record))
    log = tmp_path/'reference.log'; log.write_text(dft_log(record))
    assert qm.main(['inspect', str(path), str(log), '--method', DFT]) == 0
    assert json.loads(capsys.readouterr().out)['qm_method'] == DFT


def test_pipeline_slurm_worker_preserves_method(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout='12345\n'))
    output = tmp_path/'run'
    pipeline.main(['--model', 'source.json', '--log', 'source.log', '--output-dir', str(output), '--slurm', '--method', DFT])
    worker = tmp_path/'run_worker.slurm'
    assert "--method 'wB97XD/6-311++G(d,p)'" in worker.read_text()


def test_resp_cli_propagates_explicit_method(monkeypatch, capsys):
    seen = []
    def fit(*args, **kwargs):
        seen.append(kwargs['method'])
        return {'status': 'CLI transport fixture only'}
    monkeypatch.setattr(resp, 'fit_atomic', fit)
    resp.main(['--model', 'source.json', '--logs', 'source.log', '--output-dir', 'unused', '--method', DFT])
    assert seen == [DFT]
    assert json.loads(capsys.readouterr().out)['status'] == 'CLI transport fixture only'


def test_dft_completed_pipeline_marks_method_in_training_comparisons(tmp_path, monkeypatch, source):
    from test_pipeline import install_transport_stage_mocks
    args = install_transport_stage_mocks(tmp_path, monkeypatch, source)
    result = pipeline.run(*args, method=DFT)
    assert result['status'] == 'DFT_RESP_CONTROL_AND_TRAINING_ESP_CONTROLS_COMPLETE'
    assert result['scope']['final_parameter_validation'] is False
    scores = json.loads((args[2]/'charge_comparison.json').read_text())
    assert all(score['qm_method'] == DFT for score in scores)
    assert (args[2]/'POSTPROCESS_COMPLETE.json').is_file()


@pytest.mark.parametrize('method,label', [(HF, 'RwB97XD'), (DFT, 'RHF'), (DFT, 'RB3LYP')])
def test_mk_gate_rejects_scf_method_that_disagrees_with_route(tmp_path, formate, method, label):
    _, log, _ = write_case(tmp_path, formate)
    text = log.read_text().replace(HF, method).replace('E(RHF)', f'E({label})')
    with pytest.raises(ValueError, match='SCF.*method|energy.*method'):
        resp._gaussian_geometry(text, formate, method=method)


@pytest.mark.parametrize('method', [HF, DFT])
@pytest.mark.parametrize('energy', [None, 'NaN', 'Infinity', '-Infinity', 'malformed'])
def test_mk_gate_requires_finite_scf_energy(tmp_path, formate, method, energy):
    _, log, _ = write_case(tmp_path, formate)
    text = log.read_text().replace(HF, method)
    if method == DFT:
        text = text.replace('E(RHF)', 'E(RwB97XD)')
    text = ('\n'.join(line for line in text.splitlines() if 'SCF Done:' not in line)+'\n'
            if energy is None else text.replace('-188.123456789', energy))
    with pytest.raises(ValueError, match='SCF.*energy|energy.*SCF'):
        resp._gaussian_geometry(text, formate, method=method)


@pytest.mark.parametrize('method,label', [(HF, 'RHF'), (DFT, 'RwB97XD')])
def test_mk_gate_accepts_finite_consistent_scf_energy(tmp_path, formate, method, label):
    _, log, _ = write_case(tmp_path, formate)
    text = log.read_text().replace(HF, method).replace('E(RHF)', f'E({label})')
    text = text.replace('-188.123456789', '-1.88123456789D+02')
    coords = resp._gaussian_geometry(text, formate, method=method)
    assert coords.tolist() == [a['position'] for a in formate['atoms']]


@pytest.mark.parametrize('earlier', ['E(RHF) = -100.0', 'E(RwB97XD) = NaN'])
def test_mk_gate_checks_every_scf_record_before_final_matching_energy(tmp_path, formate, earlier):
    _, log, _ = write_case(tmp_path, formate)
    text = log.read_text().replace(HF, DFT).replace('E(RHF)', 'E(RwB97XD)')
    text = text.replace(' SCF Done:', f' SCF Done:  {earlier} A.U. after 1 cycles\n SCF Done:', 1)
    with pytest.raises(ValueError, match='SCF.*(?:method|energy)|energy.*SCF'):
        resp._gaussian_geometry(text, formate, method=DFT)


def test_mixed_hf_esp_on_accepted_dft_geometry_preserves_provenance_and_gates(tmp_path, formate):
    record = model()
    manifest = qm.prepare_orientations(record, dft_log(record), tmp_path/'mixed_prep',
                                       method=DFT, esp_method=HF)
    assert manifest['qm_method'] == f'{HF}//{DFT}'
    assert manifest['geometry_qm_method'] == DFT
    assert manifest['esp_qm_method'] == HF
    assert manifest['scope'] == 'conventional_rhf_mk_esp_on_accepted_dft_geometry'
    prepared_model = json.loads(Path(manifest['optimized_model_path']).read_text())
    assert prepared_model['metadata']['qm_method'] == DFT
    assert prepared_model['metadata']['geometry_qm_method'] == DFT
    assert prepared_model['metadata']['esp_qm_method'] == HF
    assert prepared_model['metadata']['wavefunction_stability_validated'] is False
    assert prepared_model['metadata']['dft_geometry_wavefunction_stability_validated'] is True
    for entry in manifest['orientations']:
        gjf = Path(entry['input_path']).read_text()
        assert f'#P {HF} Pop=MK ' in gjf
        assert 'Integral=UltraFine' not in gjf

    with pytest.raises(ValueError, match='Unsupported mixed'):
        qm.prepare_orientations(record, synthetic_log(record), tmp_path/'bad_mix',
                                method=HF, esp_method=DFT)

    formate['metadata'] = {'qm_method': DFT, 'coordinate_units': 'angstrom'}
    path, log, _ = write_case(tmp_path, formate)
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    result = resp.fit_atomic(path, [log], tmp_path/'mixed_fits', espgen=eg, resp=rp,
                             method=HF, geometry_method=DFT)
    assert result['protocol']['qm_method'] == f'{HF}//{DFT}'
    assert result['protocol']['geometry_qm_method'] == DFT
    assert result['protocol']['esp_qm_method'] == HF
    assert 'conventional_rhf_mk_esp_on_accepted_dft_geometry' in result['protocol']['scope']
    assert result['final_parameter_validation'] is False


def test_pipeline_mixed_hf_on_dft_propagates_separate_methods(tmp_path, monkeypatch, source):
    from test_pipeline import install_transport_stage_mocks
    args = install_transport_stage_mocks(tmp_path, monkeypatch, source)
    result = pipeline.run(*args, method=DFT, esp_method=HF)
    assert result['status'] == 'CONVENTIONAL_HF_ESP_ON_DFT_GEOMETRY_ATOMIC_RESP_AND_TRAINING_ESP_CONTROLS_COMPLETE'
    assert result['scope']['qm_method'] == f'{HF}//{DFT}'
    assert result['scope']['geometry_qm_method'] == DFT
    assert result['scope']['esp_qm_method'] == HF
    scores = json.loads((args[2]/'charge_comparison.json').read_text())
    assert all(score['qm_method'] == f'{HF}//{DFT}' for score in scores)
