"""Controller safety tests using explicit stage mocks; no QM/RESP results are claimed."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from go_param import pipeline

INSTABILITY = 'The wavefunction has an RHF -> UHF instability.'


def forbid(*args, **kwargs):
    pytest.fail('Rejected evidence must not reach a downstream scientific stage')


@pytest.fixture
def source(tmp_path):
    model = {'model_id': 'controller_fixture', 'formal_charge': 0, 'multiplicity': 1,
             'atoms': [{'id': i, 'element': e, 'position': xyz, 'symmetry_class': i,
                        'formal_charge': 0, 'role': 'test_fixture', 'roles': ['test_fixture']}
                       for i, (e, xyz) in enumerate([('O', [0., 0., 0.]), ('H', [.95, 0., 0.]),
                                                     ('H', [-.3, .9, 0.])], 1)],
             'bonds': [{'a': 1, 'b': 2, 'order': 1}, {'a': 1, 'b': 3, 'order': 1}]}
    path = tmp_path/'model.json'; path.write_text(json.dumps(model))
    log = tmp_path/'reference.log'; log.write_text('CONTROLLER TEST FIXTURE ONLY\n'+INSTABILITY+'\n')
    return path, log, model


@pytest.mark.parametrize('issues', [['stability_not_confirmed'], ['frequency_count_mismatch'],
                                   ['stability_not_confirmed', 'small_negative_frequencies_manual_review']])
def test_rejected_reference_never_runs_gaussian_or_resp(tmp_path, monkeypatch, source, issues):
    model, log, _ = source
    gate = {'accepted': False, 'issues': issues}
    monkeypatch.setattr(pipeline, 'inspect_reference', lambda *a: copy.deepcopy(gate))
    monkeypatch.setattr(pipeline, 'prepare_orientations', forbid)
    monkeypatch.setattr(pipeline, 'fit_atomic', forbid)
    monkeypatch.setattr(pipeline.subprocess, 'run', forbid)
    out = tmp_path/'attempt'
    with pytest.raises(ValueError, match='Reference rejected'):
        pipeline.run(model, log, out)
    assert json.loads((out/'reference_gate.json').read_text())['accepted'] is False
    assert (out/'POSTPROCESS_FAILED.json').is_file()
    assert not (out/'POSTPROCESS_COMPLETE.json').exists()


def test_nonfinite_rejected_gate_is_preserved_as_portable_json(tmp_path, monkeypatch, source):
    model, log, _ = source
    monkeypatch.setattr(pipeline, 'inspect_reference', lambda *a: {
        'accepted': False, 'issues': ['nonfinite_frequencies'], 'frequencies_cm1': [float('nan'), float('inf')]})
    monkeypatch.setattr(pipeline, 'prepare_orientations', forbid)
    with pytest.raises(ValueError, match='Reference rejected'):
        pipeline.run(model, log, tmp_path/'attempt')
    saved = (tmp_path/'attempt/reference_gate.json').read_text()
    assert 'NaN' not in saved and 'Infinity' not in saved
    assert json.loads(saved)['frequencies_cm1'] == [None, None]


def prepared_stub(monkeypatch, source):
    model_path, _, model = source
    def prepare(*args, **kwargs):
        out = Path(args[2]); out.mkdir()
        manifest = {'optimized_model_path': str(model_path), 'orientations': []}
        for i in range(3):
            folder = out/f'orientation_{i+1:02d}'; folder.mkdir()
            path = folder/'mk.gjf'; path.write_text('STAGE MOCK: never sent to Gaussian\n')
            manifest['orientations'].append({'input_path': str(path), 'input_sha256': f'fixture_{i}'})
        return manifest
    monkeypatch.setattr(pipeline, 'prepare_orientations', prepare)
    return model


def test_gaussian_nonzero_exit_blocks_resp_and_complete_marker(tmp_path, monkeypatch, source):
    model, log, _ = source
    monkeypatch.setattr(pipeline, 'inspect_reference', lambda *a: {'accepted': True, 'issues': []})
    prepared_stub(monkeypatch, source)
    monkeypatch.setattr(pipeline.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=1))
    monkeypatch.setattr(pipeline, 'fit_atomic', forbid)
    out = tmp_path/'attempt'
    with pytest.raises(RuntimeError, match='Gaussian orientation failed'):
        pipeline.run(model, log, out)
    assert (out/'POSTPROCESS_FAILED.json').is_file()
    assert not (out/'POSTPROCESS_COMPLETE.json').exists()


def test_zeroexit_incomplete_gaussian_logs_fail_real_resp_input_gate(tmp_path, monkeypatch, source):
    model, log, _ = source
    monkeypatch.setattr(pipeline, 'inspect_reference', lambda *a: {'accepted': True, 'issues': []})
    prepared_stub(monkeypatch, source)
    def incomplete(*args, **kw):
        # The process returncode is explicitly mocked; actual RESP validation runs.
        kw['stdout'].write(b'INCOMPLETE GAUSSIAN PROCESS FIXTURE\n')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(pipeline.subprocess, 'run', incomplete)
    out = tmp_path/'attempt'
    with pytest.raises(ValueError, match='normal termination'):
        pipeline.run(model, log, out)
    assert list((out/'fits').glob('*/FIT_FAILED.json'))
    assert not list(out.rglob('FIT_COMPLETE.json'))
    assert not (out/'POSTPROCESS_COMPLETE.json').exists()


def test_existing_output_is_never_overwritten(tmp_path, monkeypatch, source):
    model, log, _ = source
    out = tmp_path/'attempt'; out.mkdir()
    original = out/'reference_gate.json'; original.write_text('PRESERVED ORIGINAL')
    monkeypatch.setattr(pipeline, 'inspect_reference', forbid)
    with pytest.raises(FileExistsError):
        pipeline.run(model, log, out)
    assert original.read_text() == 'PRESERVED ORIGINAL'
    assert not (out/'POSTPROCESS_FAILED.json').exists()


def install_transport_stage_mocks(tmp_path, monkeypatch, source, *, gate=None, scores=None):
    model_path, log, model = source
    monkeypatch.setattr(pipeline, 'inspect_reference', lambda *a: copy.deepcopy(gate or {'accepted': True, 'issues': []}))
    prepared_stub(monkeypatch, source)
    monkeypatch.setattr(pipeline.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=0))
    # Predeclared arithmetic fixtures test controller data transport only.
    q = [-.834001, .417002, .417002]  # +3 micro-e: acceptable six-decimal RESP rounding.
    atoms = [{'atom_id': a['id'], 'charge_e': v, 'element': a['element']} for a, v in zip(model['atoms'], q)]
    charges = tmp_path/'PREDECLARED_TEST_CHARGES.json'; charges.write_text(json.dumps({'atoms': atoms}))
    fit = {'attempt_dir': str(tmp_path), 'charges_json': str(charges), 'status': 'TEST_TRANSPORT_MOCK_NOT_A_FIT'}
    monkeypatch.setattr(pipeline, 'fit_atomic', lambda *a, **kw: fit)
    xyz = np.array([a['position'] for a in model['atoms']])/pipeline.BOHR_ANGSTROM
    points = np.array([[8., 4., 3.], [-8., 4., 3.]])
    potential = (1/np.linalg.norm(points[:, None, :]-xyz[None, :, :], axis=2))@np.array(q)
    monkeypatch.setattr(pipeline, 'read_esp', lambda *a: (xyz, points, potential))
    monkeypatch.setattr(pipeline, 'CANDIDATES', ())
    if scores is not None:
        monkeypatch.setattr(pipeline, 'compare_esp', lambda *a, **kw: scores)
    return model_path, log, tmp_path/'attempt'


def test_resp_rounding_tolerance_is_recorded_without_rescaling_charges(tmp_path, monkeypatch, source):
    args = install_transport_stage_mocks(tmp_path, monkeypatch, source)
    result = pipeline.run(*args)
    score = json.loads((args[2]/'charge_comparison.json').read_text())[0]
    assert score['charge_tolerance_e'] == pytest.approx(3*5.1e-7+2e-6)
    row = score['results'][0]
    assert row['status'] == 'matched_charge_control'
    assert row['charge_error_e'] == pytest.approx(3e-6)
    assert row['rmse_au'] < 1e-14
    assert result['status'] == 'ATOMIC_REFERENCE_AND_TRAINING_ESP_CONTROLS_COMPLETE'


def test_skipped_atomic_candidate_cannot_produce_complete_marker(tmp_path, monkeypatch, source):
    scores = {'results': [{'candidate': 'new_atomic_resp_reference', 'status': 'skipped', 'reason': 'charge_mismatch'}]}
    args = install_transport_stage_mocks(tmp_path, monkeypatch, source, scores=scores)
    with pytest.raises(ValueError, match='fitted|atomic|RESP'):
        pipeline.run(*args)
    assert not (args[2]/'POSTPROCESS_COMPLETE.json').exists()


def test_explicit_rhf_instability_diagnostic_keeps_rejection_and_distinct_status(tmp_path, monkeypatch, source):
    gate = {'accepted': False, 'issues': ['stability_not_confirmed']}
    args = install_transport_stage_mocks(tmp_path, monkeypatch, source, gate=gate)
    result = pipeline.run(*args, rhf_diagnostic=True)
    saved = json.loads((args[2]/'reference_gate.json').read_text())
    assert saved['accepted'] is False and saved['issues'] == ['stability_not_confirmed']
    assert 'DIAGNOSTIC' in result['status']
    assert (args[2]/'POSTPROCESS_DIAGNOSTIC_COMPLETE.json').is_file()
    assert not (args[2]/'POSTPROCESS_COMPLETE.json').exists()


@pytest.mark.parametrize('issues,log_text', [
    (['stability_not_confirmed', 'significant_imaginary_frequencies'], INSTABILITY),
    (['stability_not_confirmed'], 'The wavefunction is unstable.'),
    (['stability_not_confirmed'], 'The wavefunction has an internal instability.'),
])
def test_rhf_diagnostic_cannot_bypass_other_evidence(tmp_path, monkeypatch, source, issues, log_text):
    model, log, _ = source; log.write_text(log_text+'\n')
    monkeypatch.setattr(pipeline, 'inspect_reference', lambda *a: {'accepted': False, 'issues': issues})
    monkeypatch.setattr(pipeline, 'prepare_orientations', forbid)
    monkeypatch.setattr(pipeline, 'fit_atomic', forbid)
    monkeypatch.setattr(pipeline.subprocess, 'run', forbid)
    with pytest.raises(ValueError, match='Reference rejected'):
        pipeline.run(model, log, tmp_path/'attempt', rhf_diagnostic=True)
