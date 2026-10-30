"""RESP identity/constraint tests; subprocess fixtures are explicitly NOT QM/RESP fits."""
import importlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pytest

BOHR = 0.529177210903


def api():
    return importlib.import_module('go_param.resp')


@pytest.fixture
def formate():
    elements = ['C', 'O', 'O', 'H']
    positions = [[0., 0., 0.], [1.2, 0., 0.], [-1.2, 0., 0.], [0., 1.1, 0.]]
    return {'model_id': 'explicit_formate_fixture', 'formal_charge': -1, 'multiplicity': 1,
            'atoms': [{'id': i, 'element': e, 'position': xyz,
                       'formal_charge': -1 if i == 3 else 0,
                       'role': 'fixture', 'roles': ['fixture'],
                       'symmetry_class': 2 if e == 'O' else i + 5}
                      for i, (e, xyz) in enumerate(zip(elements, positions), 1)],
            'bonds': [{'a': 1, 'b': i, 'order': o} for i, o in [(2, 2), (3, 1), (4, 1)]]}


def esp_text(coords, points, charges):
    potential = (1 / np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)) @ charges
    return (f'{len(coords):5d}{len(points):5d}{0:5d}\n'
            + ''.join(' ' * 17 + ''.join(f'{v:16.7E}' for v in x) + '\n' for x in coords)
            + ''.join(' ' + ''.join(f'{v:16.7E}' for v in [p, *x]) + '\n'
                      for p, x in zip(potential, points)))


def gaussian_text(model, coords, esp):
    numbers = {'H': 1, 'C': 6, 'O': 8}
    return (' FIXTURE ONLY: NOT A QUANTUM CALCULATION\n'
            ' #p HF/6-31G* NoSymm Pop=MK SCF=Tight\n ----------------\n'
            f" Charge = {model['formal_charge']} Multiplicity = {model['multiplicity']}\n"
            ' Input orientation:\n ---------------------------------------------------------------------\n'
            ' Center Atomic Atomic Coordinates (Angstroms)\n Number Number Type X Y Z\n'
            ' ---------------------------------------------------------------------\n'
            + ''.join(f" {i} {numbers[a['element']]} 0 {x[0]:.8f} {x[1]:.8f} {x[2]:.8f}\n"
                      for i, (a, x) in enumerate(zip(model['atoms'], coords), 1))
            + ' ---------------------------------------------------------------------\n'
            + ' SCF Done:  E(RHF) =  -188.123456789 A.U. after 12 cycles\n'
            + ''.join(f' Atomic Center {i} is at {x[0]:.8f} {x[1]:.8f} {x[2]:.8f}\n'
                      for i, x in enumerate(coords, 1))
            + f' TEST_ESP_SOURCE: {esp}\n Normal termination of Gaussian 16\n')


def write_case(tmp_path, model, coords=None, name='a', charges=None, nrepeat=1):
    model_path = tmp_path / 'model.json'
    model_path.write_text(json.dumps(model))
    if coords is None:
        coords = np.array([a['position'] for a in model['atoms']])
    if charges is None:
        charges = [.5, -.8, -.8, .1]
    points = np.tile(np.array([[8., 5., 4.], [-8., 4., -4.], [0., -8., 4.],
                               [4., 4., -8.], [8., -6., 6.], [-4., 10., 8.]]), (nrepeat, 1))
    esp = tmp_path / f'{name}.esp'
    esp.write_text(esp_text(np.array(coords) / BOHR, points, charges))
    log = tmp_path / f'{name}.log'
    log.write_text(gaussian_text(model, coords, esp))
    return model_path, log, esp


def external_fixtures(tmp_path, stage1, stage2=None, failure=False):
    """Executables only transport explicit fixtures; no numerical fitting is simulated."""
    eg = tmp_path / 'espgen_fixture'
    eg.write_text(f'#!{sys.executable}\n' + '''import pathlib,re,shutil,sys
args=sys.argv
source=re.search(r'TEST_ESP_SOURCE: (.+)',pathlib.Path(args[args.index('-i')+1]).read_text())[1]
shutil.copyfile(source,pathlib.Path(args[args.index('-o')+1]))
print('TEST STUB: explicit ESP fixture copied; no espgen computation')
''')
    rp = tmp_path / 'resp_fixture'
    rp.write_text(f'#!{sys.executable}\n' + f'STAGE1={stage1!r}\nSTAGE2={stage2!r}\nFAIL={failure!r}\n' + '''import pathlib,re,sys
args=sys.argv
inp=pathlib.Path(args[args.index('-i')+1]).read_text()
n=int(re.search(r'nmol\\s*=\\s*(\\d+)',inp)[1])
q=STAGE2 if STAGE2 is not None and 'stage2' in args[args.index('-i')+1] else STAGE1
q=q*n
for key in ['-o','-p','-s']:
 pathlib.Path(args[args.index(key)+1]).write_text('TEST FIXTURE ONLY\\n'+('after 24 iterations, no convergence!' if FAIL else 'Convergence in 4 iterations\\nESP relative RMS'))
pathlib.Path(args[args.index('-t')+1]).write_text(''.join(''.join(f'{v:10.6f}' for v in q[i:i+8])+'\\n' for i in range(0,len(q),8)))
print('TEST STUB: predeclared charges copied; no RESP fitting')
''')
    eg.chmod(0o755)
    rp.chmod(0o755)
    return str(eg), str(rp)


def test_fixed_width_header_can_have_adjacent_atom_and_large_grid_counts():
    assert api().esp_header(f'{28:5d}{25000:5d}{0:5d}') == (28, 25000)
    with pytest.raises(ValueError, match='header'):
        api().esp_header('   28   10    1')


def test_esp_units_and_strict_record_count(tmp_path, formate):
    _, _, esp = write_case(tmp_path, formate)
    xyz, points, potential = api().read_esp(esp)
    np.testing.assert_allclose(xyz * BOHR, [a['position'] for a in formate['atoms']], atol=1e-7)
    assert points.shape == (6, 3) and potential.shape == (6,)
    esp.write_text(esp.read_text() + ' 0.0 0.0 0.0 0.0\n')
    with pytest.raises(ValueError, match='record|trailing'):
        api().read_esp(esp)


@pytest.mark.parametrize('corruption', ['nan', 'truncated'])
def test_malformed_esp_rejected(tmp_path, formate, corruption):
    _, _, esp = write_case(tmp_path, formate)
    text = esp.read_text()
    esp.write_text(text.replace('0.0000000E+00', '            NaN', 1) if corruption == 'nan'
                   else '\n'.join(text.splitlines()[:-1]) + '\n')
    with pytest.raises(ValueError):
        api().read_esp(esp)


def molecule_blocks(text):
    """Read actual fixed-width RESP decks independently of the generator."""
    lines = text.split('&end', 1)[1].splitlines()
    lines = lines[1:]
    blocks = []
    while lines and lines[0].strip():
        weight = float(lines[0])
        charge, n = int(lines[2][:5]), int(lines[2][5:10])
        rows = [(int(x[:5]), int(x[5:10])) for x in lines[3:3+n]]
        blocks.append((weight, charge, rows))
        lines = lines[4+n:]
    return blocks, lines


def test_formate_equivalence_and_one_independent_total_constraint(formate):
    inp = api().make_resp_input(formate, stage=1, weights=[1., 1.])
    blocks, tail = molecule_blocks(inp)
    assert [b[1] for b in blocks] == [-1, -99]
    assert [b[2] for b in blocks] == [[(6, 0), (8, 0), (8, 2), (1, 0)]] * 2
    # Empty intermolecular charge section, then shared representatives across orientations.
    assert tail[0].strip() == ''
    pairs = [list(map(int, line.split())) for line in tail if line.strip() and len(line.split()) == 4]
    assert pairs == [[1, 1, 2, 1], [1, 2, 2, 2], [1, 4, 2, 4]]
    assert 'qwt = 0.0005' in inp and 'ihfree = 1' in inp


def test_stage2_only_aliphatic_methyl_methylene_groups(formate):
    assert api().stage2_groups(formate) == []
    model = {'atoms': [{'id': 1, 'element': 'C', 'symmetry_class': 1}]
                     + [{'id': i, 'element': 'H', 'symmetry_class': 2} for i in range(2, 5)]
                     + [{'id': 5, 'element': 'O', 'symmetry_class': 3},
                        {'id': 6, 'element': 'H', 'symmetry_class': 4}],
             'bonds': [{'a': 1, 'b': i, 'order': 1} for i in range(2, 6)]
                      + [{'a': 5, 'b': 6, 'order': 1}], 'formal_charge': 0}
    assert api().stage2_groups(model) == [[1, 2, 3, 4]]
    blocks, _ = molecule_blocks(api().make_resp_input(model, stage=1, weights=[1.]))
    assert [x[1] for x in blocks[0][2]][:4] == [0, 0, 0, 0]
    blocks, _ = molecule_blocks(api().make_resp_input(model, stage=2, weights=[1., 1.]))
    assert [x[1] for x in blocks[0][2]] == [0, 0, 2, 2, -1, -1]


def test_atomic_pipeline_uses_fixtures_and_preserves_every_attempt(tmp_path, formate):
    model, log, _ = write_case(tmp_path, formate)
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    result = api().fit_atomic(model, [log], tmp_path/'fits', espgen=eg, resp=rp)
    out = Path(result['attempt_dir'])
    assert result['stage2']['status'] == 'noop_no_methyl_methylene'
    assert not (out/'stage2.out').exists()
    assert result['charge_sum_e'] == pytest.approx(-1.)
    assert result['metrics'][0]['esp_rmse_au'] < 1e-8
    assert result['metrics'][0]['esp_rrms'] < 1e-6
    assert result['metrics'][0]['dipole']['origin_angstrom'] == [0., 0., 0.]
    assert (out/'charges.csv').is_file() and (out/'charges.json').is_file()
    assert (out/'FIT_COMPLETE.json').is_file()
    assert result['provenance']['model_sha256']
    assert result['provenance']['gaussian_logs'][0]['sha256']
    assert len((out/'commands.jsonl').read_text().splitlines()) == 2
    old = (out/'charges.csv').read_bytes()
    again = api().fit_atomic(model, [log], tmp_path/'fits', espgen=eg, resp=rp)
    assert again['attempt_dir'] != result['attempt_dir']
    assert (out/'charges.csv').read_bytes() == old


def test_rotated_same_conformer_and_point_count_weighting(tmp_path, formate):
    model, log1, _ = write_case(tmp_path, formate)
    coords = np.array([a['position'] for a in formate['atoms']])
    coords = coords @ np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) + [1., 2., 3.]
    _, log2, _ = write_case(tmp_path, formate, coords=coords, name='b', nrepeat=2)
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    result = api().fit_atomic(model, [log1, log2], tmp_path/'fits', espgen=eg, resp=rp)
    weights = result['orientation_weights']
    assert weights[0]**2 * 6 == pytest.approx(weights[1]**2 * 12, rel=3e-5)
    assert len(result['metrics']) == 2


@pytest.mark.parametrize('problem', ['wrong_element', 'nonrigid', 'wrong_charge', 'wrong_method', 'wrong_esp_units'])
def test_identity_and_protocol_gates_reject_before_resp(tmp_path, formate, problem):
    model, log, esp = write_case(tmp_path, formate)
    logs = [log]
    if problem == 'wrong_element':
        log.write_text(log.read_text().replace(' 2 8 0 ', ' 2 6 0 '))
    elif problem == 'wrong_charge':
        log.write_text(log.read_text().replace('Charge = -1', 'Charge = 0'))
    elif problem == 'wrong_method':
        log.write_text(log.read_text().replace('HF/6-31G*', 'B3LYP/6-31G*'))
    elif problem == 'wrong_esp_units':
        esp.write_text(esp_text(np.array([a['position'] for a in formate['atoms']]),
                                np.array([[8., 5., 4.]]), [.5, -.8, -.8, .1]))
    else:
        coords = np.array([a['position'] for a in formate['atoms']]); coords[1, 0] += .05
        _, second, _ = write_case(tmp_path, formate, coords=coords, name='b')
        logs.append(second)
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    with pytest.raises(ValueError):
        api().fit_atomic(model, logs, tmp_path/'fits', espgen=eg, resp=rp)
    attempt = next((tmp_path/'fits').glob('attempt_*'))
    assert (attempt/'FIT_FAILED.json').exists()
    assert not (attempt/'stage1.charges').exists()
    assert not (attempt/'FIT_COMPLETE.json').exists()


@pytest.mark.parametrize('charges', [[.5, -.7, -.9, .1], [.5, -.8, -.8, .2]])
def test_charge_equivalence_and_integer_sum_checked(tmp_path, formate, charges):
    model, log, _ = write_case(tmp_path, formate)
    eg, rp = external_fixtures(tmp_path, charges)
    with pytest.raises(ValueError, match='equivalence|charge sum'):
        api().fit_atomic(model, [log], tmp_path/'fits', espgen=eg, resp=rp)
    assert not list((tmp_path/'fits').glob('*/FIT_COMPLETE.json'))


def test_nonconvergence_not_accepted_even_with_finite_charge_file(tmp_path, formate):
    model, log, _ = write_case(tmp_path, formate)
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1], failure=True)
    with pytest.raises(ValueError, match='converg'):
        api().fit_atomic(model, [log], tmp_path/'fits', espgen=eg, resp=rp)


def test_duplicate_logs_and_invalid_stable_ids_rejected(tmp_path, formate):
    model, log, _ = write_case(tmp_path, formate)
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    with pytest.raises(ValueError, match='duplicate|Duplicate'):
        api().fit_atomic(model, [log, log], tmp_path/'fits', espgen=eg, resp=rp)
    formate['atoms'][1]['id'] = 1
    model.write_text(json.dumps(formate))
    with pytest.raises(ValueError, match='atom IDs|atom ids'):
        api().fit_atomic(model, [log], tmp_path/'fits', espgen=eg, resp=rp)


def test_amber_espgen_64_column_records_are_read_as_written(tmp_path, formate):
    # Amber espgen.c uses %32.7E%16.7E%16.7E and four %16.7E fields,
    # without the additional prefix blank in the RESP manual.
    _, _, esp = write_case(tmp_path, formate)
    lines = esp.read_text().splitlines()
    esp.write_text(lines[0] + '\n' + '\n'.join(x[1:] for x in lines[1:]) + '\n')
    xyz, points, potential = api().read_esp(esp)
    assert xyz.shape == (4, 3) and points.shape == (6, 3)
    assert potential[0] != 0


def methyl_model():
    xyz = [[0., 0., 0.], [.629, .629, .629], [-.629, -.629, .629],
           [-.629, .629, -.629], [.82, -.82, -.82], [1.38, -1.38, -1.38]]
    elements = ['C', 'H', 'H', 'H', 'O', 'H']
    return {'model_id': 'methanol_fixture', 'formal_charge': 0, 'multiplicity': 1,
            'atoms': [{'id': i, 'element': e, 'position': x, 'formal_charge': 0,
                       'symmetry_class': 2 if 2 <= i <= 4 else i + 10}
                      for i, (e, x) in enumerate(zip(elements, xyz), 1)],
            'bonds': [{'a': 1, 'b': i, 'order': 1} for i in range(2, 6)]
                     + [{'a': 5, 'b': 6, 'order': 1}]}


def test_stage2_executes_only_declared_groups_and_checks_frozen_atoms(tmp_path):
    model = methyl_model()
    q1 = [-.1, .05, .1, .15, -.6, .4]
    q2 = [-.13, .11, .11, .11, -.6, .4]
    path, log, _ = write_case(tmp_path, model, charges=q2)
    eg, rp = external_fixtures(tmp_path, q1, q2)
    result = api().fit_atomic(path, [log], tmp_path/'fits', espgen=eg, resp=rp)
    assert result['stage2']['status'] == 'completed'
    assert result['stage2']['max_abs_change_e'] == pytest.approx(.06)
    assert result['metrics'][0]['esp_rrms'] < 1e-6
    # Preserve net charge while violating the frozen oxygen; the gate must see it.
    eg, rp = external_fixtures(tmp_path, q1, [-.12, .11, .11, .11, -.61, .4])
    with pytest.raises(ValueError, match='frozen'):
        api().fit_atomic(path, [log], tmp_path/'fits', espgen=eg, resp=rp)


def test_dipole_uses_explicit_origin_and_translates_correctly_for_an_ion(tmp_path, formate):
    path, log1, _ = write_case(tmp_path, formate)
    coords = np.array([a['position'] for a in formate['atoms']])
    translation = np.array([2., 3., -1.])
    _, log2, _ = write_case(tmp_path, formate, coords=coords+translation, name='translated')
    eg, rp = external_fixtures(tmp_path, [.5, -.8, -.8, .1])
    result = api().fit_atomic(path, [log1, log2], tmp_path/'fits', espgen=eg, resp=rp)
    d1, d2 = [np.array(m['dipole']['vector_debye']) for m in result['metrics']]
    np.testing.assert_allclose(d2, d1-translation*4.80320471257, atol=2e-6)
    assert all(m['dipole']['origin_dependent_for_ions'] for m in result['metrics'])


@pytest.mark.parametrize('keyword', ['NoSymm', 'Symmetry=None', 'Symmetry = None'])
def test_mk_gate_accepts_symmetry_none_without_changing_coordinates(tmp_path, formate, keyword):
    _, log, _ = write_case(tmp_path, formate)
    text = log.read_text().replace('NoSymm', keyword)
    coords = api()._gaussian_geometry(text, formate)
    np.testing.assert_allclose(coords, [a['position'] for a in formate['atoms']], atol=1e-12)


def test_mk_gate_requires_explicit_symmetry_disabling(tmp_path, formate):
    _, log, _ = write_case(tmp_path, formate)
    with pytest.raises(ValueError, match='Symmetry|symmetry|NoSymm'):
        api()._gaussian_geometry(log.read_text().replace('NoSymm', 'Symmetry=On'), formate)
