"""Atomic-reference RESP fits using system AmberTools espgen and resp only.

This module performs no QM calculation, substitutes no optimizer for RESP, and
fits one physical conformer supplied in one or more rigid orientations. Parent
QM optimization, frequency and wavefunction-stability approval is a separate gate.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

BOHR_ANGSTROM = 0.529177210903
E_ANGSTROM_DEBYE = 4.80320471257
# Deliberately limited to the campaign's explicit C/H/O chemistry.
ATOMIC_NUMBERS = {'H': 1, 'C': 6, 'O': 8}
COVALENT_RADII = {'H': .31, 'C': .76, 'O': .66}
FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][-+]?\d+)?'


def _float(value):
    return float(value.replace('D', 'E').replace('d', 'e'))


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def esp_header(line):
    """Read espgen's adjacent I5 fields, including a fully occupied point field."""
    if len(line) < 10 or line[15:].strip():
        raise ValueError('Malformed ESP header: expected 2I5 or 3I5')
    try:
        nat, npoints = int(line[:5]), int(line[5:10])
        extra = int(line[10:15]) if line[10:15].strip() else 0
    except ValueError as exc:
        raise ValueError('Malformed ESP header') from exc
    if nat < 1 or npoints < 1 or extra != 0:
        raise ValueError('Unsupported ESP header count or multipole data')
    return nat, npoints


def _esp_row(line, offset, count):
    if len(line) < offset + 16 * count or line[:offset].strip() or line[offset + 16 * count:].strip():
        raise ValueError('Malformed ESP fixed-width record')
    try:
        row = [_float(line[offset + 16*i:offset + 16*(i+1)]) for i in range(count)]
    except ValueError as exc:
        raise ValueError('Malformed ESP numeric record') from exc
    if not np.isfinite(row).all():
        raise ValueError('Nonfinite ESP record')
    return row


def read_esp(path):
    """Return atomic coordinates, grid coordinates (bohr), and potential (a.u.)."""
    lines = Path(path).read_text().splitlines()
    if not lines:
        raise ValueError('Empty ESP file')
    nat, npoints = esp_header(lines[0])
    if len(lines) < 1 + nat + npoints or any(x.strip() for x in lines[1+nat+npoints:]):
        raise ValueError('ESP record count differs from header or trailing records present')
    # espgen uses 64 columns; the RESP manual adds one prefix blank (65).
    # Parse both explicit layouts without losing field boundaries.
    xyz = np.array([_esp_row(x, 16 if len(x) == 64 else 17, 3) for x in lines[1:1+nat]])
    grid = np.array([_esp_row(x, 0 if len(x) == 64 else 1, 4) for x in lines[1+nat:1+nat+npoints]])
    if np.any(np.linalg.norm(grid[:, None, 1:] - xyz[None, :, :], axis=2) < 1e-8):
        raise ValueError('ESP grid coincides with atomic center')
    if not np.any(grid[:, 0]):
        raise ValueError('ESP has zero potential everywhere; no fitting signal')
    return xyz, grid[:, 1:], grid[:, 0]


def _validate_model(model):
    atoms = model.get('atoms', [])
    if not atoms or [a.get('id') for a in atoms] != list(range(1, len(atoms)+1)):
        raise ValueError('Stable atom IDs must be contiguous 1-based and in order')
    if any(type(a['id']) is not int for a in atoms):
        raise ValueError('Stable atom IDs must be integers')
    if type(model.get('formal_charge')) is not int or not -98 <= model['formal_charge'] <= 98:
        raise ValueError('Whole-model formal charge must be an integer in RESP range')
    if type(model.get('multiplicity')) is not int or model['multiplicity'] < 1:
        raise ValueError('Multiplicity must be a positive integer')
    if any(a.get('element') not in ATOMIC_NUMBERS for a in atoms):
        raise ValueError('Atomic-reference implementation supports campaign C/H/O models only')
    coords = np.array([a.get('position', []) for a in atoms], dtype=float)
    if coords.shape != (len(atoms), 3) or not np.isfinite(coords).all():
        raise ValueError('Model positions must be finite Angstrom coordinates')
    if model.get('metadata', {}).get('coordinate_units', 'angstrom') != 'angstrom':
        raise ValueError('Model coordinates must be Angstrom')
    if any(type(a.get('formal_charge')) is not int for a in atoms):
        raise ValueError('Each atom needs its integer formal charge')
    if sum(a['formal_charge'] for a in atoms) != model['formal_charge']:
        raise ValueError('Atomic formal charges disagree with whole-model charge')
    classes = {}
    for atom in atoms:
        cls = atom.get('symmetry_class')
        if type(cls) is not int:
            raise ValueError('Every atom needs an integer symmetry_class')
        if cls in classes and classes[cls] != atom['element']:
            raise ValueError('A symmetry class cannot mix elements')
        classes[cls] = atom['element']
    seen = set()
    for bond in model.get('bonds', []):
        a, b, order = bond['a'], bond['b'], bond['order']
        if type(a) is not int or type(b) is not int or not 1 <= a <= len(atoms) or not 1 <= b <= len(atoms) or a == b:
            raise ValueError('Invalid bond atom IDs')
        pair = tuple(sorted((a, b)))
        if pair in seen or order not in (1, 1.0, 1.5, 2, 2.0, 3, 3.0):
            raise ValueError('Invalid/duplicate model bond')
        seen.add(pair)
    if len(atoms) > 1 and not seen:
        raise ValueError('Model bond graph is missing')


def stage2_groups(model):
    """Return carbon + H IDs of tetrahedral CH2/CH3 groups (classic stage 2)."""
    atoms = model['atoms']
    neighbors = {a['id']: [] for a in atoms}
    for bond in model['bonds']:
        neighbors[bond['a']].append((bond['b'], bond['order']))
        neighbors[bond['b']].append((bond['a'], bond['order']))
    groups = []
    for atom in atoms:
        connected = neighbors[atom['id']]
        hydrogens = [i for i, _ in connected if atoms[i-1]['element'] == 'H']
        if (atom['element'] == 'C' and len(connected) == 4 and len(hydrogens) in (2, 3)
                and all(order == 1 for _, order in connected)):
            groups.append([atom['id'], *sorted(hydrogens)])
    return groups


def _variation(model, stage):
    active = {i for group in stage2_groups(model) for i in group}
    # Stage 1 does not equivalence atoms reserved for the second-stage refit.
    first = {}
    values = []
    for atom in model['atoms']:
        idx, cls = atom['id'], atom['symmetry_class']
        if stage == 2 and idx not in active:
            values.append(-1)
        elif stage == 1 and idx in active:
            values.append(0)
        else:
            values.append(first.get(cls, 0))
            first.setdefault(cls, idx)
    return values


def make_resp_input(model, *, stage, weights):
    """Build fixed-width Amber RESP cards with one physical charge constraint.

    All orientation copies of each independent charge center are equivalenced.
    -99 in later molecular charge headers disables duplicate total constraints;
    shared charges still impose the physical integer total in every orientation.
    """
    if stage not in (1, 2) or not weights or any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError('Invalid RESP stage/weights')
    vary = _variation(model, stage)
    lines = [f'Atomic reference RESP stage {stage}', ' &cntrl',
             f' nmol = {len(weights)},', f' iqopt = {1 if stage == 1 else 2},',
             ' ihfree = 1,', ' irstrnt = 1,', ' ioutopt = 1,',
             f' qwt = {0.0005 if stage == 1 else 0.001},', ' &end']
    for j, weight in enumerate(weights):
        lines += [f'{weight:10.5f}', f'orientation {j+1}',
                  f'{model["formal_charge"] if j == 0 else -99:5d}{len(vary):5d}']
        lines += [f'{ATOMIC_NUMBERS[a["element"]]:5d}{v:5d}' for a, v in zip(model['atoms'], vary)]
        lines.append('')  # End intra-molecular constraints for this orientation.
    if len(weights) > 1:
        lines.append('')  # End inter-molecular charge-constraint section.
        # Equivalence only independent UNfrozen centers. Local ivary pointers
        # propagate those roots; frozen copies already inherit identical stage-1 q.
        for idx, value in enumerate(vary, 1):
            if value == 0:
                for j in range(2, len(weights)+1):
                    lines += [f'{2:5d}', f'{1:5d}{idx:5d}{j:5d}{idx:5d}']
        lines.append('')
    return '\n'.join(lines) + '\n'


def _orientation_table(text, before):
    table = list(re.finditer(r'(?:Input|Standard|Z-Matrix) orientation:', text[:before]))
    if not table:
        raise ValueError('Gaussian atomic identity/orientation table missing')
    rows = []
    begun = False
    for line in text[table[-1].end():before].splitlines():
        fields = line.split()
        if len(fields) == 6 and all(re.fullmatch(r'\d+', s) for s in fields[:3]):
            rows.append((int(fields[0]), int(fields[1]), [_float(s) for s in fields[3:]]))
            begun = True
        elif begun:
            break
    return rows


def _gaussian_geometry(text, model):
    if re.search(r'Error termination|Convergence failure|SCF has not converged', text, re.I):
        raise ValueError('Gaussian error/nonconvergence in supplied log')
    if not re.search(r'Normal termination of Gaussian[^\n]*\s*\Z', text):
        raise ValueError('Gaussian log lacks final normal termination')
    routes = re.findall(r'^\s*#[^\n]*(?:\n(?!\s*-{5,})[^\n]*)*', text, re.M)
    if not routes:
        raise ValueError('Gaussian route missing')
    # Collecting until the separator tolerates Gaussian route line wrapping.
    route = ' '.join(routes[-1].split())
    if not re.search(r'\b(?:R?HF)/6-31G(?:\*|\(d\))', route, re.I):
        raise ValueError('Atomic RESP protocol requires HF/6-31G*')
    if not re.search(r'\b(?:NoSymm(?:etry)?|Symmetry\s*=\s*None)\b', route, re.I) or not re.search(r'\b(?:Pop(?:ulation)?\s*=\s*\(?\s*MK|Pop\s*\(\s*MK)', route, re.I):
        raise ValueError('Gaussian route must request NoSymm or Symmetry=None, and Pop=MK')
    charges = re.findall(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)', text)
    if not charges or tuple(map(int, charges[-1])) != (model['formal_charge'], model['multiplicity']):
        raise ValueError('Gaussian charge/multiplicity differs from model')
    atom_re = re.compile(r'Atomic Center\s+(\d+)\s+is at\s+(' + FLOAT + r')\s+(' + FLOAT + r')\s+(' + FLOAT + r')')
    matches = list(atom_re.finditer(text))
    starts = [i for i, match in enumerate(matches) if int(match[1]) == 1]
    if not starts:
        raise ValueError('Gaussian MK atomic centers missing')
    matches = matches[starts[-1]:]
    nat = len(model['atoms'])
    if [int(m[1]) for m in matches] != list(range(1, nat+1)):
        raise ValueError('Gaussian MK atom order/count differs from model')
    coords = np.array([[_float(m[j]) for j in (2, 3, 4)] for m in matches])
    if not np.isfinite(coords).all():
        raise ValueError('Nonfinite Gaussian MK coordinates')
    table = _orientation_table(text, matches[0].start())
    if ([row[0] for row in table] != list(range(1, nat+1))
            or [row[1] for row in table] != [ATOMIC_NUMBERS[a['element']] for a in model['atoms']]):
        raise ValueError('Gaussian atomic identity/order differs from stable model IDs')
    if not np.allclose([row[2] for row in table], coords, atol=2e-5, rtol=0):
        raise ValueError('Gaussian MK coordinates differ from final atomic identity table')
    _check_geometry(model, coords)
    return coords


def _check_geometry(model, coords):
    atoms = model['atoms']
    declared = {tuple(sorted((b['a'], b['b']))) for b in model['bonds']}
    for i, atom in enumerate(atoms):
        for j in range(i+1, len(atoms)):
            radius = COVALENT_RADII[atom['element']] + COVALENT_RADII[atoms[j]['element']]
            distance = np.linalg.norm(coords[i] - coords[j])
            if (i+1, j+1) in declared:
                if not .55 * radius < distance < 1.30 * radius:
                    raise ValueError(f'QM bond graph/atom identity changed at {i+1}-{j+1}: {distance:.5f} Angstrom')
            elif distance < 1.15 * radius:
                raise ValueError(f'QM new bond/atom clash at {i+1}-{j+1}')


def _executable(name):
    found = shutil.which(str(name))
    if found is None:
        raise FileNotFoundError(f'System AmberTools executable unavailable: {name}')
    path = Path(found).resolve()
    return {'path': str(path), 'sha256': _digest(path)}


def _command(args, folder, label):
    start = datetime.now(timezone.utc).isoformat()
    with (folder / f'{label}.log').open('x') as handle:
        result = subprocess.run([str(x) for x in args], cwd=folder, stdout=handle,
                                stderr=subprocess.STDOUT, check=False)
    with (folder/'commands.jsonl').open('a') as handle:
        handle.write(json.dumps({'argv': [str(x) for x in args], 'started_utc': start,
                                 'finished_utc': datetime.now(timezone.utc).isoformat(),
                                 'returncode': result.returncode, 'log': f'{label}.log'}) + '\n')
    if result.returncode:
        raise RuntimeError(f'{label} failed with exit code {result.returncode}; see {folder}')


def _read_charges(path, expected):
    values = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if len(line) % 10 or len(line) > 80:
            raise ValueError('Malformed RESP charge file: expected 8F10.6')
        values += [_float(line[i:i+10]) for i in range(0, len(line), 10)]
    q = np.array(values)
    if len(q) != expected or not np.isfinite(q).all():
        raise ValueError('RESP charge count/nonfinite charge mismatch')
    return q


def _check_charges(q, model, stage, nmol, stage1=None):
    q = q.reshape(nmol, len(model['atoms']))
    if not np.allclose(q, q[0], atol=2e-6, rtol=0):
        raise ValueError('RESP orientation charge equivalence failed')
    # Each atom is written to 6 decimal places by RESP. Sum tolerance propagates
    # at most half a last decimal per atom, with 2 micro-e numerical slack.
    tolerance = len(model['atoms']) * 5.1e-7 + 2e-6
    if np.any(np.abs(q.sum(axis=1) - model['formal_charge']) > tolerance):
        raise ValueError('RESP whole-model integer charge sum failed')
    vary = _variation(model, stage)
    for i, value in enumerate(vary):
        if value > 0 and not np.allclose(q[:, i], q[:, value-1], atol=2e-6, rtol=0):
            raise ValueError('RESP symmetry charge equivalence failed')
        if stage == 2 and value < 0 and not np.allclose(q[:, i], stage1[:, i], atol=1.1e-6, rtol=0):
            raise ValueError('RESP second stage changed frozen atoms')
    return q


def _run_resp(exe, folder, model, stage, weights, stage1=None):
    prefix = f'stage{stage}'
    (folder/f'{prefix}.in').write_text(make_resp_input(model, stage=stage, weights=weights))
    args = [exe, '-i', f'{prefix}.in', '-o', f'{prefix}.out', '-p', f'{prefix}.punch',
            '-e', 'all.esp', '-t', f'{prefix}.charges', '-s', f'{prefix}.esout']
    if stage == 2:
        args += ['-q', 'stage1.charges']
    _command(args, folder, f'resp{stage}')
    for suffix in ('out', 'punch', 'charges'):
        path = folder/f'{prefix}.{suffix}'
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f'RESP required output absent: {path.name}')
    text = (folder/f'{prefix}.out').read_text() + (folder/f'resp{stage}.log').read_text()
    if re.search(r'\bNaN\b|\bInfinity\b|singular|no convergence|not converg|error|illegal|premature|unexpected eof|not enough.*q0', text, re.I):
        raise ValueError(f'RESP stage {stage} error/nonconvergence; see {folder}')
    if not re.search(r'Convergence in\s+\d+\s+iterations', text, re.I):
        raise ValueError(f'RESP stage {stage} convergence evidence absent')
    return _check_charges(_read_charges(folder/f'{prefix}.charges', len(weights)*len(model['atoms'])),
                          model, stage, len(weights), stage1)


def _metrics(coords, points, potential, q):
    squared_error = 0.
    for start in range(0, len(points), 8192):
        matrix = 1 / np.linalg.norm(points[start:start+8192, None, :] - coords[None, :, :], axis=2)
        residual = matrix @ q - potential[start:start+8192]
        squared_error += float(residual @ residual)
    signal = float(potential @ potential)
    vector = (q[:, None] * coords * BOHR_ANGSTROM).sum(axis=0) * E_ANGSTROM_DEBYE
    return {'grid_points': len(points), 'esp_rmse_au': math.sqrt(squared_error/len(points)),
            'esp_rrms': math.sqrt(squared_error/signal),
            'dipole': {'vector_debye': vector.tolist(), 'magnitude_debye': float(np.linalg.norm(vector)),
                       'origin_angstrom': [0., 0., 0.],
                       'frame': 'this Gaussian MK orientation; coordinate origin, not center of mass',
                       'origin_dependent_for_ions': bool(abs(q.sum()) > 1e-4)}}


def _new_attempt(output_dir):
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    i = 1
    while True:
        path = root/f'attempt_{i:03d}'
        try:
            path.mkdir()
        except FileExistsError:
            i += 1
        else:
            return path


def fit_atomic(model_path, gaussian_logs: list, output_dir, *, espgen='espgen', resp='resp') -> dict:
    """Fit one atomic reference; supplied logs are rigid orientations, not conformers.

    A fresh attempt is always retained, including failure evidence. Success writes
    charges.csv, charges.json, FIT_COMPLETE.json inside result['attempt_dir'].
    No parameter library or topology is published by this function.
    """
    out = _new_attempt(output_dir)
    try:
        model_path = Path(model_path).resolve()
        model_bytes = model_path.read_bytes()
        (out/'model.json').write_bytes(model_bytes)
        model = json.loads(model_bytes)
        _validate_model(model)
        if not gaussian_logs:
            raise ValueError('At least one Gaussian MK log is required')
        logs = [Path(p).resolve() for p in gaussian_logs]
        if len(set(logs)) != len(logs):
            raise ValueError('Duplicate Gaussian log paths would bias orientation weighting')
        inputs = out/'gaussian_logs'
        inputs.mkdir()
        evidence, coords, datasets = [], [], []
        digests = set()
        for i, path in enumerate(logs, 1):
            snapshot = inputs/f'orientation_{i:03d}.log'
            snapshot.write_bytes(path.read_bytes())
            sha = _digest(snapshot)
            if sha in digests:
                raise ValueError('Duplicate Gaussian log content would bias orientation weighting')
            digests.add(sha)
            geometry = _gaussian_geometry(snapshot.read_text(), model)
            if coords:
                d1 = np.linalg.norm(coords[0][:, None, :] - coords[0][None, :, :], axis=2)
                d2 = np.linalg.norm(geometry[:, None, :] - geometry[None, :, :], axis=2)
                if not np.allclose(d1, d2, atol=5e-5, rtol=0):
                    raise ValueError('Logs are not rigid orientations of one physical conformer')
                # Exclude reflected chiral geometries even when all distances agree.
                a, b = coords[0] - coords[0].mean(axis=0), geometry - geometry.mean(axis=0)
                u, _, vh = np.linalg.svd(a.T @ b)
                rotation = u @ np.diag([1, 1, np.linalg.det(u @ vh)]) @ vh
                if not np.allclose(a @ rotation, b, atol=5e-5, rtol=0):
                    raise ValueError('Orientation is a reflection or atom mapping changed')
            coords.append(geometry)
            evidence.append({'source_path': str(path), 'snapshot': str(snapshot.relative_to(out)), 'sha256': sha})
        programs = {'espgen': _executable(espgen), 'resp': _executable(resp)}
        for i, (geometry, entry) in enumerate(zip(coords, evidence), 1):
            esp = out/f'orientation_{i:03d}.esp'
            _command([programs['espgen']['path'], '-i', entry['snapshot'], '-o', esp.name], out, f'espgen_{i:03d}')
            data = read_esp(esp)
            if data[0].shape != geometry.shape or not np.allclose(data[0]*BOHR_ANGSTROM, geometry, atol=2e-5, rtol=0):
                raise ValueError('ESP atomic coordinates/order/bohr units disagree with Gaussian MK centers')
            datasets.append(data)
            entry['esp_sha256'] = _digest(esp)
            entry['grid_points'] = len(data[1])
        nref = float(np.mean([len(data[1]) for data in datasets]))
        # RESP reads F10.5 and squares this weight. Restraints occur once per
        # atom COPY; their sum and the ESP sum both scale with orientation count.
        weights = [float(f'{math.sqrt(nref / len(data[1])):.5f}') for data in datasets]
        with (out/'all.esp').open('x') as handle:
            for i in range(1, len(datasets)+1):
                handle.write((out/f'orientation_{i:03d}.esp').read_text().rstrip() + '\n')
        q1 = _run_resp(programs['resp']['path'], out, model, 1, weights)
        groups = stage2_groups(model)
        if groups:
            q = _run_resp(programs['resp']['path'], out, model, 2, weights, stage1=q1)
            second = {'status': 'completed', 'groups_atom_ids': groups,
                      'max_abs_change_e': float(np.max(abs(q-q1)))}
        else:
            (out/'stage2.in').write_text(make_resp_input(model, stage=2, weights=weights))
            q = q1
            second = {'status': 'noop_no_methyl_methylene', 'groups_atom_ids': [],
                      'max_abs_change_e': 0., 'explanation': 'All charges retained from stage 1; stage2.in is audit only and was not executed.'}
        metrics = [dict(_metrics(*data, q[0]), orientation=i) for i, data in enumerate(datasets, 1)]
        if any(not math.isfinite(m['esp_rmse_au']) or not math.isfinite(m['esp_rrms']) for m in metrics):
            raise ValueError('Nonfinite independent ESP metrics')
        records = [dict(atom_id=a['id'], element=a['element'], name=a.get('name', f'A{a["id"]}'),
                        role=a.get('role'), roles=a.get('roles'), symmetry_class=a['symmetry_class'],
                        charge_e=float(charge)) for a, charge in zip(model['atoms'], q[0])]
        _json(out/'charges.json', {'model_id': model.get('model_id'), 'formal_charge': model['formal_charge'], 'atoms': records})
        with (out/'charges.csv').open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['atom_id', 'element', 'name', 'role', 'symmetry_class', 'charge_e'])
            writer.writeheader()
            writer.writerows({key: atom[key] for key in writer.fieldnames} for atom in records)
        result = {'schema_version': 1, 'status': 'atomic_resp_fit_completed', 'model_id': model.get('model_id'),
                  'attempt_dir': str(out), 'charges_csv': str(out/'charges.csv'), 'charges_json': str(out/'charges.json'),
                  'charge_sum_e': float(q[0].sum()), 'formal_charge': model['formal_charge'],
                  'stage2': second, 'orientation_weights': weights, 'metrics': metrics,
                  'protocol': {'method': 'HF/6-31G* MK ESP; Amber hyperbolic RESP', 'stage1_qwt': .0005,
                               'stage2_qwt': .001, 'ihfree': 1, 'physical_conformers': 1,
                               'orientations': len(datasets), 'mean_grid_points': nref,
                               'parent_qm_optimization_frequencies_stability': 'must be approved upstream; not inferred from MK logs',
                               'scope': 'finite whole-model atomic reference; no bulk or modular-transfer validation'},
                  'provenance': {'model_source_path': str(model_path), 'model_sha256': _digest(out/'model.json'),
                                 'gaussian_logs': evidence, 'executables': programs,
                                 'fit_code_sha256': _digest(Path(__file__))}}
        result['artifacts_sha256'] = {str(p.relative_to(out)): _digest(p) for p in sorted(out.rglob('*'))
                                      if p.is_file()}
        _json(out/'FIT_COMPLETE.json', result)
        return result
    except Exception as exc:
        _json(out/'FIT_FAILED.json', {'status': 'failed', 'exception': type(exc).__name__, 'message': str(exc)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True)
    parser.add_argument('--logs', nargs='+', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--espgen', default='espgen')
    parser.add_argument('--resp', default='resp')
    args = parser.parse_args()
    result = fit_atomic(args.model, args.logs, args.output_dir, espgen=args.espgen, resp=args.resp)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
