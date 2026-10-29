"""Pure-stdlib Gaussian input, evidence gates and rigid MK orientation preparation.

No subprocess, scheduler submission, charge fit or scientific result generation.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re

ATOMIC_NUMBERS = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16, 'Cl': 17}
COVALENT_RADII = {'H': .31, 'C': .76, 'N': .71, 'O': .66, 'F': .57, 'S': 1.05, 'Cl': 1.02}
ROTATIONS_DEGREES = ((17., 31., 47.), (-29., 13., 71.), (43., -37., 19.))


def _model(record):
    return json.loads(Path(record).read_text()) if isinstance(record, (str, Path)) else copy.deepcopy(record)


def _log_text(log):
    if isinstance(log, Path):
        return log.read_text(errors='replace')
    if '\n' not in log and len(log) < 4096 and Path(log).is_file():
        return Path(log).read_text(errors='replace')
    return log


def _number(value):
    try:
        return float(value.replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return float('nan')


def _hash_bytes(value):
    return hashlib.sha256(value).hexdigest()


def _source_log_hash(log):
    if isinstance(log, Path) or ('\n' not in log and len(log) < 4096 and Path(log).is_file()):
        return _hash_bytes(Path(log).read_bytes())
    return _hash_bytes(log.encode())


def _json_text(value):
    return json.dumps(value, indent=2, allow_nan=False)+'\n'


def _check_model(record):
    atoms = record['atoms']
    if not atoms or len({a['id'] for a in atoms}) != len(atoms):
        raise ValueError('Model atoms require unique nonempty IDs')
    if not isinstance(record['formal_charge'], int) or not isinstance(record['multiplicity'], int) or record['multiplicity'] < 1:
        raise ValueError('Model requires integer formal charge and positive multiplicity')
    for a in atoms:
        if a['element'] not in ATOMIC_NUMBERS:
            raise ValueError(f'Unsupported element: {a["element"]}')
        if len(a['position']) != 3 or not all(math.isfinite(x) for x in a['position']):
            raise ValueError('Model coordinates must be finite three-vectors')
    electron_count = sum(ATOMIC_NUMBERS[a['element']] for a in atoms)-record['formal_charge']
    if electron_count < 1 or (electron_count-(record['multiplicity']-1)) % 2:
        raise ValueError('Charge/electron count and multiplicity parity mismatch')


def gaussian_input(model_record, cpus=8, mem_gb=6):
    """Two-link HF/6-31G* Opt/Freq/MK then wavefunction stability check."""
    record = _model(model_record)
    _check_model(record)
    if not isinstance(cpus, int) or cpus < 1 or not isinstance(mem_gb, int) or mem_gb < 1:
        raise ValueError('cpus and mem_gb must be positive integers')
    lines = [f'%NProcShared={cpus}', f'%Mem={mem_gb}GB', '%Chk=reference.chk',
             '#P HF/6-31G* Opt=(Tight,CalcFC,MaxCycles=200) Freq Pop=MK SCF=(Tight,XQC) Symmetry=None IOp(6/33=2,6/41=10,6/42=17)',
             '', f'{record["model_id"]}: HF reference optimization, frequency and MK evidence', '',
             f'{record["formal_charge"]} {record["multiplicity"]}']
    for a in record['atoms']:
        lines.append(f'{a["element"]:2s} '+ ' '.join(f'{x:.12f}' for x in a['position']))
    lines += ['', '--Link1--', f'%NProcShared={cpus}', f'%Mem={mem_gb}GB', '%Chk=reference.chk',
              '#P HF/6-31G* Geom=AllCheck Guess=Read Stable Symmetry=None SCF=(Tight,XQC)', '', '']
    return '\n'.join(lines)


def parse_gaussian_log(log):
    """Extract evidence, never infer missing evidence from normal termination."""
    text = _log_text(log)
    orientations = []
    lines = text.splitlines()
    routes = []
    for index, line in enumerate(lines):
        if not re.match(r'^\s*#', line):
            continue
        route_lines = [line.strip()]
        cursor = index+1
        while cursor < len(lines) and lines[cursor].strip() and not re.match(r'^\s*-{5,}', lines[cursor]):
            route_lines.append(lines[cursor].strip())
            cursor += 1
        routes.append(' '.join(route_lines))
    for index, line in enumerate(lines):
        match = re.search(r'(Input|Standard) orientation:', line)
        if not match:
            continue
        rows, parse_error = [], False
        cursor = index + 1
        # Two divider lines surround the two-line heading.
        for _ in range(2):
            while cursor < len(lines) and not re.match(r'^\s*-{5,}\s*$', lines[cursor]):
                cursor += 1
            cursor += 1
        while cursor < len(lines) and not re.match(r'^\s*-{5,}\s*$', lines[cursor]):
            fields = lines[cursor].split()
            try:
                if len(fields) != 6:
                    raise ValueError('malformed orientation row')
                rows.append({'center_number': int(fields[0]), 'atomic_number': int(fields[1]),
                             'position': [_number(x) for x in fields[3:6]]})
            except ValueError:
                parse_error = True
            cursor += 1
        orientations.append({'kind': match[1], 'line_number': index+1,
                             'atoms': rows, 'parse_error': parse_error or cursor >= len(lines)})
    harmonic_headers = list(re.finditer(r'Harmonic frequencies', text, flags=re.I))
    frequency_text = text[harmonic_headers[-1].start():] if harmonic_headers else text
    frequencies = [_number(x) for row in re.findall(r'^\s*Frequencies\s+--\s*(.*)$', frequency_text, flags=re.M)
                   for x in row.split()]
    hf_energies = [_number(x) for x in re.findall(r'SCF Done:\s+E\((?:RHF|UHF|ROHF)\)\s*=\s*(\S+)', text, flags=re.I)]
    charge_spin = [(int(c), int(m)) for c, m in re.findall(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)', text)]
    error_lines = [line.strip() for line in lines if re.search(
        r'Error termination|Convergence failure|Optimization stopped|SCF has not converged|Error in internal coordinate|Erroneous write', line, flags=re.I)]
    stability_messages = [line for line in lines if re.search(
        r'The wavefunction (?:is stable under the perturbations considered\.|has an? .*instability\.|is unstable)', line, flags=re.I)]
    return {'normal_terminations': len(re.findall(r'Normal termination of Gaussian', text)),
            'error_lines': error_lines, 'optimization_completed': 'Optimization completed.' in text,
            'final_hf_energy_hartree': hf_energies[-1] if hf_energies else None,
            'frequencies_cm1': frequencies,
            'wavefunction_stable': bool(stability_messages and re.search(
                r'The wavefunction is stable under the perturbations considered\.', stability_messages[-1], flags=re.I)),
            'charge_multiplicity_records': charge_spin, 'orientations': orientations,
            'routes': routes,
            'final_orientation': orientations[-1] if orientations else None}


def _topology(record, positions):
    ids = {a['id']: i for i, a in enumerate(record['atoms'])}
    intended = {frozenset((b['a'], b['b'])) for b in record['bonds']}
    lost, short, novel = [], [], []
    for i, a in enumerate(record['atoms']):
        for j in range(i+1, len(record['atoms'])):
            b = record['atoms'][j]
            radii_sum = COVALENT_RADII[a['element']]+COVALENT_RADII[b['element']]
            distance = math.dist(positions[i], positions[j])
            pair = {'atom_ids': [a['id'], b['id']], 'distance_angstrom': distance}
            if frozenset(pair['atom_ids']) in intended:
                if distance > 1.30*radii_sum:
                    lost.append(pair)
                if distance < .65*radii_sum:
                    short.append(pair)
            elif distance < 1.20*radii_sum:
                novel.append(pair)
    if any(b['a'] not in ids or b['b'] not in ids for b in record['bonds']):
        raise ValueError('Intended bond contains unknown atom ID')
    return {'passed': not (lost or short or novel), 'lost_bonds': lost,
            'too_short_bonds': short, 'novel_close_pairs': novel,
            'method': 'covalent-radius distance screen; does not validate bond order',
            'bond_min_radius_factor': .65, 'bond_max_radius_factor': 1.30,
            'novel_pair_radius_factor': 1.20, 'radii_angstrom': COVALENT_RADII}


def validate_gaussian_log(model_record, log):
    """Accept only complete optimization/frequency/stability/identity evidence.

    Negative modes in [−20, 0) cm⁻¹ require manual review and do not auto-pass.
    The fragment catalog is nonlinear: exactly 3N−6 frequencies are required.
    """
    record = _model(model_record)
    _check_model(record)
    parsed = parse_gaussian_log(log)
    issues = []
    routes = parsed['routes']
    relevant = [r for r in routes if re.search(r'\b(?:Opt|Freq|Stable)\b', r, re.I)]
    method = re.compile(r'(?<![A-Za-z])R?HF\s*/\s*6-\s*31G(?:\*|\(d\))', re.I)
    symmetry_off = re.compile(r'\b(?:NoSymm(?:etry)?|Symmetry\s*=\s*None)\b', re.I)
    opt_freq = [r for r in relevant if re.search(r'\bOpt\s*=', r, re.I) and re.search(r'\bFreq\b', r, re.I)]
    stable_routes = [r for r in relevant if re.search(r'\bStable\b', r, re.I)]
    bad_route = False
    preceding_explicit = None
    for route in relevant:
        if not method.search(route):
            bad_route = True
        if symmetry_off.search(route):
            preceding_explicit = route
            continue
        # Gaussian splits an explicit Opt+Freq into an internal checkpoint-based
        # frequency job. Its printed route omits inherited symmetry settings.
        # Exempt only the observed internal signature, backed by a preceding
        # explicit HF Opt+Freq route whose method and symmetry pass independently.
        internal_frequency = (re.search(r'\bFreq\b', route, re.I)
                              and not re.search(r'\b(?:Opt|Stable)\b', route, re.I)
                              and all(re.search(pattern, route, re.I) for pattern in (
                                  r'\bGeom\s*=\s*AllCheck\b', r'\bGuess\s*=\s*TCheck\b',
                                  r'\bSCRF\s*=\s*Check\b', r'\bTest\b', r'\bGenChk\b')))
        parent_opt_freq = (preceding_explicit is not None and method.search(preceding_explicit)
                           and re.search(r'\bOpt\s*=', preceding_explicit, re.I)
                           and re.search(r'\bFreq\b', preceding_explicit, re.I))
        if not internal_frequency or not parent_opt_freq:
            bad_route = True
    if not opt_freq or not stable_routes or bad_route or any(
            re.search(r'\bStable\s*=', r, re.I) for r in stable_routes):
        issues.append('route_provenance_mismatch')
    if parsed['normal_terminations'] < 2:
        issues.append('normal_terminations_missing')
    if parsed['error_lines']:
        issues.append('gaussian_errors_present')
    if not parsed['optimization_completed']:
        issues.append('optimization_not_completed')
    energy = parsed['final_hf_energy_hartree']
    if energy is None or not math.isfinite(energy):
        issues.append('finite_final_hf_energy_missing')
    frequencies = parsed['frequencies_cm1']
    if len(frequencies) != 3*len(record['atoms'])-6:
        issues.append('frequency_count_mismatch')
    if not all(math.isfinite(x) for x in frequencies):
        issues.append('nonfinite_frequencies')
    if any(x < -20 for x in frequencies):
        issues.append('significant_imaginary_frequencies')
    manual_review = any(-20 <= x < 0 for x in frequencies)
    if manual_review:
        issues.append('small_negative_frequencies_manual_review')
    if not parsed['wavefunction_stable']:
        issues.append('stability_not_confirmed')
    expected_charge_spin = (record['formal_charge'], record['multiplicity'])
    if not parsed['charge_multiplicity_records']:
        issues.append('charge_and_multiplicity_evidence_missing')
    elif any(pair != expected_charge_spin for pair in parsed['charge_multiplicity_records']):
        issues.append('charge_or_multiplicity_mismatch')
    final = parsed['final_orientation']
    optimized_atoms, topology = [], None
    if final is None or final['parse_error']:
        issues.append('final_orientation_missing_or_malformed')
    elif len(final['atoms']) != len(record['atoms']) or [a['center_number'] for a in final['atoms']] != list(range(1, len(record['atoms'])+1)) or [a['atomic_number'] for a in final['atoms']] != [ATOMIC_NUMBERS[a['element']] for a in record['atoms']]:
        issues.append('atom_identity_mismatch')
    elif not all(math.isfinite(x) for a in final['atoms'] for x in a['position']):
        issues.append('nonfinite_final_geometry')
    else:
        for original, atom in zip(record['atoms'], final['atoms']):
            optimized = copy.deepcopy(original)
            optimized['position'] = atom['position']
            optimized_atoms.append(optimized)
        topology = _topology(record, [a['position'] for a in optimized_atoms])
        if not topology['passed']:
            issues.append('topology_distance_check_failed')
    return {'accepted': not issues, 'issues': issues, 'manual_review_required': manual_review,
            'model_id': record['model_id'], 'normal_terminations': parsed['normal_terminations'],
            'error_lines': parsed['error_lines'], 'optimization_completed': parsed['optimization_completed'],
            'final_hf_energy_hartree': energy, 'frequencies_cm1': frequencies,
            'wavefunction_stable': parsed['wavefunction_stable'],
            'charge_multiplicity_records': parsed['charge_multiplicity_records'],
            'routes': routes,
            'final_orientation': final, 'optimized_atoms': optimized_atoms, 'topology': topology,
            'limits': ['Distance screen does not validate bond order or stereochemistry.',
                       'Acceptance requires a nonlinear fragment and does not establish global-minimum sampling.']}


def _rotate(position, center, angles):
    x,y,z = [p-c for p,c in zip(position, center)]
    ax,ay,az = [math.radians(a) for a in angles]
    y,z = y*math.cos(ax)-z*math.sin(ax), y*math.sin(ax)+z*math.cos(ax)
    x,z = x*math.cos(ay)+z*math.sin(ay), -x*math.sin(ay)+z*math.cos(ay)
    x,y = x*math.cos(az)-y*math.sin(az), x*math.sin(az)+y*math.cos(az)
    return [x+center[0], y+center[1], z+center[2]]


def inspect_reference(model_record, log):
    """Controller API: validated evidence plus optimized model only if accepted."""
    record = _model(model_record)
    report = validate_gaussian_log(record, log)
    report['optimized_model'] = None
    if report['accepted']:
        optimized = copy.deepcopy(record)
        optimized['atoms'] = report['optimized_atoms']
        optimized.setdefault('metadata', {}).update(
            geometry_status='HF_6-31Gstar_optimized_frequency_and_stability_validated',
            wavefunction_stability_validated=True,
            coordinate_units='angstrom', source_log_sha256=_source_log_hash(log),
            topology_validation='distance_screen_only_not_bond_order_or_stereochemical_validation')
        report['optimized_model'] = optimized
    return report


def _same_element_cloud(a, b, tolerance=1e-7):
    return all(min(math.dist(x['position'], y['position']) for y in b if x['element'] == y['element']) < tolerance for x in a)


def prepare_orientations(model_record, log, output_dir, *, cpus=8, mem_gb=6,
                         allow_rhf_instability_diagnostic=False):
    """Gate, preserve optimized geometry, then prepare three generic MK rotations.

    An existing output directory is rejected. This writes inputs only, never runs
    Gaussian. Non-equivalence checks element-resolved coordinate clouds without
    rotational alignment, since these are orientation controls, not conformers.
    """
    record = _model(model_record)
    validation = inspect_reference(record, log)
    diagnostic = (allow_rhf_instability_diagnostic and
                  validation['issues'] == ['stability_not_confirmed'] and
                  bool(re.search(r'^\s*The wavefunction has an RHF -> UHF instability\.\s*$',
                                 _log_text(log), re.M)))
    if not validation['accepted'] and not diagnostic:
        raise ValueError('Gaussian reference not accepted: '+', '.join(validation['issues']))
    # Validate resources before creating any output directory.
    gaussian_input(record, cpus=cpus, mem_gb=mem_gb)
    if diagnostic:
        optimized = copy.deepcopy(record)
        optimized['atoms'] = validation['optimized_atoms']
        optimized.setdefault('metadata', {}).update(
            geometry_status='RHF_6-31Gstar_constrained_geometry_DIAGNOSTIC_external_instability',
            wavefunction_stability_validated=False,
            scope='rhf_instability_diagnostic_control',
            scientific_status='diagnostic_control_only_not_validated_final_parameters',
            coordinate_units='angstrom', source_log_sha256=_source_log_hash(log),
            topology_validation='distance_screen_only_not_bond_order_or_stereochemical_validation')
    else:
        optimized = validation['optimized_model']
    center = [sum(a['position'][k] for a in optimized['atoms'])/len(optimized['atoms']) for k in range(3)]
    orientations = []
    for index, angles in enumerate(ROTATIONS_DEGREES, 1):
        rotated = copy.deepcopy(optimized)
        for atom in rotated['atoms']:
            atom['position'] = _rotate(atom['position'], center, angles)
        if _same_element_cloud(rotated['atoms'], optimized['atoms']) or any(_same_element_cloud(rotated['atoms'], prior['atoms']) for prior in orientations):
            raise ValueError('Requested rotations are symmetry-equivalent for this geometry')
        rotated['metadata'].update(mk_orientation_index=index, euler_rotation_xyz_degrees=list(angles),
                                   rotation_center_angstrom=center)
        orientations.append(rotated)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=False)
    optimized_text = _json_text(optimized)
    (output/'model.json').write_text(optimized_text)
    # Validation may contain NaN for rejected evidence, but accepted evidence must serialize strictly.
    (output/'validation.json').write_text(_json_text(validation))
    original_bytes = Path(model_record).read_bytes() if isinstance(model_record, (str, Path)) else _json_text(record).encode()
    manifest = {'schema_version': 1, 'status': 'prepared_inputs_only_not_executed',
                'scope': 'rhf_instability_diagnostic_control' if diagnostic else 'stable_reference',
                'reference_accepted': validation['accepted'],
                'final_parameter_validation': False,
                'model_id': record['model_id'], 'source_model_sha256': _hash_bytes(original_bytes),
                'source_log_sha256': _source_log_hash(log), 'optimized_model_sha256': _hash_bytes(optimized_text.encode()),
                'optimized_model_path': str((output/'model.json').resolve()),
                'method': 'HF/6-31G* MK single points; fixed geometry, Symmetry=None', 'orientations': []}
    for index, rotated in enumerate(orientations, 1):
        folder = output/f'orientation_{index:02d}'
        folder.mkdir()
        model_text = _json_text(rotated)
        (folder/'model.json').write_text(model_text)
        lines = [f'%NProcShared={cpus}', f'%Mem={mem_gb}GB', '%Chk=mk.chk',
                 '#P HF/6-31G* Pop=MK SCF=(Tight,XQC) Symmetry=None IOp(6/33=2,6/41=10,6/42=17)', '',
                 f'{record["model_id"]}: rigid MK orientation {index}', '',
                 f'{record["formal_charge"]} {record["multiplicity"]}']
        lines += [f'{a["element"]:2s} '+' '.join(f'{x:.12f}' for x in a['position']) for a in rotated['atoms']]
        input_text = '\n'.join(lines+['', ''])
        (folder/'mk.gjf').write_text(input_text)
        manifest['orientations'].append({'directory': folder.name, 'euler_rotation_xyz_degrees': list(ROTATIONS_DEGREES[index-1]),
                                         'input_path': str((folder/'mk.gjf').resolve()), 'model_path': str((folder/'model.json').resolve()),
                                         'input_sha256': _hash_bytes(input_text.encode()), 'model_sha256': _hash_bytes(model_text.encode())})
    (output/'manifest.json').write_text(_json_text(manifest))
    return manifest


def _safe_report(value):
    """Represent nonfinite parser evidence as null in portable JSON reports."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _safe_report(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_report(item) for item in value]
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    inspect = commands.add_parser('inspect')
    inspect.add_argument('model', type=Path)
    inspect.add_argument('log', type=Path)
    prepare = commands.add_parser('prepare-orientations')
    prepare.add_argument('model', type=Path)
    prepare.add_argument('log', type=Path)
    prepare.add_argument('output_dir', type=Path)
    prepare.add_argument('--cpus', type=int, default=8)
    prepare.add_argument('--mem-gb', type=int, default=6)
    prepare.add_argument('--rhf-diagnostic', action='store_true',
                         help='Explicit constrained-RHF control only; never validates an unstable reference')
    args = parser.parse_args(argv)
    if args.command == 'inspect':
        report = inspect_reference(args.model, args.log)
        print(_json_text(_safe_report(report)), end='')
        return 0 if report['accepted'] else 2
    report = prepare_orientations(args.model, args.log, args.output_dir, cpus=args.cpus, mem_gb=args.mem_gb,
                                  allow_rhf_instability_diagnostic=args.rhf_diagnostic)
    print(_json_text(report), end='')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
