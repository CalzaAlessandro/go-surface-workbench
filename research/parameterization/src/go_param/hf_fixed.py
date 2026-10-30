"""HF electronic-state probe at an independently accepted DFT geometry.

Preparation and evidence validation only. No optimization, RESP fitting or
claim that the HF wavefunction was used to obtain the DFT geometry.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
import re

from .protocols import DFT_METHOD, HF_METHOD, route_matches_method
from .qm import (ATOMIC_NUMBERS, _check_model, _json_text, _log_text, _model,
                 _source_log_hash, inspect_reference, parse_gaussian_log)


def prepare_probe(model_record, dft_log, output_dir, *, cpus=4, mem_gb=6):
    """Write a fixed-geometry closed-shell HF probe after the strict DFT gate."""
    record = _model(model_record)
    _check_model(record)
    if record['multiplicity'] != 1:
        raise ValueError('This pilot supports closed-shell singlet references only')
    if any(type(value) is not int or value < 1 for value in (cpus, mem_gb)):
        raise ValueError('CPU and memory requests must be positive integers')
    gate = inspect_reference(record, dft_log, method=DFT_METHOD)
    if not gate['accepted']:
        raise ValueError('DFT parent rejected: '+', '.join(gate['issues']))
    geometry = gate['optimized_model']
    lines = [f'%NProcShared={cpus}', f'%Mem={mem_gb}GB', '%Chk=hf_probe.chk',
             '#P HF/6-31G* Stable SCF=(Tight,XQC) Symmetry=None', '',
             geometry['model_id']+': HF electronic stability at fixed accepted DFT geometry', '',
             f"{geometry['formal_charge']} {geometry['multiplicity']}"]
    lines += [f"{a['element']:2s} "+' '.join(f'{x:.12f}' for x in a['position'])
              for a in geometry['atoms']]
    input_text = '\n'.join(lines+['', ''])
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=False)
    (output/'geometry_model.json').write_text(_json_text(geometry))
    (output/'DFT_GEOMETRY_GATE.json').write_text(_json_text(gate))
    (output/'hf_probe.gjf').write_text(input_text)
    result = {'status': 'PREPARED_INPUT_ONLY', 'geometry_method': DFT_METHOD,
              'probe_method': HF_METHOD, 'hf_reference_accepted': False,
              'geometry_qm_log_sha256': _source_log_hash(dft_log),
              'input_sha256': hashlib.sha256(input_text.encode()).hexdigest(),
              'final_parameter_validation': False}
    (output/'PROBE_INPUT.json').write_text(_json_text(result))
    return result


def validate_probe(model_record, hf_log):
    """Validate only fixed-geometry RHF evidence, independently of DFT gates."""
    model = _model(model_record)
    _check_model(model)
    text = _log_text(hf_log)
    parsed = parse_gaussian_log(text)
    issues = []
    routes = parsed['routes']
    if (model['multiplicity'] != 1 or not routes
            or any(not route_matches_method(r, HF_METHOD) for r in routes)
            or any(re.search(r'\b(?:Opt|Freq)\b|\bStable\s*=', r, re.I) for r in routes)
            or not any(re.search(r'\bStable\b', r, re.I)
                       and re.search(r'\b(?:NoSymm(?:etry)?|Symmetry\s*=\s*None)\b', r, re.I)
                       for r in routes)):
        issues.append('hf_route_mismatch')
    if (parsed['error_lines'] or parsed['normal_terminations'] < 1
            or not re.search(r'Normal termination of Gaussian[^\n]*\s*\Z', text)):
        issues.append('hf_incomplete_or_error')
    energies = parsed['scf_energies']
    if not energies or len(energies) != len(re.findall(r'SCF Done:', text, re.I)):
        issues.append('hf_energy_missing_or_malformed')
    if any(e['method'].upper() != 'RHF' for e in energies):
        issues.append('hf_energy_method_mismatch')
    if any(not math.isfinite(e['energy_hartree']) for e in energies):
        issues.append('nonfinite_hf_energy')
    if not parsed['wavefunction_stable']:
        issues.append('hf_stability_not_confirmed')
    pairs = parsed['charge_multiplicity_records']
    if not pairs or any(pair != (model['formal_charge'], model['multiplicity']) for pair in pairs):
        issues.append('charge_or_multiplicity_mismatch')
    final = parsed['final_orientation']
    if (final is None or final['parse_error'] or len(final['atoms']) != len(model['atoms'])
            or [a['center_number'] for a in final['atoms']] != [a['id'] for a in model['atoms']]
            or [a['atomic_number'] for a in final['atoms']] != [ATOMIC_NUMBERS[a['element']] for a in model['atoms']]):
        issues.append('atom_identity_mismatch')
    elif any(not all(math.isfinite(v) for v in actual['position'])
             or math.dist(expected['position'], actual['position']) > 2e-5
             for expected, actual in zip(model['atoms'], final['atoms'])):
        issues.append('fixed_geometry_changed')
    return {'accepted': not issues, 'issues': issues, 'model_id': model['model_id'],
            'probe_method': HF_METHOD, 'hf_wavefunction_stable': parsed['wavefunction_stable'],
            'hf_geometry_optimization_validated': False, 'scf_energies': energies,
            'geometry_source_method': model.get('metadata', {}).get('qm_method'),
            'source_hf_log_sha256': _source_log_hash(hf_log),
            'scope': 'fixed-geometry closed-shell HF electronic-state probe only',
            'final_parameter_validation': False}
