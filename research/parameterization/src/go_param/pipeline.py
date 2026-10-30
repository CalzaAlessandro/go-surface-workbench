"""Gated Gaussian orientations, atomic RESP and matched-grid charge controls.

Runs locally by default; --slurm submits the identical worker. Gaussian and
AmberTools remain external software. Every invocation uses a new output path.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys

from .compare import CANDIDATES, compare_esp
from .protocols import HF_METHOD, DFT_METHOD, RESP_METHODS, require_method
from .qm import inspect_reference, prepare_orientations
from .resp import BOHR_ANGSTROM, fit_atomic, read_esp


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def _portable_evidence(value):
    # Preserve the rejection reasons while encoding malformed numeric evidence as
    # portable nulls. This sanitizer is used only for gate evidence, never scores.
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _portable_evidence(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_portable_evidence(item) for item in value]
    return value


def run(model_path, log_path, output_dir, *, cpus=8, mem_gb=6, gaussian='g16', rhf_diagnostic=False, method=HF_METHOD):
    require_method(method)
    if rhf_diagnostic and method != HF_METHOD:
        raise ValueError('RHF instability diagnostic is HF-only (HF/6-31G*)')
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    try:
        gate = inspect_reference(model_path, log_path, method)
        write_json(output/'reference_gate.json', _portable_evidence(gate))
        diagnostic_bypass = (rhf_diagnostic and gate['issues'] == ['stability_not_confirmed']
                             and 'The wavefunction has an RHF -> UHF instability.' in
                             [line.strip() for line in Path(log_path).read_text().splitlines()])
        if not gate['accepted'] and not diagnostic_bypass:
            raise ValueError('Reference rejected: '+', '.join(gate['issues']))
        diagnostic = bool(rhf_diagnostic)
        scope = {'rhf_diagnostic_requested': bool(rhf_diagnostic),
                 'reference_accepted': gate['accepted'], 'diagnostic_rejection_bypass': bool(diagnostic_bypass),
                 'scope': ('rhf_instability_diagnostic_control' if diagnostic else
                           'dft_resp_control' if method == DFT_METHOD else 'strict_atomic_reference'),
                 'qm_method': method,
                 'final_parameter_validation': False}
        write_json(output/'execution_scope.json', scope)
        manifest = prepare_orientations(model_path, log_path, output/'orientations', cpus=cpus, mem_gb=mem_gb,
                                        allow_rhf_instability_diagnostic=bool(diagnostic_bypass), method=method)
        logs = []
        for entry in manifest['orientations']:
            folder = Path(entry['input_path']).parent
            log = folder/'gaussian.log'
            with Path(entry['input_path']).open('rb') as inp, log.open('xb') as stdout, (folder/'gaussian.stderr').open('xb') as stderr:
                result = subprocess.run([gaussian], stdin=inp, stdout=stdout, stderr=stderr, cwd=folder, check=False)
            write_json(folder/'execution.json', {'exit_code': result.returncode, 'input_sha256': entry['input_sha256']})
            if result.returncode:
                raise RuntimeError(f'Gaussian orientation failed ({result.returncode}): {folder}')
            logs.append(log)
        fit = fit_atomic(manifest['optimized_model_path'], logs, output/'fits', method=method)
        fitted = json.loads(Path(fit['charges_json']).read_text())
        candidate = {'candidate': 'new_atomic_resp_reference', 'atoms': fitted['atoms'], 'coverage': {'complete': True}}
        original_model = json.loads(Path(manifest['optimized_model_path']).read_text())
        scores = []
        # Exact charges are retained; tolerance propagates RESP's F10.6 output rounding.
        rounding_tolerance = len(original_model['atoms'])*5.1e-7+2e-6
        for index in range(1, len(logs)+1):
            xyz, points, potential = read_esp(Path(fit['attempt_dir'])/f'orientation_{index:03d}.esp')
            # Use the same espgen-rounded coordinates for every scored model.
            model = copy.deepcopy(original_model)
            for atom, position in zip(model['atoms'], xyz*BOHR_ANGSTROM):
                atom['position'] = position.tolist()
            ref = {'atom_ids': [a['id'] for a in model['atoms']],
                   'atom_positions_angstrom': (xyz*BOHR_ANGSTROM).tolist(),
                   'formal_charge': model['formal_charge'],
                   'points_angstrom': (points*BOHR_ANGSTROM).tolist(),
                   'potential_au': potential.tolist(), 'dipole_origin_angstrom': [0., 0., 0.]}
            result = compare_esp(model, ref, [candidate, *CANDIDATES], charge_tolerance_e=rounding_tolerance)
            fitted_rows = [row for row in result['results'] if row['candidate'] == candidate['candidate']]
            if len(fitted_rows) != 1 or fitted_rows[0]['status'] != 'matched_charge_control':
                raise ValueError('The fitted atomic RESP candidate was not scored on the matched ESP grid')
            result['orientation'] = index
            result['qm_method'] = method
            scores.append(result)
        write_json(output/'charge_comparison.json', scores)
        fields = ['model_id', 'qm_method', 'orientation', 'candidate', 'status', 'reason', 'rmse_au', 'relative_rmse', 'mae_au', 'charge_error_e']
        with (output/'charge_comparison.csv').open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            for score in scores:
                for row in score['results']:
                    writer.writerow(dict(row, model_id=model['model_id'], qm_method=method, orientation=score['orientation']))
        result = {'status': ('RHF_INSTABILITY_DIAGNOSTIC_ATOMIC_RESP_AND_TRAINING_ESP_COMPLETE' if diagnostic
                             else 'DFT_RESP_CONTROL_AND_TRAINING_ESP_CONTROLS_COMPLETE' if method == DFT_METHOD
                             else 'ATOMIC_REFERENCE_AND_TRAINING_ESP_CONTROLS_COMPLETE'),
                  'scope': scope, 'charge_tolerance_e': rounding_tolerance,
                  'model_id': model['model_id'], 'fit': fit, 'comparison_json': str(output/'charge_comparison.json'),
                  'limitations': ['Finite-model atomic fit; transferable block fit not completed.',
                                  'Training ESP comparisons are not held-out adsorption or full force-field validation.']}
        if diagnostic:
            result['limitations'].append('Explicit RHF diagnostic control only; RHF -> UHF instability is not a validated stable reference. No final parameters are published.')
        marker = 'POSTPROCESS_DIAGNOSTIC_COMPLETE.json' if diagnostic else 'POSTPROCESS_COMPLETE.json'
        write_json(output/marker, result)
        return result
    except Exception as exc:
        write_json(output/'POSTPROCESS_FAILED.json', {'exception': type(exc).__name__, 'message': str(exc)})
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=RESP_METHODS, default=HF_METHOD)
    parser.add_argument('--model', required=True)
    parser.add_argument('--log', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--cpus', type=int, default=8)
    parser.add_argument('--mem-gb', type=int, default=6)
    parser.add_argument('--gaussian', default='g16')
    parser.add_argument('--slurm', action='store_true')
    parser.add_argument('--rhf-diagnostic', action='store_true',
                        help='Explicit unstable-RHF diagnostic control; only the exact RHF -> UHF instability alone may bypass rejection')
    parser.add_argument('--partition', default='ciamician_cpu')
    parser.add_argument('--dependency')
    parser.add_argument('--modules', nargs='*', default=['gaussian/16', 'amber/25'])
    args = parser.parse_args(argv)
    if args.rhf_diagnostic and args.method != HF_METHOD:
        parser.error('--rhf-diagnostic is HF-only (HF/6-31G*)')
    if args.cpus < 1 or args.mem_gb < 1:
        parser.error('CPU and memory requests must be positive')
    if args.slurm:
        root = Path(__file__).resolve().parents[2]
        output = Path(args.output_dir).resolve()
        worker = output.with_name(output.name+'_worker.slurm')
        command = ['python', '-m', 'go_param.pipeline', '--model', str(Path(args.model).resolve()),
                   '--log', str(Path(args.log).resolve()), '--output-dir', str(output),
                   '--cpus', str(args.cpus), '--mem-gb', str(args.mem_gb), '--gaussian', args.gaussian,
                   '--method', args.method]
        if args.rhf_diagnostic:
            command += ['--rhf-diagnostic']
        body = ['#!/bin/bash -l', 'set -eo pipefail']
        body += ['module load '+shlex.quote(m) for m in args.modules]
        body += ['set -u', 'export OMP_NUM_THREADS='+str(args.cpus),
                 'export PYTHONPATH='+shlex.quote(str(root/'src')),
                 'export GAUSS_SCRDIR',
                 'GAUSS_SCRDIR=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/go-resp-${SLURM_JOB_ID}.XXXXXX")',
                 shlex.join(command)]
        worker.parent.mkdir(parents=True, exist_ok=True)
        with worker.open('x') as handle:
            handle.write('\n'.join(body)+'\n')
        cmd = ['sbatch', '--parsable', '--partition', args.partition, '--nodes=1', '--ntasks=1',
               '--cpus-per-task', str(args.cpus), '--mem', f'{args.mem_gb+2}G', '--time=02:00:00',
               '--job-name=GO-RESP', '--output', str(worker)+'.out', '--error', str(worker)+'.err']
        if args.dependency:
            cmd += ['--dependency', args.dependency]
        receipt = subprocess.run([*cmd, str(worker)], capture_output=True, text=True, check=True)
        write_json(worker.with_suffix('.receipt.json'), {'job_id': receipt.stdout.strip(), 'command': cmd+[str(worker)], 'status': 'SUBMITTED'})
        print(receipt.stdout.strip())
    else:
        result = run(args.model, args.log, args.output_dir, cpus=args.cpus, mem_gb=args.mem_gb, gaussian=args.gaussian, rhf_diagnostic=args.rhf_diagnostic, method=args.method)
        print(json.dumps({'status': result['status'], 'model_id': result['model_id']}, indent=2))


if __name__ == '__main__':
    main()
