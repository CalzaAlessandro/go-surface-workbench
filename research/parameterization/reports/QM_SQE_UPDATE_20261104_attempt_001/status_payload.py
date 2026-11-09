"""Read-only remote snapshot; run with SSH stdin, never a queue watcher."""
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

ROOT = Path('/export/scratch/alessandro.calza3/GO_PARAMETERIZATION_20261029_fqOFQv/calculations')


def command(args):
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(result.stderr)
    return result.stdout


snapshot = {
    'verified_at_utc': datetime.now(timezone.utc).isoformat(),
    'queue': command(['squeue', '-j', '1468880,1468993', '-h', '-o', '%i|%T|%M|%R']),
    'accounting': command(['sacct', '-j', '1468875,1468879,1468880,1468993', '-X', '-P',
                           '--format=JobID,State,ExitCode,Elapsed']),
    'whole_sheet': [],
    'expanded_RESP': {},
}
for folder in sorted((ROOT/'whole_sheet_qm_20261104_attempt_002/cases').iterdir()):
    log = folder/'attempt_001/orca.out'
    text = log.read_text(errors='replace') if log.exists() else ''
    iterations = re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+[-+\deE.]+\s+[-+\deE.]+\s+[-+\deE.]+\s+[-+\deE.]+', text, re.M)
    snapshot['whole_sheet'].append({
        'case': folder.name, 'log_exists': log.exists(),
        'last_SCF_iteration': int(iterations[-1][0]) if iterations else None,
        'last_SCF_energy_Eh': float(iterations[-1][1]) if iterations else None,
        'SCF_convergence_marker_present': 'SCF CONVERGED AFTER' in text,
        'normal_termination_marker_present': 'ORCA TERMINATED NORMALLY' in text,
        'tail': '\n'.join(text.splitlines()[-14:]),
        'scientific_acceptance_not_inferred_from_job_state': True,
    })
base = ROOT/'expanded_resp_panel_20261104_attempt_001'
cases = json.loads((base/'cases.json').read_text())
counts = {'prepared': len(cases), 'reference_logs_started': 0,
          'reference_executions_finished': 0, 'reference_executions_failed': 0, 'HF_RESP_complete': 0,
          'HF_RESP_failed_or_rejected': 0}
details = []
for case in cases:
    source = base/case['source_directory']
    log = source/'reference_attempt_001/gaussian.log'
    counts['reference_logs_started'] += int(log.exists())
    counts['reference_executions_finished'] += int((log.parent/'EXECUTION.json').exists())
    if (log.parent/'EXECUTION.json').exists():
        counts['reference_executions_failed'] += int(json.loads((log.parent/'EXECUTION.json').read_text())['exit_code'] != 0)
    counts['HF_RESP_complete'] += int((source/'HF_RESP_attempt_001/POSTPROCESS_COMPLETE.json').exists())
    counts['HF_RESP_failed_or_rejected'] += int((source/'HF_RESP_attempt_001/POSTPROCESS_FAILED.json').exists())
    if log.exists():
        text = log.read_text(errors='replace')
        details.append({'case': case['model_id'], 'log_bytes': log.stat().st_size,
                        'last_lines': '\n'.join(text.splitlines()[-5:])})
snapshot['expanded_RESP'] = {'counts': counts, 'started_reference_logs': details}
print(json.dumps(snapshot, indent=2, allow_nan=False))
