"""Literal-argv Slurm submission; no scientific behavior changes between backends."""

from __future__ import annotations

import hashlib
import json
import re
import shlex
import subprocess
import uuid
from datetime import UTC, datetime
from pathlib import Path


def submit_slurm(
    command: list[str],
    root: str | Path,
    output_dir: str | Path,
    *,
    cpus: int = 1,
    memory: str = "2G",
    time: str = "01:00:00",
    partition: str | None = None,
    account: str | None = None,
    setup: str | Path | None = None,
    dry_run: bool = False,
) -> dict:
    """Write one worker and submit it, or return a concrete unsubmitted plan."""
    if isinstance(cpus, bool) or not isinstance(cpus, int) or cpus < 1:
        raise ValueError("cpus must be a positive integer")
    if not re.fullmatch(r"[1-9][0-9]*[KMGTP]?", memory, re.IGNORECASE):
        raise ValueError("memory must be a positive Slurm size, e.g. 2G")
    if not re.fullmatch(r"(?:[0-9]+-)?[0-9]{1,3}:[0-5][0-9]:[0-5][0-9]", time):
        raise ValueError("time must be HH:MM:SS or D-HH:MM:SS")
    if not command or any(not isinstance(arg, str) or "\0" in arg for arg in command):
        raise ValueError("command must be a nonempty list of literal arguments")
    for value in (partition, account):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_.:,/-]+", value):
            raise ValueError("Invalid account or partition")
    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("Working directory must be a directory")
    if setup is not None:
        setup = Path(setup).resolve(strict=True)
    directory = Path(output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    token = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    worker = directory / f"go_{token}.sh"
    receipt = directory / f"go_{token}.json"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"cd -- {shlex.quote(str(root))}"]
    if setup:
        lines.append(f"source {shlex.quote(str(setup))}")
    lines += [
        f"export OMP_NUM_THREADS={cpus}",
        f"export OPENBLAS_NUM_THREADS={cpus}",
        "exec " + shlex.join(command),
    ]
    script = "\n".join(lines) + "\n"
    worker.write_text(script)
    worker.chmod(0o700)
    scheduler = [
        "sbatch",
        "--parsable",
        "--job-name=GO-builder",
        f"--cpus-per-task={cpus}",
        f"--mem={memory}",
        f"--time={time}",
        f"--output={directory / ('go_' + token + '.out')}",
        f"--error={directory / ('go_' + token + '.err')}",
    ]
    if partition:
        scheduler.append(f"--partition={partition}")
    if account:
        scheduler.append(f"--account={account}")
    scheduler.append(str(worker))
    result = {
        "status": "PLANNED",
        "worker": str(worker),
        "receipt": str(receipt),
        "command": command,
        "scheduler_command": scheduler,
        "worker_sha256": hashlib.sha256(script.encode()).hexdigest(),
        "job_id": None,
    }
    if not dry_run:
        completed = subprocess.run(scheduler, check=True, capture_output=True, text=True)
        job_id = completed.stdout.strip().split(";", 1)[0]
        if not job_id.isdigit():
            raise RuntimeError(f"Unrecognized sbatch receipt: {completed.stdout!r}")
        result.update(status="SUBMITTED", job_id=job_id)
    receipt.write_text(json.dumps(result, indent=2) + "\n")
    return result
