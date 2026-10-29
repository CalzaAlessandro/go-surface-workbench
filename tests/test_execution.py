"""Slurm plans execute the same argv literally and do not submit during dry runs."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from go_builder.execution import submit_slurm


def test_dry_run_worker_preserves_literal_arguments(tmp_path):
    destination = tmp_path / "result with spaces.json"
    payload = "literal $(touch SHOULD_NOT_EXIST) `echo injected` ; value"
    command = [
        sys.executable,
        "-c",
        "import json,sys;open(sys.argv[1],'w').write(json.dumps(sys.argv[2]))",
        str(destination),
        payload,
    ]
    result = submit_slurm(command, tmp_path, tmp_path / "jobs", dry_run=True)
    assert result["status"] == "PLANNED"
    assert not destination.exists()
    subprocess.run(["bash", result["worker"]], check=True, cwd=tmp_path)
    assert json.loads(destination.read_text()) == payload
    assert not (tmp_path / "SHOULD_NOT_EXIST").exists()


def test_submission_receipt_records_actual_scheduler_id(tmp_path, monkeypatch):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    sbatch = bindir / "sbatch"
    sbatch.write_text("#!/bin/sh\nprintf '12345;testcluster\\n'\n")
    sbatch.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ["PATH"])
    result = submit_slurm([sys.executable, "-c", "pass"], tmp_path, tmp_path / "jobs")
    assert result["status"] == "SUBMITTED"
    assert result["job_id"] == "12345"
    assert json.loads(Path(result["receipt"]).read_text())["job_id"] == "12345"


@pytest.mark.parametrize("kwargs", [{"cpus": 0}, {"memory": "2G\nmalicious"}, {"time": "bad"}])
def test_invalid_resources_fail_before_writing_worker(tmp_path, kwargs):
    with pytest.raises(ValueError):
        submit_slurm(["python3", "-V"], tmp_path, tmp_path / "jobs", dry_run=True, **kwargs)
    assert not (tmp_path / "jobs").exists()
