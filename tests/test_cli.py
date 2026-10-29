"""Exercise entry points as a user does, without importing scientific engines."""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_build_help_available_without_external_amber():
    result = subprocess.run(
        [sys.executable, str(ROOT / "01_build_GO.py"), "--help"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "--recipe" in result.stdout
    assert "--slurm" not in result.stdout


def test_missing_recipe_returns_nonzero_without_creating_output(tmp_path):
    output = tmp_path / "go.json"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "01_build_GO.py"),
            "--recipe",
            str(tmp_path / "absent.json"),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert not output.exists()


def test_generic_fragment_cli_preserves_explicit_selected_amine(tmp_path):
    output = tmp_path / "glycine.json"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "fragment",
            "--smiles",
            "NCC(=O)O",
            "--attachment-atom",
            "1",
            "--output",
            str(output),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    fragment = json.loads(output.read_text())["fragment"]
    assert fragment["attachment_atom"] == 1
    assert fragment["atoms"][0]["element"] == "N"
    assert len(fragment["atoms"]) == 10


def test_gnina_cli_plan_keeps_selected_engine_and_formats(tmp_path):
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "dock",
            "--engine",
            "gnina",
            "--gnina",
            "/system/gnina",
            "--boxes",
            str(tmp_path / "boxes.json"),
            "--receptor",
            str(tmp_path / "receptor.pdb"),
            "--ligand",
            str(tmp_path / "ligand.sdf"),
            "--output-dir",
            str(tmp_path / "docking"),
            "--dry-run",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    command = json.loads(result.stdout)["command"]
    assert command[command.index("--engine") + 1] == "gnina"
    assert command[command.index("--gnina") + 1] == "/system/gnina"


def test_go_creation_and_opening_report_actual_chemistry(tmp_path):
    def invoke(*args):
        result = subprocess.run(
            [sys.executable, str(ROOT / "go.py"), *map(str, args), "--json"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return json.loads(result.stdout)

    native = tmp_path / "GO.json"
    grafted = tmp_path / "GO_PEG.json"
    built = invoke(
        "build", "--recipe", ROOT / "examples/GO_for_random_grafting.json", "--output", native
    )
    opened = invoke(
        "functionalize",
        "--model",
        native,
        "--recipe",
        ROOT / "examples/open_epoxide_PEG2_amine.json",
        "--count",
        "2",
        "--seed",
        "17",
        "--output",
        grafted,
    )
    assert built["group_counts"]["epoxide"] == 4
    assert opened["group_counts"]["epoxide"] == 2
    assert opened["group_counts"]["opened_epoxide"] == 2
    assert opened["atoms"] == built["atoms"] + 36
    assert opened["composition"]["N"] == 2
    assert opened["formal_charge_e"] == built["formal_charge_e"] == 0
    assert opened["requires_parameterization"] is True
    assert opened["geometry_status"] == "constructed_unrelaxed"
    assert opened["functionalization"]["mode"] == "random_count"
    assert opened["functionalization"]["actual_fraction"] == 0.5
    assert opened["functionalization"]["seed"] == 17
    checked = invoke("check", "--model", grafted)
    assert checked["graph"]["valid"] is True
    assert checked["gbsa"]["ready"] is False

    sites = opened["functionalization"]["selected_epoxide_ids"]
    selected = invoke(
        "functionalize",
        "--model",
        native,
        "--recipe",
        ROOT / "examples/open_epoxide_PEG2_amine.json",
        "--id",
        *sites,
        "--seed",
        "17",
        "--output",
        tmp_path / "GO_selected.json",
    )
    assert selected["functionalization"]["mode"] == "explicit_sites"
    assert selected["functionalization"]["placed_count"] == 2
    assert selected["functionalization"]["selected_epoxide_ids"] == sites
    assert selected["functionalization"]["actual_fraction"] == 0.5


def test_functionalization_help_exposes_count_and_optional_id_only():
    result = subprocess.run(
        [sys.executable, str(ROOT / "02_open_epoxides.py"), "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--count" in result.stdout and "--id" in result.stdout
    assert "--fraction-epoxides" not in result.stdout
    assert "--epoxide-ids" not in result.stdout
