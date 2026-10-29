"""Amber import must survive the actual JSON/CLI boundary, not only ParmEd."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from go_builder.build import build_sheet
from go_builder.model import save_model, validate_model

ROOT = Path(__file__).resolve().parents[1]


def test_real_old_amber_import_can_be_saved_and_checked(tmp_path):
    output = tmp_path / "imported.json"
    command = [
        sys.executable,
        str(ROOT / "import_Amber.py"),
        "--topology",
        str(ROOT / "resources/reference/old/gox.parm7"),
        "--coordinates",
        str(ROOT / "resources/reference/old/gox.rst7"),
        "--output",
        str(output),
        "--json",
    ]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    model = json.loads(output.read_text())
    assert model["cell"] is None  # No invented vacuum or solvent cell.
    report = validate_model(model)
    assert report["atom_count"] == 940
    assert report["electron_count"] is None  # Amber does not encode formal charges.


def test_generic_imported_box_allows_three_dimensional_pbc_and_ion():
    model = {
        "schema_version": 1,
        "name": "Na",
        "atoms": [{"id": 1, "element": "Na", "position": [0, 0, 0], "formal_charge": 1}],
        "bonds": [],
        "cell": [[20, 0, 0], [0, 20, 0], [0, 0, 20]],
        "pbc": [True, True, True],
        "metadata": {"representation": "amber_import"},
    }
    assert validate_model(model)["valid"]


def test_hydrogen_double_bond_is_not_a_valid_constructed_graph():
    model = build_sheet({"nx": 2, "ny": 2})
    hydrogen = next(a["id"] for a in model["atoms"] if a["element"] == "H")
    next(b for b in model["bonds"] if hydrogen in (b["a"], b["b"]))["order"] = 2
    with pytest.raises(ValueError, match="Hydrogen"):
        validate_model(model)


def test_save_creates_missing_output_parent_without_overwriting(tmp_path):
    destination = tmp_path / "new" / "nested" / "graph.json"
    model = build_sheet({"nx": 2, "ny": 2})
    save_model(model, destination)
    assert destination.is_file()
    with pytest.raises(FileExistsError):
        save_model(model, destination)
