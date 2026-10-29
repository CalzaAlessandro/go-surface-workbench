"""An unmet loading is clearly reported and never saved as a complete system."""

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from go_builder.build import build_sheet
from go_builder.errors import PlacementError
from go_builder.functionalize import open_epoxides

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("machine", [False, True])
def test_classic_group_capacity_warning_has_nonzero_exit_and_no_partial_model(tmp_path, machine):
    recipe = tmp_path / "too_many.json"
    recipe.write_text(json.dumps({"nx": 3, "ny": 3, "groups": {"OH": 9999}}))
    output = tmp_path / "GO.json"
    command = [
        sys.executable,
        str(ROOT / "01_build_GO.py"),
        "--recipe",
        str(recipe),
        "--output",
        str(output),
    ]
    if machine:
        command.append("--json")
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 1
    assert not output.exists()
    if machine:
        report = json.loads(result.stdout)
        assert report["status"] == "INCOMPLETE"
        assert report["requested"]["OH"] == 9999
        assert report["attempts"] == 0
        assert report["capacity"]["eligible_basal_carbons"] < 9999
        assert report["suggestions"]
    else:
        assert "[WARNING]" in result.stderr
        assert "INCOMPLETE" in result.stderr
        assert "suggestions" in result.stderr


def test_fragment_capacity_failure_reports_the_available_rings_without_mutating_input():
    model = build_sheet({"nx": 6, "ny": 6, "groups": {"epoxide": 2}, "seed": 17})
    before = copy.deepcopy(model)
    recipe = json.loads((ROOT / "examples/open_epoxide_PEG2_amine.json").read_text())
    recipe["count"] = 3
    with pytest.raises(PlacementError) as caught:
        open_epoxides(model, recipe)
    report = caught.value.as_report()
    assert report["status"] == "INCOMPLETE"
    assert report["requested"] == {"grafts": 3}
    assert report["capacity"] == {"available_epoxides": 2}
    assert report["attempts"] == 0
    assert model == before


def test_fragment_steric_failure_reports_best_tried_loading_and_search_effort():
    model = build_sheet(json.loads((ROOT / "examples/finite_GO.json").read_text()))
    recipe = json.loads((ROOT / "examples/open_epoxide_PEG2_amine.json").read_text())
    recipe["count"] = 2
    with pytest.raises(PlacementError) as caught:
        open_epoxides(model, recipe)
    report = caught.value.as_report()
    assert report["requested"] == {"grafts": 2}
    assert report["achieved"] == {"grafts": 1}
    assert report["attempts"] > 0
    assert any("conform" in tip.lower() for tip in report["suggestions"])
