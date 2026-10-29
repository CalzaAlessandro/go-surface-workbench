"""Built-in recipes remain portable; GO creation is local-only and reproducible."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PRESETS = ("GOX8-P935", "REGOX8-P980", "REGOX8-P995", "REGOX8-P999")


def invoke(*args, cwd=None, no_site=False):
    return subprocess.run(
        [
            sys.executable,
            *(["-S"] if no_site else []),
            str(ROOT / "01_build_GO.py"),
            *map(str, args),
        ],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def test_listing_requires_no_output_or_scientific_dependencies(tmp_path):
    result = invoke("--list-presets", "--json", cwd=tmp_path, no_site=True)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert [p["name"] for p in report["presets"]] == list(PRESETS)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "args",
    [
        ["--recipe", "custom.json", "--preset", PRESETS[0], "--output", "model.json"],
        ["--preset", PRESETS[0]],
        ["--preset", "unknown", "--output", "model.json"],
        ["--list-presets", "--preset", PRESETS[0]],
        ["--preset", PRESETS[0], "--output", "model.json", "--slurm"],
        ["--list-presets", "--slurm"],
    ],
)
def test_invalid_modes_fail_without_output_or_submission(tmp_path, args):
    result = invoke(*args, cwd=tmp_path)
    assert result.returncode != 0
    assert list(tmp_path.iterdir()) == []


def test_help_has_no_scheduler_or_worker_environment_flags():
    result = invoke("--help", no_site=True)
    assert result.returncode == 0
    for option in (
        "--slurm",
        "--submission-dir",
        "--partition",
        "--account",
        "--time",
        "--memory",
        "--cpus",
        "--setup",
        "--python",
    ):
        assert option not in result.stdout
    for option in ("--preset", "--list-presets", "--seed", "--recipe", "--dry-run"):
        assert option in result.stdout


def test_packaged_recipe_load_is_fresh_and_matches_canonical_examples():
    from go_builder.presets import load_preset

    for name in PRESETS:
        recipe = load_preset(name)
        assert recipe == json.loads((ROOT / "examples/layerone" / f"{name}.json").read_text())
        original_seed = recipe["seed"]
        recipe["seed"] = 999999
        assert load_preset(name)["seed"] == original_seed


def test_custom_recipe_seed_is_preserved_unless_explicitly_overridden(tmp_path):
    recipe = tmp_path / "custom.json"
    recipe.write_text(json.dumps({"nx": 3, "ny": 3, "seed": 17, "groups": {"epoxide": 1}}))
    for suffix, extra, seed in (("default", [], 17), ("override", ["--seed", "23"], 23)):
        destination = tmp_path / f"{suffix}.json"
        result = invoke("--recipe", recipe, "--output", destination, "--json", *extra)
        assert result.returncode == 0, result.stdout + result.stderr
        model = json.loads(destination.read_text())
        assert model["metadata"]["seed"] == seed
        assert model["metadata"]["recipe"]["seed"] == seed
        assert "preset" not in model["metadata"]


def test_preset_default_and_override_resolution_preserves_source_provenance():
    from go_builder.presets import load_preset, resolve_preset

    default = load_preset(PRESETS[1])
    recipe, provenance = resolve_preset(PRESETS[1])
    assert recipe == default
    assert provenance["name"] == PRESETS[1]
    assert provenance["source_recipe"]["groups"] == default["groups"]
    assert provenance["resource"] == f"go_builder:data/layerone/{PRESETS[1]}.json"
    recipe, provenance = resolve_preset(PRESETS[1], seed=23)
    assert recipe["seed"] == 23
    assert provenance["source_recipe"]["seed"] == default["seed"]
    assert provenance["seed_override"] == 23
    assert len(provenance["recipe_sha256"]) == 64


def test_local_dry_run_does_not_build_or_create_scheduler_files(tmp_path):
    output = tmp_path / "model.json"
    result = invoke(
        "--preset",
        PRESETS[0],
        "--seed",
        "23",
        "--output",
        output,
        "--dry-run",
        "--json",
        cwd=tmp_path,
        no_site=True,
    )
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report["backend"] == "local"
    assert "--preset" in report["command"] and "--seed" in report["command"]
    assert list(tmp_path.iterdir()) == []


def test_provisional_preset_build_records_separate_source_and_constructed_data(tmp_path):
    output = tmp_path / "rego-p995.json"
    result = invoke(
        "--preset", "REGOX8-P995", "--output", output, "--seed", "23", "--json", cwd=tmp_path
    )
    assert result.returncode == 0, result.stdout + result.stderr
    model = json.loads(output.read_text())
    provenance = model["metadata"]["preset"]
    assert provenance["source_info"]["readiness"] == "provisional_unrelaxed"
    assert provenance["source_info"]["build_enabled"] is True
    assert provenance["source_composition"] == provenance["source_info"]["source_composition"]
    assert provenance["constructed_composition"]["nx"] == 12
    assert provenance["source_recipe"]["seed"] == 17
    assert provenance["seed_override"] == 23
    assert model["metadata"]["seed"] == 23


def test_provisional_unbuildable_preset_fails_clearly_before_outputs(tmp_path):
    result = invoke("--preset", "GOX8-P935", "--output", tmp_path/"GO.json", "--json", cwd=tmp_path)
    assert result.returncode != 0
    report = json.loads(result.stdout)
    assert "placement" in report["message"].lower()
    assert "parameterization" in report["message"].lower()
    assert list(tmp_path.iterdir()) == []


def test_public_source_registry_and_readiness_are_preserved():
    from go_builder.presets import list_presets, resolve_preset

    available = list_presets()
    assert available[0]["build_enabled"] is False
    assert available[0]["readiness"] == "placement_validation_pending"
    _, metadata = resolve_preset("REGOX8-P980")
    assert metadata["source_info"]["build_enabled"] is True
    assert metadata["source_info"]["readiness"] == "provisional_unrelaxed"
    assert len(metadata["source_info_sha256"]) == 64
    assert metadata["source_composition"] == metadata["source_info"]["source_composition"]
    assert "/Users/" not in json.dumps(metadata["source_info"])
