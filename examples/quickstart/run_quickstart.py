"""Run a minimal, self-contained GO Workbench build, graph check, and Amber template preparation.

Does not require AmberTools, Gaussian, or external engines.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run minimal GO Workbench quickstart build, check, and Amber template generation."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("work_quickstart"),
        help="Target output directory (default: work_quickstart)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    recipe_path = Path("examples/finite_GO.json").resolve()
    output_dir = args.output_dir.resolve()

    if not recipe_path.exists():
        print(f"Error: recipe file not found: {recipe_path}", file=sys.stderr)
        return 1

    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"Output directory must be empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    model_json = output_dir / "finite_GO.json"
    model_pdb = output_dir / "finite_GO.pdb"
    amber_dir = output_dir / "amber"

    print(f"[Quickstart] Building sheet from {recipe_path.name}...")
    build_cmd = [
        sys.executable,
        "01_build_GO.py",
        "--recipe",
        str(recipe_path),
        "--output",
        str(model_json),
        "--pdb",
        str(model_pdb),
        "--json",
    ]
    res_build = subprocess.run(build_cmd, capture_output=True, text=True, check=False)
    if res_build.returncode != 0:
        print(f"Build failed:\n{res_build.stderr}", file=sys.stderr)
        return res_build.returncode

    build_data = json.loads(res_build.stdout)
    assert build_data.get("atoms") == 179, f"Expected 179 atoms, got {build_data.get('atoms')}"
    assert build_data.get("bonds") == 232, f"Expected 232 bonds, got {build_data.get('bonds')}"
    assert build_data.get("formal_charge_e") == -1, f"Expected -1 formal charge, got {build_data.get('formal_charge_e')}"
    print(f"[Quickstart] Sheet constructed: {build_data['atoms']} atoms, {build_data['bonds']} bonds, formal charge {build_data['formal_charge_e']} e")

    print("[Quickstart] Checking molecular graph topology...")
    check_cmd = [
        sys.executable,
        "check_model.py",
        "--model",
        str(model_json),
        "--json",
    ]
    res_check = subprocess.run(check_cmd, capture_output=True, text=True, check=False)
    if res_check.returncode != 0:
        print(f"Check failed:\n{res_check.stderr}", file=sys.stderr)
        return res_check.returncode

    check_data = json.loads(res_check.stdout)
    assert check_data.get("graph", {}).get("valid") is True, "Graph validation failed"
    graph = check_data["graph"]
    assert graph["connected_components"] == 1, graph
    assert graph["undercoordinated_atom_ids"] == [], graph
    print("[Quickstart] Graph check passed: 1 connected component, 0 undercoordinated atoms")

    print("[Quickstart] Generating Amber input templates (without external tleap)...")
    amber_cmd = [
        sys.executable,
        "03_prepare_amber.py",
        "--model",
        str(model_json),
        "--output-dir",
        str(amber_dir),
        "--radii",
        "mbondi2",
    ]
    res_amber = subprocess.run(amber_cmd, capture_output=True, text=True, check=False)
    if res_amber.returncode != 0:
        print(f"Amber preparation failed:\n{res_amber.stderr}", file=sys.stderr)
        return res_amber.returncode

    expected_files = ["tleap.in", "GO_gaff2.frcmod", "assigned_model.json", "manifest.json"]
    for fname in expected_files:
        assert (amber_dir / fname).exists(), f"Missing expected template file: {fname}"

    print("[Quickstart] Success: all minimal outputs verified in", output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
