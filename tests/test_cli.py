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


def test_build_cli_exports_xyz_mol2_pdb_and_check_reports_geometry(tmp_path):
    model_json = tmp_path / "go.json"
    pdb_path = tmp_path / "go.pdb"
    xyz_path = tmp_path / "go.xyz"
    mol2_path = tmp_path / "go.mol2"

    built = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--recipe",
            str(ROOT / "examples/finite_GO.json"),
            "--output",
            str(model_json),
            "--pdb",
            str(pdb_path),
            "--xyz",
            str(xyz_path),
            "--mol2",
            str(mol2_path),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    report = json.loads(built.stdout)
    assert report["status"] == "BUILT"
    assert model_json.is_file()
    assert pdb_path.is_file()
    assert xyz_path.is_file()
    assert mol2_path.is_file()
    assert int(xyz_path.read_text().splitlines()[0]) == report["atoms"]
    assert "@<TRIPOS>MOLECULE" in mol2_path.read_text()

    checked = subprocess.run(
        [sys.executable, str(ROOT / "go.py"), "check", "--model", str(model_json), "--json"],
        capture_output=True,
        text=True,
    )
    assert checked.returncode == 0, checked.stdout + checked.stderr
    chk_report = json.loads(checked.stdout)
    assert chk_report["geometry"]["passed"] is True
    assert chk_report["geometry"]["clash_count"] == 0
    assert chk_report["geometry"]["undercoordinated_atom_ids"] == []
    assert isinstance(chk_report["geometry"]["electron_count"], int)
    assert chk_report["geometry"]["electron_parity"] in ("even", "odd")


def test_build_cli_with_o_over_c_placement_and_group_priority(tmp_path):
    recipe_path = tmp_path / "base_holes.json"
    recipe_path.write_text(
        json.dumps(
            {
                "nx": 8,
                "ny": 8,
                "seed": 17,
                "holes": {"distribution": "small", "percent_vacancies": 5.0},
            }
        )
    )
    out_path = tmp_path / "go_oc.json"
    built = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--recipe",
            str(recipe_path),
            "--output",
            str(out_path),
            "--o-over-c",
            "0.35",
            "--placement",
            "island",
            "--group-priority",
            "OH",
            "epoxide",
            "COOH",
            "COO",
            "--ph",
            "7.0",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    cli_report = json.loads(built.stdout)
    assert cli_report["status"] == "BUILT"
    assert cli_report["hole_summary"]["distribution"] == "small"
    assert abs(cli_report["composition_report"]["actual_O_over_C"] - 0.35) < 0.02
    model = json.loads(out_path.read_text())
    assert model["metadata"]["placement"]["mode"] == "island"
    assert model["metadata"]["composition"]["group_priority"] == ["OH", "epoxide", "COOH", "COO"]
    assert cli_report["group_counts"]["OH"] > cli_report["group_counts"]["epoxide"]
    assert cli_report["group_counts"]["COO"] > cli_report["group_counts"].get("COOH", 0)

    # Random composition mode via CLI
    out_rand = tmp_path / "go_rand.json"
    rand_res = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--recipe",
            str(recipe_path),
            "--output",
            str(out_rand),
            "--c-over-o",
            "3.0",
            "--composition-mode",
            "random",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert rand_res.returncode == 0, rand_res.stdout + rand_res.stderr
    rand_report = json.loads(rand_res.stdout)
    assert rand_report["composition_report"]["composition_mode"] == "random"
    assert abs(rand_report["composition_report"]["actual_C_over_O"] - 3.0) < 0.15

    # Conflict validation fails cleanly without creating output
    conflict_out = tmp_path / "conflict.json"
    conflict = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--recipe",
            str(recipe_path),
            "--output",
            str(conflict_out),
            "--c-over-o",
            "3.0",
            "--o-over-c",
            "0.3",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert conflict.returncode != 0
    assert not conflict_out.exists()


def test_config_toml_json_roundtrip_and_init_templates_loadable(tmp_path):
    from go_builder.config import (
        dump_config,
        load_config,
        render_annotated_json_template,
        render_annotated_toml_template,
    )

    recipe = {
        "name": "cfg_test",
        "lx_nm": 3.0,
        "ly_nm": 3.0,
        "shape": "rectangular",
        "pbc": [False, False, False],
        "seed": 7,
        "holes": {"distribution": "small", "percent_vacancies": 5.0},
        "composition": {
            "o_over_c": 0.3,
            "priority": ["OH", "epoxide", "COOH", "COO"],
            "mode": "balanced",
            "ph": 7.4,
        },
        "placement": {"mode": "island", "n_islands": 2},
    }
    toml_path = tmp_path / "r.toml"
    json_path = tmp_path / "r.json"
    dump_config(recipe, toml_path)
    dump_config(recipe, json_path)
    assert load_config(toml_path) == recipe
    assert load_config(json_path) == recipe

    # Annotated templates must parse and normalize cleanly without _documentation leakage
    tpl_toml = tmp_path / "tpl.toml"
    tpl_json = tmp_path / "tpl.json"
    tpl_toml.write_text(render_annotated_toml_template())
    tpl_json.write_text(render_annotated_json_template())
    loaded_toml = load_config(tpl_toml)
    loaded_json = load_config(tpl_json)
    assert "_documentation" not in loaded_json
    assert loaded_toml["lx_nm"] == 4.0 and loaded_json["lx_nm"] == 4.0
    assert loaded_toml["composition"]["o_over_c"] == 0.35


def test_direct_cli_build_without_recipe_exports_all_formats_and_human_report_and_log(tmp_path):
    out_json = tmp_path / "go.json"
    out_extxyz = tmp_path / "go.extxyz"
    out_mol2 = tmp_path / "go.mol2"
    out_sdf = tmp_path / "go.sdf"
    out_cif = tmp_path / "go.cif"
    out_log = tmp_path / "go.log"
    out_cfg = tmp_path / "used.toml"

    res = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--size",
            "3",
            "3",
            "--oc",
            "0.30",
            "--ph",
            "7.4",
            "--place",
            "island",
            "--holes",
            "small:5",
            "--out",
            str(out_json),
            "--extxyz",
            str(out_extxyz),
            "--mol2",
            str(out_mol2),
            "--sdf",
            str(out_sdf),
            "--cif",
            str(out_cif),
            "--write-config",
            str(out_cfg),
            "--log",
            str(out_log),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    for p in (out_json, out_extxyz, out_mol2, out_sdf, out_cif, out_log, out_cfg):
        assert p.is_file()

    # Human-readable screen output has key sections and no raw JSON dump
    for section in (
        "GO Builder — Build Summary",
        "Geometry & Dimensions",
        "Pores / Holes",
        "Composition & Chemistry",
        "structural graph environments, not XPS intensities",
        "Integrity Gate & QM Diagnostics",
        "Sublattice",
        "Placement",
        "Written Files",
    ):
        assert section in res.stdout

    log_text = out_log.read_text()
    assert "GO Builder — Build Summary" in log_text
    assert "Timestamp (ISO):" in log_text
    assert "Effective Configuration (JSON):" in log_text
    assert '"o_over_c": 0.3' in log_text

    # Direct CLI with --groups and --holes COUNTx works too
    grp_json = tmp_path / "grp.json"
    res_grp = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--nx",
            "4",
            "--ny",
            "4",
            "--groups",
            "OH=3,epoxide=2,COO=1",
            "--holes",
            "small:1x",
            "--out",
            str(grp_json),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert res_grp.returncode == 0, res_grp.stdout + res_grp.stderr
    grp_rep = json.loads(res_grp.stdout)
    assert grp_rep["group_counts"]["OH"] == 3
    assert grp_rep["group_counts"]["epoxide"] == 2
    assert grp_rep["group_counts"]["COO"] == 1
    assert grp_rep["formal_charge_e"] == -1


def test_build_cli_safety_distinct_paths_and_dry_run_co_oc_conflict(tmp_path):
    dup = tmp_path / "same.out"
    res_dup = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--nx",
            "3",
            "--ny",
            "3",
            "--out",
            str(dup),
            "--sdf",
            str(dup),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert res_dup.returncode != 0
    assert "distinct" in res_dup.stdout.lower()
    assert not dup.exists()

    dry_conflict = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--nx",
            "3",
            "--ny",
            "3",
            "--co",
            "3.0",
            "--oc",
            "0.3",
            "--out",
            str(tmp_path / "out.json"),
            "--dry-run",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert dry_conflict.returncode != 0
    assert "mutually exclusive" in dry_conflict.stdout.lower()


def test_example_toml_configs_load_and_build(tmp_path):
    for name in ("go_standard_ph7.toml", "go_island_porous.toml", "go_hex_qm.toml"):
        out = tmp_path / f"{name}.json"
        res = subprocess.run(
            [
                sys.executable,
                str(ROOT / "go.py"),
                "build",
                "--config",
                str(ROOT / "examples/configs" / name),
                "--size",
                "2.5",
                "2.5",
                "--out",
                str(out),
                "--json",
            ],
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0, res.stdout + res.stderr
        assert json.loads(res.stdout)["status"] == "BUILT"


def test_check_strict_exits_nonzero_on_bad_geometry_and_zero_on_healthy(tmp_path):
    healthy = tmp_path / "healthy.json"
    res_build = subprocess.run(
        [
            sys.executable,
            str(ROOT / "go.py"),
            "build",
            "--nx",
            "3",
            "--ny",
            "3",
            "--groups",
            "OH=1",
            "--out",
            str(healthy),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert res_build.returncode == 0, res_build.stdout + res_build.stderr

    # Strict check on healthy model exits 0 and includes chemistry & human summary
    chk_ok_json = subprocess.run(
        [sys.executable, str(ROOT / "go.py"), "check", "--model", str(healthy), "--strict", "--json"],
        capture_output=True,
        text=True,
    )
    assert chk_ok_json.returncode == 0, chk_ok_json.stdout + chk_ok_json.stderr
    ok_rep = json.loads(chk_ok_json.stdout)
    assert ok_rep["geometry"]["passed"] is True
    assert "chemistry" in ok_rep
    assert ok_rep["chemistry"]["qm_diagnostics"]["electron_count"] > 0

    chk_ok_human = subprocess.run(
        [sys.executable, str(ROOT / "go.py"), "check", "--model", str(healthy), "--strict"],
        capture_output=True,
        text=True,
    )
    assert chk_ok_human.returncode == 0
    assert "GO Builder — Check Summary" in chk_ok_human.stdout
    assert "Integrity Gate : PASS" in chk_ok_human.stdout

    # Mutate one C-H bond to 0.80 A -> strict must exit 1
    model = json.loads(healthy.read_text())
    h_atom = next(
        a
        for a in model["atoms"]
        if a["element"] == "H" and a.get("metadata", {}).get("role") == "edge_hydrogen"
    )
    anchor = next(a for a in model["atoms"] if a["id"] == h_atom["metadata"]["anchor_id"])
    dx, dy, dz = (h_atom["position"][k] - anchor["position"][k] for k in range(3))
    cur = (dx * dx + dy * dy + dz * dz) ** 0.5
    scale = 0.80 / cur
    h_atom["position"] = [
        anchor["position"][0] + dx * scale,
        anchor["position"][1] + dy * scale,
        anchor["position"][2] + dz * scale,
    ]
    bad_ch_path = tmp_path / "bad_ch.json"
    bad_ch_path.write_text(json.dumps(model))

    chk_bad_strict = subprocess.run(
        [sys.executable, str(ROOT / "go.py"), "check", "--model", str(bad_ch_path), "--strict", "--json"],
        capture_output=True,
        text=True,
    )
    assert chk_bad_strict.returncode == 1
    bad_rep = json.loads(chk_bad_strict.stdout)
    assert bad_rep["geometry"]["passed"] is False
    assert bad_rep["geometry"]["bond_lengths_valid"] is False

    # Non-strict check still exits 0 while reporting passed=False
    chk_bad_nonstrict = subprocess.run(
        [sys.executable, str(ROOT / "go.py"), "check", "--model", str(bad_ch_path), "--json"],
        capture_output=True,
        text=True,
    )
    assert chk_bad_nonstrict.returncode == 0
    assert json.loads(chk_bad_nonstrict.stdout)["geometry"]["passed"] is False
