"""Regression tests for Blocco 1 — Passo 1e joint audit items and 24x24 performance."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
import tomllib
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from go_builder.build import _orthonormal_edge_frame, build_sheet
from go_builder.config import dumps_toml, load_config
from go_builder.model import (
    _minimum_image_delta,
    check_geometry_integrity,
    compute_chemistry_summary,
    write_cif,
    write_extxyz,
    write_pdb,
    write_sdf,
)

WORKBENCH_ROOT = Path(__file__).resolve().parent.parent


def _parse_cif_and_check(cif_path: Path, model: dict) -> tuple[float, float]:
    txt = cif_path.read_text().splitlines()
    vals: dict[str, float] = {}
    fracs: list[list[float]] = []
    in_loop = False
    for line in txt:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("_cell_"):
            k, v = s.split()
            vals[k] = float(v)
        elif s == "loop_":
            in_loop = True
        elif in_loop and not s.startswith("_atom_site_"):
            parts = s.split()
            if len(parts) >= 5:
                fracs.append([float(parts[2]), float(parts[3]), float(parts[4])])

    a, b, c = vals["_cell_length_a"], vals["_cell_length_b"], vals["_cell_length_c"]
    al, be, ga = map(
        math.radians,
        [vals["_cell_angle_alpha"], vals["_cell_angle_beta"], vals["_cell_angle_gamma"]],
    )
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(ga), b * math.sin(ga), 0.0])
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    vc = np.array([cx, cy, cz])
    cell_mat = np.vstack([va, vb, vc])

    cart = np.array(fracs, dtype=float) @ cell_mat
    orig = np.array([at["position"] for at in model["atoms"]], dtype=float)
    pbc = model["pbc"]

    if any(pbc):
        orig_cell = np.array(model["cell"], dtype=float)
        r0 = (-1, 0, 1) if pbc[0] else (0,)
        r1 = (-1, 0, 1) if pbc[1] else (0,)
        r2 = (-1, 0, 1) if pbc[2] else (0,)
        shifts = [(sx, sy, sz) for sx in r0 for sy in r1 for sz in r2]
        max_err = 0.0
        for i in range(len(orig)):
            for j in range(i + 1, len(orig)):
                best_c = float("inf")
                best_o = float("inf")
                for sx, sy, sz in shifts:
                    sh_c = sx * va + sy * vb + sz * vc
                    best_c = min(best_c, float(np.linalg.norm(cart[j] - cart[i] + sh_c)))
                    sh_o = sx * orig_cell[0] + sy * orig_cell[1] + sz * orig_cell[2]
                    best_o = min(best_o, float(np.linalg.norm(orig[j] - orig[i] + sh_o)))
                max_err = max(max_err, abs(best_c - best_o))
    else:
        d_cif = np.linalg.norm(cart[:, None, :] - cart[None, :, :], axis=-1)
        d_orig = np.linalg.norm(orig[:, None, :] - orig[None, :, :], axis=-1)
        max_err = float(np.max(np.abs(d_cif - d_orig)))

    min_rep = float("inf")
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                if (sx, sy, sz) == (0, 0, 0):
                    continue
                if pbc == [True, False, False] and sy == 0 and sz == 0:
                    continue
                shift = sx * va + sy * vb + sz * vc
                d = np.linalg.norm(cart[:, None, :] - (cart[None, :, :] + shift[None, None, :]), axis=-1)
                min_rep = min(min_rep, float(np.min(d)))
    return max_err, min_rep


def test_prepare_composition_flat_wrapped_ambiguous_and_preflight(tmp_path: Path) -> None:
    out_flat = tmp_path / "flat_comp.json"
    rep_flat = tmp_path / "flat_rep.json"
    res = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            "examples/configs/go_hex_qm.toml",
            "--output",
            str(out_flat),
            "--report",
            str(rep_flat),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    rep_data = json.loads(rep_flat.read_text())
    assert rep_data["scaffold_carbon_atoms"] == 150

    wrapped = tmp_path / "wrapped.json"
    wrapped.write_text(
        json.dumps({"recipe": {"nx": 4, "ny": 4, "pbc": [False, False, False]}, "O_over_C": 0.2})
    )
    out_wrap = tmp_path / "wrap_comp.json"
    res_w = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            str(wrapped),
            "--output",
            str(out_wrap),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res_w.returncode == 0, res_w.stderr

    ambig = tmp_path / "ambig.json"
    ambig.write_text(json.dumps({"recipe": {"nx": 4, "ny": 4}, "nx": 4}))
    res_a = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            str(ambig),
            "--output",
            str(tmp_path / "ambig_out.json"),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res_a.returncode != 0
    assert not (tmp_path / "ambig_out.json").exists()

    # Preflight non-existence and distinctness before any write
    rep_unused = tmp_path / "unused_rep.json"
    res_exists = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            "examples/configs/go_hex_qm.toml",
            "--output",
            str(out_flat),
            "--report",
            str(rep_unused),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res_exists.returncode != 0
    assert not rep_unused.exists()

    same = tmp_path / "same.json"
    res_dup = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            "examples/configs/go_hex_qm.toml",
            "--output",
            str(same),
            "--report",
            str(same),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res_dup.returncode != 0
    assert not same.exists()


def test_write_sdf_rdkit_sanitize_all_cases(tmp_path: Path) -> None:
    for cfg in (
        "examples/configs/go_hex_qm.toml",
        "examples/configs/go_standard_ph7.toml",
        "examples/configs/go_island_porous.toml",
    ):
        model = build_sheet(load_config(WORKBENCH_ROOT / cfg))
        sdf_path = tmp_path / f"{Path(cfg).stem}.sdf"
        write_sdf(model, sdf_path)
        mol = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)[0]
        assert mol is not None
        assert mol.GetNumAtoms() == len(model["atoms"])
        assert mol.GetNumBonds() == len(model["bonds"])

    m_oh = build_sheet({"nx": 4, "ny": 4, "pbc": [False, False, False], "groups": {"OH": 1}, "seed": 1})
    sdf_oh = tmp_path / "oh1.sdf"
    write_sdf(m_oh, sdf_oh)
    mol_oh = Chem.SDMolSupplier(str(sdf_oh), sanitize=True, removeHs=False)[0]
    assert mol_oh is not None
    assert int(mol_oh.GetProp("unpaired_pi_sites")) >= 1

    m_coo = build_sheet({"nx": 4, "ny": 4, "pbc": [False, False, False], "groups": {"COO": 1}, "seed": 1})
    sdf_coo = tmp_path / "coo1.sdf"
    write_sdf(m_coo, sdf_coo)
    mol_coo = Chem.SDMolSupplier(str(sdf_coo), sanitize=True, removeHs=False)[0]
    assert mol_coo is not None
    assert Chem.GetFormalCharge(mol_coo) == -1
    assert mol_coo.GetNumAtoms() == len(m_coo["atoms"])
    assert mol_coo.GetNumBonds() == len(m_coo["bonds"])

    m_big = build_sheet({"nx": 25, "ny": 25, "pbc": [False, False, False], "seed": 1})
    assert len(m_big["atoms"]) > 999
    sdf_big = tmp_path / "big.sdf"
    write_sdf(m_big, sdf_big)
    mol_big = Chem.SDMolSupplier(str(sdf_big), sanitize=True, removeHs=False)[0]
    assert mol_big is not None
    assert mol_big.GetNumAtoms() == len(m_big["atoms"])

    # Empty/whitespace name fallback and newline rejection
    m_empty_name = build_sheet({"nx": 3, "ny": 3, "pbc": [False, False, False], "seed": 1})
    m_empty_name["name"] = "   "
    sdf_empty = tmp_path / "empty_name.sdf"
    write_sdf(m_empty_name, sdf_empty)
    assert sdf_empty.read_text().splitlines()[0] == "GO"

    m_nl = build_sheet({"nx": 3, "ny": 3, "pbc": [False, False, False], "seed": 1})
    m_nl["name"] = "bad\nname"
    with pytest.raises(ValueError):
        write_sdf(m_nl, tmp_path / "bad_nl.sdf")


def test_write_cif_finite_and_ribbon_clearance(tmp_path: Path) -> None:
    specs = [
        ("4x4", {"nx": 4, "ny": 4, "pbc": [False, False, False], "seed": 1}),
        ("4x4_cbx", {"nx": 4, "ny": 4, "pbc": [False, False, False], "groups": {"COOH": 1, "COO": 1}, "seed": 1}),
        ("ribbon_TFF", {"nx": 4, "ny": 4, "pbc": [True, False, False], "groups": {"OH": 1}, "seed": 1}),
    ]
    for label, spec in specs:
        model = build_sheet(spec)
        cif_path = tmp_path / f"{label}.cif"
        write_cif(model, cif_path)
        max_err, min_rep = _parse_cif_and_check(cif_path, model)
        assert max_err < 1e-5
        assert min_rep >= 15.0


def test_chemistry_summary_and_cli_no_oco_doubling_or_spin_overstatement(tmp_path: Path) -> None:
    model = build_sheet({"nx": 4, "ny": 4, "pbc": [False, False, False], "groups": {"COOH": 1, "COO": 1}, "seed": 1})
    chem = compute_chemistry_summary(model)
    envs = chem["chemical_environments"]
    assert envs["counts"]["O-C=O_carboxyl"] == 2
    assert envs["carboxyl_breakdown"] == {"carboxyl_C4": 2, "anchor_CY": 2}
    assert sum(envs["counts"].values()) == chem["atom_counts"]["C"]
    assert "min_spin_multiplicity_hint" not in chem["qm_diagnostics"]

    res = subprocess.run(
        [sys.executable, "01_build_GO.py", "--nx", "4", "--ny", "4", "--groups", "COOH=1,COO=1", "--out", str(tmp_path / "m4.json")],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    assert "O-C=O=2 " in res.stdout
    assert "O-C=O=4" not in res.stdout
    assert "min_mult>=" not in res.stdout


def test_motif_bond_gate_and_skew_cell_minimum_image() -> None:
    cell_skew = np.array([[10.0, 0.0, 0.0], [9.0, 2.0, 0.0], [0.0, 0.0, 40.0]], dtype=float)
    delta = np.array([[9.31, 0.98, 0.0]])
    mi = _minimum_image_delta(delta, [True, True, False], cell_skew, np.linalg.inv(cell_skew))
    assert abs(float(np.linalg.norm(mi)) - math.hypot(0.31, 1.02)) < 1e-12

    m_skew = {
        "schema_version": 1,
        "name": "skew",
        "atoms": [
            {"id": 1, "element": "C", "resname": "MOL", "name": "C1", "position": [0.0, 0.0, 20.0], "formal_charge": 0, "metadata": {}},
            {"id": 2, "element": "C", "resname": "MOL", "name": "C2", "position": [9.31, 0.98, 20.0], "formal_charge": 0, "metadata": {}},
        ],
        "bonds": [],
        "cell": cell_skew.tolist(),
        "pbc": [True, True, False],
        "metadata": {"imported": True},
    }
    rep_skew = check_geometry_integrity(m_skew)
    assert not rep_skew["passed"]
    assert abs(rep_skew["min_nonbonded_heavy_distance"] - math.hypot(0.31, 1.02)) < 1e-6

    m_co = build_sheet({"nx": 4, "ny": 4, "pbc": [False, False, False], "groups": {"carbonyl": 1}, "seed": 1})
    id_to_idx = {a["id"]: idx for idx, a in enumerate(m_co["atoms"])}
    for bond in m_co["bonds"]:
        i, j = id_to_idx[bond["a"]], id_to_idx[bond["b"]]
        if {m_co["atoms"][i]["element"], m_co["atoms"][j]["element"]} == {"C", "O"} and bond.get("order") == 2.0:
            ci = np.array(m_co["atoms"][i]["position"])
            cj = np.array(m_co["atoms"][j]["position"])
            u = (cj - ci) / np.linalg.norm(cj - ci)
            if m_co["atoms"][j]["element"] == "O":
                m_co["atoms"][j]["position"] = (ci + 1.510 * u).tolist()
            else:
                m_co["atoms"][i]["position"] = (cj - 1.510 * u).tolist()
            break
    rep_co = check_geometry_integrity(m_co)
    assert not rep_co["passed"]
    assert rep_co["bond_length_violations"][0]["motif"] == "C=O"

    with pytest.raises(ValueError):
        build_sheet({"nx": 4, "ny": 4, "bond_length": 1.25})


def test_hole_summary_physical_openings() -> None:
    m7 = build_sheet(
        {
            "nx": 7,
            "ny": 7,
            "pbc": [False, False, False],
            "seed": 2,
            "holes": [{"distribution": "small", "count": 7, "min_size": 1, "max_size": 1}],
        }
    )
    hs = m7["metadata"]["hole_summary"]
    assert hs["excavation_events"] == 7
    assert hs["deleted_atom_clusters"] == 7
    assert hs["final_pore_count"] == 6
    assert hs["internal_pores"] == 6
    assert hs["outer_edge_notches"] == 0
    assert sorted(hs["opening_ring_sizes"]) == [12, 12, 12, 12, 12, 18]


def test_exact_epoxide_matching_and_composition_ph_defaults() -> None:
    m3 = build_sheet(
        {
            "nx": 3,
            "ny": 3,
            "pbc": [False, False, False],
            "seed": 1,
            "holes": [{"distribution": "small", "count": 1}],
            "groups": {"epoxide": {"percent": 100, "basis": "eligible"}},
        }
    )
    assert m3["metadata"]["resolved_group_counts"]["epoxide"] == 2

    m_ph = build_sheet(
        {"nx": 4, "ny": 4, "pbc": [False, False, False], "composition": {"o_over_c": 0.35, "ph": 7.4}, "seed": 1}
    )
    gc_ph = m_ph["metadata"]["composition"]["group_counts"]
    assert gc_ph["COOH"] + gc_ph["COO"] > 0

    m_en = build_sheet(
        {
            "nx": 4,
            "ny": 4,
            "pbc": [False, False, False],
            "composition": {
                "o_over_c": 0.30,
                "ph": 7.4,
                "priority": ["OH", "carboxyl"],
                "enabled_groups": ["OH", "COOH", "COO"],
            },
            "seed": 1,
        }
    )
    gc_en = m_en["metadata"]["composition"]["group_counts"]
    assert gc_en["COOH"] + gc_en["COO"] > 0


def test_export_config_cli_provenance_polish(tmp_path: Path) -> None:
    f1, f2, f3 = _orthonormal_edge_frame([1e-10, 0.0, 1.0])
    mat = np.vstack([f1, f2, f3])
    assert float(np.max(np.abs(mat @ mat.T - np.eye(3)))) < 1e-14

    m_ex = build_sheet({"nx": 3, "ny": 3, "pbc": [False, False, False], "groups": {"OH": 1, "COOH": 1}, "seed": 1})
    m_ex["name"] = 'quote "test" slash\\'
    ep = tmp_path / "test.extxyz"
    write_extxyz(m_ex, ep)
    comment = ep.read_text().splitlines()[1]
    assert 'quote \\"test\\" slash\\\\' in comment

    pp = tmp_path / "test.pdb"
    write_pdb(m_ex, pp)
    for line in pp.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name_field = line[12:16]
            el = line[76:78].strip()
            if len(el) == 1 and len(name_field.strip()) < 4:
                assert name_field[0] == " "

    extless = tmp_path / "legacy_recipe"
    extless.write_text(json.dumps({"nx": 4, "ny": 4}))
    assert load_config(extless)["nx"] == 4
    with pytest.raises(ValueError):
        dumps_toml({"a": [1, None, 2]})
    assert tomllib.loads(dumps_toml({"ctrl": "x\x7fy"}))["ctrl"] == "x\x7fy"

    out_m = tmp_path / "m.json"
    r_bad = subprocess.run(
        [sys.executable, "01_build_GO.py", "--preset", "template-go-ph7", "--out", str(out_m), "--write-config", str(tmp_path / "bad.ext")],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_bad.returncode != 0
    assert not out_m.exists()

    r_dry = subprocess.run(
        [sys.executable, "01_build_GO.py", "--preset", "template-go-ph7", "--dry-run", "--write-config", str(tmp_path / "dry.toml")],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_dry.returncode == 0
    assert not (tmp_path / "dry.toml").exists()

    r_prov = subprocess.run(
        [
            sys.executable,
            "01_build_GO.py",
            "--preset",
            "template-go-ph7",
            "--size",
            "2.5",
            "2.5",
            "--oc",
            "0.15",
            "--out",
            str(tmp_path / "prov.json"),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_prov.returncode == 0, r_prov.stderr
    prov_m = json.loads((tmp_path / "prov.json").read_text())
    prov = prov_m["metadata"]["preset"]
    assert prov["constructed_composition"]["lx_nm"] == 2.5
    assert prov["constructed_composition"]["composition"]["o_over_c"] == 0.15
    assert prov["cli_overrides"]["lx_nm"] == 2.5
    assert "scaffold_diameter_nm" in prov["realized_composition"]

    r_s = subprocess.run(
        [sys.executable, "-S", "01_build_GO.py", "--list-presets"],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_s.returncode == 0


def test_finite_24x24_high_density_performance() -> None:
    t0 = time.perf_counter()
    model = build_sheet(
        {
            "nx": 24,
            "ny": 24,
            "pbc": [False, False, False],
            "composition": {"o_over_c": 0.35, "priority": ["OH", "epoxide", "COOH", "COO"], "ph": 7.4},
            "seed": 1,
        }
    )
    dt = time.perf_counter() - t0
    report = check_geometry_integrity(model)
    assert report["passed"], report
    assert dt < 3.0


def test_bond_length_endpoints_and_rejection() -> None:
    for bl in (1.30, 1.55):
        model = build_sheet({"nx": 4, "ny": 4, "pbc": [False, False, False], "bond_length": bl, "seed": 1})
        rep = check_geometry_integrity(model)
        assert rep["passed"], f"bond_length={bl} failed: {rep['bond_length_violations']}"
    for bad_bl in (1.299, 1.551):
        with pytest.raises(ValueError):
            build_sheet({"nx": 4, "ny": 4, "bond_length": bad_bl})


def test_dangling_symlink_preflight_and_periodic_pdb_preflight(tmp_path: Path) -> None:
    dangling_cfg = tmp_path / "dangling_cfg.toml"
    dangling_cfg.symlink_to(tmp_path / "missing_target.toml")
    out_model = tmp_path / "model.json"
    r_cli = subprocess.run(
        [
            sys.executable,
            "01_build_GO.py",
            "--nx",
            "3",
            "--ny",
            "3",
            "--out",
            str(out_model),
            "--write-config",
            str(dangling_cfg),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_cli.returncode != 0
    assert not out_model.exists()

    dangling_rep = tmp_path / "dangling_rep.json"
    dangling_rep.symlink_to(tmp_path / "missing_rep.json")
    out_rec = tmp_path / "rec.json"
    r_prep = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            "examples/configs/go_hex_qm.toml",
            "--output",
            str(out_rec),
            "--report",
            str(dangling_rep),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_prep.returncode != 0
    assert not out_rec.exists()

    # Periodic --pdb preflight in prepare_composition must not leave partial JSON outputs
    per_cfg = tmp_path / "periodic_comp.json"
    per_cfg.write_text(
        json.dumps(
            {
                "nx": 4,
                "ny": 4,
                "pbc": [True, True, False],
                "seed": 0,
                "composition": {"o_over_c": 0.25, "priority": ["OH", "epoxide"]},
            }
        )
    )
    per_out = tmp_path / "per_out.json"
    per_rep = tmp_path / "per_rep.json"
    per_pdb = tmp_path / "per_out.pdb"
    r_per = subprocess.run(
        [
            sys.executable,
            "prepare_composition.py",
            "--input",
            str(per_cfg),
            "--output",
            str(per_out),
            "--report",
            str(per_rep),
            "--pdb",
            str(per_pdb),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert r_per.returncode != 0
    assert not per_out.exists()
    assert not per_rep.exists()
    assert not per_pdb.exists()


def test_rotated_ribbon_cif_cartn_matches_fractional_frame(tmp_path: Path) -> None:
    model = build_sheet({"nx": 4, "ny": 4, "pbc": [True, False, False], "groups": {"OH": 1}, "seed": 1})
    theta = 0.42
    rot = np.array(
        [
            [math.cos(theta), -math.sin(theta), 0.0],
            [math.sin(theta), math.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    model["cell"] = (np.asarray(model["cell"], dtype=float) @ rot.T).tolist()
    for atom in model["atoms"]:
        atom["position"] = (np.asarray(atom["position"], dtype=float) @ rot.T).tolist()

    cif_path = tmp_path / "rotated_ribbon.cif"
    write_cif(model, cif_path)
    lines = cif_path.read_text().splitlines()
    params = {}
    frac_rows = []
    cart_rows = []
    in_loop = False
    for line in lines:
        if line.startswith("_cell_"):
            k, v = line.split()[:2]
            params[k] = float(v)
        elif line.strip() == "_atom_site_Cartn_z":
            in_loop = True
        elif in_loop and line.strip():
            parts = line.split()
            frac_rows.append([float(parts[2]), float(parts[3]), float(parts[4])])
            cart_rows.append([float(parts[5]), float(parts[6]), float(parts[7])])

    a, b, c = params["_cell_length_a"], params["_cell_length_b"], params["_cell_length_c"]
    al, be, ga = (
        math.radians(params["_cell_angle_alpha"]),
        math.radians(params["_cell_angle_beta"]),
        math.radians(params["_cell_angle_gamma"]),
    )
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(ga), b * math.sin(ga), 0.0])
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    std_cell = np.vstack([va, vb, np.array([cx, cy, cz])])
    reconstructed = np.asarray(frac_rows, dtype=float) @ std_cell
    written_cart = np.asarray(cart_rows, dtype=float)
    assert float(np.max(np.abs(reconstructed - written_cart))) < 1e-6


def test_mixed_holes_order_invariance() -> None:
    prof = {"distribution": "small", "count": 3, "min_size": 1, "max_size": 1}
    leg = {"size": 1, "seed_atom": 17}
    m_a = build_sheet({"nx": 9, "ny": 9, "pbc": [True, True, False], "seed": 14, "holes": [prof, leg]})
    m_b = build_sheet({"nx": 9, "ny": 9, "pbc": [True, True, False], "seed": 14, "holes": [leg, prof]})
    assert m_a["metadata"]["removed_carbon_ids"] == m_b["metadata"]["removed_carbon_ids"]
    assert m_a["metadata"]["hole_summary"] == m_b["metadata"]["hole_summary"]


def test_cli_shape_override_clears_inherited_edge_type(tmp_path: Path) -> None:
    out_circ = tmp_path / "circ.json"
    res = subprocess.run(
        [
            sys.executable,
            "01_build_GO.py",
            "--config",
            "examples/configs/go_hex_qm.toml",
            "--shape",
            "circular",
            "--out",
            str(out_circ),
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    m = json.loads(out_circ.read_text())
    assert m["metadata"]["morphology"]["shape"] == "circular"


def test_scaffold_diameter_convex_hull_matches_bruteforce() -> None:
    from go_builder.model import _compute_scaffold_diameter_nm

    for recipe in (
        {"nx": 8, "ny": 8, "pbc": [False, False, False], "seed": 0},
        {"nx": 8, "ny": 8, "pbc": [False, False, False], "groups": {"OH": 15, "epoxide": 10}, "seed": 0},
    ):
        model = build_sheet(recipe)
        sc_ids = {a["id"] for a in model["atoms"] if a.get("metadata", {}).get("role") == "sheet_carbon"}
        hull_d = _compute_scaffold_diameter_nm(model["atoms"], sc_ids)
        sc = np.asarray([a["position"] for a in model["atoms"] if a["id"] in sc_ids], dtype=float)
        diff = sc[:, None, :] - sc[None, :, :]
        brute_d = round(float(np.max(np.linalg.norm(diff, axis=-1))) / 10.0, 6)
        assert abs(hull_d - brute_d) < 1e-6

        # Rigid translation by 1e8 A on all axes must leave diameter invariant
        shifted_atoms = [
            {**a, "position": [c + 1e8 for c in a["position"]]} for a in model["atoms"]
        ]
        shifted_d = _compute_scaffold_diameter_nm(shifted_atoms, sc_ids)
        assert abs(shifted_d - hull_d) < 1e-6

    # 512-point planar circle (all 512 points on convex hull, exact rotating-calipers O(N log N))
    angles = np.linspace(0.0, 2.0 * math.pi, 512, endpoint=False)
    circle_atoms = [
        {
            "id": i + 1,
            "element": "C",
            "position": [25.0 * math.cos(float(theta)), 25.0 * math.sin(float(theta)), 0.0],
        }
        for i, theta in enumerate(angles)
    ]
    circle_d = _compute_scaffold_diameter_nm(circle_atoms, {a["id"] for a in circle_atoms})
    assert abs(circle_d - 5.0) < 1e-6
