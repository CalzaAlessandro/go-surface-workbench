"""Tests for Block 2 GAFF2 v3 + L0_v2 force field, provenance, OpenMM/GROMACS export, and CLI."""

from __future__ import annotations

import hashlib
import json
import random
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest

from go_builder.build import build_sheet
from go_builder.ff import DEFAULT_REFERENCE
from go_builder.ff_v3 import (
    DEFAULT_V3_REFERENCE,
    V3_RESIDUE_NAMES,
    assign_gaff2_v3,
    classify_go_v3_roles,
    color_epoxides_v3,
    compute_current_lib_role_charges,
    compute_l0_v2_charges,
    list_force_field_presets,
    load_custom_charges,
    resolve_force_field_preset,
    verify_v3_bundle,
)
from go_builder.topology import export_gromacs, export_openmm, prepare_amber

WORKBENCH_ROOT = Path(__file__).resolve().parents[1]


def _make_graph(elements: list[str], edges: list[tuple[int, int]], *, orders: list[float] | None = None, metadata: dict | None = None) -> dict:
    atoms = [
        {
            "id": i + 1,
            "element": el,
            "name": f"{el}{i + 1}",
            "resname": "UNK",
            "resid": 1,
            "position": [float(i) * 1.4, 0.0, 0.0],
        }
        for i, el in enumerate(elements)
    ]
    bonds = [
        {
            "a": a + 1,
            "b": b + 1,
            "order": 1.0 if orders is None else orders[k],
            "image": [0, 0, 0],
        }
        for k, (a, b) in enumerate(edges)
    ]
    return {"atoms": atoms, "bonds": bonds, "metadata": metadata or {}}


def _scramble_model(model: dict, *, seed: int = 123, kekule_alt: bool = True) -> tuple[dict, list[int]]:
    """Permute atom order, reassign non-contiguous atom IDs, clear labels, and vary Kekule orders."""
    rng = random.Random(seed)
    n = len(model["atoms"])
    perm = list(range(n))
    rng.shuffle(perm)
    new_ids = [100 + 7 * k for k in range(n)]
    rng.shuffle(new_ids)
    orig_id_to_new_idx = {model["atoms"][orig_i]["id"]: new_i for new_i, orig_i in enumerate(perm)}

    new_atoms = []
    for new_i, orig_i in enumerate(perm):
        atom = deepcopy(model["atoms"][orig_i])
        atom["id"] = new_ids[new_i]
        atom["resid"] = 99
        atom["resname"] = "UNK"
        atom["name"] = f"X{new_i}"
        new_atoms.append(atom)

    new_bonds = []
    for k, bond in enumerate(model["bonds"]):
        nb = deepcopy(bond)
        ia = new_atoms[orig_id_to_new_idx[bond["a"]]]["id"]
        ib = new_atoms[orig_id_to_new_idx[bond["b"]]]["id"]
        if k % 2 == 0:
            nb["a"], nb["b"] = ia, ib
        else:
            nb["a"], nb["b"] = ib, ia
        if kekule_alt and nb.get("order") in (1.0, 1.5, 2.0):
            el_a = next(a["element"] for a in new_atoms if a["id"] == nb["a"])
            el_b = next(a["element"] for a in new_atoms if a["id"] == nb["b"])
            if el_a == "C" and el_b == "C":
                nb["order"] = 1.5 if (k % 3 == 0) else (2.0 if k % 3 == 1 else 1.0)
        new_bonds.append(nb)
    rng.shuffle(new_bonds)

    scrambled = deepcopy(model)
    scrambled["atoms"] = new_atoms
    scrambled["bonds"] = new_bonds
    return scrambled, perm


def test_v3_bundle_integrity_and_tamper_detection(tmp_path):
    info = verify_v3_bundle()
    assert len(info["provenance"]["files"]) == 20
    assert Path(info["frcmod_path"]).name == "GO_gaff2_v3.frcmod"
    assert Path(info["l0_v2_path"]).name == "L0_v2_charges.json"
    assert len(V3_RESIDUE_NAMES) == 18

    copied = tmp_path / "gaff2_v3"
    shutil.copytree(DEFAULT_V3_REFERENCE, copied)
    verify_v3_bundle(copied)

    frcmod_copy = copied / "GO_gaff2_v3.frcmod"
    frcmod_copy.write_text(frcmod_copy.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_v3_bundle(copied)

    # Removing both a .lib file and its manifest entry must fail completeness check
    copied2 = tmp_path / "gaff2_v3_incomplete"
    shutil.copytree(DEFAULT_V3_REFERENCE, copied2)
    (copied2 / "EPT_gaff2_v3.lib").unlink()
    prov_file = copied2 / "PROVENANCE.json"
    prov_data = json.loads(prov_file.read_text(encoding="utf-8"))
    prov_data["files"].pop("EPT_gaff2_v3.lib")
    prov_file.write_text(json.dumps(prov_data, indent=2), encoding="utf-8")
    with pytest.raises(ValueError, match="Incomplete gaff2_v3 bundle manifest"):
        verify_v3_bundle(copied2)


def test_force_field_preset_resolution_and_listing():
    presets = list_force_field_presets()
    names = [p["name"] for p in presets]
    assert names == ["gaff2_v3", "gaff2_current"]
    assert presets[0]["default"] is True
    assert presets[1]["default"] is False

    for alias in ("gaff2_v3", "gaff2-v3", "v3", "recommended", "resp-l0-v2", None):
        assert resolve_force_field_preset(alias) == "gaff2_v3"
    for alias in ("gaff2_current", "gaff2-current", "current", "resp-current", "legacy"):
        assert resolve_force_field_preset(alias) == "gaff2_current"
    with pytest.raises(ValueError, match="Unknown force-field preset"):
        resolve_force_field_preset("nonexistent_preset")


def test_priority12_perception_25_roles_and_invariance():
    neutral_sheet = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "shape": "rectangular",
            "seed": 11,
            "groups": {
                "OH": 3,
                "epoxide": 3,
                "epoxide_eph": 1,
                "epoxide_ept": 1,
                "carbonyl": 2,
                "COOH": 2,
            },
        }
    )
    anionic_sheet = build_sheet(
        {
            "nx": 5,
            "ny": 5,
            "shape": "rectangular",
            "seed": 19,
            "groups": {"OH": 2, "epoxide": 2, "COO": 2, "COOH": 1},
        }
    )

    roles_neu = classify_go_v3_roles(neutral_sheet)
    roles_ani = classify_go_v3_roles(anionic_sheet)
    seen_roles = set(roles_neu) | set(roles_ani)
    expected_25 = {
        "scaffold_C",
        "edge_C",
        "edge_H",
        "hydroxyl_C",
        "hydroxyl_O",
        "hydroxyl_H",
        "epoxide_C",
        "epoxide_O",
        "EPH_C",
        "EPH_CH",
        "EPH_H",
        "EPH_O",
        "EPT_C",
        "EPT_H",
        "EPT_O",
        "carbonyl_C",
        "carbonyl_O",
        "COOH_anchor",
        "COOH_C",
        "COOH_O",
        "COOH_OH",
        "COOH_H",
        "COO_anchor",
        "COO_C",
        "COO_O",
    }
    assert seen_roles == expected_25

    assigned_neu = assign_gaff2_v3(neutral_sheet)
    assigned_ani = assign_gaff2_v3(anionic_sheet)
    seen_types = {a["atom_type"] for a in assigned_neu["atoms"] + assigned_ani["atoms"]}
    assert {
        "ga",
        "gb",
        "gh",
        "g6",
        "gm",
        "gw",
        "gp",
        "g1",
        "gk",
        "gq",
        "gc",
        "go",
        "gd",
        "gu",
    }.issubset(seen_types)
    assert any(t in seen_types for t in ("gx", "gy", "gz", "ge"))

    # Permutation + label-clearing + Kekule bond-order invariance
    for sheet, orig_roles in ((neutral_sheet, roles_neu), (anionic_sheet, roles_ani)):
        scrambled, perm = _scramble_model(sheet, seed=77, kekule_alt=True)
        scr_roles = classify_go_v3_roles(scrambled)
        assert scr_roles == [orig_roles[orig_i] for orig_i in perm]


def test_priority3_and_junctions_rejected():
    # 1. Phenolic sp2 OH (benzene + OH)
    phenol = _make_graph(
        ["C", "C", "C", "C", "C", "C", "O", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
            (0, 6), (6, 7), (1, 8), (2, 9), (3, 10), (4, 11), (5, 12),
        ],
    )
    with pytest.raises(ValueError, match="Priority-3"):
        classify_go_v3_roles(phenol)

    # 2. Planar non-epoxide ether (C-O-C without C-C bond)
    ether = _make_graph(
        ["C", "C", "C", "C", "O"],
        [(0, 4), (4, 1), (0, 2), (1, 3)],
    )
    with pytest.raises(ValueError, match="Priority-3"):
        classify_go_v3_roles(ether)

    # 3. Ester / lactone (bare sp2 anchor C - C(=O)-O-C)
    lactone = _make_graph(
        ["C", "C", "C", "C", "O", "O", "C"],
        [(0, 1), (0, 2), (0, 3), (3, 4), (3, 5), (5, 6)],
    )
    with pytest.raises(ValueError, match="Priority-3"):
        classify_go_v3_roles(lactone)

    # 4. Lactol (sp3 C connected to C, H, -OH, -O-C)
    lactol = _make_graph(
        ["C", "C", "H", "O", "H", "O", "C"],
        [(0, 1), (0, 2), (0, 3), (3, 4), (0, 5), (5, 6)],
    )
    with pytest.raises(ValueError, match="Priority-3"):
        classify_go_v3_roles(lactol)

    # 5. Aldehyde (C-CH=O)
    aldehyde = _make_graph(
        ["C", "C", "H", "O"],
        [(0, 1), (1, 2), (1, 3)],
    )
    with pytest.raises(ValueError, match="Priority-3"):
        classify_go_v3_roles(aldehyde)

    # 6. Unparameterized covalent junction
    sheet = build_sheet({"nx": 4, "ny": 4, "seed": 1, "groups": {"OH": 1}})
    flagged = deepcopy(sheet)
    flagged["metadata"]["requires_parameterization"] = True
    with pytest.raises(ValueError, match="junction"):
        assign_gaff2_v3(flagged)


def test_dsatur_four_color_adjacent_epoxides():
    # Build a synthetic coronene-like core with 3 mutually adjacent epoxide pairs (K3 -> needs 3 colors)
    # 6-ring core C0..C5 where (0,1), (2,3), (4,5) carry O6, O7, O8, plus outer C9..C14 carrying third bonds
    # Simpler: two stacked C2 pairs on a honeycomb patch via build_sheet dense epoxides
    dense = build_sheet({"nx": 5, "ny": 5, "seed": 5, "groups": {"epoxide": 10}})
    roles = classify_go_v3_roles(dense)
    atom_color, epx_units = color_epoxides_v3(dense, roles)
    assert len(epx_units) == 10
    assert all(c in (0, 1, 2, 3) for c in atom_color.values())

    # Also build an explicit K3 mutually adjacent epoxide triangle patch (coronene + 3 inner-ring epoxides):
    coronene_atoms = ["C"] * 24 + ["H"] * 12 + ["O"] * 3
    # Inner hexagon 0..5, middle 6 bonds 0-6,1-7,2-8,3-9,4-10,5-11? Standard coronene:
    # Inner 6: 0-5; outer 18: 6..23; H on outer CH 24..35; epoxides on (0,1)=O36, (2,3)=O37, (4,5)=O38.
    cor_edges = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (0, 6), (1, 9), (2, 12), (3, 15), (4, 18), (5, 21),
        (6, 7), (7, 8), (8, 9),
        (9, 10), (10, 11), (11, 12),
        (12, 13), (13, 14), (14, 15),
        (15, 16), (16, 17), (17, 18),
        (18, 19), (19, 20), (20, 21),
        (21, 22), (22, 23), (23, 6),
    ]
    # Add Hydrogens to outer 2-coord carbons: 7,8, 10,11, 13,14, 16,17, 19,20, 22,23
    h_idx = 24
    for oc in (7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23):
        cor_edges.append((oc, h_idx))
        h_idx += 1
    cor_edges.extend([(0, 36), (1, 36), (2, 37), (3, 37), (4, 38), (5, 38)])
    tri_epx = _make_graph(coronene_atoms, cor_edges)
    tri_roles = classify_go_v3_roles(tri_epx)
    tri_colors, _ = color_epoxides_v3(tri_epx, tri_roles)
    # 3 mutually adjacent epoxides MUST receive 3 distinct colors {0, 1, 2}
    assert len({tri_colors[0], tri_colors[2], tri_colors[4]}) == 3
    assigned_tri = assign_gaff2_v3(tri_epx)
    epx_resnames = {a["resname"] for a in assigned_tri["atoms"] if a["resname"].startswith("EP")}
    assert epx_resnames == {"EPX", "EP2", "EP3"}


def test_l0_v2_dual_regime_conservation_and_invariance():
    neutral = build_sheet({"nx": 5, "ny": 5, "seed": 23, "groups": {"OH": 2, "epoxide": 1, "COOH": 1}})
    anionic = build_sheet({"nx": 5, "ny": 5, "seed": 29, "groups": {"OH": 2, "COO": 2, "COOH": 1}})

    q_neu, meta_neu = compute_l0_v2_charges(neutral)
    roles_neu = classify_go_v3_roles(neutral)
    assert meta_neu["regime"] == "neutral_fixed_roles"
    assert abs(float(q_neu.sum())) < 1e-12
    for i, r in enumerate(roles_neu):
        if r == "scaffold_C":
            assert q_neu[i] == 0.0

    q_ani, meta_ani = compute_l0_v2_charges(anionic)
    roles_ani = classify_go_v3_roles(anionic)
    assert meta_ani["regime"] == "anionic_sp2_delocalization"
    assert meta_ani["n_coo"] == 2
    assert abs(float(q_ani.sum()) - (-2.0)) < 1e-12
    for i, r in enumerate(roles_ani):
        if r == "COO_C":
            assert q_ani[i] == pytest.approx(0.835809, abs=1e-6)
        elif r == "COO_O":
            assert q_ani[i] == pytest.approx(-0.799783, abs=1e-6)
    # At least some sp2 acceptor carbons near COO- receive negative delocalized charge
    assert any(q_ani[i] < -1e-4 for i, r in enumerate(roles_ani) if r == "scaffold_C")

    # Atom permutation & Kekule invariance < 1e-12 e
    for sheet, q_orig in ((neutral, q_neu), (anionic, q_ani)):
        scr, perm = _scramble_model(sheet, seed=314, kekule_alt=True)
        q_scr, _ = compute_l0_v2_charges(scr)
        expected = np.array([q_orig[orig_i] for orig_i in perm])
        assert np.allclose(q_scr, expected, atol=1e-12)


def test_custom_charges_and_combinations(tmp_path):
    anionic = build_sheet({"nx": 5, "ny": 5, "seed": 31, "groups": {"OH": 2, "COO": 1}})

    # Combo: gaff2_v3 + current_lib
    q_cur, meta_cur = compute_current_lib_role_charges(anionic)
    assert abs(float(q_cur.sum()) - (-1.0)) < 1e-6
    assert meta_cur["charge_model"] == "current_lib"
    assigned_combo = assign_gaff2_v3(anionic, charge_model="current_lib")
    assert assigned_combo["metadata"]["force_field"]["custom_combination"] is True
    assert assigned_combo["metadata"]["force_field"]["tier"] == "custom_combination"

    # Custom JSON charges_e + provenance
    n = len(anionic["atoms"])
    custom_vec = [0.0] * n
    custom_vec[0] = -1.0
    cj = tmp_path / "custom_q.json"
    cj.write_text(
        json.dumps({"charges_e": custom_vec, "provenance": {"method": "test_custom"}}),
        encoding="utf-8",
    )
    q_c, meta_c = load_custom_charges(cj, anionic)
    assert meta_c["source_sha256"]
    assert meta_c["external_provenance"] == {"method": "test_custom"}
    assert abs(float(q_c.sum()) - (-1.0)) < 1e-12

    # Custom JSON charges_by_id
    by_id = {str(a["id"]): (-1.0 if i == 0 else 0.0) for i, a in enumerate(anionic["atoms"])}
    cj_id = tmp_path / "custom_by_id.json"
    cj_id.write_text(json.dumps({"charges_by_id": by_id}), encoding="utf-8")
    assigned_c = assign_gaff2_v3(anionic, charge_model="custom", custom_charges=cj_id)
    assert assigned_c["metadata"]["force_field"]["charge_model"] == "custom"
    assert abs(sum(a["charge"] for a in assigned_c["atoms"]) - (-1.0)) < 1e-12

    # Bad net charge fails
    bad_cj = tmp_path / "bad_q.json"
    bad_cj.write_text(json.dumps({"charges_e": [0.0] * n}), encoding="utf-8")
    with pytest.raises(ValueError, match="differs from formal charge"):
        load_custom_charges(bad_cj, anionic)

    # Malformed custom charges JSON raises ValueError (not TypeError), and booleans are rejected
    for bad_payload in (
        {"charges_e": None},
        {"charges_e": 1.0},
        {"charges_e": ["bad"] * n},
        {"charges_e": [True] + [0.0] * (n - 2) + [-2.0]},
        {"charges_by_id": None},
        [None],
        [True] + [0.0] * (n - 2) + [-2.0],
    ):
        mal_cj = tmp_path / "malformed_q.json"
        mal_cj.write_text(json.dumps(bad_payload), encoding="utf-8")
        with pytest.raises(ValueError):
            load_custom_charges(mal_cj, anionic)

    # Supplying custom_charges / charges_json without charge_model='custom' must error cleanly
    with pytest.raises(ValueError, match="requires charge_model='custom'"):
        assign_gaff2_v3(anionic, charge_model="l0_v2", custom_charges=cj)
    with pytest.raises(ValueError, match="requires charge_model='custom'"):
        prepare_amber(anionic, tmp_path / "ignored_err_v3", charges_json=cj, run=False)
    with pytest.raises(ValueError, match="requires charge_model='custom'"):
        prepare_amber(
            anionic,
            tmp_path / "ignored_err_cur",
            preset="gaff2_current",
            charge_model="current_lib",
            charges_json=cj,
            run=False,
        )

    # parameter_templates conflicts and clean provenance
    frc_tpl = tmp_path / "junc.frcmod"
    frc_tpl.write_text("remark\n", encoding="utf-8")
    tpl = {
        "atoms": [
            {"id": a["id"], "element": a["element"], "atom_type": "ca", "charge": (-1.0 if i == 0 else 0.0)}
            for i, a in enumerate(anionic["atoms"])
        ],
        "bonds": deepcopy(anionic["bonds"]),
        "frcmod": str(frc_tpl),
        "provenance": {"source": "test_tpl"},
    }
    with pytest.raises(ValueError, match="parameter_templates"):
        prepare_amber(
            anionic,
            tmp_path / "tpl_conflict",
            charge_model="custom",
            charges_json=tmp_path / "nonexistent.json",
            parameter_templates=tpl,
            run=False,
        )
    with pytest.raises(ValueError, match="parameter_templates"):
        assign_gaff2_v3(
            anionic,
            charge_model="l0_v2",
            parameter_templates=tpl,
        )
    res_tpl = prepare_amber(
        anionic,
        tmp_path / "tpl_ok",
        preset="gaff2_current",
        parameter_templates=tpl,
        run=False,
    )
    tpl_prov = json.loads((tmp_path / "tpl_ok/PROVENANCE.json").read_text(encoding="utf-8"))
    assert res_tpl["charge_model"] == "dedicated_parameter_templates"
    assert tpl_prov["charge_model"] == "dedicated_parameter_templates"
    assert tpl_prov["charge_info"]["charge_model"] == "dedicated_parameter_templates"

    # Reapplying parameter_templates on an already-L0_v2-assigned model must cleanly overwrite charge metadata
    pre_assigned = assign_gaff2_v3(anionic)
    res_re_v3 = prepare_amber(
        pre_assigned,
        tmp_path / "tpl_reapply_v3",
        preset="gaff2_v3",
        parameter_templates=tpl,
        run=False,
    )
    re_prov = json.loads((tmp_path / "tpl_reapply_v3/PROVENANCE.json").read_text(encoding="utf-8"))
    re_assigned = json.loads(
        (tmp_path / "tpl_reapply_v3/assigned_model.json").read_text(encoding="utf-8")
    )
    assert res_re_v3["charge_model"] == "dedicated_parameter_templates"
    assert re_prov["charge_model"] == "dedicated_parameter_templates"
    assert re_prov["charge_info"]["charge_model"] == "dedicated_parameter_templates"
    assert all(
        a["metadata"]["charge_source"] == "dedicated_parameter_templates"
        and "v3_role" not in a["metadata"]
        for a in re_assigned["atoms"]
    )

    # Duplicate normalized atom IDs ("1" vs "01") must be rejected instead of silently picking the last value
    dup_by_id = {str(a["id"]): (-1.0 if i == 0 else 0.0) for i, a in enumerate(anionic["atoms"])}
    first_id = str(anionic["atoms"][0]["id"])
    dup_by_id[f"0{first_id}"] = 0.0
    dup_cj = tmp_path / "dup_id.json"
    dup_cj.write_text(json.dumps({"charges_by_id": dup_by_id}), encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate atom ID"):
        load_custom_charges(dup_cj, anionic)
    dup_top = tmp_path / "dup_top.json"
    dup_top.write_text(json.dumps(dup_by_id), encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate atom ID"):
        load_custom_charges(dup_top, anionic)

    # Repeated literal key inside JSON text ("1": 0.9, "1": 0.0) must be rejected at JSON parse time
    raw_dup_json = tmp_path / "raw_dup.json"
    raw_dup_json.write_text(
        f'{{"charges_by_id": {{"{first_id}": 0.9, "{first_id}": -1.0}}}}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate JSON object key"):
        load_custom_charges(raw_dup_json, anionic)


def test_l0_v2_rejects_sp3_only_anion_fallback():
    # Minimal 12-atom valid graph: C0(COO_anchor)-C1(O2,O3) + sp3 epoxides (C4-C5 bridged by O6)
    # + EPT carbons (C7-C8 bridged by O9 + H10, H11) -> zero sp2 carbons in 1 <= d <= d_max.
    sp3_only = _make_graph(
        ["C", "C", "O", "O", "C", "C", "O", "C", "C", "O", "H", "H"],
        [
            (0, 1), (1, 2), (1, 3),
            (0, 4), (0, 5), (4, 5), (4, 6), (5, 6),
            (4, 7), (5, 8), (7, 8), (7, 9), (8, 9), (7, 10), (8, 11),
        ],
    )
    roles = classify_go_v3_roles(sp3_only)
    assert "COO_anchor" in roles
    assert not any(r in ("scaffold_C", "edge_C", "COOH_anchor", "carbonyl_C") for r in roles)
    with pytest.raises(ValueError, match="No eligible sp2 acceptor carbon"):
        compute_l0_v2_charges(sp3_only, roles)


def test_prepare_amber_defaults_and_backward_compatibility(tmp_path):
    sheet = build_sheet({"nx": 5, "ny": 5, "seed": 41, "groups": {"OH": 2, "epoxide": 1, "COOH": 1}})

    # Default preset -> gaff2_v3
    out_v3 = tmp_path / "amber_v3"
    res_v3 = prepare_amber(sheet, out_v3, run=False)
    assert res_v3["preset"] == "gaff2_v3"
    assert res_v3["charge_model"] == "L0_v2c_a030"
    assert res_v3["custom_combination"] is False
    assert (out_v3 / "GO_gaff2_v3.frcmod").is_file()
    assert len(list((out_v3 / "v3_libraries").glob("*_gaff2_v3.lib"))) == 18
    assert (out_v3 / "PROVENANCE.json").is_file()

    # Passing reference_dir without explicit preset preserves legacy gaff2_current behavior
    out_cur = tmp_path / "amber_cur"
    res_cur = prepare_amber(sheet, out_cur, reference_dir=DEFAULT_REFERENCE, run=False)
    assert res_cur["preset"] == "gaff2_current"
    assert res_cur["charge_model"] == "current_lib"
    assert res_cur["custom_combination"] is False
    assert (out_cur / "GO_gaff2.frcmod").is_file()
    assert (out_cur / "PROVENANCE.json").is_file()

    # Cross combo: gaff2_current + l0_v2
    out_combo = tmp_path / "amber_combo"
    res_combo = prepare_amber(sheet, out_combo, preset="gaff2_current", charge_model="l0_v2", run=False)
    assert res_combo["preset"] == "gaff2_current"
    assert res_combo["charge_model"] == "L0_v2c_a030"
    assert res_combo["custom_combination"] is True


def test_gromacs_bond_precision_and_openmm_pure_python_xml(tmp_path):
    amber_dir = tmp_path / "amber_src"
    amber_dir.mkdir()
    prmtop = amber_dir / "go.parm7"
    inpcrd = amber_dir / "go.rst7"
    shutil.copyfile(WORKBENCH_ROOT / "resources/reference/old/gox.parm7", prmtop)
    shutil.copyfile(WORKBENCH_ROOT / "resources/reference/old/gox.rst7", inpcrd)
    top_sha = hashlib.sha256(prmtop.read_bytes()).hexdigest()
    crd_sha = hashlib.sha256(inpcrd.read_bytes()).hexdigest()
    (amber_dir / "PROVENANCE.json").write_text(
        json.dumps(
            {
                "preset": "gaff2_v3",
                "charge_model": "L0_v2c_a030",
                "files": {"prmtop_sha256": top_sha, "inpcrd_sha256": crd_sha},
            }
        )
        + "\n",
        encoding="utf-8",
    )

    gmx = export_gromacs(prmtop, inpcrd, tmp_path / "gmx/go")
    assert gmx["status"] == "exported"
    assert gmx["itp_bond_precision_upgraded"] > 0
    assert Path(gmx["provenance"]).is_file()
    itp_text = Path(gmx["itp"]).read_text(encoding="utf-8")
    # High precision scientific notation in upgraded [ bonds ]
    assert re.search(r"\d\.\d{12,}e[+-]\d{2}", itp_text)
    gmx_prov = json.loads(Path(gmx["provenance"]).read_text(encoding="utf-8"))
    assert gmx_prov["preset"] == "gaff2_v3"
    assert gmx_prov["amber_provenance_verified"] is True
    assert gmx_prov["amber_provenance_sha256"]

    # Stale/mismatched sidecar PROVENANCE.json must NOT be trusted or inherited
    stale_dir = tmp_path / "amber_stale"
    stale_dir.mkdir()
    prm_stale = stale_dir / "go.parm7"
    crd_stale = stale_dir / "go.rst7"
    shutil.copyfile(prmtop, prm_stale)
    shutil.copyfile(inpcrd, crd_stale)
    (stale_dir / "PROVENANCE.json").write_text(
        json.dumps(
            {
                "preset": "gaff2_v3",
                "charge_model": "L0_v2c_a030",
                "files": {"prmtop_sha256": "0" * 64, "inpcrd_sha256": "0" * 64},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    gmx_stale = export_gromacs(prm_stale, crd_stale, tmp_path / "gmx_stale/go")
    gmx_stale_prov = json.loads(Path(gmx_stale["provenance"]).read_text(encoding="utf-8"))
    assert gmx_stale_prov["amber_provenance_verified"] is False
    assert gmx_stale_prov["preset"] is None
    omm_stale = export_openmm(
        prm_stale, crd_stale, tmp_path / "openmm/stale.xml", _force_pure_python=True
    )
    omm_stale_prov = json.loads(Path(omm_stale["provenance"]).read_text(encoding="utf-8"))
    assert omm_stale_prov["amber_provenance_verified"] is False
    assert omm_stale_prov["preset"] is None

    omm_out = tmp_path / "openmm/go.xml"
    omm = export_openmm(
        prmtop,
        inpcrd,
        omm_out,
        nonbonded_method="NoCutoff",
        _force_pure_python=True,
    )
    assert omm["status"] == "exported"
    assert omm["serializer"] == "pure_python_amber_xml"
    assert Path(omm["xml"]).is_file()
    assert Path(omm["manifest"]).is_file()
    assert Path(omm["provenance"]).is_file()

    root = ET.fromstring(Path(omm["xml"]).read_text(encoding="utf-8"))
    assert root.tag == "System"
    assert root.find("PeriodicBoxVectors") is not None
    assert len(root.findall("./Particles/Particle")) == omm["atom_count"]
    forces = {f.attrib["name"]: f for f in root.findall("./Forces/Force")}
    assert {"HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "NonbondedForce"} <= set(forces)
    assert forces["NonbondedForce"].attrib["method"] == "0"
    assert len(forces["NonbondedForce"].findall("./Exceptions/Exception")) == omm["exception_count"]

    # Overwrite rejection
    with pytest.raises(FileExistsError):
        export_openmm(prmtop, inpcrd, omm_out, _force_pure_python=True)

    # CutoffNonPeriodic -> method 1
    omm_cut = export_openmm(
        prmtop,
        inpcrd,
        tmp_path / "openmm/go_cut.xml",
        nonbonded_method="CutoffNonPeriodic",
        cutoff_nm=1.2,
        _force_pure_python=True,
    )
    root_cut = ET.fromstring(Path(omm_cut["xml"]).read_text(encoding="utf-8"))
    nb_cut = next(f for f in root_cut.findall("./Forces/Force") if f.attrib["name"] == "NonbondedForce")
    assert nb_cut.attrib["method"] == "1"
    assert nb_cut.attrib["cutoff"] == "1.2"
    assert nb_cut.attrib["rfDielectric"] == "78.3"

    # PME rejected without periodic box
    with pytest.raises(ValueError, match="PME requires periodic box"):
        export_openmm(
            prmtop,
            inpcrd,
            tmp_path / "openmm/go_pme.xml",
            nonbonded_method="PME",
            _force_pure_python=True,
        )

    # Boxed PME uses method="4" and matches OpenMM AmberPrmtopFile usesPeriodic="0"/exceptionsUsePeriodic="0"
    boxed = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    boxed.box = [40.0, 40.0, 40.0, 90.0, 90.0, 90.0]
    pme_prm = tmp_path / "amber_src/boxed.prmtop"
    pme_crd = tmp_path / "amber_src/boxed.inpcrd"
    boxed.write_parm(str(pme_prm))
    boxed.write_rst7(str(pme_crd))
    omm_pme = export_openmm(
        pme_prm,
        pme_crd,
        tmp_path / "openmm/go_boxed_pme.xml",
        nonbonded_method="PME",
        cutoff_nm=1.0,
        _force_pure_python=True,
    )
    root_pme = ET.fromstring(Path(omm_pme["xml"]).read_text(encoding="utf-8"))
    forces_pme = {f.attrib["name"]: f for f in root_pme.findall("./Forces/Force")}
    assert forces_pme["NonbondedForce"].attrib["method"] == "4"
    assert forces_pme["HarmonicBondForce"].attrib["usesPeriodic"] == "0"
    assert forces_pme["HarmonicAngleForce"].attrib["usesPeriodic"] == "0"
    assert forces_pme["PeriodicTorsionForce"].attrib["usesPeriodic"] == "0"
    assert forces_pme["NonbondedForce"].attrib["exceptionsUsePeriodic"] == "0"

    # NBFIX topologies must be rejected cleanly rather than dropping pair-specific LJ overrides
    nb_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    nb_parm.parm_data["LENNARD_JONES_ACOEF"][1] *= 1.8
    nbfix_prm = tmp_path / "amber_src/nbfix.prmtop"
    nb_parm.write_parm(str(nbfix_prm))
    with pytest.raises(ValueError, match="NBFIX"):
        export_openmm(
            nbfix_prm,
            inpcrd,
            tmp_path / "openmm/nbfix.xml",
            _force_pure_python=True,
        )

    # Geometric combining rule / pair mismatch must be rejected even when has_NBFIX() is False
    geo_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    geo_parm.combining_rule = "geometric"
    geo_prm = tmp_path / "amber_src/geo.prmtop"
    geo_parm.write_parm(str(geo_prm))
    # Also test reloaded prmtop where combining_rule attribute defaults to lorentz but pairs use geometric mixing
    import math as _math

    geo2 = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    ntypes = int(geo2.parm_data["POINTERS"][1])
    nb_map = {int(a.nb_idx): (float(a.rmin), float(a.epsilon)) for a in geo2.atoms if int(a.nb_idx) > 0}
    nb_idx_arr = geo2.parm_data["NONBONDED_PARM_INDEX"]
    for i in sorted(nb_map):
        for j in sorted(nb_map):
            if i >= j:
                continue
            ptr = int(nb_idx_arr[ntypes * (i - 1) + (j - 1)])
            if ptr > 0:
                ri, ei = nb_map[i]
                rj, ej = nb_map[j]
                r_geo = 2.0 * _math.sqrt(ri * rj)
                e_geo = _math.sqrt(ei * ej)
                geo2.parm_data["LENNARD_JONES_ACOEF"][ptr - 1] = e_geo * (r_geo**12)
                geo2.parm_data["LENNARD_JONES_BCOEF"][ptr - 1] = 2.0 * e_geo * (r_geo**6)
    geo2.combining_rule = "geometric"
    assert geo2.has_NBFIX() is False
    geo2_prm = tmp_path / "amber_src/geo2.prmtop"
    geo2.write_parm(str(geo2_prm))
    with pytest.raises(ValueError, match="Lorentz-Berthelot"):
        export_openmm(
            geo2_prm,
            inpcrd,
            tmp_path / "openmm/geo2.xml",
            _force_pure_python=True,
        )

    # Nonzero LENNARD_JONES_CCOEF (12-6-4 terms) must be rejected cleanly even when has_NBFIX() is False
    c4_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    c4_parm.add_flag("LENNARD_JONES_CCOEF", "5E16.8", num_items=len(c4_parm.parm_data["LENNARD_JONES_ACOEF"]))
    c4_parm.parm_data["LENNARD_JONES_CCOEF"][0] = 12.5
    c4_prm = tmp_path / "amber_src/c4.prmtop"
    c4_parm.write_parm(str(c4_prm))
    c4_chk = pmd.load_file(str(c4_prm), xyz=str(inpcrd))
    assert c4_chk.has_NBFIX() is False
    with pytest.raises(ValueError, match="LENNARD_JONES_CCOEF"):
        export_openmm(
            c4_prm,
            inpcrd,
            tmp_path / "openmm/c4.xml",
            _force_pure_python=True,
        )

    # Active LJ pair pointing to zero-valued 10-12 HBOND entry (has_NBFIX() == False) must be rejected
    hb_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    nt_hb = int(hb_parm.parm_data["POINTERS"][1])
    nb_i0 = int(hb_parm.atoms[0].nb_idx)
    nb_j0 = int(hb_parm.atoms[1].nb_idx)
    hb_parm.parm_data["NONBONDED_PARM_INDEX"][nt_hb * (nb_i0 - 1) + (nb_j0 - 1)] = -1
    hb_parm.parm_data["NONBONDED_PARM_INDEX"][nt_hb * (nb_j0 - 1) + (nb_i0 - 1)] = -1
    hb_parm.parm_data["HBOND_ACOEF"] = [0.0]
    hb_parm.parm_data["HBOND_BCOEF"] = [0.0]
    hb_parm.parm_data["POINTERS"][19] = 1
    hb_prm = tmp_path / "amber_src/hb_zero.prmtop"
    hb_parm.write_parm(str(hb_prm))
    hb_chk = pmd.load_file(str(hb_prm), xyz=str(inpcrd))
    assert hb_chk.has_NBFIX() is False
    with pytest.raises(ValueError, match="10-12 HBOND"):
        export_openmm(
            hb_prm,
            inpcrd,
            tmp_path / "openmm/hb_zero.xml",
            _force_pure_python=True,
        )

    # Truncated/corrupted sidecar PROVENANCE.json must be treated safely as unverified without partial crash
    trunc_dir = tmp_path / "amber_trunc"
    trunc_dir.mkdir()
    prm_tr = trunc_dir / "go.parm7"
    crd_tr = trunc_dir / "go.rst7"
    shutil.copyfile(prmtop, prm_tr)
    shutil.copyfile(inpcrd, crd_tr)
    (trunc_dir / "PROVENANCE.json").write_text("{", encoding="utf-8")
    gmx_tr = export_gromacs(prm_tr, crd_tr, tmp_path / "gmx_trunc/go")
    gmx_tr_prov = json.loads(Path(gmx_tr["provenance"]).read_text(encoding="utf-8"))
    assert gmx_tr_prov["amber_provenance_verified"] is False
    assert gmx_tr_prov["preset"] is None
    omm_tr = export_openmm(
        prm_tr, crd_tr, tmp_path / "openmm/trunc.xml", _force_pure_python=True
    )
    omm_tr_prov = json.loads(Path(omm_tr["provenance"]).read_text(encoding="utf-8"))
    assert omm_tr_prov["amber_provenance_verified"] is False
    assert omm_tr_prov["preset"] is None


def test_cli_list_presets_and_openmm_script(tmp_path):
    res_human = subprocess.run(
        [sys.executable, "03_prepare_amber.py", "--list-presets"],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert res_human.returncode == 0, res_human.stderr
    assert "gaff2_v3" in res_human.stdout
    assert "gaff2_current" in res_human.stdout

    res_json = subprocess.run(
        [sys.executable, "03_prepare_amber.py", "--list-presets", "--json"],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert res_json.returncode == 0, res_json.stderr
    payload = json.loads(res_json.stdout)
    assert payload["status"] == "AVAILABLE_FORCE_FIELD_PRESETS"
    assert [p["name"] for p in payload["presets"]] == ["gaff2_v3", "gaff2_current"]

    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    xml_out = tmp_path / "cli_go.xml"
    res_omm = subprocess.run(
        [
            sys.executable,
            "08_export_openmm.py",
            "--topology",
            str(prmtop),
            "--coordinates",
            str(inpcrd),
            "--output",
            str(xml_out),
            "--json",
        ],
        cwd=WORKBENCH_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert res_omm.returncode == 0, res_omm.stderr
    omm_data = json.loads(res_omm.stdout)
    assert omm_data["status"] == "exported"
    assert xml_out.is_file()
