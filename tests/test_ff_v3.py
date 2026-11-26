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
    decompose_go_v3,
    list_force_field_presets,
    load_custom_charges,
    resolve_force_field_preset,
    verify_v3_bundle,
)
from go_builder.topology import export_gromacs, export_openmm, import_amber, prepare_amber

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
    assert len(info["provenance"]["files"]) == 25
    assert Path(info["frcmod_path"]).name == "GO_gaff2_v3.frcmod"
    assert Path(info["l0_v2_path"]).name == "L0_v2_charges.json"
    assert len(V3_RESIDUE_NAMES) == 18

    # Standalone multi-engine files present and well-formed
    ref_dir = Path(info["reference_dir"])
    leaprc_txt = (ref_dir / "leaprc.go_gaff2_v3").read_text(encoding="utf-8")
    assert "addAtomTypes {" in leaprc_txt and "loadAmberParams GO_gaff2_v3.frcmod" in leaprc_txt
    for res in V3_RESIDUE_NAMES:
        assert f"loadoff {res}_gaff2_v3.lib" in leaprc_txt
    nb_itp = (ref_dir / "GO_gaff2_v3_ffnonbonded.itp").read_text(encoding="utf-8")
    assert "[ atomtypes ]" in nb_itp and "ga          6 " in nb_itp and "gw          1 " in nb_itp
    b_itp = (ref_dir / "GO_gaff2_v3_ffbonded.itp").read_text(encoding="utf-8")
    assert "[ bondtypes ]" in b_itp and "[ angletypes ]" in b_itp and "[ dihedraltypes ]" in b_itp
    omm_xml = (ref_dir / "GO_gaff2_v3_openmm.xml").read_text(encoding="utf-8")
    assert "<ForceField>" in omm_xml and '<PeriodicTorsionForce ordering="amber">' in omm_xml
    lmp_inc = (ref_dir / "GO_gaff2_v3_lammps.inc").read_text(encoding="utf-8")
    assert "improper_style  cvff" in lmp_inc and "dihedral_style  fourier" in lmp_inc

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


def test_openmm_and_gromacs_reject_cmap_and_gromacs_rejects_c4_nbfix_and_explicit_exclusions(tmp_path):
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"

    # 1. CMAP torsions must be rejected by both OpenMM XML and GROMACS exporters without writing partial files
    cm_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    ct = pmd.CmapType(4, [float(i) * 0.1 for i in range(16)])
    cm_parm.cmap_types.append(ct)
    cm_parm.cmaps.append(
        pmd.Cmap(cm_parm.atoms[0], cm_parm.atoms[1], cm_parm.atoms[2], cm_parm.atoms[3], cm_parm.atoms[4], type=ct)
    )
    cm_parm.add_flag("CMAP_COUNT", "2I8", data=[1, 1])
    cm_parm.add_flag("CMAP_RESOLUTION", "20I4", data=[4])
    cm_parm.add_flag("CMAP_PARAMETER_01", "8F9.5", data=[float(i) * 0.1 for i in range(16)])
    cm_parm.add_flag("CMAP_INDEX", "6I8", data=[1, 2, 3, 4, 5, 1])
    cm_parm.pointers["CMAP"] = 1
    cm_parm.pointers["CMAP_TYPES"] = 1
    cmap_prm = tmp_path / "cmap.prmtop"
    cm_parm.save(str(cmap_prm))
    assert len(pmd.load_file(str(cmap_prm), xyz=str(inpcrd)).cmaps) == 1

    with pytest.raises(ValueError, match="CMAP"):
        export_openmm(cmap_prm, inpcrd, tmp_path / "openmm_cmap/go.xml", _force_pure_python=True)
    assert not (tmp_path / "openmm_cmap").exists()

    with pytest.raises(ValueError, match="CMAP"):
        export_gromacs(cmap_prm, inpcrd, tmp_path / "gmx_cmap/go")
    assert not (tmp_path / "gmx_cmap").exists()

    # 2. Explicit Amber nonbonded exclusion: OpenMM XML emits a zeroed <Exception>, GROMACS rejects cleanly
    excl_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    excl_parm.atoms[0].exclude(excl_parm.atoms[939])
    excl_parm.remake_parm()
    excl_prm = tmp_path / "excl.prmtop"
    excl_parm.save(str(excl_prm))
    reloaded_excl = pmd.load_file(str(excl_prm), xyz=str(inpcrd))
    assert 939 in [p.idx for p in reloaded_excl.atoms[0].exclusion_partners]

    omm_excl = export_openmm(
        excl_prm,
        inpcrd,
        tmp_path / "openmm_excl/go.xml",
        _force_pure_python=True,
    )
    omm_excl_prov = json.loads(Path(omm_excl["provenance"]).read_text(encoding="utf-8"))
    assert omm_excl_prov["counts"]["explicit_exclusions"] == 1
    root_excl = ET.fromstring(Path(omm_excl["xml"]).read_text(encoding="utf-8"))
    exceptions = root_excl.findall("./Forces/Force[@name='NonbondedForce']/Exceptions/Exception")
    matched = [e for e in exceptions if e.attrib["p1"] == "0" and e.attrib["p2"] == "939"]
    assert len(matched) == 1
    assert matched[0].attrib["q"] == "0"
    assert matched[0].attrib["eps"] == "0"

    with pytest.raises(ValueError, match="explicit nonbonded exclusion"):
        export_gromacs(excl_prm, inpcrd, tmp_path / "gmx_excl/go")
    assert not (tmp_path / "gmx_excl").exists()

    # 3. GROMACS export must reject 12-6-4 LENNARD_JONES_CCOEF and NBFIX overrides cleanly
    c4_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    c4_parm.add_flag("LENNARD_JONES_CCOEF", "5E16.8", num_items=len(c4_parm.parm_data["LENNARD_JONES_ACOEF"]))
    c4_parm.parm_data["LENNARD_JONES_CCOEF"][0] = 12.5
    c4_prm = tmp_path / "gmx_c4.prmtop"
    c4_parm.write_parm(str(c4_prm))
    with pytest.raises(ValueError, match="LENNARD_JONES_CCOEF"):
        export_gromacs(c4_prm, inpcrd, tmp_path / "gmx_c4/go")
    assert not (tmp_path / "gmx_c4").exists()

    nb_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    nb_parm.parm_data["LENNARD_JONES_ACOEF"][1] *= 1.8
    nbfix_prm = tmp_path / "gmx_nbfix.prmtop"
    nb_parm.write_parm(str(nbfix_prm))
    with pytest.raises(ValueError, match="NBFIX"):
        export_gromacs(nbfix_prm, inpcrd, tmp_path / "gmx_nbfix/go")
    assert not (tmp_path / "gmx_nbfix").exists()


def test_v3_hydrogen_valence_order_independence(tmp_path):
    base = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    roles = classify_go_v3_roles(base)
    edge_h_ids = [a["id"] for a, r in zip(base["atoms"], roles, strict=True) if r == "edge_H"]
    assert len(edge_h_ids) >= 2

    hh_orig = deepcopy(base)
    hh_orig["bonds"].append({"a": edge_h_ids[0], "b": edge_h_ids[1], "order": 1, "image": [0, 0, 0]})

    hh_h_first = deepcopy(hh_orig)
    hh_h_first["atoms"] = [a for a in hh_h_first["atoms"] if a["element"] == "H"] + [
        a for a in hh_h_first["atoms"] if a["element"] != "H"
    ]

    for idx, broken in enumerate((hh_orig, hh_h_first)):
        with pytest.raises(ValueError, match="Unsupported hydrogen valence"):
            classify_go_v3_roles(broken)
        with pytest.raises(ValueError, match="Unsupported hydrogen valence"):
            prepare_amber(broken, tmp_path / f"hh_fail_{idx}", run=False)


def test_dedicated_templates_validation_and_junction_frcmod_provenance(tmp_path):
    anionic = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    frc1 = tmp_path / "junc1.frcmod"
    frc1.write_text("remark junc1\nMASS\n", encoding="utf-8")
    frc2 = tmp_path / "junc2.frcmod"
    frc2.write_text("remark junc2 different\nMASS\n", encoding="utf-8")

    def make_tpl(frc_path: Path) -> dict:
        return {
            "atoms": [
                {"id": a["id"], "element": a["element"], "atom_type": "ca", "charge": (-1.0 if i == 0 else 0.0)}
                for i, a in enumerate(anionic["atoms"])
            ],
            "bonds": deepcopy(anionic["bonds"]),
            "frcmod": str(frc_path),
            "provenance": {"source": "same_label"},
        }

    # Duplicate atom record must fail rather than last-write-wins
    tpl_dup = make_tpl(frc1)
    tpl_dup["atoms"].append(dict(tpl_dup["atoms"][0], charge=0.5))
    with pytest.raises(ValueError, match="Duplicate atom ID"):
        prepare_amber(anionic, tmp_path / "tpl_dup", parameter_templates=tpl_dup, run=False)

    # Boolean charge must be rejected
    tpl_bool = make_tpl(frc1)
    tpl_bool["atoms"][0]["charge"] = True
    with pytest.raises(ValueError, match="finite numeric"):
        prepare_amber(anionic, tmp_path / "tpl_bool", parameter_templates=tpl_bool, run=False)

    # Duplicate bond must be rejected
    tpl_dup_bond = make_tpl(frc1)
    tpl_dup_bond["bonds"].append(dict(tpl_dup_bond["bonds"][0]))
    with pytest.raises(ValueError, match="Duplicate bond"):
        prepare_amber(anionic, tmp_path / "tpl_dup_bond", parameter_templates=tpl_dup_bond, run=False)

    # Byte-distinct junction frcmod files with identical provenance labels must record distinct SHA-256 digests
    res1 = prepare_amber(anionic, tmp_path / "tpl_1", parameter_templates=make_tpl(frc1), run=False)
    res2 = prepare_amber(anionic, tmp_path / "tpl_2", parameter_templates=make_tpl(frc2), run=False)
    prov1 = json.loads(Path(res1["provenance"]).read_text(encoding="utf-8"))
    prov2 = json.loads(Path(res2["provenance"]).read_text(encoding="utf-8"))
    assert prov1["files"]["junction_frcmod"] == "junction.frcmod"
    assert prov1["files"]["junction_frcmod_sha256"] == hashlib.sha256(frc1.read_bytes()).hexdigest()
    assert prov2["files"]["junction_frcmod_sha256"] == hashlib.sha256(frc2.read_bytes()).hexdigest()
    assert prov1["files"]["junction_frcmod_sha256"] != prov2["files"]["junction_frcmod_sha256"]


def test_v3_bundle_rejects_unmanifested_files_and_input_hardening(tmp_path):
    anionic = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))

    # Unmanifested file inside gaff2_v3 reference directory must be rejected
    bundle_copy = tmp_path / "gaff2_v3_extra"
    shutil.copytree(DEFAULT_V3_REFERENCE, bundle_copy)
    (bundle_copy / "EXTRA_gaff2_v3.lib").write_text("bogus\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unmanifested"):
        verify_v3_bundle(bundle_copy)
    with pytest.raises(ValueError, match="Unmanifested"):
        prepare_amber(anionic, tmp_path / "amber_extra", preset="gaff2_v3", reference_dir=bundle_copy, run=False)

    # Custom charges must reject specifying both charges_e and charges_by_id
    q_vec, _ = compute_l0_v2_charges(anionic)
    dual_json = tmp_path / "dual_custom.json"
    dual_json.write_text(
        json.dumps({"charges_e": [float(x) for x in q_vec], "charges_by_id": {"bogus": None}}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="cannot specify both"):
        load_custom_charges(dual_json, anionic)

    # Custom charges must reject mixing structured charges_e with direct atom-ID keys
    mixed_json = tmp_path / "mixed_custom.json"
    mixed_json.write_text(
        json.dumps({"charges_e": [float(x) for x in q_vec], str(anionic["atoms"][0]["id"]): 0.0}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="cannot mix"):
        load_custom_charges(mixed_json, anionic)

    # Python custom charge maps must reject float/fractional atom IDs without silently truncating
    frac_id_map = {float(a["id"]) + 0.75: float(x) for a, x in zip(anionic["atoms"], q_vec, strict=True)}
    with pytest.raises(ValueError, match="exact integer"):
        load_custom_charges(frac_id_map, anionic)
    int_float_map = {float(a["id"]): float(x) for a, x in zip(anionic["atoms"], q_vec, strict=True)}
    with pytest.raises(ValueError, match="exact integer"):
        load_custom_charges(int_float_map, anionic)

    # Fractional or boolean model/per-atom formal_charge must be rejected
    bad_fc = deepcopy(anionic)
    bad_fc["formal_charge"] = -1.9
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_fc)
    with pytest.raises(ValueError, match="exact integer formal charge"):
        load_custom_charges([float(x) for x in q_vec], bad_fc)
    bad_fc["formal_charge"] = True
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_fc)
    bad_fc_atom = deepcopy(anionic)
    bad_fc_atom.pop("formal_charge", None)
    bad_fc_atom["atoms"][0]["formal_charge"] = -1.5
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_fc_atom)


def test_exporters_reject_distinct_14_lj_chamber_tables(tmp_path):
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    s = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    cp = pmd.amber.ChamberParm.from_structure(s)
    cp.parm_data["LENNARD_JONES_14_ACOEF"] = [x * 2.0 for x in cp.parm_data["LENNARD_JONES_14_ACOEF"]]
    cp.parm_data["LENNARD_JONES_14_BCOEF"] = [x * 2.0 for x in cp.parm_data["LENNARD_JONES_14_BCOEF"]]
    ch_prm = tmp_path / "chamber_14.parm7"
    cp.write_parm(str(ch_prm))

    with pytest.raises(ValueError, match="1-4 Lennard-Jones"):
        export_openmm(ch_prm, inpcrd, tmp_path / "omm_ch/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_ch").exists()

    with pytest.raises(ValueError, match="1-4 Lennard-Jones"):
        export_gromacs(ch_prm, inpcrd, tmp_path / "gmx_ch/go")
    assert not (tmp_path / "gmx_ch").exists()


def test_exporters_hash_parsed_bytes_and_prepare_amber_verifies_bundle_once(tmp_path, monkeypatch):
    import go_builder.ff_v3 as ff_v3
    import go_builder.topology as topo_mod

    prmtop_src = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd_src = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"

    # Prepare two distinct valid topologies A and B with matching sidecars
    dir_work = tmp_path / "swap_amber"
    dir_work.mkdir()
    prm = dir_work / "go.parm7"
    crd = dir_work / "go.rst7"
    shutil.copyfile(prmtop_src, prm)
    shutil.copyfile(inpcrd_src, crd)
    sha_a = hashlib.sha256(prm.read_bytes()).hexdigest()
    sha_crd = hashlib.sha256(crd.read_bytes()).hexdigest()
    parm_a = pmd.load_file(str(prm), xyz=str(crd))
    q_a0 = float(parm_a.atoms[0].charge)

    parm_b = pmd.load_file(str(prm), xyz=str(crd))
    parm_b.atoms[0].charge = q_a0 + 0.125
    parm_b.atoms[1].charge = float(parm_b.atoms[1].charge) - 0.125
    prm_b = tmp_path / "go_b.parm7"
    parm_b.write_parm(str(prm_b))
    sha_b = hashlib.sha256(prm_b.read_bytes()).hexdigest()
    assert sha_b != sha_a
    prov_b = json.dumps(
        {
            "preset": "gaff2_v3",
            "charge_model": "L0_v2c_a030",
            "files": {"prmtop_sha256": sha_b, "inpcrd_sha256": sha_crd},
        }
    )

    orig_link = topo_mod._link_upstream_provenance

    def swapping_link(topology, coordinates, **kwargs):
        shutil.copyfile(prm_b, prm)
        (dir_work / "PROVENANCE.json").write_text(prov_b + "\n", encoding="utf-8")
        return orig_link(topology, coordinates, **kwargs)

    monkeypatch.setattr(topo_mod, "_link_upstream_provenance", swapping_link)

    omm = export_openmm(prm, crd, tmp_path / "omm_swap/go.xml", _force_pure_python=True)
    assert omm["source_topology_sha256"] == sha_a
    omm_prov = json.loads(Path(omm["provenance"]).read_text(encoding="utf-8"))
    assert omm_prov["amber_provenance_verified"] is False
    assert omm_prov["preset"] is None
    root = ET.fromstring(Path(omm["xml"]).read_text(encoding="utf-8"))
    q_emitted = float(root.findall("./Forces/Force[@name='NonbondedForce']/Particles/Particle")[0].attrib["q"])
    assert q_emitted == pytest.approx(q_a0, abs=1e-12)

    shutil.copyfile(prmtop_src, prm)
    (dir_work / "PROVENANCE.json").unlink()
    gmx = export_gromacs(prm, crd, tmp_path / "gmx_swap/go")
    assert gmx["source_topology_sha256"] == sha_a
    gmx_prov = json.loads(Path(gmx["provenance"]).read_text(encoding="utf-8"))
    assert gmx_prov["amber_provenance_verified"] is False
    assert gmx_prov["preset"] is None
    monkeypatch.undo()

    # Default prepare_amber must verify/snapshot the v3 bundle strictly once
    sheet = build_sheet({"nx": 5, "ny": 5, "seed": 41, "groups": {"OH": 2, "epoxide": 1, "COOH": 1}})
    calls = 0
    orig_verify = ff_v3.verify_v3_bundle

    def counted_verify(reference_dir=None):
        nonlocal calls
        calls += 1
        return orig_verify(reference_dir)

    monkeypatch.setattr(ff_v3, "verify_v3_bundle", counted_verify)
    monkeypatch.setattr(topo_mod, "verify_v3_bundle", counted_verify)
    res = prepare_amber(sheet, tmp_path / "single_snap", run=False)
    assert calls == 1
    assert Path(res["provenance"]).is_file()


def test_import_amber_annotates_all_v3_basal_epoxide_colors(tmp_path):
    sheet = build_sheet({"nx": 5, "ny": 5, "seed": 5, "groups": {"epoxide": 10}})
    assigned = assign_gaff2_v3(sheet)
    resnames = {a["resname"] for a in assigned["atoms"]}
    assert "EPX" in resnames and "EP2" in resnames

    struct = pmd.Structure()
    atype_map: dict[str, pmd.topologyobjects.AtomType] = {}
    by_res: dict[tuple[int, str], list[dict]] = {}
    for atom in assigned["atoms"]:
        by_res.setdefault((atom["resid"], atom["resname"]), []).append(atom)
    coords = []
    id_to_pmd: dict[int, pmd.Atom] = {}
    for (resid, resname), group in by_res.items():
        for a in group:
            at_name = a["atom_type"]
            if at_name not in atype_map:
                num = {"C": 6, "H": 1, "O": 8}[a["element"]]
                mass = {"C": 12.01, "H": 1.008, "O": 16.0}[a["element"]]
                at = pmd.topologyobjects.AtomType(at_name, len(atype_map) + 1, mass, num)
                at.set_lj_params(0.1, 1.7)
                atype_map[at_name] = at
            at = atype_map[at_name]
            pa = pmd.Atom(name=a["name"], type=at_name, atomic_number=at.atomic_number, mass=at.mass, charge=a["charge"])
            pa.atom_type = at
            pa.solvent_radius = 1.5
            pa.screen = 0.8
            struct.add_atom(pa, resname, resid)
            id_to_pmd[a["id"]] = pa
            coords.append(a["position"])
    btype = pmd.BondType(300.0, 1.45)
    struct.bond_types.append(btype)
    for b in assigned["bonds"]:
        struct.bonds.append(pmd.Bond(id_to_pmd[b["a"]], id_to_pmd[b["b"]], type=btype))
    struct.coordinates = np.asarray(coords, dtype=float)
    prm = tmp_path / "epx10.prmtop"
    crd = tmp_path / "epx10.rst7"
    struct.save(str(prm))
    struct.save(str(crd))

    imported = import_amber(prm, crd)
    epx_groups = [g for g in imported["metadata"]["groups"] if g["type"] == "epoxide"]
    assert len(epx_groups) == 10


def test_formal_charge_per_atom_always_checked_and_nonfinite_sums_rejected(tmp_path):
    anionic = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    anionic["formal_charge"] = -1

    bad_atom_fc = deepcopy(anionic)
    bad_atom_fc["atoms"][0]["formal_charge"] = -1.9
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_atom_fc)

    bad_atom_bool = deepcopy(anionic)
    bad_atom_bool["atoms"][0]["formal_charge"] = True
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_atom_bool)

    n = len(anionic["atoms"])
    overflow_vec = [1e308, 1e308, -1e308, -1e308] + [0.0] * (n - 4)
    with pytest.raises(ValueError, match="differs from formal charge"):
        load_custom_charges(overflow_vec, anionic)

    bundle_ov = tmp_path / "bundle_overflow"
    shutil.copytree(DEFAULT_V3_REFERENCE, bundle_ov)
    l0_path = bundle_ov / "L0_v2_charges.json"
    l0_data = json.loads(l0_path.read_text(encoding="utf-8"))
    l0_data["role_charges_e"]["scaffold_C"] = 1e308
    l0_raw = json.dumps(l0_data, indent=2).encode("utf-8")
    l0_path.write_bytes(l0_raw)
    prov_path = bundle_ov / "PROVENANCE.json"
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    prov["files"]["L0_v2_charges.json"]["sha256"] = hashlib.sha256(l0_raw).hexdigest()
    prov_path.write_text(json.dumps(prov, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="L0_v2 charge conservation check failed"):
        compute_l0_v2_charges(anionic, reference_dir=bundle_ov)


def test_verify_v3_bundle_rejects_hidden_unmanifested_files(tmp_path):
    bundle_hidden = tmp_path / "gaff2_v3_hidden"
    shutil.copytree(DEFAULT_V3_REFERENCE, bundle_hidden)
    (bundle_hidden / ".extra.lib").write_text("hidden\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unmanifested"):
        verify_v3_bundle(bundle_hidden)


def test_exporters_reject_nonfinite_params_fractional_periodicity_negative_scee_cap_and_missing_bonded_exclusions(tmp_path):
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"

    # 1. NaN mass and NaN bond force constant must be rejected cleanly by both exporters
    for flag, val, label in (
        ("MASS", float("nan"), "nan_mass"),
        ("BOND_FORCE_CONSTANT", float("nan"), "nan_bond_k"),
    ):
        af = pmd.amber.AmberFormat(str(prmtop))
        af.parm_data[flag][0] = val
        bad_prm = tmp_path / f"{label}.parm7"
        af.write_parm(str(bad_prm))
        with pytest.raises(ValueError):
            export_openmm(bad_prm, inpcrd, tmp_path / f"omm_{label}/go.xml", _force_pure_python=True)
        assert not (tmp_path / f"omm_{label}").exists()
        with pytest.raises(ValueError):
            export_gromacs(bad_prm, inpcrd, tmp_path / f"gmx_{label}/go")
        assert not (tmp_path / f"gmx_{label}").exists()

    # 2. Fractional dihedral periodicity and negative SCEE_SCALE_FACTOR on an active dihedral must be rejected by both exporters
    ref_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    dih_ids = list(ref_parm.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"])
    active_dih_type_idx = next(
        dih_ids[i + 4] for i in range(0, len(dih_ids), 5) if dih_ids[i + 2] > 0 and dih_ids[i + 3] > 0
    ) - 1
    for flag, val, match_re, label in (
        ("DIHEDRAL_PERIODICITY", 1.5, "periodicity", "frac_per"),
        ("SCEE_SCALE_FACTOR", -1.0, "scaling factor", "neg_scee"),
    ):
        af = pmd.amber.AmberFormat(str(prmtop))
        af.parm_data[flag][active_dih_type_idx] = val
        bad_prm = tmp_path / f"{label}.parm7"
        af.write_parm(str(bad_prm))
        with pytest.raises(ValueError, match=match_re):
            export_openmm(bad_prm, inpcrd, tmp_path / f"omm_{label}/go.xml", _force_pure_python=True)
        assert not (tmp_path / f"omm_{label}").exists()
        with pytest.raises(ValueError, match=match_re):
            export_gromacs(bad_prm, inpcrd, tmp_path / f"gmx_{label}/go")
        assert not (tmp_path / f"gmx_{label}").exists()

    # 3. Solvent cap (CAP_INFO / POINTERS[29] != 0) must be rejected by both exporters
    cap_af = pmd.amber.AmberFormat(str(prmtop))
    cap_af.parm_data["POINTERS"][29] = 1
    cap_af.add_flag("CAP_INFO", "10I8", data=[1])
    cap_af.add_flag("CAP_INFO2", "5E16.8", data=[10.0, 0.0, 0.0, 0.0])
    cap_prm = tmp_path / "cap.parm7"
    cap_af.write_parm(str(cap_prm))
    with pytest.raises(ValueError, match="CAP_INFO"):
        export_openmm(cap_prm, inpcrd, tmp_path / "omm_cap/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_cap").exists()
    with pytest.raises(ValueError, match="CAP_INFO"):
        export_gromacs(cap_prm, inpcrd, tmp_path / "gmx_cap/go")
    assert not (tmp_path / "gmx_cap").exists()

    # 4. Omitting a bonded 1-2 pair from EXCLUDED_ATOMS_LIST must be rejected by both exporters
    excl_af = pmd.amber.AmberFormat(str(prmtop))
    s_ref = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    bonded_1based = {b.idx + 1 for b in s_ref.atoms[0].bond_partners}
    n0 = excl_af.parm_data["NUMBER_EXCLUDED_ATOMS"][0]
    excl_list = excl_af.parm_data["EXCLUDED_ATOMS_LIST"]
    rm_pos = next(k for k in range(n0) if excl_list[k] in bonded_1based)
    del excl_list[rm_pos]
    excl_af.parm_data["NUMBER_EXCLUDED_ATOMS"][0] -= 1
    excl_af.parm_data["POINTERS"][10] -= 1
    missing_excl_prm = tmp_path / "missing_bonded_excl.parm7"
    excl_af.write_parm(str(missing_excl_prm))
    with pytest.raises(ValueError, match="EXCLUDED_ATOMS_LIST"):
        export_openmm(missing_excl_prm, inpcrd, tmp_path / "omm_miss_excl/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_miss_excl").exists()
    with pytest.raises(ValueError, match="EXCLUDED_ATOMS_LIST"):
        export_gromacs(missing_excl_prm, inpcrd, tmp_path / "gmx_miss_excl/go")
    assert not (tmp_path / "gmx_miss_excl").exists()


def test_charge_calculators_reject_mismatched_roles_argument():
    model = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    n = len(model["atoms"])
    custom_vec = [0.0] * n
    custom_vec[0] = -1.0
    for bad_roles in ([], ["scaffold_C"] * (n - 1), ["not_a_real_role"] * n):
        with pytest.raises(ValueError, match="roles"):
            compute_l0_v2_charges(model, roles=bad_roles)
        with pytest.raises(ValueError, match="roles"):
            compute_current_lib_role_charges(model, roles=bad_roles)
        with pytest.raises(ValueError, match="roles"):
            load_custom_charges(custom_vec, model, roles=bad_roles)


def test_per_atom_formal_charge_type_validated_even_when_formal_charges_unknown():
    model = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    model["formal_charge"] = -1
    model.setdefault("metadata", {})["formal_charges_unknown"] = True

    bad_bool = deepcopy(model)
    bad_bool["atoms"][0]["formal_charge"] = True
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_bool)

    bad_frac = deepcopy(model)
    bad_frac["atoms"][0]["formal_charge"] = -1.5
    with pytest.raises(ValueError, match="exact integer formal charge"):
        compute_l0_v2_charges(bad_frac)


def test_exporters_reject_invalid_box_and_zero_periodicity_nonzero_k_and_decompose_roles(tmp_path):
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"

    # 1. Invalid / nonfinite / nonpositive periodic box must be rejected cleanly before writing files
    for box_vals, label in (
        ([float("nan"), 40.0, 40.0, 90.0, 90.0, 90.0], "nan_box"),
        ([0.0, 40.0, 40.0, 90.0, 90.0, 90.0], "zero_box"),
        ([40.0, 40.0, 40.0, 180.0, 90.0, 90.0], "deg_angle"),
    ):
        bx = pmd.load_file(str(prmtop), xyz=str(inpcrd))
        bx.box = box_vals
        bx_prm = tmp_path / f"{label}.parm7"
        bx_crd = tmp_path / f"{label}.rst7"
        bx.write_parm(str(bx_prm))
        bx.write_rst7(str(bx_crd))
        with pytest.raises(ValueError, match="periodic box"):
            export_openmm(
                bx_prm,
                bx_crd,
                tmp_path / f"omm_{label}/go.xml",
                nonbonded_method="PME",
                _force_pure_python=True,
            )
        assert not (tmp_path / f"omm_{label}").exists()
        with pytest.raises(ValueError, match="periodic box"):
            export_gromacs(bx_prm, bx_crd, tmp_path / f"gmx_{label}/go")
        assert not (tmp_path / f"gmx_{label}").exists()

    # 2. Zero dihedral periodicity with nonzero force constant must be rejected by both exporters
    ref_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    dih_ids = list(ref_parm.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"])
    active_dih_idx = next(
        dih_ids[i + 4] for i in range(0, len(dih_ids), 5) if dih_ids[i + 2] > 0 and dih_ids[i + 3] > 0
    ) - 1
    zp = pmd.amber.AmberFormat(str(prmtop))
    zp.parm_data["DIHEDRAL_PERIODICITY"][active_dih_idx] = 0.0
    zp.parm_data["DIHEDRAL_FORCE_CONSTANT"][active_dih_idx] = 2.5
    zp_prm = tmp_path / "zero_per_nonzero_k.parm7"
    zp.write_parm(str(zp_prm))
    with pytest.raises(ValueError, match="Zero dihedral periodicity"):
        export_openmm(zp_prm, inpcrd, tmp_path / "omm_zp/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_zp").exists()
    with pytest.raises(ValueError, match="Zero dihedral periodicity"):
        export_gromacs(zp_prm, inpcrd, tmp_path / "gmx_zp/go")
    assert not (tmp_path / "gmx_zp").exists()

    # 3. Direct decompose_go_v3 and color_epoxides_v3 must validate explicit roles lengths/names cleanly
    model = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    for fn in (decompose_go_v3, color_epoxides_v3):
        with pytest.raises(ValueError, match="roles"):
            fn(model, [])
        with pytest.raises(ValueError, match="roles"):
            fn(model, ["unknown_role"] * len(model["atoms"]))


def test_astra_v8_exporter_and_box_hardening_guards(tmp_path):
    from go_builder.topology import (
        _amber_parm_to_openmm_xml,
        _verify_parameter_roundtrip,
        _verify_standard_lorentz_berthelot_lj,
    )

    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    ref_parm = pmd.load_file(str(prmtop), xyz=str(inpcrd))

    dih_ids = list(ref_parm.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"])
    dih_type_idx = next(
        dih_ids[k + 4] for k in range(0, len(dih_ids), 5) if dih_ids[k + 2] > 0 and dih_ids[k + 3] > 0
    ) - 1
    first_14_pair = next(
        tuple(sorted((d.atom1.idx, d.atom4.idx)))
        for d in ref_parm.dihedrals
        if not d.improper and not d.ignore_end
    )

    # BLOCKER #1: Suppressed 1-4 dihedral end + removed EXCLUDED_ATOMS_LIST pair -> rejected before files exist,
    # and _verify_parameter_roundtrip catches mismatched 1-4 adjusts directly
    b1_af = pmd.amber.AmberFormat(str(prmtop))
    a_lo, a_hi = first_14_pair
    for flag in ("DIHEDRALS_WITHOUT_HYDROGEN", "DIHEDRALS_INC_HYDROGEN"):
        arr = b1_af.parm_data[flag]
        for k in range(0, len(arr), 5):
            at1 = abs(arr[k]) // 3
            at4 = abs(arr[k + 3]) // 3
            if tuple(sorted((at1, at4))) == first_14_pair and arr[k + 2] > 0:
                arr[k + 2] = -arr[k + 2]
    num_excl = list(b1_af.parm_data["NUMBER_EXCLUDED_ATOMS"])
    excl_list = list(b1_af.parm_data["EXCLUDED_ATOMS_LIST"])
    offset = sum(num_excl[:a_lo])
    span = num_excl[a_lo]
    target_1based = a_hi + 1
    assert target_1based in excl_list[offset : offset + span]
    del excl_list[offset + excl_list[offset : offset + span].index(target_1based)]
    num_excl[a_lo] -= 1
    b1_af.parm_data["NUMBER_EXCLUDED_ATOMS"] = num_excl
    b1_af.parm_data["EXCLUDED_ATOMS_LIST"] = excl_list
    b1_af.parm_data["POINTERS"][10] -= 1
    b1_prm = tmp_path / "b1_unexcl14.parm7"
    b1_af.write_parm(str(b1_prm))
    with pytest.raises(ValueError, match="unexcluded|EXCLUDED_ATOMS_LIST"):
        export_gromacs(b1_prm, inpcrd, tmp_path / "gmx_b1/go")
    assert not (tmp_path / "gmx_b1").exists()

    gmx_ok = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_ref/go")
    assert gmx_ok["roundtrip_verified"]["pairs_14"] > 0
    rt_bad = pmd.load_file(gmx_ok["top"], xyz=gmx_ok["gro"])
    rt_bad.adjusts.pop()
    with pytest.raises(ValueError, match="1-4 pair identity mismatch"):
        _verify_parameter_roundtrip(ref_parm, rt_bad)

    # MAJOR #1: Nonstandard (1.6) and zero (0.0) SCEE / SCNB on active 1-4 dihedral -> rejected by both exporters
    for flag, val, label in (
        ("SCEE_SCALE_FACTOR", 1.6, "scee_16"),
        ("SCEE_SCALE_FACTOR", 0.0, "scee_0"),
        ("SCNB_SCALE_FACTOR", 0.0, "scnb_0"),
    ):
        af = pmd.amber.AmberFormat(str(prmtop))
        af.parm_data[flag][dih_type_idx] = val
        bad_prm = tmp_path / f"{label}.parm7"
        af.write_parm(str(bad_prm))
        with pytest.raises(ValueError, match="1-4 scaling factor"):
            export_openmm(bad_prm, inpcrd, tmp_path / f"omm_{label}/go.xml", _force_pure_python=True)
        assert not (tmp_path / f"omm_{label}").exists()
        with pytest.raises(ValueError, match="1-4 scaling factor"):
            export_gromacs(bad_prm, inpcrd, tmp_path / f"gmx_{label}/go")
        assert not (tmp_path / f"gmx_{label}").exists()

    # MAJOR #2 & MINOR #2: Triclinic (40,40,40, 90,90,45) emits OpenMM reduced box vectors; oversized PME cutoff rejected
    tri = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    tri.box = [40.0, 40.0, 40.0, 90.0, 90.0, 45.0]
    tri_prm = tmp_path / "tri.parm7"
    tri_crd = tmp_path / "tri.rst7"
    tri.write_parm(str(tri_prm))
    tri.write_rst7(str(tri_crd))
    tri_out = export_openmm(
        tri_prm,
        tri_crd,
        tmp_path / "omm_tri/go.xml",
        nonbonded_method="PME",
        cutoff_nm=1.0,
        _force_pure_python=True,
    )
    tri_root = ET.fromstring(Path(tri_out["xml"]).read_text(encoding="utf-8"))
    vecs = tri_root.find("PeriodicBoxVectors")
    assert vecs is not None
    ax = float(vecs.find("A").attrib["x"])
    bx = float(vecs.find("B").attrib["x"])
    by = float(vecs.find("B").attrib["y"])
    cx = float(vecs.find("C").attrib["x"])
    cy = float(vecs.find("C").attrib["y"])
    cz = float(vecs.find("C").attrib["z"])
    assert ax > 0 and by > 0 and cz > 0
    assert ax >= 2.0 * abs(bx) - 1e-9
    assert ax >= 2.0 * abs(cx) - 1e-9
    assert by >= 2.0 * abs(cy) - 1e-9

    with pytest.raises(ValueError, match="cutoff_nm"):
        export_openmm(
            tri_prm,
            tri_crd,
            tmp_path / "omm_bad_cut/go.xml",
            nonbonded_method="PME",
            cutoff_nm=1.5,
            _force_pure_python=True,
        )
    assert not (tmp_path / "omm_bad_cut").exists()

    # MAJOR #3: Perturbation flag IFPERT != 0 rejected by both exporters
    pert_af = pmd.amber.AmberFormat(str(prmtop))
    pert_af.parm_data["POINTERS"][20] = 1
    pert_prm = tmp_path / "pert.parm7"
    pert_af.write_parm(str(pert_prm))
    with pytest.raises(ValueError, match="IFPERT"):
        export_openmm(pert_prm, inpcrd, tmp_path / "omm_pert/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_pert").exists()
    with pytest.raises(ValueError, match="IFPERT"):
        export_gromacs(pert_prm, inpcrd, tmp_path / "gmx_pert/go")
    assert not (tmp_path / "gmx_pert").exists()

    # MAJOR #4: Negative torsion periodicity rejected instead of abs-normalized
    neg_af = pmd.amber.AmberFormat(str(prmtop))
    neg_af.parm_data["DIHEDRAL_PERIODICITY"][dih_type_idx] = -3.0
    neg_prm = tmp_path / "neg_per.parm7"
    neg_af.write_parm(str(neg_prm))
    with pytest.raises(ValueError, match="Negative.*periodicity"):
        export_openmm(neg_prm, inpcrd, tmp_path / "omm_negper/go.xml", _force_pure_python=True)
    with pytest.raises(ValueError, match="Negative.*periodicity"):
        export_gromacs(neg_prm, inpcrd, tmp_path / "gmx_negper/go")
    mem_neg = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    mem_neg.dihedrals[0].type.per = -3.0
    with pytest.raises(ValueError, match="periodicity"):
        _amber_parm_to_openmm_xml(mem_neg)

    # MAJOR #5: NaN in LENNARD_JONES_CCOEF / LENNARD_JONES_14_CCOEF rejected
    for flag, label in (("LENNARD_JONES_CCOEF", "nan_c4"), ("LENNARD_JONES_14_CCOEF", "nan_c4_14")):
        c4_af = pmd.amber.AmberFormat(str(prmtop))
        ntypes = int(c4_af.parm_data["POINTERS"][1])
        c4_af.add_flag(flag, "5E16.8", data=[float("nan")] * (ntypes * (ntypes + 1) // 2))
        c4_prm = tmp_path / f"{label}.parm7"
        c4_af.write_parm(str(c4_prm))
        with pytest.raises(ValueError, match=flag):
            export_openmm(c4_prm, inpcrd, tmp_path / f"omm_{label}/go.xml", _force_pure_python=True)
        with pytest.raises(ValueError, match=flag):
            export_gromacs(c4_prm, inpcrd, tmp_path / f"gmx_{label}/go")

    # MINOR #1: Explicit zero per-atom rmin_14/epsilon_14 override rejected on nonzero LJ atom
    zero14_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    nonzero_atom = next(a for a in zero14_s.atoms if a.epsilon > 0.0 and a.rmin > 0.0)
    nonzero_atom.rmin_14 = 0.0
    nonzero_atom.epsilon_14 = 0.0
    with pytest.raises(ValueError, match="rmin_14/epsilon_14"):
        _verify_standard_lorentz_berthelot_lj(zero14_s, target="OpenMM XML")


def test_astra_v8_charge_role_dedicated_and_packaging_guards(tmp_path):
    from go_builder.ff_v3 import _resolve_reference_dir

    # MAJOR #6 & MINOR #3: Swapped valid roles & complex NumPy custom charges rejected
    model = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    roles = classify_go_v3_roles(model)
    idx_h = roles.index("edge_H")
    idx_o = roles.index("hydroxyl_O")
    swapped_roles = list(roles)
    swapped_roles[idx_h], swapped_roles[idx_o] = swapped_roles[idx_o], swapped_roles[idx_h]
    q_ok, _ = compute_l0_v2_charges(model)
    with pytest.raises(ValueError, match="Explicit role mismatch"):
        compute_l0_v2_charges(model, roles=swapped_roles)
    with pytest.raises(ValueError, match="Explicit role mismatch"):
        compute_current_lib_role_charges(model, roles=swapped_roles)
    with pytest.raises(ValueError, match="Explicit role mismatch"):
        load_custom_charges(q_ok, model, roles=swapped_roles)
    with pytest.raises(ValueError, match="Explicit role mismatch"):
        color_epoxides_v3(model, swapped_roles)
    with pytest.raises(ValueError, match="Explicit role mismatch"):
        decompose_go_v3(model, swapped_roles)
    with pytest.raises(ValueError, match="finite real number"):
        load_custom_charges(q_ok.astype(complex) + 12j, model)

    # MAJOR #7: Dedicated parameter_templates cannot silently alter known molecular formal charge
    assigned = assign_gaff2_v3(model)
    frc_j = tmp_path / "junc.frcmod"
    frc_j.write_text("remark\n", encoding="utf-8")
    tpl_bad_q = {
        "atoms": [
            {
                "id": a["id"],
                "element": a["element"],
                "atom_type": a["atom_type"],
                "charge": float(a["charge"]) + (1.0 if idx == 0 else 0.0),
            }
            for idx, a in enumerate(assigned["atoms"])
        ],
        "bonds": deepcopy(assigned["bonds"]),
        "frcmod": str(frc_j),
        "provenance": {"source": "bad_total_q_test"},
    }
    with pytest.raises(ValueError, match="total charge"):
        prepare_amber(model, tmp_path / "bad_tpl_q", parameter_templates=tpl_bad_q, run=False)

    # MAJOR #8 & MINOR #4: Packaged wheel resources resolver + VALIDATION.md exact atom types check
    assert (_resolve_reference_dir("gaff2_v3") / "PROVENANCE.json").is_file()
    val_text = (WORKBENCH_ROOT / "VALIDATION.md").read_text(encoding="utf-8")
    assert "ga..gu" not in val_text
    assert "ga, gb, gc" in val_text


def test_astra_v9_exclusions_cleanup_pert_formal_charge_and_direct_xml_guards(tmp_path, monkeypatch):
    from go_builder.topology import _amber_parm_to_openmm_xml, _verify_parameter_roundtrip

    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"

    # BLOCKER 1: Active proper dihedral whose endpoints are 1-2 bonded must fail in BOTH exporters before mkdir
    ref_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    bonded_pair = (ref_s.bonds[0].atom1.idx, ref_s.bonds[0].atom2.idx)
    other_atoms = [a.idx for a in ref_s.atoms if a.idx not in bonded_pair and a.idx > 0]
    b1_af = pmd.amber.AmberFormat(str(prmtop))
    d_arr = list(b1_af.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"])
    # Overwrite first heavy-atom proper dihedral with 4 distinct atoms whose endpoints (1, 4) are 1-2 bonded and at3 > 0
    d_arr[0] = bonded_pair[0] * 3
    d_arr[1] = other_atoms[0] * 3
    d_arr[2] = other_atoms[1] * 3
    d_arr[3] = bonded_pair[1] * 3
    b1_af.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"] = d_arr
    b1_prm = tmp_path / "active_14_bonded.parm7"
    b1_af.write_parm(str(b1_prm))
    with pytest.raises(ValueError, match="active 1-4 dihedral on bonded/excluded pair"):
        export_openmm(b1_prm, inpcrd, tmp_path / "omm_b1/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_b1").exists()
    with pytest.raises(ValueError, match="active 1-4 dihedral on bonded/excluded pair"):
        export_gromacs(b1_prm, inpcrd, tmp_path / "gmx_b1/go")
    assert not (tmp_path / "gmx_b1").exists()

    # BLOCKER 2: 3-atom chain with two bonds, no angle term, and (0, 2) unexcluded in Amber
    chain = pmd.Structure()
    at_type = pmd.AtomType("C", None, 12.011, 6, charge=0.4)
    at_type.set_lj_params(0.086, 1.908)
    for idx in range(3):
        atom = pmd.Atom(name=f"C{idx + 1}", type="C", charge=0.4 if idx != 1 else -0.8, mass=12.011, atomic_number=6)
        atom.atom_type = at_type
        chain.add_atom(atom, "CH3", 1)
    bt = pmd.BondType(300.0, 1.40)
    chain.bond_types.append(bt)
    chain.bonds.append(pmd.Bond(chain.atoms[0], chain.atoms[1], type=bt))
    chain.bonds.append(pmd.Bond(chain.atoms[1], chain.atoms[2], type=bt))
    chain.positions = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.4, 0.9, 0.0]]
    amber_chain = pmd.amber.AmberParm.from_structure(chain)
    b2_raw = tmp_path / "chain3_raw.parm7"
    b2_prm = tmp_path / "chain3_unexcl13.parm7"
    b2_crd = tmp_path / "chain3_unexcl13.rst7"
    amber_chain.write_parm(str(b2_raw))
    amber_chain.write_rst7(str(b2_crd))
    # Strip auto-generated angle term and omit (0, 2) from Amber EXCLUDED_ATOMS_LIST via raw AmberFormat
    b2_af = pmd.amber.AmberFormat(str(b2_raw))
    b2_af.parm_data["ANGLES_WITHOUT_HYDROGEN"] = []
    b2_af.parm_data["POINTERS"][5] = 0
    b2_af.parm_data["POINTERS"][13] = 0
    b2_af.parm_data["NUMBER_EXCLUDED_ATOMS"] = [1, 1, 1]
    b2_af.parm_data["EXCLUDED_ATOMS_LIST"] = [2, 3, 0]
    b2_af.parm_data["POINTERS"][10] = 3
    b2_af.write_parm(str(b2_prm))
    b2_loaded = pmd.load_file(str(b2_prm), xyz=str(b2_crd))
    assert len(b2_loaded.angles) == 0
    with pytest.raises(ValueError, match="nrexcl=3"):
        export_gromacs(b2_prm, b2_crd, tmp_path / "gmx_b2/go")
    assert not (tmp_path / "gmx_b2").exists()
    # Direct roundtrip check also catches effective exclusion mismatch
    gmx_chain = pmd.gromacs.GromacsTopologyFile.from_structure(b2_loaded)
    gmx_chain.nrexcl = 3
    gmx_chain.molecules = {"system": (gmx_chain, 3)}
    with pytest.raises(ValueError, match="effective nonbonded exclusion mismatch"):
        _verify_parameter_roundtrip(b2_loaded, gmx_chain)

    # MAJOR 1: Mid-write OSError during GROMACS .top write AND sidecar write leaves zero files so retry succeeds
    orig_write = pmd.gromacs.GromacsTopologyFile.write

    def failing_write(self, top, molfile=None, **kw):
        Path(top).write_text("partial top\n", encoding="utf-8")
        raise OSError("simulated mid-write disk error")

    monkeypatch.setattr(pmd.gromacs.GromacsTopologyFile, "write", failing_write)
    with pytest.raises(OSError, match="simulated mid-write"):
        export_gromacs(prmtop, inpcrd, tmp_path / "gmx_retry/go")
    assert not (tmp_path / "gmx_retry").exists()
    monkeypatch.setattr(pmd.gromacs.GromacsTopologyFile, "write", orig_write)

    import go_builder.topology as topo_mod

    orig_json_fn = topo_mod._json
    calls = {"n": 0}

    def failing_manifest_json(path, payload):
        calls["n"] += 1
        if str(path).endswith(".manifest.json"):
            raise OSError("simulated sidecar write error")
        return orig_json_fn(path, payload)

    monkeypatch.setattr(topo_mod, "_json", failing_manifest_json)
    with pytest.raises(OSError, match="simulated sidecar"):
        export_gromacs(prmtop, inpcrd, tmp_path / "gmx_retry/go")
    assert not (tmp_path / "gmx_retry").exists()
    monkeypatch.setattr(topo_mod, "_json", orig_json_fn)
    ok_retry = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_retry/go")
    assert ok_retry["status"] == "exported"

    # MAJOR 2: Populated PERT_CHARGE bypass rejected even when IFPERT pointers are 0
    pc_af = pmd.amber.AmberFormat(str(prmtop))
    pc_af.add_flag("PERT_CHARGE", "5E16.8", data=[0.0] * len(pc_af.parm_data["CHARGE"]))
    pc_prm = tmp_path / "pert_charge.parm7"
    pc_af.write_parm(str(pc_prm))
    with pytest.raises(ValueError, match="perturbation|IFPERT"):
        export_openmm(pc_prm, inpcrd, tmp_path / "omm_pc/go.xml", _force_pure_python=True)
    with pytest.raises(ValueError, match="perturbation|IFPERT"):
        export_gromacs(pc_prm, inpcrd, tmp_path / "gmx_pc/go")

    # MAJOR 3: Dedicated template cannot contradict explicit molecular formal_charge even if unknown=True or all-zero atomic fcs
    model = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    assigned = assign_gaff2_v3(model)
    frc_j = tmp_path / "junc.frcmod"
    frc_j.write_text("remark\n", encoding="utf-8")
    tpl_zero = {
        "atoms": [
            {"id": a["id"], "element": a["element"], "atom_type": a["atom_type"], "charge": 0.0}
            for a in assigned["atoms"]
        ],
        "bonds": deepcopy(assigned["bonds"]),
        "frcmod": str(frc_j),
        "provenance": {"source": "zero_tpl"},
    }
    m_unk = deepcopy(model)
    m_unk["formal_charge"] = -1
    m_unk.setdefault("metadata", {})["formal_charges_unknown"] = True
    with pytest.raises(ValueError, match="total charge"):
        prepare_amber(m_unk, tmp_path / "bad_unk_tpl", parameter_templates=tpl_zero, run=False)

    m_zero_fc = deepcopy(model)
    m_zero_fc["formal_charge"] = -1
    for atom in m_zero_fc["atoms"]:
        atom["formal_charge"] = 0
    m_zero_fc["metadata"].pop("formal_charges_unknown", None)
    with pytest.raises(ValueError, match="Per-atom formal_charge sum"):
        prepare_amber(m_zero_fc, tmp_path / "bad_zero_atom_fc", parameter_templates=tpl_zero, run=False)

    # MINOR 1: Direct XML serializer rejects periodicity == 0 with nonzero phi_k
    mem_zero = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    mem_zero.dihedrals[0].type.per = 0.0
    mem_zero.dihedrals[0].type.phi_k = 1.0
    with pytest.raises(ValueError, match="periodicity 0 with nonzero force constant"):
        _amber_parm_to_openmm_xml(mem_zero)


def test_astra_v10_reversed_dihedral_nrexcl_concurrency_dup14_and_v3_fc_guards(tmp_path, monkeypatch):
    import re

    import go_builder.topology as topo_mod
    from go_builder.topology import _amber_parm_to_openmm_xml, _verify_parameter_roundtrip

    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    ref_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))

    # BLOCKER: Reversed proper dihedral placing atom 0 in slot 4 must fail in BOTH exporters pre-mkdir
    rev_af = pmd.amber.AmberFormat(str(prmtop))
    d_arr = list(rev_af.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"])
    k0 = next(
        k
        for k in range(0, len(d_arr), 5)
        if d_arr[k] == 0 and d_arr[k + 2] > 0 and d_arr[k + 3] > 0
    )
    d_arr[k0 : k0 + 4] = [d_arr[k0 + 3], abs(d_arr[k0 + 2]), d_arr[k0 + 1], d_arr[k0]]
    rev_af.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"] = d_arr
    rev_prm = tmp_path / "reversed_zero_slot4.parm7"
    rev_af.write_parm(str(rev_prm))
    with pytest.raises(ValueError, match="atom index 0"):
        export_openmm(rev_prm, inpcrd, tmp_path / "omm_rev/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_rev").exists()
    with pytest.raises(ValueError, match="atom index 0"):
        export_gromacs(rev_prm, inpcrd, tmp_path / "gmx_rev/go")
    assert not (tmp_path / "gmx_rev").exists()

    # MAJOR 1: Real parsed GROMACS .top/.itp edited from nrexcl 3 to 2 or 1 must be caught by _verify_parameter_roundtrip
    gmx_ok = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_nrexcl/go")
    itp_path = Path(gmx_ok["itp"])
    orig_itp = itp_path.read_text(encoding="utf-8")
    for bad_n in ("2", "1"):
        itp_path.write_text(
            re.sub(r"(system\s+)3\b", rf"\g<1>{bad_n}", orig_itp, count=1),
            encoding="utf-8",
        )
        rt_bad_nrexcl = pmd.load_file(gmx_ok["top"], xyz=gmx_ok["gro"])
        with pytest.raises(ValueError, match="nrexcl=|exclusion mismatch"):
            _verify_parameter_roundtrip(ref_s, rt_bad_nrexcl)
    itp_path.write_text(orig_itp, encoding="utf-8")

    # MAJOR 2: Overlapping exports to same prefix must NEVER unlink winning concurrent outputs
    orig_json_fn = topo_mod._json
    winner_dir = tmp_path / "concurrent_gmx"
    winner_dir.mkdir()
    winner_files = {
        k: winner_dir / f"go{ext}"
        for k, ext in (
            ("top", ".top"),
            ("itp", ".itp"),
            ("gro", ".gro"),
            ("manifest", ".manifest.json"),
            ("provenance", ".provenance.json"),
        )
    }

    def racer_json_gmx(path, payload):
        orig_json_fn(path, payload)
        if str(path).endswith(".manifest.json"):
            for fp in winner_files.values():
                fp.write_text(f"WINNER_{fp.name}", encoding="utf-8")

    monkeypatch.setattr(topo_mod, "_json", racer_json_gmx)
    with pytest.raises(FileExistsError):
        export_gromacs(prmtop, inpcrd, winner_dir / "go")
    monkeypatch.setattr(topo_mod, "_json", orig_json_fn)
    for fp in winner_files.values():
        assert fp.read_text(encoding="utf-8") == f"WINNER_{fp.name}"

    omm_dir = tmp_path / "concurrent_omm"
    omm_dir.mkdir()
    omm_winner = {
        "xml": omm_dir / "go.xml",
        "manifest": omm_dir / "go.manifest.json",
        "provenance": omm_dir / "go.provenance.json",
    }

    def racer_json_omm(path, payload):
        orig_json_fn(path, payload)
        if str(path).endswith(".manifest.json"):
            for fp in omm_winner.values():
                fp.write_text(f"WINNER_{fp.name}", encoding="utf-8")

    monkeypatch.setattr(topo_mod, "_json", racer_json_omm)
    with pytest.raises(FileExistsError):
        export_openmm(prmtop, inpcrd, omm_dir / "go.xml", _force_pure_python=True)
    monkeypatch.setattr(topo_mod, "_json", orig_json_fn)
    for fp in omm_winner.values():
        assert fp.read_text(encoding="utf-8") == f"WINNER_{fp.name}"

    # MAJOR 3: Duplicate active 1-4 dihedral records on same pair rejected before mkdir
    dup_af = pmd.amber.AmberFormat(str(prmtop))
    d_arr2 = list(dup_af.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"])
    active_by_pair: dict[tuple[int, int], int] = {}
    for k in range(0, len(d_arr2), 5):
        if d_arr2[k + 2] > 0 and d_arr2[k + 3] > 0:
            active_by_pair[tuple(sorted((abs(d_arr2[k]) // 3, abs(d_arr2[k + 3]) // 3)))] = k
    # Activate a second suppressed multi-term torsion sharing the exact same endpoint pair (or append + bump POINTERS[7])
    toggled = False
    for k in range(0, len(d_arr2), 5):
        p = tuple(sorted((abs(d_arr2[k]) // 3, abs(d_arr2[k + 3]) // 3)))
        if d_arr2[k + 2] < 0 and d_arr2[k + 3] > 0 and p in active_by_pair:
            d_arr2[k + 2] = -d_arr2[k + 2]
            toggled = True
            break
    if not toggled:
        k0_act = next(iter(active_by_pair.values()))
        d_arr2.extend(d_arr2[k0_act : k0_act + 5])
        dup_af.parm_data["POINTERS"][7] += 1
        dup_af.parm_data["POINTERS"][14] += 1
    dup_af.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"] = d_arr2
    dup_prm = tmp_path / "dup_active14.parm7"
    dup_af.write_parm(str(dup_prm))
    with pytest.raises(ValueError, match="duplicate active 1-4"):
        export_openmm(dup_prm, inpcrd, tmp_path / "omm_dup/go.xml", _force_pure_python=True)
    assert not (tmp_path / "omm_dup").exists()
    with pytest.raises(ValueError, match="duplicate active 1-4"):
        export_gromacs(dup_prm, inpcrd, tmp_path / "gmx_dup/go")
    assert not (tmp_path / "gmx_dup").exists()

    # MAJOR 4: assign_gaff2_v3 rejects explicit all-zero per-atom formal_charge contradicting formal_charge=-1
    model = build_sheet(json.loads((WORKBENCH_ROOT / "examples/finite_GO.json").read_text(encoding="utf-8")))
    m_v3_bad = deepcopy(model)
    m_v3_bad["formal_charge"] = -1
    for atom in m_v3_bad["atoms"]:
        atom["formal_charge"] = 0
    m_v3_bad.get("metadata", {}).pop("formal_charges_unknown", None)
    with pytest.raises(ValueError, match="Per-atom formal_charge sum"):
        assign_gaff2_v3(m_v3_bad)

    # MINOR: Direct XML serializer rejects nonfinite scee / scnb even on improper dihedrals
    mem_imp = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    imp = next(d for d in mem_imp.dihedrals if d.improper)
    imp.type.scee = float("nan")
    with pytest.raises(ValueError, match="Nonfinite torsion scee"):
        _amber_parm_to_openmm_xml(mem_imp)
    imp.type.scee = 1.2
    imp.type.scnb = float("inf")
    with pytest.raises(ValueError, match="Nonfinite torsion scnb"):
        _amber_parm_to_openmm_xml(mem_imp)


def test_astra_v11_roundtrip_exclusions_combrule_box_and_pkg_resource_priority(tmp_path, monkeypatch):
    from go_builder.charges import _resolve_reference_root
    from go_builder.ff import _resolve_current_reference
    from go_builder.ff_v3 import _resolve_reference_dir
    from go_builder.topology import _verify_parameter_roundtrip

    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    ref_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))

    # MAJOR 1: Appending [ exclusions ] to exported .itp must be caught by molecule-level roundtrip verification
    gmx_excl = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_mol_excl/go")
    itp_excl = Path(gmx_excl["itp"])
    itp_excl.write_text(itp_excl.read_text(encoding="utf-8") + "\n[ exclusions ]\n2 12\n", encoding="utf-8")
    rt_excl = pmd.load_file(gmx_excl["top"], xyz=gmx_excl["gro"])
    with pytest.raises(ValueError, match="explicit nonbonded exclusion mismatch"):
        _verify_parameter_roundtrip(ref_s, rt_excl)

    # MAJOR 2: Unbonded 2-atom system with different LJ radii catches [ defaults ] comb-rule 2 -> 3
    two = pmd.Structure()
    t1 = pmd.AtomType("C", None, 12.011, 6, charge=0.0)
    t1.set_lj_params(0.10, 1.90)
    t2 = pmd.AtomType("O", None, 15.999, 8, charge=0.0)
    t2.set_lj_params(0.20, 1.60)
    a1 = pmd.Atom(name="C1", type="C", charge=0.0, mass=12.011, atomic_number=6)
    a1.atom_type = t1
    a2 = pmd.Atom(name="O1", type="O", charge=0.0, mass=15.999, atomic_number=8)
    a2.atom_type = t2
    two.add_atom(a1, "UN1", 1)
    two.add_atom(a2, "UN2", 2)
    two.positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]]
    two_amb = pmd.amber.AmberParm.from_structure(two)
    two_prm = tmp_path / "two_unbonded.parm7"
    two_crd = tmp_path / "two_unbonded.rst7"
    two_amb.write_parm(str(two_prm))
    two_amb.write_rst7(str(two_crd))
    two_loaded = pmd.load_file(str(two_prm), xyz=str(two_crd))
    gmx_two = export_gromacs(two_prm, two_crd, tmp_path / "gmx_comb/go")
    top_p = Path(gmx_two["top"])
    lines = top_p.read_text(encoding="utf-8").splitlines()
    for idx, line in enumerate(lines):
        if line.strip().startswith("[ defaults ]"):
            for j in range(idx + 1, len(lines)):
                st = lines[j].strip()
                if st and not st.startswith(";"):
                    toks = st.split()
                    assert toks[1] == "2"
                    toks[1] = "3"
                    lines[j] = "  ".join(toks)
                    break
            break
    top_p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    rt_comb = pmd.load_file(gmx_two["top"], xyz=gmx_two["gro"])
    with pytest.raises(ValueError, match="combining_rule|comb_rule"):
        _verify_parameter_roundtrip(two_loaded, rt_comb)

    # MAJOR 3: Boxed GO validates periodic box vectors & rejects halved GRO box dimensions
    boxed = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    boxed.box = [80.0, 81.0, 82.0, 90.0, 90.0, 90.0]
    b_prm = tmp_path / "boxed.parm7"
    b_crd = tmp_path / "boxed.rst7"
    boxed.write_parm(str(b_prm))
    boxed.write_rst7(str(b_crd))
    gmx_box = export_gromacs(b_prm, b_crd, tmp_path / "gmx_box/go")
    assert gmx_box["roundtrip_verified"]["periodic_box_verified"] is True
    gro_p = Path(gmx_box["gro"])
    gro_lines = gro_p.read_text(encoding="utf-8").splitlines()
    gro_lines[-1] = "   4.00000   4.05000   4.10000"
    gro_p.write_text("\n".join(gro_lines) + "\n", encoding="utf-8")
    rt_half_box = pmd.load_file(gmx_box["top"], xyz=gmx_box["gro"])
    with pytest.raises(ValueError, match="periodic box mismatch"):
        _verify_parameter_roundtrip(pmd.load_file(str(b_prm), xyz=str(b_crd)), rt_half_box)

    # MINOR 1: Package-local _resources/reference takes precedence over ancestor resources/reference
    import go_builder.charges as chg_mod
    import go_builder.ff as ff_mod
    import go_builder.ff_v3 as v3_mod

    pkg_dir = tmp_path / "a/b/site/go_builder"
    pkg_dir.mkdir(parents=True)
    for sub in ("gaff2_v3", "current", "old", "opls"):
        (pkg_dir / "_resources/reference" / sub).mkdir(parents=True, exist_ok=True)
        (tmp_path / "a/resources/reference" / sub).mkdir(parents=True, exist_ok=True)
    fake_init = pkg_dir / "dummy.py"
    fake_init.write_text("", encoding="utf-8")
    monkeypatch.setattr(v3_mod, "__file__", str(fake_init))
    monkeypatch.setattr(ff_mod, "__file__", str(fake_init))
    monkeypatch.setattr(chg_mod, "__file__", str(fake_init))
    assert "_resources" in str(_resolve_reference_dir("gaff2_v3"))
    assert "_resources" in str(_resolve_current_reference())
    assert "_resources" in str(_resolve_reference_root())


def test_astra_v12_sdist_only_include_and_gromacs_directive_roundtrip_guards(tmp_path):
    import tomllib

    from go_builder.topology import _verify_parameter_roundtrip

    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    ref_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    gmx = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_dir_guard/go")
    itp_p = Path(gmx["itp"])
    orig_itp = itp_p.read_text(encoding="utf-8")

    rt_ok = pmd.load_file(gmx["top"], xyz=gmx["gro"])
    summary = _verify_parameter_roundtrip(ref_s, rt_ok, top_path=gmx["top"])
    assert summary["atoms"] == len(ref_s.atoms)

    for section_name, body in [
        ("constraints", "1 2 2 0.3"),
        ("position_restraints", "1 1 1000 1000 1000"),
        ("pairs_nb", "1 2 1 0.0 0.0 0.3 0.1"),
        ("settles", "; empty non-water [ settles ] directive must still be rejected"),
    ]:
        itp_p.write_text(f"{orig_itp}\n[ {section_name} ]\n{body}\n", encoding="utf-8")
        rt_bad = pmd.load_file(gmx["top"], xyz=gmx["gro"])
        with pytest.raises(
            ValueError,
            match=rf"unsupported/unverified directive \[{section_name}\]",
        ):
            _verify_parameter_roundtrip(ref_s, rt_bad, top_path=gmx["top"])

    itp_p.write_text(orig_itp, encoding="utf-8")
    rt_unknown = pmd.load_file(gmx["top"], xyz=gmx["gro"])
    rt_unknown.unknown_functional = True
    with pytest.raises(ValueError, match="unknown_functional"):
        _verify_parameter_roundtrip(ref_s, rt_unknown)

    rt_rb = pmd.load_file(gmx["top"], xyz=gmx["gro"])
    rt_rb.rb_torsions.append(object())
    with pytest.raises(ValueError, match="rb_torsions count mismatch"):
        _verify_parameter_roundtrip(ref_s, rt_rb)

    pyproject_data = tomllib.loads((WORKBENCH_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    sdist_cfg = pyproject_data["tool"]["hatch"]["build"]["targets"]["sdist"]
    assert "only-include" in sdist_cfg
    assert "include" not in sdist_cfg
    assert "/src" in sdist_cfg["only-include"]
    assert "/resources/reference" in sdist_cfg["only-include"]

    from types import SimpleNamespace

    from go_builder.topology import _verify_gromacs_topology_sections

    wat_mol = pmd.Structure()
    wat_mol.add_atom(pmd.Atom(name="O", type="OW", atomic_number=8), "WAT", 1)
    wat_mol.add_atom(pmd.Atom(name="H1", type="HW", atomic_number=1), "WAT", 1)
    wat_mol.add_atom(pmd.Atom(name="H2", type="HW", atomic_number=1), "WAT", 1)
    bt_wat = pmd.BondType(553.0, 0.9572)
    at_wat = pmd.AngleType(100.0, 104.52)
    wat_mol.bond_types.append(bt_wat)
    wat_mol.angle_types.append(at_wat)
    wat_mol.bonds.append(pmd.Bond(wat_mol.atoms[0], wat_mol.atoms[1], type=bt_wat))
    wat_mol.bonds.append(pmd.Bond(wat_mol.atoms[0], wat_mol.atoms[2], type=bt_wat))
    wat_mol.angles.append(pmd.Angle(wat_mol.atoms[1], wat_mol.atoms[0], wat_mol.atoms[2], type=at_wat))
    go_mol = pmd.Structure()
    go_mol.add_atom(pmd.Atom(name="C1", type="ca", atomic_number=6), "GOX", 1)
    mixed_conv = SimpleNamespace(molecules={"WAT": (wat_mol, 3), "GOX": (go_mol, 3)})

    mix_dir = tmp_path / "mix_settles"
    mix_dir.mkdir()
    (mix_dir / "wat.itp").write_text(
        "[ moleculetype ]\nWAT 3\n[ atoms ]\n; ok\n[ settles ]\n1 1 0.09572 0.15139\n",
        encoding="utf-8",
    )
    (mix_dir / "gox.itp").write_text(
        "[ moleculetype ]\nGOX 3\n[ atoms ]\n; no settles\n",
        encoding="utf-8",
    )
    (mix_dir / "ok.top").write_text(
        '[ defaults ]\n1 2 yes 0.5 0.833333\n#include "wat.itp"\n#include "gox.itp"\n'
        "[ system ]\nmix\n[ molecules ]\nWAT 1\nGOX 1\n",
        encoding="utf-8",
    )
    _verify_gromacs_topology_sections(mix_dir / "ok.top", converted=mixed_conv)

    (mix_dir / "gox_bad.itp").write_text(
        "[ moleculetype ]\nGOX 3\n[ settles ]\n1 1 0.09572 0.15139\n",
        encoding="utf-8",
    )
    (mix_dir / "bad.top").write_text(
        '[ defaults ]\n1 2 yes 0.5 0.833333\n#include "wat.itp"\n#include "gox_bad.itp"\n'
        "[ system ]\nmix\n[ molecules ]\nWAT 1\nGOX 1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match=r"unsupported/unverified directive \[settles\]"):
        _verify_gromacs_topology_sections(mix_dir / "bad.top", converted=mixed_conv)


def test_astra_v13_settle_include_bstate_and_wheel_cli_worker_guards(tmp_path, monkeypatch, capsys):
    import sys
    from types import SimpleNamespace

    import go_builder.cli as cli_mod
    from go_builder.topology import _verify_gromacs_topology_sections, _verify_parameter_roundtrip

    # BLOCKER 1: H1-O-H2 ordering and asymmetric O-H lengths must be rejected by export_gromacs
    def _build_water(order: tuple[str, str, str], r1: float, r2: float, tag: str) -> tuple[Path, Path]:
        w = pmd.Structure()
        to = pmd.AtomType("OW", None, 16.0, 8)
        to.set_lj_params(0.1521, 1.7682)
        th1 = pmd.AtomType("HW1", None, 1.0, 1)
        th1.set_lj_params(0.0, 0.0)
        th2 = pmd.AtomType("HW2", None, 1.0, 1)
        th2.set_lj_params(0.0, 0.0)
        specs = {
            "O": ("O", "OW", -0.834, 16.0, 8, to, [0.0, 0.0, 0.0]),
            "H1": ("H1", "HW1", 0.417, 1.0, 1, th1, [r1, 0.0, 0.0]),
            "H2": ("H2", "HW2", 0.417, 1.0, 1, th2, [-0.24, r2, 0.0]),
        }
        amap = {}
        for k in order:
            n, t, q, m, z, at, _pos = specs[k]
            atom = pmd.Atom(name=n, type=t, charge=q, mass=m, atomic_number=z)
            atom.atom_type = at
            w.add_atom(atom, "WAT", 1)
            amap[k] = atom
        b1 = pmd.BondType(553.0, r1)
        b2 = pmd.BondType(553.0, r2)
        w.bond_types.extend([b1, b2])
        w.bonds.append(pmd.Bond(amap["O"], amap["H1"], type=b1))
        w.bonds.append(pmd.Bond(amap["O"], amap["H2"], type=b2))
        at = pmd.AngleType(100.0, 104.52)
        w.angle_types.append(at)
        w.angles.append(pmd.Angle(amap["H1"], amap["O"], amap["H2"], type=at))
        w.positions = [specs[k][6] for k in order]
        w_amb = pmd.amber.AmberParm.from_structure(w)
        if abs(r1 - r2) > 1e-6:
            w_amb.bond_types[0].req = r1
            bt2 = pmd.BondType(553.0, r2, list=w_amb.bond_types)
            w_amb.bond_types.append(bt2)
            w_amb.bonds[1].type = bt2
            w_amb.remake_parm()
        wp = tmp_path / f"{tag}.parm7"
        wc = tmp_path / f"{tag}.rst7"
        w_amb.write_parm(str(wp))
        w_amb.write_rst7(str(wc))
        return wp, wc

    hoh_p, hoh_c = _build_water(("H1", "O", "H2"), 0.96, 0.96, "wat_hoh")
    with pytest.raises(ValueError, match="SETTLE water geometry/ordering mismatch"):
        export_gromacs(hoh_p, hoh_c, tmp_path / "out_hoh/wat")
    assert not (tmp_path / "out_hoh/wat.top").exists()

    asym_p, asym_c = _build_water(("O", "H1", "H2"), 0.96, 1.02, "wat_asym")
    with pytest.raises(ValueError, match="SETTLE water geometry/ordering mismatch"):
        export_gromacs(asym_p, asym_c, tmp_path / "out_asym/wat")
    assert not (tmp_path / "out_asym/wat.top").exists()

    # MAJOR 1 & 2: shared #include reuse across moltypes, malformed #include, compact #include
    wat_mol = pmd.Structure()
    wat_mol.add_atom(pmd.Atom(name="O", type="OW", atomic_number=8), "WAT", 1)
    wat_mol.add_atom(pmd.Atom(name="H1", type="HW", atomic_number=1), "WAT", 1)
    wat_mol.add_atom(pmd.Atom(name="H2", type="HW", atomic_number=1), "WAT", 1)
    bw = pmd.BondType(553.0, 0.9572)
    aw = pmd.AngleType(100.0, 104.52)
    wat_mol.bond_types.append(bw)
    wat_mol.angle_types.append(aw)
    wat_mol.bonds.append(pmd.Bond(wat_mol.atoms[0], wat_mol.atoms[1], type=bw))
    wat_mol.bonds.append(pmd.Bond(wat_mol.atoms[0], wat_mol.atoms[2], type=bw))
    wat_mol.angles.append(pmd.Angle(wat_mol.atoms[1], wat_mol.atoms[0], wat_mol.atoms[2], type=aw))
    go_mol = pmd.Structure()
    go_mol.add_atom(pmd.Atom(name="C1", type="ca", atomic_number=6), "GOX", 1)
    mixed_conv = SimpleNamespace(molecules={"WAT": (wat_mol, 3), "GOX": (go_mol, 3)})

    inc_dir = tmp_path / "inc_guards"
    inc_dir.mkdir()
    (inc_dir / "rigidity.inc").write_text("[ settles ]\n1 1 0.09572 0.15139\n", encoding="utf-8")
    (inc_dir / "shared.top").write_text(
        '[ defaults ]\n1 2 yes 0.5 0.833333\n'
        '[ moleculetype ]\nWAT 3\n#ifndef FLEXIBLE\n#include "rigidity.inc"\n#endif\n'
        '[ moleculetype ]\nGOX 3\n#ifndef FLEXIBLE\n#include "rigidity.inc"\n#endif\n'
        '[ system ]\nmix\n[ molecules ]\nWAT 1\nGOX 1\n',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match=r"unsupported/unverified directive \[settles\]"):
        _verify_gromacs_topology_sections(inc_dir / "shared.top", converted=mixed_conv)

    (inc_dir / "malformed.top").write_text('#include "oops\n', encoding="utf-8")
    with pytest.raises(ValueError, match="malformed #include directive"):
        _verify_gromacs_topology_sections(inc_dir / "malformed.top", converted=mixed_conv)

    (inc_dir / "restraint.inc").write_text("[ position_restraints ]\n1 1 1000 1000 1000\n", encoding="utf-8")
    (inc_dir / "compact.top").write_text('#include"restraint.inc"\n', encoding="utf-8")
    with pytest.raises(ValueError, match=r"unsupported/unverified directive \[position_restraints\]"):
        _verify_gromacs_topology_sections(inc_dir / "compact.top", converted=mixed_conv)

    # MAJOR 3: GROMACS B-state columns on staged [ atoms ] must be rejected
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    ref_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    gmx = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_bstate/go")
    itp_p = Path(gmx["itp"])
    itp_lines = itp_p.read_text(encoding="utf-8").splitlines()
    for i, ln in enumerate(itp_lines):
        if ln.strip().startswith("1 "):
            itp_lines[i] = f"{ln.split(';')[0].rstrip()}  C  9.0  99.0"
            break
    itp_p.write_text("\n".join(itp_lines) + "\n", encoding="utf-8")
    rt_b = pmd.load_file(gmx["top"], xyz=gmx["gro"])
    with pytest.raises(ValueError, match="unsupported B-state columns"):
        _verify_parameter_roundtrip(ref_s, rt_b, top_path=gmx["top"])

    # MAJOR 4: Wheel-mode CLI without checkout go.py falls back to `-m go_builder`
    fake_root = tmp_path / "no_go_py_root"
    fake_root.mkdir()
    monkeypatch.setattr(cli_mod, "ROOT", fake_root)
    capsys.readouterr()
    planned = cli_mod.main(
        "openmm",
        [
            "--topology",
            str(prmtop),
            "--coordinates",
            str(inpcrd),
            "--output",
            str(tmp_path / "omm_out/go.xml"),
            "--python",
            sys.executable,
            "--dry-run",
            "--json",
        ],
    )
    assert planned == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["command"][:4] == [sys.executable, "-m", "go_builder", "openmm"]
    assert all("go.py" not in tok for tok in payload["command"])


def test_neon_v1_settle_dup_continuation_define_ignore_end_and_decoy_gopy_guards(tmp_path):
    import os
    from types import SimpleNamespace

    from go_builder.topology import _verify_gromacs_topology_sections, _verify_parameter_roundtrip

    # Helper water moleculetype (TIP3P-like O-H1-H2)
    def _make_water_mol(resname: str) -> pmd.Structure:
        w = pmd.Structure()
        w.add_atom(pmd.Atom(name="O", type="OW", atomic_number=8), resname, 1)
        w.add_atom(pmd.Atom(name="H1", type="HW", atomic_number=1), resname, 1)
        w.add_atom(pmd.Atom(name="H2", type="HW", atomic_number=1), resname, 1)
        bw = pmd.BondType(553.0, 0.9572)
        aw = pmd.AngleType(100.0, 104.52)
        w.bond_types.append(bw)
        w.angle_types.append(aw)
        w.bonds.append(pmd.Bond(w.atoms[0], w.atoms[1], type=bw))
        w.bonds.append(pmd.Bond(w.atoms[0], w.atoms[2], type=bw))
        w.angles.append(pmd.Angle(w.atoms[1], w.atoms[0], w.atoms[2], type=aw))
        return w

    # 1. MAJOR 1: Duplicate [ settles ] blocks in one WAT moleculetype (inline & repeated #include)
    #    must be rejected, while two distinct water moleculetypes (WAT1 + WAT2) pass.
    two_wat_conv = SimpleNamespace(
        molecules={"WAT1": (_make_water_mol("WAT1"), 3), "WAT2": (_make_water_mol("WAT2"), 3)}
    )
    settle_dir = tmp_path / "settle_dup"
    settle_dir.mkdir()
    (settle_dir / "wat1.itp").write_text(
        "[ moleculetype ]\nWAT1 3\n[ atoms ]\n; ok\n[ settles ]\n1 1 0.09572 0.15139\n",
        encoding="utf-8",
    )
    (settle_dir / "wat2.itp").write_text(
        "[ moleculetype ]\nWAT2 3\n[ atoms ]\n; ok\n[ settles ]\n1 1 0.09572 0.15139\n",
        encoding="utf-8",
    )
    (settle_dir / "two_waters_ok.top").write_text(
        '[ defaults ]\n1 2 yes 0.5 0.833333\n'
        '#include "wat1.itp"\n#include "wat2.itp"\n'
        "[ system ]\nwaters\n[ molecules ]\nWAT1 1\nWAT2 1\n",
        encoding="utf-8",
    )
    _verify_gromacs_topology_sections(settle_dir / "two_waters_ok.top", converted=two_wat_conv)

    (settle_dir / "wat1_inline_dup.itp").write_text(
        "[ moleculetype ]\nWAT1 3\n[ atoms ]\n; ok\n"
        "[ settles ]\n1 1 0.09572 0.15139\n"
        "[ settles ]\n1 1 0.09572 0.15139\n",
        encoding="utf-8",
    )
    (settle_dir / "inline_dup.top").write_text(
        '[ defaults ]\n1 2 yes 0.5 0.833333\n#include "wat1_inline_dup.itp"\n'
        "[ system ]\nwaters\n[ molecules ]\nWAT1 1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match=r"duplicate \[settles\] directive in moleculetype 'WAT1'"):
        _verify_gromacs_topology_sections(settle_dir / "inline_dup.top", converted=two_wat_conv)

    (settle_dir / "settle_row.inc").write_text(
        "[ settles ]\n1 1 0.09572 0.15139\n", encoding="utf-8"
    )
    (settle_dir / "wat1_inc_dup.itp").write_text(
        "[ moleculetype ]\nWAT1 3\n[ atoms ]\n; ok\n"
        '#include "settle_row.inc"\n#include "settle_row.inc"\n',
        encoding="utf-8",
    )
    (settle_dir / "inc_dup.top").write_text(
        '[ defaults ]\n1 2 yes 0.5 0.833333\n#include "wat1_inc_dup.itp"\n'
        "[ system ]\nwaters\n[ molecules ]\nWAT1 1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match=r"duplicate \[settles\] directive in moleculetype 'WAT1'"):
        _verify_gromacs_topology_sections(settle_dir / "inc_dup.top", converted=two_wat_conv)

    # 2, 3, 5: MAJOR 2 (backslash continuations), MAJOR 3 (#define/#undef), MINOR 1 (malformed [ ... header)
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    ref_s = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    gmx = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_scan_base/go")
    top_p = Path(gmx["top"])
    itp_p = Path(gmx["itp"])
    orig_itp = itp_p.read_text(encoding="utf-8")
    rt_gox = pmd.load_file(str(top_p), xyz=gmx["gro"])

    # 2a. [ atoms ] B-state split via backslash continuation
    atoms_split_lines = orig_itp.splitlines()
    for i, ln in enumerate(atoms_split_lines):
        if ln.strip().startswith("1 "):
            atoms_split_lines[i] = f"{ln.split(';')[0].rstrip()} \\\n  C  9.0  99.0"
            break
    itp_p.write_text("\n".join(atoms_split_lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported line continuation"):
        _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    # 2b. [ position_restraints \\\n] header split via backslash continuation
    itp_p.write_text(
        f"{orig_itp}\n[ position_restraints \\\n]\n1 1 1000 1000 1000\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unsupported line continuation"):
        _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    # 2c. ; comment \ continuation
    itp_p.write_text(f"{orig_itp}\n; innocent comment \\\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported line continuation"):
        _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    # 3. #define (B-state macro and directive macro) and #undef rejected
    for pp_line in (
        "#define B_STATE C 9.0 99.0",
        "#define POSRES",
        "#undef FLEXIBLE",
    ):
        itp_p.write_text(f"{orig_itp}\n{pp_line}\n", encoding="utf-8")
        with pytest.raises(ValueError, match="unsupported preprocessor directive"):
            _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    # 5. Malformed [ position_restraints ] trailing section header rejected
    itp_p.write_text(
        f"{orig_itp}\n[ position_restraints ] trailing\n1 1 1000 1000 1000\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="malformed section header"):
        _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    # Restore clean .itp for tampering checks below
    itp_p.write_text(orig_itp, encoding="utf-8")

    # 4. MAJOR 4: ignore_end=True proper dihedral with SCEE=0.0, SCNB=0.0 on gox.parm7 succeeds
    #    in BOTH export_openmm and export_gromacs; valid 5-membered ring Amber topology also succeeds
    #    in BOTH; tampering with an active 1-4 pair's fudgeQQ / [ pairs ] is still rejected.
    gox_ign0 = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    ign_d = next(d for d in gox_ign0.dihedrals if (not d.improper) and d.ignore_end)
    zero_dt = pmd.DihedralType(
        ign_d.type.phi_k,
        ign_d.type.per,
        ign_d.type.phase,
        scee=0.0,
        scnb=0.0,
        list=gox_ign0.dihedral_types,
    )
    gox_ign0.dihedral_types.append(zero_dt)
    ign_d.type = zero_dt
    gox_ign0.remake_parm()
    gox_ign0_prm = tmp_path / "gox_ign0.parm7"
    gox_ign0.write_parm(str(gox_ign0_prm))
    gox_ign0_loaded = pmd.load_file(str(gox_ign0_prm), xyz=str(inpcrd))
    assert any(
        (not d.improper) and d.ignore_end and d.type.scee == 0.0 and d.type.scnb == 0.0
        for d in gox_ign0_loaded.dihedrals
    )
    omm_ign0 = export_openmm(
        gox_ign0_prm, inpcrd, tmp_path / "omm_ign0/go.xml", _force_pure_python=True
    )
    assert omm_ign0["status"] == "exported"
    gmx_ign0 = export_gromacs(gox_ign0_prm, inpcrd, tmp_path / "gmx_ign0/go")
    assert gmx_ign0["status"] == "exported"

    # Valid 5-membered ring Amber topology (all dihedrals ignore_end=True, scee=1.2, scnb=2.0, 0 active 1-4 pairs)
    ring = pmd.Structure()
    c_type = pmd.AtomType("C", None, 12.011, 6, charge=0.0)
    c_type.set_lj_params(0.1094, 1.908)
    for idx in range(5):
        a = pmd.Atom(name=f"C{idx + 1}", type="C", charge=0.0, mass=12.011, atomic_number=6)
        a.atom_type = c_type
        ring.add_atom(a, "RNG", 1)
    r_bt = pmd.BondType(300.0, 1.52)
    r_at = pmd.AngleType(60.0, 108.0)
    r_dt = pmd.DihedralType(0.5, 3, 0.0, scee=1.2, scnb=2.0)
    ring.bond_types.append(r_bt)
    ring.angle_types.append(r_at)
    ring.dihedral_types.append(r_dt)
    for idx in range(5):
        ring.bonds.append(pmd.Bond(ring.atoms[idx], ring.atoms[(idx + 1) % 5], type=r_bt))
        ring.angles.append(
            pmd.Angle(ring.atoms[idx], ring.atoms[(idx + 1) % 5], ring.atoms[(idx + 2) % 5], type=r_at)
        )
        ring.dihedrals.append(
            pmd.Dihedral(
                ring.atoms[idx],
                ring.atoms[(idx + 1) % 5],
                ring.atoms[(idx + 2) % 5],
                ring.atoms[(idx + 3) % 5],
                type=r_dt,
                ignore_end=True,
            )
        )
    import math as _math

    rad = 1.52 / (2.0 * _math.sin(_math.pi / 5.0))
    ring.positions = [
        [rad * _math.cos(2.0 * _math.pi * i / 5.0), rad * _math.sin(2.0 * _math.pi * i / 5.0), 0.0]
        for i in range(5)
    ]
    ring_amber = pmd.amber.AmberParm.from_structure(ring)
    ring_prm = tmp_path / "ring5.parm7"
    ring_crd = tmp_path / "ring5.rst7"
    ring_amber.write_parm(str(ring_prm))
    ring_amber.write_rst7(str(ring_crd))
    ring_reloaded = pmd.load_file(str(ring_prm), xyz=str(ring_crd))
    assert len(ring_reloaded.dihedrals) == 5
    assert all((not d.improper) and d.ignore_end for d in ring_reloaded.dihedrals)
    assert all(d.type.scee == pytest.approx(1.2) and d.type.scnb == pytest.approx(2.0) for d in ring_reloaded.dihedrals)
    omm_ring = export_openmm(ring_prm, ring_crd, tmp_path / "omm_ring/go.xml", _force_pure_python=True)
    assert omm_ring["status"] == "exported"
    gmx_ring = export_gromacs(ring_prm, ring_crd, tmp_path / "gmx_ring/go")
    assert gmx_ring["status"] == "exported"
    assert gmx_ring["roundtrip_verified"]["pairs_14"] == 0

    # Tampering with fudgeQQ in [ defaults ] on active 1-4 pairs must still be rejected
    orig_top = top_p.read_text(encoding="utf-8")
    tampered_top_lines = orig_top.splitlines()
    in_defaults = False
    for idx, ln in enumerate(tampered_top_lines):
        s = ln.strip()
        if s.startswith("["):
            in_defaults = s.lower().strip("[] ") == "defaults"
            continue
        if in_defaults and s and not s.startswith(";"):
            toks = s.split()
            if len(toks) >= 5:
                toks[4] = "0.5"
                tampered_top_lines[idx] = " ".join(toks)
                break
    top_p.write_text("\n".join(tampered_top_lines) + "\n", encoding="utf-8")
    rt_bad_fudge = pmd.load_file(str(top_p), xyz=gmx["gro"])
    with pytest.raises(ValueError, match="1-4 pair LJ/electrostatic parameter mismatch"):
        _verify_parameter_roundtrip(ref_s, rt_bad_fudge, top_path=top_p)
    top_p.write_text(orig_top, encoding="utf-8")

    # Tampering with [ pairs ] in .itp must also still be rejected
    tampered_pairs_lines = orig_itp.splitlines()
    in_pairs = False
    removed_pair = False
    filtered_lines = []
    for ln in tampered_pairs_lines:
        s = ln.strip()
        if s.startswith("["):
            in_pairs = s.lower().strip("[] ") == "pairs"
        elif in_pairs and not removed_pair and s and not s.startswith(";"):
            removed_pair = True
            continue
        filtered_lines.append(ln)
    assert removed_pair
    itp_p.write_text("\n".join(filtered_lines) + "\n", encoding="utf-8")
    rt_bad_pairs = pmd.load_file(str(top_p), xyz=gmx["gro"])
    with pytest.raises(ValueError, match="1-4 pair LJ/electrostatic parameter mismatch"):
        _verify_parameter_roundtrip(ref_s, rt_bad_pairs, top_path=top_p)

    # 6. MINOR 2: Packaged wheel layout (<tmp>/decoy_root/site/go_builder/_resources + decoy <tmp>/decoy_root/go.py)
    #    must invoke `-m go_builder`, NEVER decoy go.py.
    decoy_root = tmp_path / "decoy_root"
    site_dir = decoy_root / "site"
    pkg_dir = site_dir / "go_builder"
    shutil.copytree(
        WORKBENCH_ROOT / "src/go_builder",
        pkg_dir,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    shutil.copytree(WORKBENCH_ROOT / "resources/reference", pkg_dir / "_resources/reference")
    (decoy_root / "go.py").write_text(
        'raise SystemExit("DECOY_GO_PY_INVOKED")\n', encoding="utf-8"
    )
    dummy_model = tmp_path / "dummy_model.json"
    dummy_model.write_text("{}", encoding="utf-8")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(site_dir)
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "go_builder",
            "check",
            "--model",
            str(dummy_model),
            "--dry-run",
            "--json",
        ],
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    decoy_payload = json.loads(proc.stdout)
    assert decoy_payload["command"][:4] == [sys.executable, "-m", "go_builder", "check"]
    assert all("go.py" not in tok for tok in decoy_payload["command"])


def test_neon_v2_preprocessor_branch_isolation_and_complex_nonstandard_14_scaling(tmp_path):
    import math as _math
    from types import SimpleNamespace

    from parmed.tools.changeradii import ChRad

    from go_builder.complexes import assemble_complex
    from go_builder.topology import _verify_gromacs_topology_sections

    def _make_water_mol(resname: str, r1: float = 0.9572, r2: float = 0.9572) -> pmd.Structure:
        w = pmd.Structure()
        to = pmd.AtomType("OW", None, 16.0, 8)
        to.set_lj_params(0.1521, 1.7682)
        th = pmd.AtomType("HW", None, 1.008, 1)
        th.set_lj_params(0.0, 0.0)
        ao = pmd.Atom(name="O", type="OW", charge=-0.834, mass=16.0, atomic_number=8)
        ao.atom_type = to
        ah1 = pmd.Atom(name="H1", type="HW", charge=0.417, mass=1.008, atomic_number=1)
        ah1.atom_type = th
        ah2 = pmd.Atom(name="H2", type="HW", charge=0.417, mass=1.008, atomic_number=1)
        ah2.atom_type = th
        w.add_atom(ao, resname, 1)
        w.add_atom(ah1, resname, 1)
        w.add_atom(ah2, resname, 1)
        b1 = pmd.BondType(553.0, r1)
        w.bond_types.append(b1)
        w.bonds.append(pmd.Bond(ao, ah1, type=b1))
        if abs(r1 - r2) > 1e-6:
            b2 = pmd.BondType(553.0, r2)
            w.bond_types.append(b2)
        else:
            b2 = b1
        w.bonds.append(pmd.Bond(ao, ah2, type=b2))
        aw = pmd.AngleType(100.0, 104.52)
        w.angle_types.append(aw)
        w.angles.append(pmd.Angle(ah1, ao, ah2, type=aw))
        ang_rad = _math.radians(104.52)
        w.positions = [
            [0.0, 0.0, 0.0],
            [r1, 0.0, 0.0],
            [r2 * _math.cos(ang_rad), r2 * _math.sin(ang_rad), 0.0],
        ]
        return w

    # 1. B-state bypass via #ifdef NEVER_DEFINED / #ifdef FLEXIBLE / #ifndef FLEXIBLE [ system ]
    #    inside [ atoms ] and [ bonds ] -> all rejected.
    prmtop = WORKBENCH_ROOT / "resources/reference/old/gox.parm7"
    inpcrd = WORKBENCH_ROOT / "resources/reference/old/gox.rst7"
    gmx = export_gromacs(prmtop, inpcrd, tmp_path / "gmx_pp_base/go")
    top_p = Path(gmx["top"])
    itp_p = Path(gmx["itp"])
    orig_itp = itp_p.read_text(encoding="utf-8")
    rt_gox = pmd.load_file(str(top_p), xyz=gmx["gro"])

    orig_lines = orig_itp.splitlines()
    atoms_idx = next(i for i, ln in enumerate(orig_lines) if ln.strip().startswith("1 "))
    bonds_hdr = next(i for i, ln in enumerate(orig_lines) if ln.strip() == "[ bonds ]")
    bonds_idx = next(
        i
        for i in range(bonds_hdr + 1, len(orig_lines))
        if orig_lines[i].strip() and not orig_lines[i].strip().startswith(";")
    )
    atoms_bstate = f"{orig_lines[atoms_idx].split(';')[0].rstrip()}  C  9.0  99.0"
    bonds_bstate = f"{orig_lines[bonds_idx].split(';')[0].rstrip()}  0.15  400000.0"

    for cond_hdr, err_match in (
        ("#ifdef NEVER_DEFINED", "unsupported preprocessor directive"),
        ("#ifdef FLEXIBLE", r"unsupported B-state columns in \[atoms\]"),
        ("#ifndef FLEXIBLE", "unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional"),
    ):
        patched = list(orig_lines)
        patched[atoms_idx] = f"{cond_hdr}\n[ system ]\n#endif\n{atoms_bstate}"
        itp_p.write_text("\n".join(patched) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    for cond_hdr, err_match in (
        ("#ifdef NEVER_DEFINED", "unsupported preprocessor directive"),
        ("#ifdef FLEXIBLE", r"unsupported B-state columns in \[bonds\]"),
        ("#ifndef FLEXIBLE", "unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional"),
    ):
        patched = list(orig_lines)
        patched[bonds_idx] = f"{cond_hdr}\n[ system ]\n#endif\n{bonds_bstate}"
        itp_p.write_text("\n".join(patched) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    # Rigid-vs-flexible conditional divergence on non-water gox must be rejected
    for idx_target, cond_block, err_match in (
        (
            bonds_hdr,
            "#ifndef FLEXIBLE\n[ bonds ]\n1 2 1 0.14 999999.0\n#endif",
            "unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional",
        ),
        (
            bonds_idx,
            "#ifndef FLEXIBLE\n1 2 1 0.14 999999.0\n#endif",
            "unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional",
        ),
        (
            bonds_hdr,
            "#ifdef FLEXIBLE\n[ bonds ]\n1 2 1 0.14 300000.0\n#endif",
            "unverified non-water/non-bonded directive/data inside FLEXIBLE conditional",
        ),
    ):
        patched = list(orig_lines)
        patched[idx_target] = f"{patched[idx_target]}\n{cond_block}"
        itp_p.write_text("\n".join(patched) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(top_p, converted=rt_gox)

    itp_p.write_text(orig_itp, encoding="utf-8")

    # 2. Duplicate [ settles ] separated by #ifdef NEVER_DEFINED / #ifdef FLEXIBLE / #ifndef FLEXIBLE
    #    [ moleculetype ] WAT 3 -> all rejected.
    wat_ok = _make_water_mol("WAT")
    wat_bad = _make_water_mol("WAT_BAD", r1=0.9572, r2=1.05)
    wat_conv = SimpleNamespace(
        molecules={"WAT": (wat_ok, 3), "WAT_BAD": (wat_bad, 3)}
    )
    pp_dir = tmp_path / "pp_cases"
    pp_dir.mkdir()

    for tag, cond_hdr, err_match in (
        ("never", "#ifdef NEVER_DEFINED", "unsupported preprocessor directive"),
        ("flex", "#ifdef FLEXIBLE", r"duplicate \[settles\] directive in moleculetype 'WAT'"),
        (
            "nflex",
            "#ifndef FLEXIBLE",
            "unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional",
        ),
    ):
        (pp_dir / f"dup_{tag}.top").write_text(
            "[ defaults ]\n1 2 yes 0.5 0.833333\n"
            "[ moleculetype ]\nWAT 3\n[ atoms ]\n; ok\n"
            "[ settles ]\n1 1 0.09572 0.15139\n"
            f"{cond_hdr}\n[ moleculetype ]\nWAT 3\n#endif\n"
            "[ settles ]\n1 1 0.09572 0.15139\n"
            "[ system ]\nwat\n[ molecules ]\nWAT 1\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(pp_dir / f"dup_{tag}.top", converted=wat_conv)

    # 3. Asymmetric WAT_BAD + valid WAT: injecting [ moleculetype ] WAT 3 inside WAT_BAD before [ settles ]
    for tag, cond_hdr, err_match in (
        ("never", "#ifdef NEVER_DEFINED", "unsupported preprocessor directive"),
        ("flex", "#ifdef FLEXIBLE", "SETTLE water geometry/ordering mismatch"),
        (
            "nflex",
            "#ifndef FLEXIBLE",
            "unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional",
        ),
    ):
        (pp_dir / f"asym_{tag}.top").write_text(
            "[ defaults ]\n1 2 yes 0.5 0.833333\n"
            "[ moleculetype ]\nWAT_BAD 3\n[ atoms ]\n; ok\n"
            f"{cond_hdr}\n[ moleculetype ]\nWAT 3\n#endif\n"
            "[ settles ]\n1 1 0.09572 0.15139\n"
            "[ system ]\nwat\n[ molecules ]\nWAT_BAD 1\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(pp_dir / f"asym_{tag}.top", converted=wat_conv)

    # 4. Runtime-empty [ settles ] with sole row wrapped in #ifdef NEVER_DEFINED / #ifdef FLEXIBLE / #ifndef FLEXIBLE
    for tag, cond_hdr, err_match in (
        ("never", "#ifdef NEVER_DEFINED", "unsupported preprocessor directive"),
        ("flex", "#ifdef FLEXIBLE", r"invalid \[settles\] row count"),
        ("nflex", "#ifndef FLEXIBLE", r"invalid \[settles\] row count"),
    ):
        (pp_dir / f"empty_settle_{tag}.top").write_text(
            "[ defaults ]\n1 2 yes 0.5 0.833333\n"
            "[ moleculetype ]\nWAT 3\n[ atoms ]\n; ok\n"
            "[ settles ]\n"
            f"{cond_hdr}\n1 1 0.09572 0.15139\n#endif\n"
            "[ system ]\nwat\n[ molecules ]\nWAT 1\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(
                pp_dir / f"empty_settle_{tag}.top", converted=wat_conv
            )

    # 5. Nested conditionals (inline & via #include), unterminated, unexpected #else/#endif,
    #    duplicate #else, and #ifdef FLEXIBLE extra -> all rejected; valid ParmEd water passes.
    (pp_dir / "nested_child.inc").write_text("#ifndef FLEXIBLE\n#endif\n", encoding="utf-8")
    for tag, body, err_match in (
        (
            "nested_inline",
            "#ifdef FLEXIBLE\n#ifndef FLEXIBLE\n#endif\n#endif\n",
            "nested preprocessor conditional",
        ),
        (
            "nested_include",
            '#ifdef FLEXIBLE\n#include "nested_child.inc"\n#endif\n',
            "nested preprocessor conditional",
        ),
        (
            "unterminated",
            "#ifdef FLEXIBLE\n[ bonds ]\n1 2 1 0.09572 462750.4\n",
            "unterminated preprocessor conditional",
        ),
        ("unexpected_else", "#else\n", "unexpected #else directive"),
        ("unexpected_endif", "#endif\n", "unexpected #endif directive"),
        ("dup_else", "#ifdef FLEXIBLE\n#else\n#else\n#endif\n", "unexpected #else directive"),
        ("extra_tok", "#ifdef FLEXIBLE extra\n#endif\n", "unsupported preprocessor directive"),
    ):
        (pp_dir / f"bad_{tag}.top").write_text(
            "[ defaults ]\n1 2 yes 0.5 0.833333\n"
            "[ moleculetype ]\nWAT 3\n[ atoms ]\n; ok\n"
            f"{body}"
            "[ system ]\nwat\n[ molecules ]\nWAT 1\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match=err_match):
            _verify_gromacs_topology_sections(pp_dir / f"bad_{tag}.top", converted=wat_conv)

    # Valid ParmEd-style #ifdef FLEXIBLE ... [ bonds ]/[ angles ] ... #else ... [ settles ] ... #endif water
    (pp_dir / "valid_parmed_water.top").write_text(
        "[ defaults ]\n1 2 yes 0.5 0.833333\n"
        "[ moleculetype ]\nWAT 3\n[ atoms ]\n; ok\n"
        "#ifdef FLEXIBLE\n"
        "[ bonds ]\n1 2 1 0.09572 462750.4\n1 3 1 0.09572 462750.4\n"
        "[ angles ]\n2 1 3 1 104.52 836.8\n"
        "#else\n"
        "[ settles ]\n1 1 0.09572 0.15139\n"
        "#endif\n"
        "[ exclusions ]\n1 2 3\n2 1 3\n3 1 2\n"
        "[ system ]\nwat\n[ molecules ]\nWAT 1\n",
        encoding="utf-8",
    )
    _verify_gromacs_topology_sections(pp_dir / "valid_parmed_water.top", converted=wat_conv)

    # Eligible water whose [ bonds ]/[ angles ] are wrapped in #ifdef FLEXIBLE ... #endif with NO #else [ settles ] -> rejected
    (pp_dir / "water_missing_rigid_settles.top").write_text(
        "[ defaults ]\n1 2 yes 0.5 0.833333\n"
        "[ moleculetype ]\nWAT 3\n[ atoms ]\n; ok\n"
        "#ifdef FLEXIBLE\n"
        "[ bonds ]\n1 2 1 0.09572 462750.4\n1 3 1 0.09572 462750.4\n"
        "[ angles ]\n2 1 3 1 104.52 836.8\n"
        "#endif\n"
        "[ exclusions ]\n1 2 3\n2 1 3\n3 1 2\n"
        "[ system ]\nwat\n[ molecules ]\nWAT 1\n",
        encoding="utf-8",
    )
    with pytest.raises(
        ValueError, match=r"FLEXIBLE-conditional water missing rigid \[settles\] for \['WAT'\]"
    ):
        _verify_gromacs_topology_sections(
            pp_dir / "water_missing_rigid_settles.top", converted=wat_conv
        )

    wat_parm = pmd.amber.AmberParm.from_structure(wat_ok)
    wat_prm = tmp_path / "wat_valid.parm7"
    wat_crd = tmp_path / "wat_valid.rst7"
    wat_parm.write_parm(str(wat_prm))
    wat_parm.write_rst7(str(wat_crd))
    wat_gmx = export_gromacs(wat_prm, wat_crd, tmp_path / "gmx_wat_ok/wat")
    assert wat_gmx["status"] == "exported"
    assert wat_gmx["roundtrip_verified"]["atoms"] == 3

    # 6. MAJOR 2: assemble_complex succeeds on uniform scee=1.0, scnb=1.0 AND scee=2.0, scnb=2.0
    #    Amber receptor+ligand topologies (with active 1-4 dihedrals), while export_gromacs still rejects them.
    def _build_chain_amber(
        tag: str, *, resname: str, scee: float, scnb: float, offset: tuple[float, float, float]
    ) -> tuple[Path, Path]:
        s = pmd.Structure()
        ct = pmd.AtomType("c3", 1, 12.01, 6)
        ct.set_lj_params(0.1094, 1.908)
        charges = [0.1, -0.1, -0.1, 0.1]
        for idx, q in enumerate(charges):
            atom = pmd.Atom(
                name=f"C{idx + 1}",
                type="c3",
                atomic_number=6,
                mass=12.01,
                charge=q,
            )
            atom.atom_type = ct
            s.add_atom(atom, resname, 1)
        bt = pmd.BondType(340.0, 1.52)
        at = pmd.AngleType(60.0, 109.5)
        dt = pmd.DihedralType(0.5, 3, 0.0, scee=scee, scnb=scnb)
        s.bond_types.append(bt)
        s.angle_types.append(at)
        s.dihedral_types.append(dt)
        for idx in range(3):
            s.bonds.append(pmd.Bond(s.atoms[idx], s.atoms[idx + 1], type=bt))
        for idx in range(2):
            s.angles.append(pmd.Angle(s.atoms[idx], s.atoms[idx + 1], s.atoms[idx + 2], type=at))
        s.dihedrals.append(
            pmd.Dihedral(
                s.atoms[0],
                s.atoms[1],
                s.atoms[2],
                s.atoms[3],
                type=dt,
                ignore_end=False,
                improper=False,
            )
        )
        s.coordinates = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.52, 0.0, 0.0],
                [2.03, 1.43, 0.0],
                [3.55, 1.43, 0.0],
            ]
        ) + np.asarray(offset)
        parm = pmd.amber.AmberParm.from_structure(s)
        ChRad(parm, "mbondi")
        prm = tmp_path / f"{tag}.parm7"
        crd = tmp_path / f"{tag}.rst7"
        parm.write_parm(str(prm))
        parm.write_rst7(str(crd))
        return prm, crd

    for label, scee_val, scnb_val in (("glycam_1_1", 1.0, 1.0), ("opls_2_2", 2.0, 2.0)):
        rec_p, rec_c = _build_chain_amber(
            f"{label}_rec", resname="REC", scee=scee_val, scnb=scnb_val, offset=(0.0, 0.0, 0.0)
        )
        lig_p, lig_c = _build_chain_amber(
            f"{label}_lig", resname="LIG", scee=scee_val, scnb=scnb_val, offset=(0.0, 0.0, 5.0)
        )
        rec_loaded = pmd.load_file(str(rec_p), xyz=str(rec_c))
        assert len(rec_loaded.dihedrals) == 1
        assert (not rec_loaded.dihedrals[0].improper) and (not rec_loaded.dihedrals[0].ignore_end)
        assert rec_loaded.dihedrals[0].type.scee == pytest.approx(scee_val)
        assert rec_loaded.dihedrals[0].type.scnb == pytest.approx(scnb_val)

        comp_out = assemble_complex(
            rec_p, rec_c, [lig_p], [lig_c], tmp_path / f"complex_{label}"
        )
        assert comp_out["status"] == "assembled"
        assert comp_out["parameterized"] is True
        assert comp_out["roundtrip_verified"]["complex"]["torsion_terms"] == 2
        for key in (
            "complex_topology",
            "complex_coordinates",
            "complex_pdb",
            "receptor_topology",
            "ligand_topology",
        ):
            assert Path(comp_out[key]).is_file()
        assert (tmp_path / f"complex_{label}/masks.json").is_file()
        assert (tmp_path / f"complex_{label}/manifest.json").is_file()
        comp_loaded = pmd.load_file(
            comp_out["complex_topology"], xyz=comp_out["complex_coordinates"]
        )
        assert all(
            d.type.scee == pytest.approx(scee_val) and d.type.scnb == pytest.approx(scnb_val)
            for d in comp_loaded.dihedrals
        )

        gmx_bad_dir = tmp_path / f"gmx_reject_{label}"
        with pytest.raises(ValueError, match="Non-standard Amber 1-4 scaling factors"):
            export_gromacs(rec_p, rec_c, gmx_bad_dir / "rec")
        assert not gmx_bad_dir.exists()
