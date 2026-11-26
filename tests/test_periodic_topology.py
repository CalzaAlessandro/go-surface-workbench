"""Tests for Block 5 direct periodic Amber, GROMACS, and OpenMM topology writers."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest

from go_builder.build import build_sheet
from go_builder.periodic_topology import (
    build_v3_parmed_structure,
    export_periodic_amber,
    export_periodic_gromacs,
    export_periodic_openmm,
    load_v3_frcmod_parameters,
    validate_periodic_model,
)
from go_builder.topology import _verify_supported_amber_terms

WORKBENCH_ROOT = Path(__file__).resolve().parents[1]


def _make_periodic_porous_anionic_model() -> dict:
    return build_sheet(
        {
            "name": "periodic_porous_anionic_test",
            "lx_nm": 3.5,
            "ly_nm": 3.5,
            "pbc": [True, True, False],
            "vacuum": 40.0,
            "composition": {"o_over_c": 0.25, "ph": 7.4},
            "holes": {"distribution": "mixed", "percent_vacancies": 5.0},
            "placement": {"mode": "irregular_island", "n_islands": 2, "branching": 0.7},
            "clar": True,
            "seed": 42,
        }
    )


def _mic_distance(p1: np.ndarray, p2: np.ndarray, cell: np.ndarray, pbc: list[bool]) -> float:
    inv_cell = np.linalg.inv(cell)
    frac = (p2 - p1) @ inv_cell
    for ax in range(3):
        if pbc[ax]:
            frac[ax] -= round(float(frac[ax]))
    return float(np.linalg.norm(frac @ cell))


def test_frcmod_loader_and_finite_reference_build() -> None:
    frc = load_v3_frcmod_parameters()
    assert ("ga", "ga") in frc["bonds"]
    assert ("ga", "gb", "gc") in frc["angles"]
    assert ("ga", "gb", "gc", "gm") in frc["dihe_specific"]
    assert ("ga", "ga") in frc["dihe_wildcard"]
    assert len(frc["impropers"]) >= 20

    finite_recipe = json.loads((WORKBENCH_ROOT / "examples" / "finite_GO.json").read_text())
    finite_model = build_sheet(finite_recipe)

    amber_parm, assigned, meta = build_v3_parmed_structure(finite_model, allow_finite=True)
    _verify_supported_amber_terms(amber_parm, target="GROMACS", allow_explicit_exclusions=False)

    expected_fc = sum(a.get("formal_charge", 0) for a in finite_model["atoms"])
    assert math.isclose(sum(a.charge for a in amber_parm.atoms), expected_fc, abs_tol=1e-8)
    assert meta["atom_count"] == len(finite_model["atoms"])
    assert meta["bond_count"] == len(finite_model["bonds"])

    # Every degree-3 trigonal sp2 center receives at most one improper and all bonds/angles/dihedrals are typed
    assert all(b.type is not None for b in amber_parm.bonds)
    assert all(a.type is not None for a in amber_parm.angles)
    assert all(d.type is not None for d in amber_parm.dihedrals)
    imp_centers = [d.atom3.idx for d in amber_parm.dihedrals if d.improper]
    assert len(imp_centers) == len(set(imp_centers)) == meta["improper_dihedral_count"]
    assert meta["improper_dihedral_count"] > 0

    gmx = pmd.gromacs.GromacsTopologyFile.from_structure(amber_parm)
    assert len(gmx.atoms) == len(amber_parm.atoms)


def test_periodic_gromacs_export_and_mic_verification(tmp_path: Path) -> None:
    per_recipe = json.loads((WORKBENCH_ROOT / "examples" / "periodic_GO.json").read_text())
    for label, model in [
        ("example_periodic", build_sheet(per_recipe)),
        ("porous_anionic", _make_periodic_porous_anionic_model()),
    ]:
        prefix = tmp_path / label / "go"
        res = export_periodic_gromacs(model, prefix)
        assert res["status"] == "exported"
        assert res["direct_periodic"] is True
        assert res["periodic_molecules"] == "yes"
        assert res["roundtrip_verified"]["periodic_box_verified"] is True
        assert res["roundtrip_verified"]["periodic_molecules"] is True
        assert res["roundtrip_verified"]["max_gro_mic_bond_err_A"] < 0.02

        top_path = Path(res["top"])
        itp_path = Path(res["itp"])
        gro_path = Path(res["gro"])
        mdp_path = Path(res["mdp"])
        assert top_path.is_file() and itp_path.is_file() and gro_path.is_file() and mdp_path.is_file()
        assert Path(res["manifest"]).is_file() and Path(res["provenance"]).is_file()

        assert "periodic-molecules = yes" in top_path.read_text()
        mdp_txt = mdp_path.read_text()
        assert "periodic-molecules       = yes" in mdp_txt
        assert "coulombtype              = PME" in mdp_txt
        assert "vdw-modifier             = None" in mdp_txt

        # Verify every cross-boundary bond maintains exact MIC distance in exported .gro
        reloaded = pmd.load_file(str(top_path), xyz=str(gro_path))
        manifest = json.loads(Path(res["manifest"]).read_text())
        assert manifest["atom_count"] == len(reloaded.atoms)
        _, _, meta = build_v3_parmed_structure(model)
        id_to_idx = meta["model_id_to_topology_index"]
        cell = np.asarray(model["cell"], dtype=float)
        pbc = [bool(x) for x in model["pbc"]]
        gro_coords = np.asarray(reloaded.coordinates, dtype=float)
        orig_coords = np.asarray([a["position"] for a in model["atoms"]], dtype=float)
        orig_idx = {int(a["id"]): i for i, a in enumerate(model["atoms"])}

        cross_bonds = [b for b in model["bonds"] if any(b.get("image", [0, 0, 0]))]
        assert len(cross_bonds) > 0
        for b in model["bonds"]:
            gi, gj = id_to_idx[int(b["a"])], id_to_idx[int(b["b"])]
            d_gro = _mic_distance(gro_coords[gi], gro_coords[gj], cell, pbc)
            oi, oj = orig_idx[int(b["a"])], orig_idx[int(b["b"])]
            img = np.asarray(b.get("image", [0, 0, 0]), dtype=float)
            d_orig = float(np.linalg.norm(orig_coords[oj] + img @ cell - orig_coords[oi]))
            assert abs(d_gro - d_orig) < 0.02
            assert 0.8 <= d_gro < 2.0


def test_periodic_openmm_xml_and_amber_exports(tmp_path: Path) -> None:
    model = _make_periodic_porous_anionic_model()
    expected_fc = sum(a.get("formal_charge", 0) for a in model["atoms"])

    omm_out = tmp_path / "openmm" / "periodic_go.xml"
    res_omm = export_periodic_openmm(model, omm_out, _force_pure_python=True)
    assert res_omm["status"] == "exported"
    assert res_omm["direct_periodic"] is True
    assert res_omm["nonbonded_method"] == "PME"
    assert res_omm["uses_periodic_bonded"] is True
    assert res_omm["exceptions_use_periodic"] is True
    assert math.isclose(res_omm["net_charge"], expected_fc, abs_tol=1e-8)

    pdb_path = Path(res_omm["pdb"])
    assert pdb_path.is_file() and Path(res_omm["inpcrd"]).is_file()
    assert "CRYST1" in pdb_path.read_text()

    root = ET.parse(omm_out).getroot()
    assert root.find("PeriodicBoxVectors") is not None
    forces = {f.attrib["name"]: f for f in root.findall("./Forces/Force")}
    assert forces["HarmonicBondForce"].attrib["usesPeriodic"] == "1"
    assert forces["HarmonicAngleForce"].attrib["usesPeriodic"] == "1"
    assert forces["PeriodicTorsionForce"].attrib["usesPeriodic"] == "1"
    nb = forces["NonbondedForce"]
    assert nb.attrib["method"] == "4"
    assert nb.attrib["exceptionsUsePeriodic"] == "1"

    # Non-periodic nonbonded method or unsafe cutoff must be rejected
    with pytest.raises(ValueError, match="nonbonded_method='PME'"):
        export_periodic_openmm(
            model, tmp_path / "bad_method" / "go.xml", nonbonded_method="NoCutoff", _force_pure_python=True
        )
    with pytest.raises(ValueError, match="strictly less than half"):
        export_periodic_openmm(
            model, tmp_path / "bad_cutoff" / "go.xml", cutoff_nm=2.5, _force_pure_python=True
        )

    # Periodic Amber prmtop/inpcrd export
    amb_dir = tmp_path / "amber"
    res_amb = export_periodic_amber(model, amb_dir)
    assert res_amb["status"] == "parameterized"
    assert res_amb["direct_periodic"] is True
    assert res_amb["ifbox"] == 1
    assert res_amb["atoms_per_molecule"] == [len(model["atoms"])]
    assert math.isclose(res_amb["net_charge"], expected_fc, abs_tol=1e-6)

    reloaded_amb = pmd.load_file(res_amb["topology"], xyz=res_amb["coordinates"])
    assert int(reloaded_amb.ptr("IFBOX")) == 1
    assert list(reloaded_amb.parm_data["ATOMS_PER_MOLECULE"]) == [len(model["atoms"])]


def test_periodic_safety_gates_custom_charges_and_overwrite(tmp_path: Path) -> None:
    # Minimum nx=2, ny=2 periodic sheet and squeezed box must fail MIC / 1-4 unambiguity gates
    tiny = build_sheet({"name": "tiny_2_2", "nx": 2, "ny": 2, "pbc": [True, True, False], "vacuum": 40.0})
    with pytest.raises(ValueError, match="MIC|half-box|winding|too small"):
        validate_periodic_model(tiny, require_periodic=True)

    per = build_sheet(json.loads((WORKBENCH_ROOT / "examples" / "periodic_GO.json").read_text()))
    squeezed = deepcopy(per)
    squeezed["cell"][1][1] = 2.5
    with pytest.raises(ValueError, match="MIC|half-box|winding|too small|disagrees"):
        validate_periodic_model(squeezed, require_periodic=True)

    # Non-lower-triangular rotated cell rejected cleanly
    rot = deepcopy(per)
    rot["cell"][0][1] = 1.0
    with pytest.raises(ValueError, match="lower-triangular"):
        validate_periodic_model(rot, require_periodic=True)

    # Finite model rejected by periodic-only exporters
    fin = build_sheet(json.loads((WORKBENCH_ROOT / "examples" / "finite_GO.json").read_text()))
    with pytest.raises(ValueError, match="non-periodic"):
        export_periodic_gromacs(fin, tmp_path / "fin_fail" / "go")

    # Priority-3 motif rejected cleanly
    p3 = deepcopy(per)
    p3["atoms"][0].setdefault("metadata", {})["role"] = "extended_phenolic_OH_O"
    p3["atoms"][0]["element"] = "O"
    with pytest.raises(ValueError):
        export_periodic_openmm(p3, tmp_path / "p3_fail" / "go.xml", _force_pure_python=True)

    # Custom charges accepted and conserved across periodic exporters
    fc = sum(a.get("formal_charge", 0) for a in per["atoms"])
    custom_q = np.zeros(len(per["atoms"]), dtype=float)
    custom_q[0] = 0.25
    custom_q[1] = -0.25
    q_path = tmp_path / "custom_q.json"
    q_path.write_text(json.dumps({"charges_e": custom_q.tolist()}), encoding="utf-8")
    res_q = export_periodic_gromacs(per, tmp_path / "custom_gmx" / "go", charge_model="custom", charges_json=q_path)
    assert math.isclose(res_q["net_charge"], fc, abs_tol=1e-9)

    # Silent overwrite refused
    with pytest.raises(FileExistsError):
        export_periodic_gromacs(per, tmp_path / "custom_gmx" / "go")


def test_cli_periodic_model_workflows_and_conflicts(tmp_path: Path) -> None:
    per = build_sheet(json.loads((WORKBENCH_ROOT / "examples" / "periodic_GO.json").read_text()))
    model_path = tmp_path / "periodic_model.json"
    model_path.write_text(json.dumps(per), encoding="utf-8")

    # GROMACS via 07_export_gromacs.py --model
    gmx_prefix = tmp_path / "cli_gmx" / "sheet"
    r_gmx = subprocess.run(
        [sys.executable, str(WORKBENCH_ROOT / "07_export_gromacs.py"), "--model", str(model_path), "--output-prefix", str(gmx_prefix), "--json"],
        capture_output=True,
        text=True,
        check=True,
    )
    out_gmx = json.loads(r_gmx.stdout)
    assert out_gmx["direct_periodic"] is True
    assert Path(out_gmx["top"]).is_file()

    # OpenMM via 08_export_openmm.py --model
    omm_xml = tmp_path / "cli_omm" / "sheet.xml"
    r_omm = subprocess.run(
        [sys.executable, str(WORKBENCH_ROOT / "08_export_openmm.py"), "--model", str(model_path), "--output", str(omm_xml), "--json"],
        capture_output=True,
        text=True,
        check=True,
    )
    out_omm = json.loads(r_omm.stdout)
    assert out_omm["direct_periodic"] is True
    assert out_omm["nonbonded_method"] == "PME"
    assert omm_xml.is_file()

    # Amber via 03_prepare_amber.py --model --direct-periodic
    amb_dir = tmp_path / "cli_amb"
    r_amb = subprocess.run(
        [sys.executable, str(WORKBENCH_ROOT / "03_prepare_amber.py"), "--model", str(model_path), "--output-dir", str(amb_dir), "--direct-periodic", "--json"],
        capture_output=True,
        text=True,
        check=True,
    )
    out_amb = json.loads(r_amb.stdout)
    assert out_amb["direct_periodic"] is True
    assert (amb_dir / "go.prmtop").is_file()

    # Conflicting --model + --topology rejected cleanly
    bad = subprocess.run(
        [
            sys.executable,
            str(WORKBENCH_ROOT / "07_export_gromacs.py"),
            "--model",
            str(model_path),
            "--topology",
            str(amb_dir / "go.prmtop"),
            "--coordinates",
            str(amb_dir / "go.inpcrd"),
            "--output-prefix",
            str(tmp_path / "bad_cli" / "go"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert bad.returncode != 0
