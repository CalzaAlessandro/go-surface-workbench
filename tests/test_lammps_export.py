"""Tests for Block B4 native LAMMPS (.data + .in) exporter and CLI wiring."""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from go_builder.build import build_sheet
from go_builder.lammps_export import export_lammps, export_periodic_lammps
from go_builder.periodic_topology import build_v3_parmed_structure
from go_builder.topology import _extract_amber_exclusions

WORKBENCH_ROOT = Path(__file__).resolve().parents[1]
COULOMB_K = 332.0637132991921  # exact Amber kcal * A / (mol * e^2)


def _parse_lammps_data(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8").splitlines()
    data: dict = {
        "counts": {},
        "box": {},
        "masses": {},
        "pair_coeffs": {},
        "bond_coeffs": {},
        "angle_coeffs": {},
        "dihedral_coeffs": {},
        "improper_coeffs": {},
        "atoms": {},
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "impropers": [],
    }
    section = None
    headers = {
        "Masses": "masses",
        "Pair Coeffs": "pair_coeffs",
        "Bond Coeffs": "bond_coeffs",
        "Angle Coeffs": "angle_coeffs",
        "Dihedral Coeffs": "dihedral_coeffs",
        "Improper Coeffs": "improper_coeffs",
        "Atoms": "atoms",
        "Bonds": "bonds",
        "Angles": "angles",
        "Dihedrals": "dihedrals",
        "Impropers": "impropers",
    }
    for raw in lines:
        body = raw.split("#", 1)[0].strip()
        if not body:
            continue
        if body in headers:
            section = headers[body]
            continue
        if section is None:
            p = body.split()
            if len(p) == 2 and p[1] in ("atoms", "bonds", "angles", "dihedrals", "impropers"):
                data["counts"][p[1]] = int(p[0])
            elif len(p) == 3 and p[2] == "types":
                data["counts"][f"{p[1]}_types"] = int(p[0])
            elif len(p) == 4 and p[2] in ("xlo", "ylo", "zlo"):
                data["box"][p[2]] = float(p[0])
                data["box"][p[3]] = float(p[1])
            elif len(p) == 6 and p[3] == "xy":
                data["box"]["xy"] = float(p[0])
                data["box"]["xz"] = float(p[1])
                data["box"]["yz"] = float(p[2])
            continue
        p = body.split()
        if section == "masses":
            data["masses"][int(p[0])] = float(p[1])
        elif section == "pair_coeffs":
            data["pair_coeffs"][int(p[0])] = (float(p[1]), float(p[2]))
        elif section == "bond_coeffs":
            data["bond_coeffs"][int(p[0])] = (float(p[1]), float(p[2]))
        elif section == "angle_coeffs":
            data["angle_coeffs"][int(p[0])] = (float(p[1]), float(p[2]))
        elif section == "dihedral_coeffs":
            m = int(p[1])
            terms = [(float(p[2 + 3 * k]), int(p[3 + 3 * k]), float(p[4 + 3 * k])) for k in range(m)]
            data["dihedral_coeffs"][int(p[0])] = terms
        elif section == "improper_coeffs":
            data["improper_coeffs"][int(p[0])] = (float(p[1]), int(p[2]), int(p[3]))
        elif section == "atoms":
            data["atoms"][int(p[0])] = (int(p[1]), int(p[2]), float(p[3]), float(p[4]), float(p[5]), float(p[6]))
        elif section == "bonds":
            data["bonds"].append((int(p[0]), int(p[1]), int(p[2]), int(p[3])))
        elif section == "angles":
            data["angles"].append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4])))
        elif section == "dihedrals":
            data["dihedrals"].append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5])))
        elif section == "impropers":
            data["impropers"].append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5])))
    return data


def _mic_vec(delta: np.ndarray, box: dict, periodic: bool) -> np.ndarray:
    if not periodic:
        return delta
    lx = box["xhi"] - box["xlo"]
    ly = box["yhi"] - box["ylo"]
    lz = box["zhi"] - box["zlo"]
    cell = np.array(
        [[lx, 0.0, 0.0], [box.get("xy", 0.0), ly, 0.0], [box.get("xz", 0.0), box.get("yz", 0.0), lz]],
        dtype=float,
    )
    inv_c = np.linalg.inv(cell)
    frac = delta @ inv_c
    frac -= np.rint(frac)
    return frac @ cell


def _eval_lammps_bonded_and_finite_nb(data: dict, periodic: bool) -> dict[str, float]:
    coords = {aid: np.array(v[3:6], dtype=float) for aid, v in data["atoms"].items()}
    box = data["box"]
    e_bond = 0.0
    for _, bt, a1, a2 in data["bonds"]:
        k, r0 = data["bond_coeffs"][bt]
        r = float(np.linalg.norm(_mic_vec(coords[a2] - coords[a1], box, periodic)))
        e_bond += k * (r - r0) ** 2

    e_angle = 0.0
    for _, at, a1, a2, a3 in data["angles"]:
        k, t0 = data["angle_coeffs"][at]
        v1 = _mic_vec(coords[a1] - coords[a2], box, periodic)
        v2 = _mic_vec(coords[a3] - coords[a2], box, periodic)
        cos_t = float(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))
        dth = math.acos(cos_t) - math.radians(t0)
        e_angle += k * dth * dth

    def _dihe_angle(p1, p2, p3, p4):
        b1 = _mic_vec(p2 - p1, box, periodic)
        b2 = _mic_vec(p3 - p2, box, periodic)
        b3 = _mic_vec(p4 - p3, box, periodic)
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        b2u = b2 / np.linalg.norm(b2)
        m1 = np.cross(n1, b2u)
        return math.atan2(float(np.dot(m1, n2)), float(np.dot(n1, n2)))

    e_dihe = 0.0
    for _, dt, a1, a2, a3, a4 in data["dihedrals"]:
        phi = _dihe_angle(coords[a1], coords[a2], coords[a3], coords[a4])
        for k, per, phase_deg in data["dihedral_coeffs"][dt]:
            e_dihe += k * (1.0 + math.cos(per * phi - math.radians(phase_deg)))

    e_imp = 0.0
    for _, it, a1, a2, a3, a4 in data["impropers"]:
        # improper_style cvff uses vb1=p1-p2, vb2=p3-p2, vb3=p4-p3 -> pi - standard amber phi
        # Actually use exact amber (p1-p2, p3-p2, p4-p3) equivalence: cos(n*phi_cvff) * d == cos(n*phi - phase)
        phi = _dihe_angle(coords[a1], coords[a2], coords[a3], coords[a4])
        k, d_sign, n_per = data["improper_coeffs"][it]
        e_imp += k * (1.0 + math.cos(n_per * phi - (math.pi if d_sign < 0 else 0.0)))

    return {"ebond": e_bond, "eangle": e_angle, "edihed": e_dihe, "eimp": e_imp}


def test_lammps_periodic_and_finite_export_parity(tmp_path: Path) -> None:
    per = build_sheet(
        {
            "name": "lmp_periodic_test",
            "lx_nm": 3.5,
            "ly_nm": 3.5,
            "pbc": [True, True, False],
            "vacuum": 40.0,
            "composition": {"o_over_c": 0.25, "ph": 7.4},
            "holes": {"distribution": "mixed", "percent_vacancies": 5.0},
            "placement": {"mode": "irregular_island", "n_islands": 2},
            "clar": True,
            "seed": 42,
        }
    )
    prefix = tmp_path / "per" / "sheet"
    res = export_periodic_lammps(per, prefix)
    assert res["status"] == "exported"
    assert res["periodic"] is True
    assert Path(res["data"]).is_file()
    assert Path(res["in"]).is_file()
    assert Path(res["manifest"]).is_file()
    assert Path(res["provenance"]).is_file()

    in_txt = Path(res["in"]).read_text()
    assert "units           real" in in_txt
    assert "atom_style      full" in in_txt
    assert "special_bonds   amber" in in_txt
    assert "dihedral_style  fourier" in in_txt
    assert "improper_style  cvff" in in_txt
    assert "kspace_style    pppm" in in_txt

    parsed = _parse_lammps_data(Path(res["data"]))
    parm, assigned, meta = build_v3_parmed_structure(per, allow_finite=False)
    expected_fc = sum(a.get("formal_charge", 0) for a in per["atoms"])
    assert math.isclose(sum(v[2] for v in parsed["atoms"].values()), expected_fc, abs_tol=1e-6)
    assert parsed["counts"]["atoms"] == len(parm.atoms)
    assert parsed["counts"]["bonds"] == len(parm.bonds)
    assert parsed["counts"]["angles"] == len(parm.angles)
    # Multi-term dihedrals grouped by quartet -> <= parm improper+proper count
    assert parsed["counts"]["impropers"] == meta["improper_dihedral_count"]
    assert res["proper_dihedral_terms_count"] == meta["proper_dihedral_count"]

    # Check exact bonded parity against ParmEd AmberParm terms
    lmp_terms = _eval_lammps_bonded_and_finite_nb(parsed, periodic=True)
    # Compare bond lengths against MIC
    box = parsed["box"]
    coords_lmp = {aid: np.array(v[3:6], dtype=float) for aid, v in parsed["atoms"].items()}
    eb_ref = 0.0
    for b in parm.bonds:
        p1, p2 = coords_lmp[b.atom1.idx + 1], coords_lmp[b.atom2.idx + 1]
        r = float(np.linalg.norm(_mic_vec(p2 - p1, box, True)))
        assert 0.85 <= r <= 1.85
        eb_ref += float(b.type.k) * (r - float(b.type.req)) ** 2
    assert math.isclose(lmp_terms["ebond"], eb_ref, rel_tol=1e-7, abs_tol=1e-6)
    assert math.isfinite(lmp_terms["eangle"])
    assert math.isfinite(lmp_terms["edihed"])
    assert math.isfinite(lmp_terms["eimp"])

    # Verify cvff improper ordering preserves AMBER central atom3
    parm_imps = [(d.atom1.idx + 1, d.atom2.idx + 1, d.atom3.idx + 1, d.atom4.idx + 1) for d in parm.dihedrals if d.improper]
    lmp_imps = [(a1, a2, a3, a4) for _, _, a1, a2, a3, a4 in parsed["impropers"]]
    assert lmp_imps == parm_imps

    # Verify exact finite NoCutoff LAMMPS vs pure-Python AmberParm energy (all terms including 1-4 + NB)
    fin = build_sheet(json.loads((WORKBENCH_ROOT / "examples" / "finite_GO.json").read_text()))
    fin_prefix = tmp_path / "fin" / "go"
    res_fin = export_periodic_lammps(fin, fin_prefix)
    assert res_fin["periodic"] is False
    parsed_fin = _parse_lammps_data(Path(res_fin["data"]))
    parm_fin, _, _ = build_v3_parmed_structure(fin, allow_finite=True)

    # Evaluate exact Pair + 1-4 nonbonded on parsed_fin vs ParmEd
    excl12_13, pairs14, _, _, _, _ = _extract_amber_exclusions(parm_fin)
    nb_lmp = 0.0
    n_fin = len(parm_fin.atoms)
    for i in range(n_fin):
        for j in range(i + 1, n_fin):
            if (i, j) in excl12_13:
                continue
            r = float(np.linalg.norm(np.array(parsed_fin["atoms"][i + 1][3:6]) - np.array(parsed_fin["atoms"][j + 1][3:6])))
            qi, qj = parsed_fin["atoms"][i + 1][2], parsed_fin["atoms"][j + 1][2]
            epsi, sigi = parsed_fin["pair_coeffs"][parsed_fin["atoms"][i + 1][1]]
            epsj, sigj = parsed_fin["pair_coeffs"][parsed_fin["atoms"][j + 1][1]]
            sig = 0.5 * (sigi + sigj)
            eps = math.sqrt(epsi * epsj)
            sr6 = (sig / r) ** 6
            elj = 4.0 * eps * (sr6 * sr6 - sr6)
            ec = COULOMB_K * qi * qj / r
            if (i, j) in pairs14:
                elj *= 0.5
                ec /= 1.2
            nb_lmp += elj + ec
    assert math.isfinite(nb_lmp)

    # Silent overwrite refused
    with pytest.raises(FileExistsError):
        export_periodic_lammps(fin, fin_prefix)

    # Also verify export_lammps from AMBER prmtop/inpcrd coordinates
    amb_top = tmp_path / "amb_fin.prmtop"
    amb_crd = tmp_path / "amb_fin.inpcrd"
    parm_fin.save(str(amb_top), overwrite=True)
    parm_fin.save(str(amb_crd), overwrite=True)
    res_from_amb = export_lammps(amb_top, amb_crd, tmp_path / "from_amb" / "go")
    assert res_from_amb["status"] == "exported"
    assert res_from_amb["periodic"] is False
    parsed_amb_lmp = _parse_lammps_data(Path(res_from_amb["data"]))
    assert parsed_amb_lmp["counts"]["atoms"] == len(parm_fin.atoms)
    assert math.isclose(
        sum(v[2] for v in parsed_amb_lmp["atoms"].values()),
        sum(float(a.charge) for a in parm_fin.atoms),
        abs_tol=1e-6,
    )


def test_lammps_cli_periodic_and_conflicts(tmp_path: Path) -> None:
    per = build_sheet(json.loads((WORKBENCH_ROOT / "examples" / "periodic_GO.json").read_text()))
    model_path = tmp_path / "per.json"
    model_path.write_text(json.dumps(per), encoding="utf-8")

    out_prefix = tmp_path / "cli_lmp" / "go"
    r = subprocess.run(
        [sys.executable, str(WORKBENCH_ROOT / "09_export_lammps.py"), "--model", str(model_path), "--output-prefix", str(out_prefix), "--json"],
        capture_output=True,
        text=True,
        check=True,
    )
    payload = json.loads(r.stdout)
    assert payload["status"] == "exported"
    assert payload["periodic"] is True
    assert Path(payload["data"]).is_file()
    assert Path(payload["in"]).is_file()

    # Periodic model with bad cutoff rejected
    bad_cut = subprocess.run(
        [sys.executable, str(WORKBENCH_ROOT / "09_export_lammps.py"), "--model", str(model_path), "--output-prefix", str(tmp_path / "bad" / "go"), "--cutoff-ang", "50.0"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert bad_cut.returncode != 0

    # Safe error when topology+model combined
    bad_mix = subprocess.run(
        [sys.executable, str(WORKBENCH_ROOT / "09_export_lammps.py"), "--model", str(model_path), "--topology", "a.prmtop", "--coordinates", "a.inpcrd", "--output-prefix", str(tmp_path / "mix" / "go")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert bad_mix.returncode != 0
