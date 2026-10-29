import copy
import json
import subprocess

import pytest

from go_builder.qm import compare_esp, espgen_to_reference, prepare_qm


def hydrogen():
    return {
        "atoms": [
            {
                "id": 7,
                "element": "H",
                "position": [0.0, 0.0, 0.0],
                "charge": 1.0,
                "formal_charge": 1,
            }
        ],
        "bonds": [],
        "pbc": [False, False, False],
        "metadata": {},
    }


def test_qm_writes_explicit_finite_single_point_and_never_runs(tmp_path):
    result = prepare_qm(hydrogen(), tmp_path / "qm", total_charge=1)
    assert result["executed"] is False
    text = (tmp_path / "qm/reference.gjf").read_text()
    assert "HF/6-31G*" in text and "Pop=MK" in text
    assert "\n1 1\n" in text
    with pytest.raises(FileExistsError):
        prepare_qm(hydrogen(), tmp_path / "qm", total_charge=1)


@pytest.mark.parametrize("mutator", ["periodic", "image", "open", "spin"])
def test_qm_rejects_invalid_or_uncapped_models_before_writing(tmp_path, mutator):
    m = hydrogen()
    charge = 1
    if mutator == "periodic":
        m["pbc"][0] = True
    elif mutator == "image":
        m["bonds"] = [{"a": 7, "b": 7, "order": 1, "image": [1, 0, 0]}]
    elif mutator == "open":
        m["atoms"][0]["element"] = "C"
        charge = 0
    else:
        charge = 0
    with pytest.raises(ValueError):
        prepare_qm(m, tmp_path / mutator, total_charge=charge)
    assert not (tmp_path / mutator).exists()


def reference_file(tmp_path, **overrides):
    # Unit +1 at origin gives V=1/2 Hartree/e at r=2 bohr.
    ref = {
        "atoms": hydrogen()["atoms"],
        "total_charge": 1,
        "position_unit": "bohr",
        "potential_unit": "hartree/e",
        "grid": [
            {"position": [2, 0, 0], "potential": 0.5},
            {"position": [0, 4, 0], "potential": 0.25},
        ],
        "dipole": [0, 0, 0],
        "dipole_unit": "debye",
    }
    ref.update(overrides)
    p = tmp_path / "ref.json"
    p.write_text(json.dumps(ref))
    return p


def test_esp_analytic_coulomb_and_units(tmp_path):
    result = compare_esp(reference_file(tmp_path), {"exact": hydrogen()}, tmp_path / "esp")
    assert result["candidates"]["exact"]["rmse_hartree_per_e"] == pytest.approx(0, abs=1e-14)
    assert result["candidates"]["exact"]["dipole_error_debye"] == pytest.approx(0)
    # Same physical coordinates/grid expressed in Angstroms.
    p = reference_file(
        tmp_path,
        position_unit="angstrom",
        grid=[{"position": [1.058354421806, 0, 0], "potential": 0.5}],
    )
    result = compare_esp(p, {"exact": hydrogen()}, tmp_path / "esp2")
    assert result["candidates"]["exact"]["rmse_hartree_per_e"] == pytest.approx(0, abs=1e-12)


def test_esp_dipole_and_nonzero_error_are_analytic(tmp_path):
    # Neutral dipole: +1 at x=1 Angstrom and -1 at x=0, V at x=2 is .5 atomic unit.
    model = {
        "atoms": [
            dict(id=1, element="H", position=[0, 0, 0], charge=-1),
            dict(id=2, element="H", position=[0.529177210903, 0, 0], charge=1),
        ]
    }
    p = reference_file(
        tmp_path,
        atoms=model["atoms"],
        total_charge=0,
        position_unit="angstrom",
        grid=[{"position": [1.058354421806, 0, 0], "potential": 0.4}],
        dipole=[2.541746473, 0, 0],
    )
    r = compare_esp(p, {"d": model}, tmp_path / "out")["candidates"]["d"]
    assert r["rmse_hartree_per_e"] == pytest.approx(0.1, abs=1e-10)
    assert r["dipole_debye"][0] == pytest.approx(2.541746473, abs=1e-8)


@pytest.mark.parametrize(
    "fault", ["geometry", "charge", "element", "ids", "nan", "units", "singular"]
)
def test_esp_rejects_unmatched_candidates_and_invalid_grids(tmp_path, fault):
    model = copy.deepcopy(hydrogen())
    p = reference_file(tmp_path)
    if fault == "geometry":
        model["atoms"][0]["position"][0] = 1.0
    if fault == "charge":
        model["atoms"][0]["charge"] = 0.0
    if fault == "element":
        model["atoms"][0]["element"] = "C"
    if fault == "ids":
        model["atoms"][0]["id"] = 8
    if fault == "nan":
        model["atoms"][0]["charge"] = float("nan")
    if fault == "units":
        p = reference_file(tmp_path, potential_unit="unknown")
    if fault == "singular":
        p = reference_file(tmp_path, grid=[{"position": [0, 0, 0], "potential": 0}])
    with pytest.raises(ValueError):
        compare_esp(p, {"bad": model}, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_espgen_bridge_preserves_atomic_units_and_checks_geometry(tmp_path):
    esp = tmp_path / "ref.esp"
    esp.write_text("    1    2\n  0.0 0.0 0.0\n  .5 2.0 0.0 0.0\n  .25 0.0 4.0 0.0\n")
    p = tmp_path / "reference.json"
    espgen_to_reference(esp, hydrogen(), p, total_charge=1)
    ref = json.loads(p.read_text())
    assert ref["grid"][0]["potential"] == 0.5
    assert ref["position_unit"] == "bohr"
    scored = compare_esp(p, {"h": hydrogen()}, tmp_path / "score")
    assert scored["candidates"]["h"]["rmse_hartree_per_e"] == pytest.approx(0)
    altered = hydrogen()
    altered["atoms"][0]["position"] = [1, 0, 0]
    with pytest.raises(ValueError):
        espgen_to_reference(esp, altered, tmp_path / "bad.json", total_charge=1)
    with pytest.raises(FileExistsError):
        espgen_to_reference(esp, hydrogen(), p, total_charge=1)


def test_qm_accepts_aromatic_capped_benzene(tmp_path):
    import math

    atoms, bonds = [], []
    for i in range(6):
        angle = i * math.pi / 3
        atoms.append(
            dict(
                id=i + 1,
                element="C",
                formal_charge=0,
                position=[1.4 * math.cos(angle), 1.4 * math.sin(angle), 0],
            )
        )
        atoms.append(
            dict(
                id=i + 7,
                element="H",
                formal_charge=0,
                position=[2.48 * math.cos(angle), 2.48 * math.sin(angle), 0],
            )
        )
        bonds.extend([dict(a=i + 1, b=(i + 1) % 6 + 1, order=1.5), dict(a=i + 1, b=i + 7, order=1)])
    model = dict(atoms=atoms, bonds=bonds)
    result = prepare_qm(model, tmp_path / "benzene", total_charge=0)
    assert result["atom_order_ids"] == [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]


def test_espgen_fixed_width_header_handles_five_digit_grid_count(tmp_path):
    esp = tmp_path / "many.esp"
    esp.write_text(f"{1:5d}{10000:5d}{0:5d}\n 0 0 0\n" + " .5 2 0 0\n" * 10000)
    ref = espgen_to_reference(esp, hydrogen(), tmp_path / "many.json", total_charge=1)
    assert len(ref["grid"]) == 10000


def test_esp_charged_dipole_respects_origin_and_candidate_hash(tmp_path):
    ref = reference_file(tmp_path, dipole_origin=[1, 0, 0], dipole=[-2.541746473, 0, 0])
    result = compare_esp(ref, {"h": hydrogen()}, tmp_path / "origin")
    score = result["candidates"]["h"]
    assert score["dipole_error_debye"] == pytest.approx(0, abs=1e-8)
    assert len(score["charge_geometry_sha256"]) == 64


def test_reference_periodic_state_is_rejected(tmp_path):
    ref = reference_file(tmp_path, pbc=[True, False, False])
    with pytest.raises(ValueError):
        compare_esp(ref, {"h": hydrogen()}, tmp_path / "periodic")


def test_resp_helper_refuses_existing_punch_before_launching_external_tools(tmp_path):
    out = tmp_path / "qm"
    prepare_qm(hydrogen(), out, total_charge=1)
    punch = out / "resp-stage1.punch"
    punch.write_text("preserve this fit\n")
    run = subprocess.run(
        ["/bin/bash", str(out / "prepare_resp.sh")], capture_output=True, text=True
    )
    assert run.returncode != 0
    assert "Refusing overwrite: resp-stage1.punch" in run.stderr
    assert punch.read_text() == "preserve this fit\n"
