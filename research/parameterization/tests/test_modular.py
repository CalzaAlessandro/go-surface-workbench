"""L0 modular constraints: synthetic charge fixtures are not fitted science."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
NAMES = [
    "pyrene_pristine",
    "pyrene_diol_cis",
    "pyrene_diol_trans",
    "pyrene_diol_epoxide",
    "pyrene_dione",
    "pyrene_COOH",
    "pyrene_COO",
]


def models(names=NAMES):
    return [
        json.loads((ROOT / "models" / name / "model.json").read_text())
        for name in names
    ]


def test_full_graph_partition_scaffold_zero_and_independent_block_constraints():
    from go_param.modular import build_problem, partition_model

    p = partition_model(models(["pyrene_pristine"])[0])
    assert len(p["groups"]) == 10
    assert len(p["scaffold_atom_ids"]) == 6
    covered = [i for group in p["groups"] for i in group["atom_ids"]] + p[
        "scaffold_atom_ids"
    ]
    assert len(covered) == len(set(covered)) == 26
    problem = build_problem(models())
    assert len(problem["roles"]) == 17
    assert problem["constraint_rank"] == 6
    assert problem["free_dimensions_after_constraints"] == 11
    assert {c["block_type"]: c["target_charge_e"] for c in problem["constraints"]} == {
        "GRT": 0.0,
        "OH": 0.0,
        "EPX": 0.0,
        "CAB": 0.0,
        "COOH": 0.0,
        "COO": -1.0,
    }
    assert problem["status"] == "PILOT_LIMITED_CONTEXT"


def test_uncovered_chemistry_infeasible_totals_and_failed_epoxide_are_rejected():
    from go_param.modular import build_problem, partition_model

    model = models(["pyrene_pristine"])[0]
    bad = copy.deepcopy(model)
    bad["atoms"][0]["role"] = "unknown_C"
    with pytest.raises(ValueError):
        partition_model(bad)
    bad = copy.deepcopy(model)
    bad["formal_charge"] = 1
    bad["atoms"][0]["formal_charge"] = 1
    with pytest.raises(ValueError, match="formal|charge"):
        partition_model(bad)
    with pytest.raises(ValueError, match="excluded"):
        build_problem(models(["pyrene_epoxide"]))


def test_resp_cards_have_zero_frozen_background_and_global_constraints_only_once():
    from go_param.modular import build_problem, make_resp_input

    problem = build_problem(models(["pyrene_pristine", "pyrene_diol_cis"]))
    records = [
        {"model_id": name, "grid_points": 100 + i * 10}
        for name in ("pyrene_pristine", "pyrene_diol_cis")
        for i in range(3)
    ]
    text, plan = make_resp_input(problem, records)
    assert "iqopt = 2" in text
    assert len(plan["molecules"]) == 6
    assert all(m["header_charge"] == -99 for m in plan["molecules"])
    assert len(plan["global_constraints"]) == 2
    assert plan["constraint_rank"] == 2
    for molecule in plan["molecules"]:
        for atom in molecule["atoms"]:
            if atom["role"] == "scaffold_C":
                assert atom["ivary"] == -1 and atom["initial_charge_e"] == 0.0
    assert plan["intermolecular_equivalences"]
    assert plan["stage2"]["status"] == "noop_no_methyl_methylene"
    with pytest.raises(ValueError, match="orientation"):
        make_resp_input(problem, records[:-1])


def test_constrained_bcc_projection_is_shared_neutral_and_keeps_raw_evidence():
    from go_param.modular import build_problem, evaluate_library, project_am1_bcc

    ms = models()
    problem = build_problem(ms)
    raw = {
        m["model_id"]: {
            "atoms": [
                {"atom_id": a["id"], "charge_e": 0.03 * (a["id"] % 7 - 3)}
                for a in m["atoms"]
            ]
        }
        for m in ms
    }
    preserved = copy.deepcopy(raw)
    result = project_am1_bcc(problem, raw)
    assert raw == preserved
    assert (
        result["method"] == "AM1_BCC_constrained_L2_block_projection_NOT_native_AM1_BCC"
    )
    assert result["raw_evidence"] == raw
    evaluation = evaluate_library(problem, result["shared_role_charges_e"])
    assert all(abs(m["charge_error_e"]) < 1e-9 for m in evaluation["models"])
    assert all(abs(g["charge_error_e"]) < 1e-9 for g in evaluation["groups"])
    assert result["shared_role_charges_e"]["scaffold_C"] == 0.0
    assert np.isfinite(list(result["shared_role_charges_e"].values())).all()


def test_external_resp_fixture_only_validates_io_not_real_fitting(tmp_path):
    from go_param.modular import build_problem, make_resp_input, validate_fitted_charges

    problem = build_problem(models(["pyrene_pristine"]))
    _text, plan = make_resp_input(
        problem, [{"model_id": "pyrene_pristine", "grid_points": 100}] * 3
    )
    q = [
        -0.15 if a["role"] == "edge_C" else 0.15 if a["role"] == "edge_H" else 0.0
        for m in plan["molecules"]
        for a in m["atoms"]
    ]
    checked = validate_fitted_charges(problem, plan, q)
    assert checked["shared_role_charges_e"] == {
        "edge_C": -0.15,
        "edge_H": 0.15,
        "scaffold_C": 0.0,
    }
    q[0] += 0.01
    with pytest.raises(ValueError, match="equivalence|constraint"):
        validate_fitted_charges(problem, plan, q)


def test_registry_covers_indexed_special_units_without_claiming_terminal_fits():
    from go_param.modular import BLOCK_REGISTRY, build_problem, current_library

    assert set(BLOCK_REGISTRY) == {
        "GGG",
        "GRT",
        "HYD",
        "EPX",
        "EPH",
        "EPT",
        "CAB",
        "CBX",
        "CBD",
    }
    assert BLOCK_REGISTRY["EPH"]["current_role_charges_e"]["EPH_C"] == 0.053
    assert BLOCK_REGISTRY["CAB"]["current_role_charges_e"]["carbonyl_C"] == 0.581
    problem = build_problem(models())
    assert problem["unfitted_terminal_blocks"] == ["EPH", "EPT"]
    control = current_library(problem)
    assert all(
        abs(m["charge_error_e"]) < 1e-12 for m in control["evaluation"]["models"]
    )


def test_synthetic_external_resp_fixture_preserves_inputs_and_reports_matched_control(
    tmp_path,
):
    import sys

    from go_param.modular import fit_modular

    model = models(["pyrene_pristine"])[0]
    mp = tmp_path / "model.json"
    mp.write_text(json.dumps(model))
    paths = []
    xyz = np.array([a["position"] for a in model["atoms"]]) / 0.529177210903
    for j, rotation in enumerate(
        (
            np.eye(3),
            np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]),
        )
    ):
        p = tmp_path / f"synthetic_{j}.esp"
        lines = [f"{len(xyz):5d}{1:5d}"]
        lines += [
            " " * 16 + "".join(f"{v:16.7E}" for v in row) for row in xyz @ rotation
        ]
        lines += ["".join(f"{v:16.7E}" for v in (0.001 * (j + 1), 20.0, 20.0, 20.0))]
        p.write_text("\n".join(lines) + "\n")
        paths.append(str(p))
    q = [
        -0.15 if a["role"] == "edge_C" else 0.15 if a["role"] == "edge_H" else 0.0
        for a in model["atoms"]
    ] * 3
    fixture = tmp_path / "synthetic_resp"
    fixture.write_text(
        "#!"
        + sys.executable
        + "\nfrom pathlib import Path\nq="
        + repr(q)
        + "\n"
        + "Path('resp.out').write_text('SYNTHETIC FIXTURE ONLY\\nConvergence in 1 iterations\\n')\n"
        + "Path('resp.punch').write_text('SYNTHETIC NO OPTIMIZATION\\n')\n"
        + "Path('resp.charges').write_text('\\n'.join(''.join(f'{v:10.6f}' for v in q[i:i+8]) for i in range(0,len(q),8))+'\\n')\n"
    )
    fixture.chmod(0o755)
    result = fit_modular(
        [
            {
                "model_path": str(mp),
                "esp_paths": paths,
                "scope": "synthetic_fixture_not_science",
            }
        ],
        tmp_path / "out",
        resp=fixture,
    )
    assert result["status"] == "PILOT_LIMITED_CONTEXT"
    assert result["final_parameter_validation"] is False
    assert len(result["metrics"]["current"]) == 3
    assert len(result["metrics"]["modular"]) == 3
    assert result["unfitted_terminal_blocks"] == ["EPH", "EPT"]
    assert (tmp_path / "out/current_library.json").exists()
