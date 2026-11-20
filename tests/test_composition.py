"""Composition preferences must preserve atom ratios and chemical capacities."""

from collections import Counter

import pytest

from go_builder.build import build_sheet
from go_builder.composition import resolve_composition


def test_C_over_O_is_carbon_over_oxygen_and_group_preference_changes_counts():
    base = {"nx": 6, "ny": 6, "seed": 7}
    hydroxyl, report = resolve_composition(base, 5.0, {"OH": 4, "epoxide": 1})
    epoxide, other = resolve_composition(base, 5.0, {"OH": 1, "epoxide": 4})
    assert hydroxyl["groups"]["OH"] > epoxide["groups"]["OH"]
    assert hydroxyl["groups"]["epoxide"] < epoxide["groups"]["epoxide"]
    for recipe, summary in [(hydroxyl, report), (epoxide, other)]:
        model = build_sheet(recipe)
        counts = Counter(a["element"] for a in model["atoms"])
        assert summary["predicted_C_over_O"] == pytest.approx(counts["C"] / counts["O"])
        assert abs(counts["O"] - counts["C"] / 5.0) <= 0.5 + 1e-8


def test_carboxyl_carbon_is_counted_and_zero_weights_disable_groups():
    recipe, report = resolve_composition(
        {"nx": 6, "ny": 6, "seed": 2},
        [65.4, 65.6],
        {"OH": 0, "epoxide": 0, "carbonyl": 0, "COOH": 1, "COO": 0},
    )
    assert recipe["groups"]["COOH"] > 0
    assert all(v == 0 for k, v in recipe["groups"].items() if k != "COOH")
    model = build_sheet(recipe)
    counts = Counter(a["element"] for a in model["atoms"])
    assert 65.4 <= counts["C"] / counts["O"] <= 65.6
    assert counts["C"] == report["scaffold_carbon_atoms"] + recipe["groups"]["COOH"]


def test_scalar_target_accounts_for_two_oxygen_atoms_per_carboxyl_group():
    recipe, report = resolve_composition(
        {"nx": 6, "ny": 6}, 20, {"OH": 0, "epoxide": 0, "carbonyl": 0, "COOH": 1, "COO": 0}
    )
    assert recipe["groups"]["COOH"] == 3
    assert abs(report["predicted_oxygen_atoms"] - report["predicted_carbon_atoms"] / 20) <= 1


def test_infeasible_interval_does_not_change_composition_request():
    with pytest.raises(ValueError, match="No integer composition"):
        resolve_composition(
            {"nx": 3, "ny": 3, "pbc": [True, True, False]},
            [1.0, 1.1],
            {"OH": 0, "epoxide": 0, "COOH": 1},
        )


@pytest.mark.parametrize(
    "preferences",
    [{"OH": -1}, {"OH": float("nan")}, {"ether": 1}, {"OH": 0, "epoxide": 0}, {"OH": True}],
)
def test_invalid_preferences_fail(preferences):
    with pytest.raises(ValueError):
        resolve_composition({"nx": 6, "ny": 6}, 4, preferences)


def test_reproducible_and_does_not_modify_input():
    base = {"nx": 6, "ny": 6, "seed": 5}
    before = dict(base)
    first = resolve_composition(base, 4, {"OH": 3})
    assert first == resolve_composition(base, 4, {"OH": 3}) and base == before


def test_existing_explicit_group_counts_are_not_silently_overwritten():
    with pytest.raises(ValueError, match="groups"):
        resolve_composition({"nx": 6, "ny": 6, "groups": {"OH": 2}}, 4)


def test_ph_henderson_hasselbalch_partitions_carboxyl_groups_and_documents_pka():
    base = {"nx": 8, "ny": 8, "seed": 3}
    acidic, rep_low = resolve_composition(base, 5.0, {"OH": 2, "epoxide": 2, "carboxyl": 2}, ph=2.0, pka=4.5)
    neutral, rep_mid = resolve_composition(base, 5.0, {"OH": 2, "epoxide": 2, "COOH": 2}, ph=4.5, pka=4.5)
    basic, rep_high = resolve_composition(base, 5.0, {"OH": 2, "epoxide": 2, "carboxyl": 2}, ph=9.0, pka=4.5)

    assert acidic["groups"]["COOH"] > acidic["groups"]["COO"]
    assert basic["groups"]["COO"] > basic["groups"]["COOH"]
    assert neutral["groups"]["COOH"] > 0 and neutral["groups"]["COO"] > 0
    assert rep_low["protonation"]["mode"] == "henderson_hasselbalch"
    assert rep_low["protonation"]["assumed_pka"] == 4.5
    assert "3-6" in rep_low["protonation"]["caveat"] or "3" in rep_low["protonation"]["caveat"]
    assert rep_mid["protonation"]["target_deprotonation_fraction"] == pytest.approx(0.5)
    assert rep_high["expected_formal_charge"] == -basic["groups"]["COO"]

    # Resulting recipe builds cleanly with exact charge
    model = build_sheet(basic)
    assert sum(a.get("formal_charge", 0) for a in model["atoms"]) == -basic["groups"]["COO"]


def test_explicit_deprotonation_fraction_and_validation():
    base = {"nx": 6, "ny": 6, "seed": 4}
    half, rep = resolve_composition(
        base, 5.0, {"OH": 1, "epoxide": 1, "COOH": 1, "COO": 1}, deprotonation_fraction=0.5
    )
    assert rep["protonation"]["mode"] == "explicit_deprotonation_fraction"
    assert rep["protonation"]["assumed_pka"] is None
    assert rep["protonation"]["target_deprotonation_fraction"] == 0.5
    assert half["groups"]["COOH"] > 0 and half["groups"]["COO"] > 0

    full_deprot, _ = resolve_composition(
        base, 5.0, {"OH": 1, "carboxyl": 1}, deprotonation_fraction=1.0
    )
    assert full_deprot["groups"]["COOH"] == 0
    assert full_deprot["groups"]["COO"] > 0

    with pytest.raises(ValueError, match="not both"):
        resolve_composition(base, 5.0, {"carboxyl": 1}, ph=7.0, deprotonation_fraction=0.5)
    with pytest.raises(ValueError, match="between 0 and 1"):
        resolve_composition(base, 5.0, {"carboxyl": 1}, deprotonation_fraction=1.2)
    with pytest.raises(ValueError, match="requires positive COOH, COO, or carboxyl"):
        resolve_composition(base, 5.0, {"OH": 1, "epoxide": 1}, ph=7.0)
    with pytest.raises(ValueError, match="requires ph or deprotonation_fraction"):
        resolve_composition(base, 5.0, {"OH": 1, "carboxyl": 1})


def test_priority2_edge_epoxides_in_resolve_composition_and_backward_compat():
    from go_builder.model import check_geometry_integrity

    base = {"nx": 6, "ny": 6, "seed": 7}
    # Legacy preferences with omitted P2 weights keep P2 counts at 0
    legacy_recipe, _ = resolve_composition(base, 5.0, {"OH": 2, "epoxide": 1})
    assert legacy_recipe["groups"]["epoxide_eph"] == 0
    assert legacy_recipe["groups"]["epoxide_ept"] == 0

    # Positive weights on epoxide_eph and epoxide_ept produce valid buildable sheets
    p2_recipe, report = resolve_composition(
        base,
        6.0,
        {"OH": 1, "epoxide": 1, "epoxide_eph": 2, "epoxide_ept": 2},
    )
    assert p2_recipe["groups"]["epoxide_eph"] > 0
    assert p2_recipe["groups"]["epoxide_ept"] > 0
    model = build_sheet(p2_recipe)
    assert check_geometry_integrity(model)["passed"] is True
    counts = Counter(a["element"] for a in model["atoms"])
    assert report["predicted_C_over_O"] == pytest.approx(counts["C"] / counts["O"])
    assert abs(counts["O"] - counts["C"] / 6.0) <= 0.5 + 1e-8


def test_o_over_c_scalar_and_interval_and_report_fields():
    base = {"nx": 8, "ny": 8, "seed": 5}
    _, rep = resolve_composition(base, o_over_c=0.35)
    assert rep["requested_O_over_C"] == 0.35
    assert rep["requested_C_over_O"] == pytest.approx(1.0 / 0.35)
    assert rep["predicted_O_over_C"] == pytest.approx(
        rep["predicted_oxygen_atoms"] / rep["predicted_carbon_atoms"]
    )
    assert rep["composition_mode"] == "balanced"
    assert rep["group_priority"] is None
    assert abs(rep["predicted_O_over_C"] - 0.35) < 0.01

    _, rep_iv = resolve_composition(base, o_over_c=[0.25, 0.35])
    assert rep_iv["requested_O_over_C"] == [0.25, 0.35]
    assert 0.25 - 1e-8 <= rep_iv["predicted_O_over_C"] <= 0.35 + 1e-8


def test_priority_list_and_carboxyl_with_ph():
    base = {"nx": 8, "ny": 8, "seed": 7}
    r_oh, rep_oh = resolve_composition(base, o_over_c=0.35, priority=["OH", "epoxide"])
    r_ep, rep_ep = resolve_composition(base, o_over_c=0.35, priority=["epoxide", "OH"])
    assert rep_oh["group_priority"] == ["OH", "epoxide"]
    assert rep_ep["group_priority"] == ["epoxide", "OH"]
    assert r_oh["groups"]["OH"] > r_ep["groups"]["OH"]
    assert r_ep["groups"]["epoxide"] > r_oh["groups"]["epoxide"]

    r_carb, rep_carb = resolve_composition(
        base, o_over_c=0.35, priority=["OH", "epoxide", "carboxyl"], ph=8.0
    )
    assert r_carb["groups"]["COO"] > r_carb["groups"]["COOH"]
    assert rep_carb["protonation"]["ph"] == 8.0


def test_random_mode_reproducible_distinct_across_seeds_and_enabled_groups():
    base = {"nx": 8, "ny": 8, "seed": 1}
    r1, rep1 = resolve_composition(base, o_over_c=0.35, mode="random", seed=10)
    r1_dup, _ = resolve_composition(base, o_over_c=0.35, mode="random", seed=10)
    assert r1 == r1_dup
    assert rep1["composition_mode"] == "random"

    seen = set()
    for s in (10, 20, 30, 40, 50):
        rec, rep = resolve_composition(base, o_over_c=0.35, mode="random", seed=s)
        assert abs(rep["predicted_O_over_C"] - 0.35) < 0.015
        seen.add(tuple(rec["groups"][g] for g in ("OH", "epoxide", "COOH", "COO")))
    assert len(seen) >= 3

    # enabled_groups restricts active groups
    r_sub, _ = resolve_composition(
        base, o_over_c=0.35, mode="random", seed=10, enabled_groups=["OH", "carbonyl"]
    )
    assert r_sub["groups"]["epoxide"] == 0
    assert r_sub["groups"]["COOH"] == 0
    assert r_sub["groups"]["COO"] == 0
    assert r_sub["groups"]["OH"] > 0 and r_sub["groups"]["carbonyl"] > 0

    # Periodic sheet in random mode defaults to basal groups only
    per_rec, per_rep = resolve_composition(
        {"nx": 6, "ny": 6, "pbc": [True, True, False], "seed": 3},
        o_over_c=0.25,
        mode="random",
    )
    assert per_rec["groups"]["COOH"] == 0 and per_rec["groups"]["COO"] == 0
    assert abs(per_rep["predicted_O_over_C"] - 0.25) < 0.015


def test_new_composition_validation_errors():
    base = {"nx": 6, "ny": 6}
    with pytest.raises(ValueError, match="not both"):
        resolve_composition(base, c_over_o=4.0, o_over_c=0.25)
    with pytest.raises(ValueError, match="must be provided"):
        resolve_composition(base)
    with pytest.raises(ValueError, match="not both"):
        resolve_composition(base, o_over_c=0.25, preferences={"OH": 1}, priority=["OH"])
    with pytest.raises(ValueError, match="mode"):
        resolve_composition(base, o_over_c=0.25, mode="stochastic")
    with pytest.raises(ValueError, match="duplicate"):
        resolve_composition(base, o_over_c=0.25, priority=["OH", "OH"])
    with pytest.raises(ValueError, match="enabled_groups"):
        resolve_composition(base, o_over_c=0.25, enabled_groups=[])


def test_enabled_groups_vs_ph_deprotonation_conflicts_and_extreme_ph_no_overflow():
    base = {"nx": 6, "ny": 6, "seed": 3}

    # Conflicting enabled_groups vs ph / deprotonation_fraction must raise ValueError
    with pytest.raises(ValueError, match="excludes 'COO'"):
        resolve_composition(
            base, o_over_c=0.25, enabled_groups=["OH", "COOH"], deprotonation_fraction=1.0
        )
    with pytest.raises(ValueError, match="excludes 'COO'"):
        resolve_composition(base, o_over_c=0.25, enabled_groups=["OH", "COOH"], ph=7.0)
    with pytest.raises(ValueError, match="excludes 'COOH'"):
        resolve_composition(
            base, o_over_c=0.25, enabled_groups=["OH", "COO"], deprotonation_fraction=0.0
        )
    with pytest.raises(ValueError, match="excludes 'COOH'"):
        resolve_composition(base, o_over_c=0.25, enabled_groups=["OH", "COO"], ph=2.0)

    # Boundary-compatible single carboxyl species are allowed without ph
    rec_cooh, _ = resolve_composition(
        base, o_over_c=0.25, enabled_groups=["OH", "COOH"], deprotonation_fraction=0.0
    )
    assert rec_cooh["groups"]["COO"] == 0
    rec_coo, _ = resolve_composition(
        base, o_over_c=0.25, enabled_groups=["OH", "COO"], deprotonation_fraction=1.0
    )
    assert rec_coo["groups"]["COOH"] == 0

    # Extreme pH values must not trigger OverflowError
    rec_low, rep_low = resolve_composition(
        base, o_over_c=0.25, enabled_groups=["OH", "carboxyl"], ph=-400.0
    )
    assert rep_low["protonation"]["target_deprotonation_fraction"] == 0.0
    assert rec_low["groups"]["COO"] == 0
    rec_high, rep_high = resolve_composition(
        base, o_over_c=0.25, enabled_groups=["OH", "carboxyl"], ph=400.0
    )
    assert rep_high["protonation"]["target_deprotonation_fraction"] == 1.0
    assert rec_high["groups"]["COOH"] == 0

    # Balanced mode zeros zero-capacity edge groups on periodic sheets without distorting basal counts
    per_bal, _ = resolve_composition(
        {"nx": 6, "ny": 6, "pbc": [True, True, False]},
        o_over_c=0.25,
        mode="balanced",
        enabled_groups=["OH", "epoxide", "COOH"],
    )
    assert per_bal["groups"]["COOH"] == 0
    assert per_bal["groups"]["OH"] > 0 and per_bal["groups"]["epoxide"] > 0


def test_prepare_composition_script_supports_o_over_c_ph_priority_and_toml(tmp_path):
    import json
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    cfg_toml = tmp_path / "comp.toml"
    cfg_toml.write_text(
        """
[recipe]
nx = 6
ny = 6
seed = 5

[composition]
o_over_c = 0.32
priority = ["OH", "epoxide", "COOH", "COO"]
ph = 7.4
mode = "balanced"
""".strip()
        + "\n"
    )
    out_recipe = tmp_path / "resolved.json"
    out_report = tmp_path / "resolved.comp.json"
    out_model = tmp_path / "model.json"

    res = subprocess.run(
        [
            sys.executable,
            str(root / "prepare_composition.py"),
            "--input",
            str(cfg_toml),
            "--output",
            str(out_recipe),
            "--report",
            str(out_report),
            "--model",
            str(out_model),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    summary = json.loads(res.stdout)
    assert summary["status"] == "BUILT_AND_COMPOSITION_VERIFIED"
    assert abs(summary["O_over_C"] - 0.32) < 0.02
    assert summary["groups"]["COO"] > summary["groups"].get("COOH", 0)
    assert out_recipe.is_file() and out_report.is_file() and out_model.is_file()
