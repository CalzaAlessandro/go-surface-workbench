"""Tests for irregular_island (Lerf-Klinowski) placement and Clar closed-shell filtering."""

from __future__ import annotations

import json

import pytest

from go_builder.build import build_sheet
from go_builder.cli import _resolve_build_recipe, main, parser_for
from go_builder.config import (
    load_config,
    render_annotated_json_template,
    render_annotated_toml_template,
)
from go_builder.electronic import (
    analyze_aromatic_electronic_structure,
    compute_oxidation_morphology_metrics,
)
from go_builder.errors import PlacementError
from go_builder.model import check_geometry_integrity


def test_placement_backwards_compatibility_random_and_compact_island():
    r1 = build_sheet(
        {
            "lx_nm": 2.5,
            "ly_nm": 2.5,
            "seed": 17,
            "groups": {"OH": 6, "epoxide": 5, "COOH": 2},
            "placement": "random",
        }
    )
    r2 = build_sheet(
        {
            "lx_nm": 2.5,
            "ly_nm": 2.5,
            "seed": 17,
            "groups": {"OH": 6, "epoxide": 5, "COOH": 2},
            "placement": "random",
        }
    )
    assert [a["position"] for a in r1["atoms"]] == [a["position"] for a in r2["atoms"]]
    assert "clar" not in r1["metadata"]

    isl = build_sheet(
        {
            "lx_nm": 3.0,
            "ly_nm": 3.0,
            "seed": 17,
            "groups": {"OH": 10, "epoxide": 8},
            "placement": {"mode": "island", "n_islands": 2},
        }
    )
    assert isl["metadata"]["placement"]["mode"] == "island"
    assert "style" not in isl["metadata"]["placement"]
    assert "branching" not in isl["metadata"]["placement"]


def test_irregular_island_api_forms_and_validation():
    m1 = build_sheet(
        {
            "lx_nm": 3.2,
            "ly_nm": 3.2,
            "seed": 101,
            "composition": {"o_over_c": 0.35, "ph": 7.4},
            "placement": {"mode": "irregular_island", "n_islands": 2, "branching": 0.7},
        }
    )
    m2 = build_sheet(
        {
            "lx_nm": 3.2,
            "ly_nm": 3.2,
            "seed": 101,
            "composition": {"o_over_c": 0.35, "ph": 7.4},
            "placement": {
                "mode": "island",
                "style": "irregular",
                "n_islands": 2,
                "branching": 0.7,
            },
        }
    )
    assert m1["metadata"]["placement"]["mode"] == "irregular_island"
    assert m1["metadata"]["placement"]["style"] == "irregular"
    assert m1["metadata"]["placement"]["branching"] == 0.7
    assert m2["metadata"]["placement"]["mode"] == "island"
    assert m2["metadata"]["placement"]["style"] == "irregular"
    assert [a["position"] for a in m1["atoms"]] == [a["position"] for a in m2["atoms"]]

    bad_specs = [
        {"mode": "unknown"},
        {"mode": "island", "style": "fractal"},
        {"mode": "island", "branching": 0.5},
        {"mode": "island", "style": "compact", "branching": 0.5},
        {"mode": "irregular_island", "style": "compact"},
        {"mode": "irregular_island", "branching": 1.2},
        {"mode": "irregular_island", "branching": -0.1},
    ]
    for bad in bad_specs:
        with pytest.raises(ValueError):
            build_sheet(
                {"lx_nm": 2.0, "ly_nm": 2.0, "groups": {"OH": 2}, "placement": bad}
            )


def test_irregular_island_separates_from_compact_island_and_random():
    metrics_by_mode: dict[str, list[dict]] = {"random": [], "island": [], "irregular_island": []}
    for seed in range(10, 20):
        base_recipe = {
            "lx_nm": 3.6,
            "ly_nm": 3.6,
            "seed": seed,
            "composition": {"o_over_c": 0.35, "ph": 7.4},
        }
        m_rnd = build_sheet({**base_recipe, "placement": "random"})
        m_cmp = build_sheet({**base_recipe, "placement": {"mode": "island", "n_islands": 2}})
        m_irr = build_sheet(
            {
                **base_recipe,
                "placement": {"mode": "irregular_island", "n_islands": 2, "branching": 0.75},
            }
        )
        metrics_by_mode["random"].append(
            compute_oxidation_morphology_metrics(m_rnd, compute_huckel=False)
        )
        metrics_by_mode["island"].append(
            compute_oxidation_morphology_metrics(m_cmp, compute_huckel=False)
        )
        metrics_by_mode["irregular_island"].append(
            compute_oxidation_morphology_metrics(m_irr, compute_huckel=False)
        )

    def _mean(mode: str, key: str) -> float:
        vals = [row[key] for row in metrics_by_mode[mode]]
        return sum(vals) / len(vals)

    # vs random: irregular_island forms few connected domains instead of scattered fragments
    assert _mean("irregular_island", "n_oxidized_clusters") < _mean(
        "random", "n_oxidized_clusters"
    ) - 3.0
    assert _mean("irregular_island", "clustered_group_fraction_ge4") > _mean(
        "random", "clustered_group_fraction_ge4"
    )
    assert _mean("irregular_island", "intact_sp2_ring_fraction") > _mean(
        "random", "intact_sp2_ring_fraction"
    ) + 0.10
    assert _mean("irregular_island", "largest_intact_sp2_ring_cluster") > _mean(
        "random", "largest_intact_sp2_ring_cluster"
    ) + 10.0
    assert _mean("irregular_island", "singlet_basal_fraction") < _mean(
        "random", "singlet_basal_fraction"
    )
    assert _mean("irregular_island", "mean_basal_group_nn_distance") < _mean(
        "random", "mean_basal_group_nn_distance"
    )

    # vs compact island: irregular_island has higher asphericity & lower convex hull fill
    assert _mean("irregular_island", "mean_domain_asphericity") > _mean(
        "island", "mean_domain_asphericity"
    ) + 0.08
    assert _mean("irregular_island", "mean_domain_hull_fill_fraction") < _mean(
        "island", "mean_domain_hull_fill_fraction"
    ) - 0.25
    assert _mean("irregular_island", "mean_domain_convex_circularity") < _mean(
        "island", "mean_domain_convex_circularity"
    ) - 0.05


def test_irregular_island_pbc_and_holes_integrity():
    m = build_sheet(
        {
            "lx_nm": 3.5,
            "ly_nm": 3.5,
            "pbc": [True, True, False],
            "seed": 33,
            "holes": {"distribution": "small", "percent_vacancies": 4.0},
            "composition": {"o_over_c": 0.30, "ph": 7.4},
            "placement": {"mode": "irregular_island", "n_islands": 2, "branching": 0.65},
        }
    )
    assert check_geometry_integrity(m)["passed"]
    morph = compute_oxidation_morphology_metrics(m, compute_huckel=False)
    assert morph["clustered_group_fraction_ge4"] > 0.70


@pytest.mark.parametrize(
    "placement",
    [
        "random",
        "dispersed",
        {"mode": "island", "n_islands": 2},
        {"mode": "irregular_island", "n_islands": 2, "branching": 0.65},
    ],
)
def test_clar_closed_shell_enforced_across_modes(placement):
    m = build_sheet(
        {
            "lx_nm": 3.0,
            "ly_nm": 3.0,
            "seed": 42,
            "composition": {"o_over_c": 0.30, "ph": 7.4},
            "placement": placement,
            "clar": {"enabled": True, "min_huckel_gap": 0.0, "max_attempts": 32},
        }
    )
    assert m["metadata"]["clar"]["enabled"] is True
    elec = analyze_aromatic_electronic_structure(m, compute_huckel=True)
    assert elec["bipartite"] is True
    assert elec["sublattice_imbalance"] == 0
    assert elec["kekule_perfect"] is True


def test_clar_min_huckel_gap_enforced_on_gapped_nanoflake():
    m = build_sheet(
        {
            "shape": "hexagonal",
            "edge_type": "armchair",
            "diameter_nm": 3.0,
            "seed": 7,
            "composition": {"o_over_c": 0.25, "ph": 7.4},
            "placement": {"mode": "irregular_island", "n_islands": 1, "branching": 0.65},
            "clar": {"enabled": True, "min_huckel_gap": 0.05, "max_attempts": 32},
        }
    )
    elec = analyze_aromatic_electronic_structure(m, compute_huckel=True)
    assert elec["sublattice_imbalance"] == 0
    assert elec["kekule_perfect"] is True
    assert elec["huckel_gap_beta"] is not None
    assert elec["huckel_gap_beta"] >= 0.05


def test_clar_explicit_counts_and_impossible_rejections():
    m_ok = build_sheet(
        {
            "lx_nm": 2.4,
            "ly_nm": 2.4,
            "seed": 5,
            "groups": {"OH": 6, "epoxide": 5, "carbonyl": 2},
            "placement": {"mode": "irregular_island", "branching": 0.6},
            "clar": True,
        }
    )
    elec = analyze_aromatic_electronic_structure(m_ok, compute_huckel=True)
    assert elec["sublattice_imbalance"] == 0
    assert elec["kekule_perfect"] is True

    # Odd number of single-sublattice groups on bipartite even-carbon scaffold => impossible
    with pytest.raises(PlacementError, match="sublattice balance"):
        build_sheet(
            {
                "nx": 6,
                "ny": 6,
                "seed": 1,
                "groups": {"OH": 3, "epoxide": 2},
                "clar": True,
            }
        )

    # Unrealistically large Huckel pi-gap must fail cleanly
    with pytest.raises(PlacementError, match="Clar closed-shell"):
        build_sheet(
            {
                "nx": 6,
                "ny": 6,
                "seed": 1,
                "groups": {"OH": 4, "epoxide": 2},
                "clar": {"enabled": True, "min_huckel_gap": 10.0, "max_attempts": 4},
            }
        )


def test_electronic_diagnostics_single_oh_radical_vs_pristine():
    pristine = build_sheet({"nx": 4, "ny": 4, "seed": 1, "groups": {}})
    ep_pristine = analyze_aromatic_electronic_structure(pristine, compute_huckel=True)
    assert ep_pristine["bipartite"] is True
    assert ep_pristine["sublattice_imbalance"] == 0
    assert ep_pristine["kekule_perfect"] is True
    assert ep_pristine["huckel_gap_beta"] > 0.0

    radical = build_sheet({"nx": 4, "ny": 4, "seed": 1, "groups": {"OH": 1}})
    ep_rad = analyze_aromatic_electronic_structure(radical, compute_huckel=True)
    assert ep_rad["sublattice_imbalance"] == 1
    assert ep_rad["kekule_perfect"] is False
    assert ep_rad["huckel_gap_beta"] == 0.0


def test_cli_irregular_island_clar_and_config_templates(tmp_path):
    parser = parser_for("build")
    args = parser.parse_args(
        [
            "--size",
            "3.0",
            "3.0",
            "--oc",
            "0.30",
            "--place",
            "irregular_island",
            "--islands",
            "2",
            "--branching",
            "0.8",
            "--clar",
            "--min-huckel-gap",
            "0.03",
            "--seed",
            "7",
        ]
    )
    recipe, _ = _resolve_build_recipe(args)
    assert recipe["placement"] == {
        "mode": "irregular_island",
        "n_islands": 2,
        "branching": 0.8,
    }
    assert recipe["clar"] == {"enabled": True, "min_huckel_gap": 0.03}

    out_json = tmp_path / "irreg_clar.json"
    rc = main(
        "build",
        [
            "--size",
            "2.8",
            "2.8",
            "--oc",
            "0.30",
            "--place",
            "island",
            "--island-style",
            "irregular",
            "--branching",
            "0.7",
            "--clar",
            "--seed",
            "9",
            "--out",
            str(out_json),
        ],
    )
    assert rc == 0
    built = json.loads(out_json.read_text(encoding="utf-8"))
    assert built["metadata"]["placement"]["style"] == "irregular"
    assert built["metadata"]["clar"]["enabled"] is True
    assert built["metadata"]["clar"]["kekule_perfect"] is True

    # Invalid CLI combinations fail cleanly
    for bad_argv in (
        ["--place", "random", "--branching", "0.5"],
        ["--place", "island", "--island-style", "compact", "--branching", "0.5"],
        ["--place", "irregular_island", "--island-style", "compact"],
        ["--place", "irregular_island", "--branching", "1.5"],
        ["--clar", "--min-huckel-gap", "-0.1"],
    ):
        with pytest.raises(ValueError):
            _resolve_build_recipe(parser.parse_args(bad_argv))

    # Annotated templates remain valid and loadable
    toml_path = tmp_path / "tpl.toml"
    json_path = tmp_path / "tpl.json"
    toml_path.write_text(render_annotated_toml_template(), encoding="utf-8")
    json_path.write_text(render_annotated_json_template(), encoding="utf-8")
    assert load_config(toml_path)["placement"]["mode"] == "random"
    assert load_config(json_path)["placement"]["mode"] == "random"


def test_geometry_integrity_ester_lactone_co_single_bond():
    # Planar/tetrahedral methyl formate H-C(=O)-O-CH3: C2-O4 single bond (1.37 A, deg(O4)=2)
    # must pass as C-O_sp3, while stretching C2=O3 (deg(O3)=1) to 1.40 A must fail C=O range.
    atoms = [
        {"id": 1, "element": "H", "position": [-0.95, -0.54, 0.0], "formal_charge": 0},
        {"id": 2, "element": "C", "position": [0.0, 0.0, 0.0], "formal_charge": 0},
        {"id": 3, "element": "O", "position": [0.0, 1.21, 0.0], "formal_charge": 0},
        {"id": 4, "element": "O", "position": [1.18, -0.69, 0.0], "formal_charge": 0},
        {"id": 5, "element": "C", "position": [2.42, 0.04, 0.0], "formal_charge": 0},
        {"id": 6, "element": "H", "position": [3.22, -0.69, 0.0], "formal_charge": 0},
        {"id": 7, "element": "H", "position": [2.50, 0.67, 0.89], "formal_charge": 0},
        {"id": 8, "element": "H", "position": [2.50, 0.67, -0.89], "formal_charge": 0},
    ]
    bonds = [
        {"a": 1, "b": 2, "order": 1},
        {"a": 2, "b": 3, "order": 2},
        {"a": 2, "b": 4, "order": 1},
        {"a": 4, "b": 5, "order": 1},
        {"a": 5, "b": 6, "order": 1},
        {"a": 5, "b": 7, "order": 1},
        {"a": 5, "b": 8, "order": 1},
    ]
    model = {
        "schema_version": 1,
        "name": "ester_bond_check",
        "pbc": [False, False, False],
        "cell": None,
        "atoms": atoms,
        "bonds": bonds,
    }
    geom_ok = check_geometry_integrity(model)
    assert geom_ok["passed"] is True, geom_ok

    stretched = {
        **model,
        "atoms": [
            {**a, "position": [0.0, 1.40, 0.0]} if a["id"] == 3 else a for a in atoms
        ],
    }
    geom_bad = check_geometry_integrity(stretched)
    assert geom_bad["passed"] is False
    assert geom_bad["bond_length_violations"][0]["motif"] == "C=O"


@pytest.mark.parametrize(
    "placement",
    ["random", "dispersed", "island", "irregular_island"],
)
def test_cli_clar_defaults_on_across_modes_and_opt_out(tmp_path, placement):
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    out_json = tmp_path / f"cli_default_{placement}.json"
    res = subprocess.run(
        [
            sys.executable,
            str(root / "01_build_GO.py"),
            "--size",
            "2.8",
            "2.8",
            "--oc",
            "0.25",
            "--place",
            placement,
            "--seed",
            "11",
            "--out",
            str(out_json),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    model = json.loads(out_json.read_text(encoding="utf-8"))
    assert model["metadata"]["clar"]["enabled"] is True
    assert model["metadata"]["clar"]["kekule_perfect"] is True

    odd_fail = subprocess.run(
        [
            sys.executable,
            str(root / "01_build_GO.py"),
            "--nx",
            "3",
            "ny_unused",
        ]
        if False
        else [
            sys.executable,
            str(root / "01_build_GO.py"),
            "--nx",
            "3",
            "--ny",
            "3",
            "--groups",
            "OH=1",
            "--out",
            str(tmp_path / "odd_fail.json"),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert odd_fail.returncode != 0
    assert "Clar" in odd_fail.stdout

    odd_ok_path = tmp_path / "odd_noclar.json"
    odd_ok = subprocess.run(
        [
            sys.executable,
            str(root / "01_build_GO.py"),
            "--nx",
            "3",
            "--ny",
            "3",
            "--groups",
            "OH=1",
            "--no-clar",
            "--out",
            str(odd_ok_path),
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    assert odd_ok.returncode == 0, odd_ok.stdout + odd_ok.stderr
    odd_model = json.loads(odd_ok_path.read_text(encoding="utf-8"))
    assert "clar" not in odd_model["metadata"]


def test_edge_carboxyl_graph_spacing_across_cooh_and_coo():
    from collections import deque
    from go_builder.model import adjacency

    model = build_sheet(
        {"nx": 6, "ny": 6, "seed": 17, "groups": {"COOH": 2, "COO": 2}}
    )
    by_id = {int(a["id"]): a for a in model["atoms"]}
    adj = adjacency(model)
    scaffold_ids = {
        int(a["id"])
        for a in model["atoms"]
        if a["element"] == "C" and a.get("metadata", {}).get("role") == "sheet_carbon"
    }
    acid_anchors = [
        int(g["site_ids"][0])
        for g in model["metadata"]["groups"]
        if g.get("type") in ("COOH", "COO")
    ]
    assert len(acid_anchors) == 4
    for i, src in enumerate(acid_anchors):
        dist = {src: 0}
        q = deque([src])
        while q:
            u = q.popleft()
            for v, _, _ in adj[u]:
                if v in scaffold_ids and v not in dist and by_id[v]["element"] == "C":
                    dist[v] = dist[u] + 1
                    q.append(v)
        for dst in acid_anchors[i + 1 :]:
            assert dist[dst] >= 3, f"Adjacent acid anchors {src} and {dst} at graph dist {dist[dst]}"

