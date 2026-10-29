"""Construction tests catch missing seam bonds, undershoot and invalid site reuse."""

from collections import Counter

import numpy as np
import pytest

from go_builder.build import build_sheet
from go_builder.model import validate_model


def test_finite_hydrogen_capped_connected_sheet_is_reproducible():
    m = build_sheet({"nx": 4, "ny": 4, "seed": 9})
    assert m == build_sheet({"nx": 4, "ny": 4, "seed": 9})
    assert Counter(a["element"] for a in m["atoms"])["H"] > 0
    assert validate_model(m)["connected_components"] == 1
    degrees = Counter(x for b in m["bonds"] for x in (b["a"], b["b"]))
    assert all(degrees[a["id"]] == 3 for a in m["atoms"] if a["element"] == "C")
    assert not any(any(b["image"]) for b in m["bonds"])
    assert m["metadata"]["geometry_status"] == "constructed_unrelaxed"


def test_periodic_tile_has_three_carbon_neighbors_and_true_seam_geometry():
    m = build_sheet({"nx": 3, "ny": 4, "pbc": [True, True, False]})
    assert len(m["atoms"]) == 48
    degrees = Counter(x for b in m["bonds"] for x in (b["a"], b["b"]))
    assert set(degrees.values()) == {3}
    atoms = {a["id"]: np.array(a["position"]) for a in m["atoms"]}
    assert any(any(b["image"]) for b in m["bonds"])
    for b in m["bonds"]:
        d = atoms[b["b"]] + np.array(b["image"]) @ np.array(m["cell"]) - atoms[b["a"]]
        assert np.linalg.norm(d) == pytest.approx(1.42)


def test_explicit_vacancy_and_hole_remove_original_ids_and_stay_connected():
    pristine = build_sheet({"nx": 8, "ny": 8})
    center = next(a for a in pristine["atoms"] if a["id"] == 145)["position"]
    m = build_sheet(
        {"nx": 8, "ny": 8, "vacancies": [145], "holes": [{"center": center[:2], "radius": 1.5}]}
    )
    assert 145 not in {a["id"] for a in m["atoms"]}
    assert len(m["metadata"]["removed_carbon_ids"]) >= 3
    assert validate_model(m)["connected_components"] == 1


def test_exact_functional_counts_charges_and_sites():
    groups = {"OH": 4, "epoxide": 3, "carbonyl": 2, "COOH": 2, "COO": 2}
    m = build_sheet({"nx": 10, "ny": 8, "groups": groups, "seed": 17})
    assert Counter(g["type"] for g in m["metadata"]["groups"]) == Counter(groups)
    assert sum(a.get("formal_charge", 0) for a in m["atoms"]) == -2
    sites = [s for g in m["metadata"]["groups"] for s in g["site_ids"]]
    assert len(sites) == len(set(sites))
    assert validate_model(m)["valid"]


def test_up_face_places_basal_oxygen_above_sheet():
    m = build_sheet({"nx": 5, "ny": 5, "groups": {"OH": 3, "epoxide": 2}, "faces": "up"})
    atoms = {a["id"]: a for a in m["atoms"]}
    for g in m["metadata"]["groups"]:
        assert atoms[g["oxygen_id"]]["position"][2] > 0


def test_group_count_objects_preserve_random_placement_and_seed():
    counts = {"OH": 4, "epoxide": 3, "carbonyl": 2, "COOH": 2, "COO": 2}
    recipe = {"nx": 10, "ny": 8, "seed": 17, "groups": counts}
    integer_model = build_sheet(recipe)
    object_model = build_sheet(
        {**recipe, "groups": {group: {"count": count} for group, count in counts.items()}}
    )
    # The submitted recipes differ, but geometry, group IDs and counts must not.
    integer_model["metadata"].pop("recipe")
    object_model["metadata"].pop("recipe")
    assert object_model == integer_model


@pytest.mark.parametrize(
    ("group", "sites", "expected_sites"),
    [
        ("OH", [58], [58]),
        ("epoxide", [[58, 59]], [58, 59]),
        ("carbonyl", [26], [26]),
        ("COOH", [26], [26]),
        ("COO", [26], [26]),
    ],
)
@pytest.mark.parametrize("include_count", [False, True])
def test_explicit_group_ids_place_requested_original_carbons(
    group, sites, expected_sites, include_count
):
    specification = {"id": sites}
    if include_count:
        specification["count"] = 1
    model = build_sheet(
        {"nx": 6, "ny": 6, "seed": 17, "faces": "down", "groups": {group: specification}}
    )
    assert len(model["metadata"]["groups"]) == 1
    placed = model["metadata"]["groups"][0]
    assert placed["type"] == group
    assert placed["site_ids"] == expected_sites
    if group in ("OH", "epoxide"):
        oxygen = next(a for a in model["atoms"] if a["id"] == placed["oxygen_id"])
        assert oxygen["position"][2] < 0
    assert validate_model(model)["valid"]


@pytest.mark.parametrize("seed", range(8))
def test_explicit_sites_are_reserved_before_other_groups_choose_random_sites(seed):
    model = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": seed,
            "groups": {
                "OH": {"id": [58]},
                "carbonyl": {"id": [26]},
                "epoxide": 3,
                "COOH": 1,
                "COO": 1,
            },
        }
    )
    groups = model["metadata"]["groups"]
    assert Counter(g["type"] for g in groups) == {
        "OH": 1,
        "carbonyl": 1,
        "epoxide": 3,
        "COOH": 1,
        "COO": 1,
    }
    assert next(g["site_ids"] for g in groups if g["type"] == "OH") == [58]
    assert next(g["site_ids"] for g in groups if g["type"] == "carbonyl") == [26]
    sites = [site for group in groups for site in group["site_ids"]]
    assert len(sites) == len(set(sites))


@pytest.mark.parametrize(
    "groups",
    [
        {"OH": {}},
        {"OH": {"count": 1.5}},
        {"OH": {"count": True}},
        {"OH": {"count": 2, "id": [58]}},
        {"OH": {"count": 0, "id": [58]}},
        {"OH": {"count": 1, "id": []}},
        {"OH": {"id": 58}},
        {"OH": {"id": [58.0]}},
        {"OH": {"id": [True]}},
        {"OH": {"id": [9999]}},
        {"OH": {"id": [1]}},  # Trimmed original carbon ID.
        {"OH": {"id": [26]}},  # Edge carbon is not a basal site.
        {"OH": {"id": [58, 58]}},
        {"OH": {"id": [58], "faces": "up"}},
        {"carbonyl": {"id": [58]}},
        {"COOH": {"id": [58]}},
        {"COO": {"id": [58]}},
        {"carbonyl": {"id": [26]}, "COOH": {"id": [26]}},
        {"OH": {"id": [58]}, "epoxide": {"id": [[58, 59]]}},
        {"epoxide": {"id": [58, 59]}},
        {"epoxide": {"id": [[58]]}},
        {"epoxide": {"id": [[58, 59, 60]]}},
        {"epoxide": {"id": [[58, 58]]}},
        {"epoxide": {"id": [[58, 59.0]]}},
        {"epoxide": {"id": [[58, 9999]]}},
        {"epoxide": {"id": [[25, 26]]}},
        {"epoxide": {"id": [[58, 60]]}},  # Basal, but not bonded.
        {"epoxide": {"id": [[58, 59], [59, 62]]}},
    ],
)
def test_invalid_explicit_group_specifications_fail(groups):
    with pytest.raises(ValueError):
        build_sheet({"nx": 6, "ny": 6, "groups": groups})


def test_explicit_site_removed_by_vacancy_fails():
    with pytest.raises(ValueError, match="(?i)(site|carbon|id)"):
        build_sheet({"nx": 6, "ny": 6, "vacancies": [58], "groups": {"OH": {"id": [58]}}})


def test_explicit_empty_site_list_requests_zero_groups():
    model = build_sheet({"nx": 6, "ny": 6, "groups": {"OH": {"id": []}}})
    assert model["metadata"]["groups"] == []
    assert model["metadata"]["events"][-1]["group_counts"]["OH"] == 0


def test_colliding_explicit_sites_fail_without_substituting_other_carbons():
    # Adjacent horizontal OH sites produce a nonbonded O-H overlap on one face.
    with pytest.raises(ValueError, match="128 attempts"):
        build_sheet({"nx": 6, "ny": 6, "faces": "up", "groups": {"OH": {"id": [57, 58]}}})


def test_infeasible_counts_report_requested_groups_and_available_capacity():
    with pytest.raises(ValueError) as caught:
        build_sheet({"nx": 6, "ny": 6, "groups": {"OH": 1000}})
    report = caught.value.as_report()
    assert report["status"] == "INCOMPLETE"
    assert report["requested"] == {"OH": 1000, "epoxide": 0, "carbonyl": 0, "COOH": 0, "COO": 0}
    assert report.get("achieved") is None
    assert report["attempts"] == 0
    assert report["capacity"]["eligible_basal_carbons"] == 98
    assert report["capacity"]["eligible_edge_carbons"] == 32


def test_exhausted_fixed_sites_report_best_tried_counts_without_returning_partial_model():
    with pytest.raises(ValueError) as caught:
        build_sheet({"nx": 6, "ny": 6, "faces": "up", "groups": {"OH": {"id": [57, 58]}}})
    report = caught.value.as_report()
    assert report["status"] == "INCOMPLETE"
    assert report["requested"] == {"OH": 2, "epoxide": 0, "carbonyl": 0, "COOH": 0, "COO": 0}
    assert report["achieved"] == {"OH": 1, "epoxide": 0, "carbonyl": 0, "COOH": 0, "COO": 0}
    assert report["attempts"] == 128
    assert "best" in str(caught.value).lower()


@pytest.mark.parametrize(
    "recipe",
    [
        {"nx": 1},
        {"pbc": [True, True, True]},
        {"groups": {"unobtainium": 1}},
        {"groups": {"epoxide": 100000}},
        {"groups": {"OH": -1}},
        {"vacancies": 99999},
        {"holes": [{"center": [0, 0], "radius": 99999}]},
    ],
)
def test_infeasible_or_bad_recipes_fail(recipe):
    with pytest.raises(ValueError):
        build_sheet(recipe)


def test_electron_parity_is_explicit_without_altering_requested_composition():
    # nx=ny=2 trims to naphthalene C10H8 (68 electrons). Each added OH
    # contributes nine electrons, so a single OH gives an odd electron count.
    odd = build_sheet({"nx": 2, "ny": 2, "groups": {"OH": 1}})
    even = build_sheet({"nx": 2, "ny": 2, "groups": {"OH": 2}})
    assert odd["metadata"]["electron_count"] == 77
    assert odd["metadata"]["electron_parity"] == "odd"
    assert even["metadata"]["electron_count"] == 86
    assert even["metadata"]["electron_parity"] == "even"
    assert odd["metadata"]["electronic_state_not_validated"] is True
    assert odd["metadata"]["pi_bonding_not_validated"] is True


def test_mixed_sheet_has_unambiguous_native_fragment_residues():
    from go_builder.ff import assign_current_fragments

    m = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": 17,
            "groups": {"OH": 4, "epoxide": 4, "carbonyl": 2, "COOH": 1, "COO": 1},
        }
    )
    original_ids = [a["id"] for a in m["atoms"]]
    native = [a for a in m["atoms"] if a["resname"] == "GGG"]
    assert len({a["resid"] for a in native}) == len(native)
    assigned = assign_current_fragments(m)
    assert [a["id"] for a in assigned["atoms"]] == original_ids
    assert all(a["charge"] is not None and a["atom_type"] for a in assigned["atoms"])
