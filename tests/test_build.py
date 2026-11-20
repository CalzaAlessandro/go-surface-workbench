"""Construction tests catch missing seam bonds, undershoot and invalid site reuse."""

from collections import Counter

import numpy as np
import pytest

from go_builder.build import build_sheet
from go_builder.model import adjacency, validate_model


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


def _min_14_plus_distances(model):
    from go_builder.build import _minimum_image
    from go_builder.model import adjacency

    atoms = model["atoms"]
    coords = np.asarray([a["position"] for a in atoms], dtype=float)
    ids = [a["id"] for a in atoms]
    elements = [a["element"] for a in atoms]
    ns = adjacency(model)
    bonded = {(b["a"], b["b"]) for b in model["bonds"]} | {(b["b"], b["a"]) for b in model["bonds"]}
    pair_13 = set()
    for aid in ids:
        nbrs = [nbr for nbr, _, _ in ns[aid]]
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                pair_13.add((nbrs[i], nbrs[j]))
                pair_13.add((nbrs[j], nbrs[i]))
    min_14_xh = float("inf")
    min_14_heavy = float("inf")
    for i in range(len(atoms)):
        dists = np.linalg.norm(_minimum_image(coords - coords[i], model), axis=1)
        for j in range(i + 1, len(atoms)):
            pair = (ids[i], ids[j])
            if pair in bonded or pair in pair_13:
                continue
            d = float(dists[j])
            if elements[i] == "H" or elements[j] == "H":
                min_14_xh = min(min_14_xh, d)
            else:
                min_14_heavy = min(min_14_heavy, d)
    return min_14_xh, min_14_heavy


@pytest.mark.parametrize(
    "label,recipe",
    [
        (
            "small_qm",
            {
                "nx": 3,
                "ny": 3,
                "seed": 17,
                "groups": {"OH": 2, "epoxide": 2, "carbonyl": 1, "COOH": 1, "COO": 1},
            },
        ),
        (
            "dense_oh_epx",
            {"nx": 6, "ny": 6, "seed": 17, "groups": {"OH": 25, "epoxide": 15}},
        ),
        (
            "medium_std",
            {
                "nx": 8,
                "ny": 8,
                "seed": 17,
                "groups": {"OH": 25, "epoxide": 20, "carbonyl": 4, "COOH": 6, "COO": 6},
            },
        ),
        (
            "dense_edges",
            {
                "nx": 6,
                "ny": 6,
                "seed": 17,
                "groups": {"OH": 15, "epoxide": 10, "carbonyl": 4, "COOH": 10, "COO": 8},
            },
        ),
        (
            "large_20x20",
            {
                "nx": 20,
                "ny": 20,
                "seed": 17,
                "groups": {"OH": 100, "epoxide": 80, "carbonyl": 10, "COOH": 20, "COO": 20},
            },
        ),
    ],
)
def test_healthy_by_construction_sheets_pass_steric_and_valence_gates(label, recipe):
    from go_builder.model import check_geometry_integrity

    model = build_sheet(recipe)
    report = check_geometry_integrity(model)
    assert report["passed"] is True, label
    assert report["clash_count"] == 0, label
    assert report["undercoordinated_atom_ids"] == [], label
    assert report["min_nonbonded_h_distance"] >= 0.90
    assert report["min_nonbonded_heavy_distance"] >= 1.25
    assert report["formal_charge"] == -recipe["groups"].get("COO", 0)
    assert model["metadata"]["geometry_status"] == "constructed_unrelaxed"
    min_14_xh, min_14_heavy = _min_14_plus_distances(model)
    assert min_14_xh >= 1.20
    assert min_14_heavy >= 1.60


def test_functional_group_internal_bond_lengths_and_angles():
    model = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": 17,
            "groups": {"OH": 1, "epoxide": 1, "carbonyl": 1, "COOH": 1, "COO": 1},
        }
    )
    by_id = {a["id"]: a for a in model["atoms"]}

    def pos(aid):
        return np.asarray(by_id[aid]["position"], dtype=float)

    def angle_deg(a, b, c):
        v1 = pos(a) - pos(b)
        v2 = pos(c) - pos(b)
        cos_v = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return float(np.degrees(np.arccos(np.clip(cos_v, -1.0, 1.0))))

    for g in model["metadata"]["groups"]:
        gtype = g["type"]
        atoms_g = {by_id[aid]["name"]: aid for aid in g["atom_ids"]}
        if gtype == "OH":
            assert np.linalg.norm(pos(atoms_g["OL"]) - pos(atoms_g["CL"])) == pytest.approx(1.43)
            assert np.linalg.norm(pos(atoms_g["HL"]) - pos(atoms_g["OL"])) == pytest.approx(0.97)
            assert angle_deg(atoms_g["CL"], atoms_g["OL"], atoms_g["HL"]) == pytest.approx(108.5)
        elif gtype == "epoxide":
            assert np.linalg.norm(pos(atoms_g["OE"]) - pos(atoms_g["CY"])) == pytest.approx(1.45)
            assert np.linalg.norm(pos(atoms_g["OE"]) - pos(atoms_g["CZ"])) == pytest.approx(1.45)
        elif gtype == "COO":
            assert np.linalg.norm(pos(atoms_g["OJ"]) - pos(atoms_g["C4"])) == pytest.approx(1.26)
            assert np.linalg.norm(pos(atoms_g["OK"]) - pos(atoms_g["C4"])) == pytest.approx(1.26)
            assert angle_deg(atoms_g["OJ"], atoms_g["C4"], atoms_g["OK"]) == pytest.approx(126.0)
        elif gtype == "COOH":
            assert np.linalg.norm(pos(atoms_g["OJ"]) - pos(atoms_g["C4"])) == pytest.approx(1.21)
            assert np.linalg.norm(pos(atoms_g["OK"]) - pos(atoms_g["C4"])) == pytest.approx(1.34)
            assert np.linalg.norm(pos(atoms_g["HK"]) - pos(atoms_g["OK"])) == pytest.approx(0.97)
            assert angle_deg(atoms_g["C4"], atoms_g["OK"], atoms_g["HK"]) == pytest.approx(108.0)


def test_epoxide_eph_and_ept_random_and_explicit_ids_with_normalization(tmp_path):
    from go_builder.model import (
        adjacency,
        check_geometry_integrity,
        write_mol2,
        write_pdb,
        write_xyz,
    )

    # Random Priority 2 placement + Priority 1
    m = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": 17,
            "groups": {"OH": 2, "epoxide": 2, "epoxide_eph": 3, "epoxide_ept": 2},
        }
    )
    report = check_geometry_integrity(m)
    assert report["passed"] is True
    assert report["clash_count"] == 0
    assert report["undercoordinated_atom_ids"] == []
    by_id = {a["id"]: a for a in m["atoms"]}
    ns = adjacency(m)

    eph_groups = [g for g in m["metadata"]["groups"] if g["type"] == "epoxide_eph"]
    ept_groups = [g for g in m["metadata"]["groups"] if g["type"] == "epoxide_ept"]
    assert len(eph_groups) == 3
    assert len(ept_groups) == 2

    for g in eph_groups:
        assert len(g["atom_ids"]) == 4
        names = {by_id[aid]["name"]: aid for aid in g["atom_ids"]}
        assert set(names) == {"CZ", "CK", "HK", "OE"}
        assert all(by_id[aid]["resname"] == "EPH" for aid in g["atom_ids"])
        # CK must not receive a second GRT edge hydrogen; degree == 4
        assert len(ns[names["CK"]]) == 4
        assert len(ns[names["CZ"]]) == 4
        for aid in (names["CZ"], names["CK"]):
            assert all(order == 1.0 for _, _, order in ns[aid])
        p_cz = np.asarray(by_id[names["CZ"]]["position"])
        p_ck = np.asarray(by_id[names["CK"]]["position"])
        p_oe = np.asarray(by_id[names["OE"]]["position"])
        p_hk = np.asarray(by_id[names["HK"]]["position"])
        assert np.linalg.norm(p_oe - p_cz) == pytest.approx(1.45)
        assert np.linalg.norm(p_oe - p_ck) == pytest.approx(1.45)
        assert np.linalg.norm(p_hk - p_ck) == pytest.approx(1.09)

    for g in ept_groups:
        assert len(g["atom_ids"]) == 5
        names = {by_id[aid]["name"]: aid for aid in g["atom_ids"]}
        assert set(names) == {"CP", "CK", "HP", "HK", "OE"}
        assert all(by_id[aid]["resname"] == "EPT" for aid in g["atom_ids"])
        assert len(ns[names["CP"]]) == 4
        assert len(ns[names["CK"]]) == 4
        for aid in (names["CP"], names["CK"]):
            assert all(order == 1.0 for _, _, order in ns[aid])
        p_cp = np.asarray(by_id[names["CP"]]["position"])
        p_ck = np.asarray(by_id[names["CK"]]["position"])
        p_oe = np.asarray(by_id[names["OE"]]["position"])
        p_hp = np.asarray(by_id[names["HP"]]["position"])
        p_hk = np.asarray(by_id[names["HK"]]["position"])
        assert np.linalg.norm(p_oe - p_cp) == pytest.approx(1.45)
        assert np.linalg.norm(p_oe - p_ck) == pytest.approx(1.45)
        assert np.linalg.norm(p_hp - p_cp) == pytest.approx(1.09)
        assert np.linalg.norm(p_hk - p_ck) == pytest.approx(1.09)

    write_xyz(m, tmp_path / "p2.xyz")
    write_mol2(m, tmp_path / "p2.mol2")
    write_pdb(m, tmp_path / "p2.pdb")

    # Explicit IDs with inverted [edge, basal] order normalized to [basal, edge]
    # On 6x6 sheet: 3 is edge (degree 2), 4 is basal (degree 3), bonded to 3;
    # 6 is edge (degree 2) bonded to 3? Let's pick two disjoint edge pairs dynamically.
    base = build_sheet({"nx": 6, "ny": 6})
    b_ns = adjacency(base)
    carbons = {a["id"] for a in base["atoms"] if a["element"] == "C"}
    c_deg = {cid: sum(1 for n, _, _ in b_ns[cid] if n in carbons) for cid in carbons}
    eph_b, eph_e = sorted(
        (a, b) for a in carbons for b, _, _ in b_ns[a] if b in carbons and c_deg[a] == 3 and c_deg[b] == 2
    )[0]
    ept_a, ept_b = sorted(
        (a, b)
        for a in carbons
        for b, _, _ in b_ns[a]
        if b in carbons and a < b and c_deg[a] == 2 and c_deg[b] == 2 and a != eph_e and b != eph_e
    )[0]
    explicit = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": 17,
            "groups": {
                "epoxide_eph": {"id": [[eph_e, eph_b]]},  # inverted [edge, basal]
                "epoxide_ept": {"id": [[ept_a, ept_b]]},
            },
        }
    )
    eph_placed = next(g for g in explicit["metadata"]["groups"] if g["type"] == "epoxide_eph")
    ept_placed = next(g for g in explicit["metadata"]["groups"] if g["type"] == "epoxide_ept")
    assert eph_placed["site_ids"] == [eph_b, eph_e]
    assert ept_placed["site_ids"] == [ept_a, ept_b]
    assert check_geometry_integrity(explicit)["passed"] is True


def test_invalid_eph_and_ept_sites_fail():
    # On 6x6: 58, 59 are basal; 26, 27? Let's inspect actual sites
    base = build_sheet({"nx": 6, "ny": 6})
    b_ns = adjacency(base)
    carbons = {a["id"] for a in base["atoms"] if a["element"] == "C"}
    c_deg = {cid: sum(1 for n, _, _ in b_ns[cid] if n in carbons) for cid in carbons}
    basal_pair = sorted(
        (a, b) for a in carbons for b, _, _ in b_ns[a] if b in carbons and a < b and c_deg[a] == 3 and c_deg[b] == 3
    )[0]
    edge_pair = sorted(
        (a, b) for a in carbons for b, _, _ in b_ns[a] if b in carbons and a < b and c_deg[a] == 2 and c_deg[b] == 2
    )[0]
    eph_pair = sorted(
        (a, b) for a in carbons for b, _, _ in b_ns[a] if b in carbons and c_deg[a] == 3 and c_deg[b] == 2
    )[0]

    for bad_groups in [
        {"epoxide_eph": {"id": [list(basal_pair)]}},
        {"epoxide_eph": {"id": [list(edge_pair)]}},
        {"epoxide_eph": {"id": [[58, 60]]}},  # not bonded
        {"epoxide_ept": {"id": [list(eph_pair)]}},
        {"epoxide_ept": {"id": [list(basal_pair)]}},
        {"epoxide_ept": {"id": [[edge_pair[0], edge_pair[0]]]}},
    ]:
        with pytest.raises(ValueError):
            build_sheet({"nx": 6, "ny": 6, "groups": bad_groups})


@pytest.mark.parametrize("n_islands", [1, 2])
def test_placement_modes_ordering_island_lt_random_lt_dispersed(n_islands):
    from go_builder.model import check_geometry_integrity

    common = {"nx": 10, "ny": 10, "seed": 17, "groups": {"OH": 15, "epoxide": 10}}
    m_island = build_sheet({**common, "placement": {"mode": "island", "n_islands": n_islands}})
    m_random = build_sheet({**common, "placement": "random"})
    m_dispersed = build_sheet({**common, "placement": "dispersed"})

    for m in (m_island, m_random, m_dispersed):
        rep = check_geometry_integrity(m)
        assert rep["passed"] is True
        assert rep["clash_count"] == 0

    d_island = m_island["metadata"]["placement"]["mean_basal_group_nn_distance"]
    d_random = m_random["metadata"]["placement"]["mean_basal_group_nn_distance"]
    d_dispersed = m_dispersed["metadata"]["placement"]["mean_basal_group_nn_distance"]
    assert d_island < d_random < d_dispersed


def test_island_circular_disk_centers_radius_and_error_on_small_radius():
    from go_builder.build import _minimum_image
    from go_builder.errors import PlacementError
    from go_builder.model import check_geometry_integrity

    base = build_sheet({"nx": 10, "ny": 10})
    xs = [a["position"][0] for a in base["atoms"]]
    ys = [a["position"][1] for a in base["atoms"]]
    center = [float(np.mean(xs)), float(np.mean(ys))]

    disk = build_sheet(
        {
            "nx": 10,
            "ny": 10,
            "seed": 17,
            "groups": {"OH": 8, "epoxide": 4, "COOH": 2},
            "placement": {
                "mode": "island",
                "centers": [center],
                "radius": 6.0,
                "edge_mode": "near_islands",
            },
        }
    )
    assert check_geometry_integrity(disk)["passed"] is True
    by_id = {a["id"]: np.asarray(a["position"][:2]) for a in base["atoms"]}
    for g in disk["metadata"]["groups"]:
        if g["type"] in ("OH", "epoxide"):
            for sid in g["site_ids"]:
                d = np.linalg.norm(
                    _minimum_image(
                        np.array([by_id[sid][0] - center[0], by_id[sid][1] - center[1], 0.0]),
                        disk,
                    )[:2]
                )
                assert d <= 6.0 + 1e-9

    with pytest.raises(PlacementError, match="radius"):
        build_sheet(
            {
                "nx": 10,
                "ny": 10,
                "seed": 17,
                "groups": {"OH": 10},
                "placement": {"mode": "island", "centers": [center], "radius": 1.0},
            }
        )


def test_island_seed_atoms_and_invalid_placement_specs():
    base = build_sheet({"nx": 8, "ny": 8})
    m = build_sheet(
        {
            "nx": 8,
            "ny": 8,
            "seed": 5,
            "groups": {"OH": 6},
            "placement": {"mode": "island", "seed_atoms": [100], "compactness": 1.0},
        }
    )
    assert m["metadata"]["placement"]["seed_atoms"] == [100]

    for bad_placement in [
        "grid",
        123,
        {"mode": "random", "radius": 5.0},
        {"mode": "island", "centers": [[1.0, 2.0]], "seed_atoms": [100]},
        {"mode": "island", "centers": [[1.0, 2.0]], "n_islands": 2},
        {"mode": "island", "radius": -1.0},
        {"mode": "island", "compactness": 1.5},
        {"mode": "island", "edge_mode": "invalid"},
        {"mode": "island", "seed_atoms": [999999]},
    ]:
        with pytest.raises(ValueError):
            build_sheet({"nx": 6, "ny": 6, "groups": {"OH": 2}, "placement": bad_placement})


def test_group_percent_and_fraction_resolution_and_validation():
    # Post-hole resolution on 8x8 with a hole
    base_holed = build_sheet(
        {"nx": 8, "ny": 8, "seed": 17, "holes": [{"size": 4, "seed_atom": 145}]}
    )
    ns = adjacency(base_holed)
    carbons = {a["id"] for a in base_holed["atoms"] if a["element"] == "C"}
    c_deg = {c: sum(1 for n, _, _ in ns[c] if n in carbons) for c in carbons}
    n_c = len(carbons)
    n_basal = sum(1 for c in carbons if c_deg[c] == 3)
    n_edge = sum(1 for c in carbons if c_deg[c] == 2)

    m = build_sheet(
        {
            "nx": 8,
            "ny": 8,
            "seed": 17,
            "holes": [{"size": 4, "seed_atom": 145}],
            "groups": {
                "OH": {"percent": 10.0},
                "epoxide": {"fraction": 0.1, "basis": "eligible"},
                "COOH": {"percent": 10.0, "basis": "eligible"},
                "carbonyl": {"fraction": 0.02, "basis": "carbons"},
            },
        }
    )
    resolved = m["metadata"]["resolved_group_counts"]
    assert resolved["OH"] == int(round(0.10 * n_basal))
    assert resolved["epoxide"] == int(round(0.10 * (n_basal // 2)))
    assert resolved["COOH"] == int(round(0.10 * n_edge))
    assert resolved["carbonyl"] == int(round(0.02 * n_c))
    assert m["metadata"]["events"][-1]["group_counts"] == resolved
    assert Counter(g["type"] for g in m["metadata"]["groups"]) == {
        k: v for k, v in resolved.items() if v > 0
    }

    for bad in [
        {"OH": {"count": 1, "percent": 10}},
        {"OH": {"percent": 10, "fraction": 0.1}},
        {"OH": {"id": [58], "percent": 10}},
        {"OH": {"count": 2, "basis": "carbons"}},
        {"OH": {"percent": -1}},
        {"OH": {"percent": 101}},
        {"OH": {"fraction": 1.2}},
        {"OH": {"fraction": 0.1, "basis": "oxygen"}},
    ]:
        with pytest.raises(ValueError):
            build_sheet({"nx": 6, "ny": 6, "groups": bad})


@pytest.mark.parametrize("distribution", ["small", "medium", "large", "mixed", "random"])
@pytest.mark.parametrize("as_list", [False, True])
def test_hiergo_hole_distributions_with_percent_vacancies_keep_sheet_connected_and_healthy(
    distribution, as_list
):
    from go_builder.model import check_geometry_integrity

    profile = {"distribution": distribution, "percent_vacancies": 10.0}
    holes_arg = [profile] if as_list else profile
    model = build_sheet({"nx": 12, "ny": 12, "seed": 17, "holes": holes_arg})
    assert validate_model(model)["connected_components"] == 1
    geom = check_geometry_integrity(model)
    assert geom["passed"] is True
    assert geom["clash_count"] == 0
    assert geom["undercoordinated_atom_ids"] == []
    summary = model["metadata"]["hole_summary"]
    assert summary["distribution"] == distribution
    assert summary["requested_percent_vacancies"] == 10.0
    assert summary["requested_hole_count"] is None
    assert summary["holes_created"] == len(summary["hole_sizes"]) >= 1
    assert summary["initial_carbons"] == 550
    assert summary["removed_carbons"] == sum(summary["hole_sizes"])
    assert summary["actual_vacancy_percent"] == pytest.approx(
        100.0 * summary["removed_carbons"] / 550
    )
    assert 8.0 <= summary["actual_vacancy_percent"] <= 12.5


def test_hiergo_hole_mean_sizes_order_small_lt_medium_lt_large_and_count_mode():
    from go_builder.model import check_geometry_integrity

    means = {}
    for dist in ("small", "medium", "large"):
        m = build_sheet(
            {
                "nx": 14,
                "ny": 14,
                "seed": 23,
                "holes": {"distribution": dist, "percent_vacancies": 10.0},
            }
        )
        assert validate_model(m)["connected_components"] == 1
        assert check_geometry_integrity(m)["passed"] is True
        sizes = m["metadata"]["hole_summary"]["hole_sizes"]
        means[dist] = sum(sizes) / len(sizes)

    assert means["small"] < means["medium"] < means["large"]

    counted = build_sheet(
        {
            "nx": 12,
            "ny": 12,
            "seed": 11,
            "holes": {"distribution": "medium", "count": 4, "min_size": 6, "max_size": 8},
        }
    )
    summary = counted["metadata"]["hole_summary"]
    assert summary["requested_hole_count"] == 4
    assert summary["holes_created"] == 4
    assert all(6 <= sz <= 10 for sz in summary["hole_sizes"])
    assert check_geometry_integrity(counted)["passed"] is True

    for bad_holes in [
        {},
        {"distribution": "huge", "count": 1},
        {"distribution": "small"},
        {"distribution": "small", "percent_vacancies": 0.0},
        {"distribution": "small", "percent_vacancies": 55.0},
        {"distribution": "small", "count": 0},
        {"distribution": "small", "count": 2, "min_size": 8, "max_size": 4},
        {"distribution": "small", "count": 1, "avoid_outer_edge": "yes"},
    ]:
        with pytest.raises(ValueError):
            build_sheet({"nx": 6, "ny": 6, "holes": bad_holes})


def test_build_sheet_direct_o_over_c_c_over_o_modes_and_priority():
    from go_builder.model import check_geometry_integrity

    # Top-level O/C and C/O shortcuts
    m_oc = build_sheet({"nx": 8, "ny": 8, "seed": 17, "o_over_c": 0.35})
    m_co = build_sheet({"nx": 8, "ny": 8, "seed": 17, "c_over_o": 2.5})
    for m in (m_oc, m_co):
        geom = check_geometry_integrity(m)
        assert geom["passed"] is True
        assert geom["clash_count"] == 0
        assert "composition" in m["metadata"]
    assert abs(m_oc["metadata"]["composition"]["actual_O_over_C"] - 0.35) < 0.02
    assert abs(m_co["metadata"]["composition"]["actual_C_over_O"] - 2.5) < 0.1

    # mode="balanced" vs mode="random" across different seeds
    m_bal = build_sheet(
        {"nx": 8, "ny": 8, "seed": 17, "composition": {"o_over_c": 0.35, "mode": "balanced"}}
    )
    mixes = []
    for s in (1, 2, 3, 4):
        mr = build_sheet(
            {"nx": 8, "ny": 8, "seed": s, "composition": {"o_over_c": 0.35, "mode": "random"}}
        )
        assert check_geometry_integrity(mr)["passed"] is True
        assert abs(mr["metadata"]["composition"]["actual_O_over_C"] - 0.35) < 0.02
        counts = Counter(g["type"] for g in mr["metadata"]["groups"])
        mixes.append((counts["OH"], counts["epoxide"], counts["COOH"], counts["COO"]))
    assert len(set(mixes)) >= 3
    assert m_bal["metadata"]["composition"]["composition_mode"] == "balanced"

    # Priority ordering flips OH vs epoxide dominance at identical O/C
    m_oh_first = build_sheet(
        {
            "nx": 8,
            "ny": 8,
            "seed": 17,
            "composition": {"o_over_c": 0.35, "priority": ["OH", "epoxide"]},
        }
    )
    m_epx_first = build_sheet(
        {
            "nx": 8,
            "ny": 8,
            "seed": 17,
            "composition": {"o_over_c": 0.35, "priority": ["epoxide", "OH"]},
        }
    )
    g_oh = Counter(g["type"] for g in m_oh_first["metadata"]["groups"])
    g_epx = Counter(g["type"] for g in m_epx_first["metadata"]["groups"])
    assert g_oh["OH"] > g_epx["OH"]
    assert g_epx["epoxide"] > g_oh["epoxide"]


def test_build_sheet_full_hiergo_combo_mixed_holes_o_over_c_priority_ph_island():
    from go_builder.model import check_geometry_integrity

    model = build_sheet(
        {
            "nx": 10,
            "ny": 10,
            "seed": 17,
            "holes": {"distribution": "mixed", "percent_vacancies": 10.0},
            "composition": {
                "o_over_c": 0.35,
                "priority": ["OH", "epoxide", "COOH", "COO"],
                "ph": 7.0,
            },
            "placement": "island",
        }
    )
    geom = check_geometry_integrity(model)
    assert geom["passed"] is True
    assert geom["clash_count"] == 0
    assert geom["undercoordinated_atom_ids"] == []
    assert validate_model(model)["connected_components"] == 1
    assert model["metadata"]["geometry_status"] == "constructed_unrelaxed"
    comp = model["metadata"]["composition"]
    assert abs(comp["actual_O_over_C"] - 0.35) < 0.02
    counts = Counter(g["type"] for g in model["metadata"]["groups"])
    assert counts["OH"] > counts["epoxide"]
    assert counts["COO"] > counts["COOH"]
    assert model["metadata"]["placement"]["mode"] == "island"
    assert model["metadata"]["hole_summary"]["distribution"] == "mixed"

    for bad_recipe in [
        {"nx": 6, "ny": 6, "groups": {"OH": 2}, "o_over_c": 0.2},
        {"nx": 6, "ny": 6, "c_over_o": 4.0, "o_over_c": 0.25},
        {"nx": 6, "ny": 6, "c_over_o": 4.0, "composition": {"c_over_o": 4.0}},
    ]:
        with pytest.raises(ValueError):
            build_sheet(bad_recipe)


def _assert_exact_edge_bond_lengths(model):
    from go_builder.model import check_geometry_integrity

    rep = check_geometry_integrity(model)
    assert rep["passed"] is True
    assert rep["clash_count"] == 0
    coords = {a["id"]: np.array(a["position"], dtype=float) for a in model["atoms"]}
    by_id = {a["id"]: a for a in model["atoms"]}
    cell = np.array(model["cell"], dtype=float)
    for b in model["bonds"]:
        a, bb = by_id[b["a"]], by_id[b["b"]]
        d = float(
            np.linalg.norm(coords[bb["id"]] + np.array(b["image"]) @ cell - coords[a["id"]])
        )
        els = tuple(sorted((a["element"], bb["element"])))
        if els == ("C", "H"):
            assert abs(d - 1.090) <= 1e-3
        elif els == ("H", "O"):
            assert abs(d - 0.970) <= 1e-3
        elif els == ("C", "O") and "CBD" in (a["resname"], bb["resname"]):
            assert abs(d - 1.260) <= 1e-3


def test_orthonormal_edge_frame_fixes_puckered_edge_bond_lengths():
    for seed in (1, 31):
        m = build_sheet(
            {
                "nx": 5,
                "ny": 5,
                "seed": seed,
                "holes": {"distribution": "mixed", "percent_vacancies": 12},
                "groups": {"OH": 4, "epoxide_eph": 2, "epoxide_ept": 1, "COO": 2},
            }
        )
        _assert_exact_edge_bond_lengths(m)

    m0 = build_sheet(
        {
            "nx": 5,
            "ny": 5,
            "seed": 0,
            "groups": {"epoxide_eph": 7, "epoxide_ept": 2, "COOH": 4, "COO": 3},
        }
    )
    _assert_exact_edge_bond_lengths(m0)

    m_isl = build_sheet(
        {
            "nx": 8,
            "ny": 8,
            "seed": 17,
            "placement": {
                "mode": "island",
                "n_islands": 2,
                "radius": 8.0,
                "edge_mode": "near_islands",
            },
            "groups": {
                "OH": 6,
                "epoxide": 4,
                "carbonyl": 2,
                "COOH": 2,
                "COO": 2,
                "epoxide_eph": 2,
                "epoxide_ept": 1,
            },
        }
    )
    _assert_exact_edge_bond_lengths(m_isl)


def test_large_sheet_48x48_scalability_and_zero_recursion_error():
    import time
    from go_builder.model import check_geometry_integrity

    t0 = time.perf_counter()
    model = build_sheet({"nx": 48, "ny": 48, "pbc": [True, True, False]})
    elapsed_build = time.perf_counter() - t0
    assert len(model["atoms"]) == 9216
    t1 = time.perf_counter()
    rep = check_geometry_integrity(model)
    elapsed_check = time.perf_counter() - t1
    assert rep["passed"] is True
    assert rep["clash_count"] == 0
    assert elapsed_build + elapsed_check < 3.0


def test_periodic_holes_island_adaptive_composition_retry():
    from go_builder.model import check_geometry_integrity

    model = build_sheet(
        {
            "nx": 10,
            "ny": 10,
            "seed": 7,
            "pbc": [True, True, False],
            "holes": {"distribution": "small", "count": 3, "min_size": 2, "max_size": 4},
            "composition": {
                "c_over_o": 3.0,
                "enabled_groups": ["OH", "epoxide", "carbonyl", "carboxyl"],
                "ph": 7.0,
                "mode": "balanced",
            },
            "placement": {"mode": "island", "n_islands": 1, "edge_mode": "near_islands"},
        }
    )
    rep = check_geometry_integrity(model)
    assert rep["passed"] is True
    assert rep["clash_count"] == 0
    comp = model["metadata"]["composition"]
    assert comp["composition_attempts"] >= 1
    assert abs(comp["actual_C_over_O"] - 3.0) < 0.15


def test_hole_constraints_avoid_outer_edge_and_pore_summary():
    # Contradictory constraints must fail explicitly
    for bad in [
        {"distribution": "small", "percent_vacancies": 1.0, "min_size": 20},
        {"distribution": "small", "percent_vacancies": 5.0, "count": 4, "min_size": 10},
        {"distribution": "large", "percent_vacancies": 30.0, "count": 2, "max_size": 5},
        {"distribution": "small", "count": 20, "min_size": 10},
    ]:
        with pytest.raises(ValueError, match="Contradictory"):
            build_sheet({"nx": 8, "ny": 8, "seed": 3, "holes": bad})

    # avoid_outer_edge=True keeps outer edge untouched (outer_edge_notch_count == 0)
    # Small 4-6 C holes are narrower than 2*(R_C + R_water)=6.2 A -> 0 accessible area;
    # large 14-20 C pores admit the 1.4 A water probe -> positive accessible area.
    m_small = build_sheet(
        {
            "nx": 12,
            "ny": 12,
            "seed": 17,
            "holes": {
                "distribution": "medium",
                "count": 3,
                "min_size": 4,
                "max_size": 6,
                "avoid_outer_edge": True,
                "probe_radius_angstrom": 1.4,
            },
        }
    )
    hs_small = m_small["metadata"]["hole_summary"]
    assert hs_small["outer_edge_notch_count"] == 0
    assert hs_small["internal_pore_count"] == hs_small["final_pore_components"]
    assert hs_small["pore_accessible_area_angstrom2"] == 0.0

    m = build_sheet(
        {
            "nx": 14,
            "ny": 14,
            "seed": 17,
            "holes": {
                "distribution": "large",
                "count": 2,
                "min_size": 14,
                "max_size": 20,
                "avoid_outer_edge": True,
                "probe_radius_angstrom": 1.4,
            },
        }
    )
    hs = m["metadata"]["hole_summary"]
    assert hs["excavation_events"] == 2
    assert len(hs["excavation_sizes"]) == 2
    assert hs["pruned_tip_carbons"] >= 0
    assert hs["final_pore_components"] >= 1
    assert len(hs["final_pore_sizes"]) == hs["final_pore_components"]
    assert hs["outer_edge_notch_count"] == 0
    assert hs["internal_pore_count"] == hs["final_pore_components"]
    assert hs["probe_radius_angstrom"] == 1.4
    assert hs["pore_accessible_area_angstrom2"] > 0.0
    assert hs["pore_accessible_area_nm2"] == pytest.approx(
        hs["pore_accessible_area_angstrom2"] / 100.0, abs=1e-5
    )


def test_nanoflake_shapes_and_nm_physical_dimensions():
    from go_builder.model import check_geometry_integrity

    hex_zz = build_sheet(
        {"diameter_nm": 2.0, "shape": "hexagonal", "edge_type": "zigzag", "seed": 1}
    )
    hex_ac = build_sheet(
        {"diameter_nm": 2.0, "shape": "hexagonal", "edge_type": "armchair", "seed": 1}
    )
    circ = build_sheet(
        {
            "radius_nm": 1.0,
            "shape": "circular",
            "seed": 1,
            "groups": {"OH": 4, "epoxide": 2, "COO": 2},
        }
    )
    rect_pbc = build_sheet(
        {"lx_nm": 2.0, "ly_nm": 2.0, "pbc": [True, True, False], "shape": "rectangular"}
    )
    rect_fin = build_sheet({"lx_nm": 2.0, "ly_nm": 2.0, "shape": "rectangular"})

    for m in (hex_zz, hex_ac, circ, rect_pbc, rect_fin):
        assert validate_model(m)["connected_components"] == 1
        assert check_geometry_integrity(m)["passed"] is True
        assert "dimensions" in m["metadata"]
        assert "morphology" in m["metadata"]

    mz = hex_zz["metadata"]["morphology"]
    ma = hex_ac["metadata"]["morphology"]
    mc = circ["metadata"]["morphology"]
    assert mz["shape"] == "hexagonal" and mz["edge_type"] == "zigzag"
    assert ma["shape"] == "hexagonal" and ma["edge_type"] == "armchair"
    assert mz["zigzag_edge_carbons"] > ma["zigzag_edge_carbons"]
    assert ma["armchair_edge_pairs"] > mz["armchair_edge_pairs"]
    assert mc["shape"] == "circular" and mc["area_nm2"] > 0.0
    assert mc["edge_to_area_ratio_nm_inv"] > 0.0

    assert rect_pbc["metadata"]["dimensions"]["rounding_policy"] == "ceil_cell_multiple_min2"
    assert rect_pbc["metadata"]["dimensions"]["cell_nm"][0] >= 2.0
    assert rect_pbc["metadata"]["dimensions"]["cell_nm"][1] >= 2.0
    assert rect_fin["metadata"]["dimensions"]["rounding_policy"] == "round_cell_multiple_min2"

    for bad in [
        {"nx": 6, "ny": 6, "shape": "hexagonal", "pbc": [True, True, False]},
        {"nx": 6, "ny": 6, "shape": "circular", "pbc": [True, True, False]},
        {"nx": 6, "ny": 6, "lx_nm": 2.0, "ly_nm": 2.0},
        {"diameter_nm": 2.0, "radius_nm": 1.0, "shape": "circular"},
        {"diameter_nm": 2.0, "shape": "rectangular"},
        {"nx": 6, "ny": 6, "shape": "circular", "edge_type": "zigzag"},
        {"lx_nm": 2.0, "shape": "rectangular"},
    ]:
        with pytest.raises(ValueError):
            build_sheet(bad)
