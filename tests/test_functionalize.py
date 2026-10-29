"""Anti ring opening conserves atoms/protons, topology and seam geometry."""

import copy
from collections import Counter

import numpy as np
import pytest

from go_builder.build import build_sheet
from go_builder.functionalize import open_epoxides
from go_builder.model import validate_model


def ammonia():
    return {
        "atoms": [
            {"id": 1, "element": "N", "position": [0, 0, 0], "formal_charge": 0},
            {"id": 2, "element": "H", "position": [0, 0, 1.01], "formal_charge": 0},
            {"id": 3, "element": "H", "position": [0.95, 0, -0.34], "formal_charge": 0},
            {"id": 4, "element": "H", "position": [-0.47, 0.82, -0.34], "formal_charge": 0},
        ],
        "bonds": [{"a": 1, "b": h, "order": 1, "image": [0, 0, 0]} for h in [2, 3, 4]],
        "attachment_atom": 1,
        "transfer_hydrogen": 2,
    }


def test_grafted_fragment_preserves_provenance_and_unique_residue_atom_names():
    m = build_sheet({"nx": 5, "ny": 5, "groups": {"epoxide": 1}})
    fragment = ammonia()
    fragment["metadata"] = {"source": "explicit_test_fragment", "site_choice": "N1"}
    fragment["atoms"][0]["metadata"] = {"source_atom_id": 1}
    fragment["bonds"][1]["metadata"] = {"source_stereo": "retained"}
    out = open_epoxides(m, {"count": 1, "fragment": fragment})
    added = [
        a
        for a in out["atoms"]
        if a["metadata"].get("role") in ("nucleophile", "transferred_hydrogen")
    ]
    assert len({a["name"] for a in added}) == len(added)
    event = out["metadata"]["events"][-1]
    assert event["fragment_preparation"]["source"] == "explicit_test_fragment"
    attachment = next(a for a in added if a["id"] == event["attachment_id"])
    assert attachment["metadata"]["fragment_atom_metadata"]["source_atom_id"] == 1
    assert any(
        b.get("metadata", {}).get("fragment_bond_metadata", {}).get("source_stereo") == "retained"
        for b in out["bonds"]
    )


def test_opening_is_anti_and_keeps_atoms_charge_and_proton_balance():
    m = build_sheet({"nx": 6, "ny": 6, "groups": {"epoxide": 1}, "faces": "up"})
    before = copy.deepcopy(m)
    out = open_epoxides(m, {"count": 1, "fragment": ammonia()})
    assert m == before
    assert Counter(a["element"] for a in out["atoms"]) == Counter(
        a["element"] for a in m["atoms"]
    ) + Counter(N=1, H=3)
    assert sum(a["formal_charge"] for a in out["atoms"]) == 0
    event = out["metadata"]["events"][-1]
    atoms = {a["id"]: a for a in out["atoms"]}
    assert atoms[event["attachment_id"]]["position"][2] < 0
    assert atoms[event["oxygen_id"]]["position"][2] > 0
    assert atoms[event["transferred_hydrogen_id"]]["position"][2] > 0
    pairs = {frozenset((b["a"], b["b"])) for b in out["bonds"]}
    assert frozenset((event["oxygen_id"], event["attacked_carbon_id"])) not in pairs
    assert frozenset((event["oxygen_id"], event["retained_carbon_id"])) in pairs
    assert frozenset((event["attachment_id"], event["attacked_carbon_id"])) in pairs
    assert validate_model(out)["valid"]


def test_opening_requires_explicit_supported_nucleophile_and_available_epoxide():
    m = build_sheet({"nx": 4, "ny": 4, "groups": {"epoxide": 1}})
    with pytest.raises(ValueError):
        open_epoxides(m, {"count": 1, "nucleophile": "PEG"})
    with pytest.raises(ValueError):
        open_epoxides(m, {"count": 2, "fragment": ammonia()})
    bad = ammonia()
    bad["transfer_hydrogen"] = 1
    with pytest.raises(ValueError):
        open_epoxides(m, {"count": 1, "fragment": bad})


def test_periodic_epoxide_opening_preserves_local_bond_lengths():
    m = build_sheet(
        {"nx": 3, "ny": 3, "pbc": [True, True, False], "groups": {"epoxide": 1}, "seed": 1}
    )
    out = open_epoxides(m, {"count": 1, "fragment": ammonia()})
    event = out["metadata"]["events"][-1]
    coords = {a["id"]: np.array(a["position"]) for a in out["atoms"]}
    for b in out["bonds"]:
        if event["attachment_id"] in [b["a"], b["b"]] or event["oxygen_id"] in [b["a"], b["b"]]:
            d = coords[b["b"]] + np.array(b["image"]) @ np.array(out["cell"]) - coords[b["a"]]
            assert 0.85 < np.linalg.norm(d) < 1.65


def test_actual_seam_epoxide_retains_images_after_opening():
    m = build_sheet(
        {"nx": 3, "ny": 3, "pbc": [True, True, False], "groups": {"epoxide": 1}, "seed": 4}
    )
    group = m["metadata"]["groups"][0]
    assert any(
        set((b["a"], b["b"])) == set(group["site_ids"]) and any(b["image"]) for b in m["bonds"]
    )
    out = open_epoxides(m, {"count": 1, "fragment": ammonia(), "seed": 4})
    event = out["metadata"]["events"][-1]
    old_atoms = {a["id"]: a for a in m["atoms"]}
    assert all(
        a["position"] == old_atoms[a["id"]]["position"]
        for a in out["atoms"]
        if a["element"] == "C" and a["id"] in old_atoms
    )
    coords = {a["id"]: np.array(a["position"]) for a in out["atoms"]}
    bond = next(
        b
        for b in out["bonds"]
        if set((b["a"], b["b"])) == {event["oxygen_id"], event["retained_carbon_id"]}
    )
    displacement = (
        coords[bond["b"]] + np.array(bond["image"]) @ np.array(out["cell"]) - coords[bond["a"]]
    )
    assert np.linalg.norm(displacement) == pytest.approx(1.43)
    assert out["metadata"]["requires_parameterization"] is True
    assert out["metadata"]["parameterization_status"] == "new_covalent_junction_unparameterized"


def test_incomplete_fragment_hydrogen_and_charged_donor_are_rejected():
    m = build_sheet({"nx": 4, "ny": 4, "groups": {"epoxide": 1}})
    incomplete = ammonia()
    incomplete["atoms"] = incomplete["atoms"][:-1]
    incomplete["bonds"] = incomplete["bonds"][:-1]
    with pytest.raises(ValueError, match="neutral amines"):
        open_epoxides(m, {"fragment": incomplete})
    charged = ammonia()
    charged["atoms"][0]["formal_charge"] = 1
    with pytest.raises(ValueError, match="neutral amines"):
        open_epoxides(m, {"fragment": charged})


def test_opening_updates_electron_count_for_added_neutral_ammonia():
    m = build_sheet({"nx": 4, "ny": 4, "groups": {"epoxide": 1}})
    out = open_epoxides(m, {"fragment": ammonia()})
    assert out["metadata"]["electron_count"] == m["metadata"]["electron_count"] + 10
    assert out["metadata"]["electron_parity"] == m["metadata"]["electron_parity"]
    assert out["metadata"]["electronic_state_not_validated"] is True


def test_anti_opening_uses_local_sheet_normal_after_rigid_rotation():
    m = build_sheet({"nx": 4, "ny": 4, "groups": {"epoxide": 1}, "faces": "up"})
    # Rotate the entire explicit graph about x: the sheet now lies in xz.
    rotation = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    for atom in m["atoms"]:
        atom["position"] = (rotation @ atom["position"]).tolist()
    m["cell"] = (np.array(m["cell"]) @ rotation.T).tolist()
    out = open_epoxides(m, {"fragment": ammonia()})
    event = out["metadata"]["events"][-1]
    coords = {a["id"]: np.array(a["position"]) for a in out["atoms"]}
    oxygen_direction = coords[event["oxygen_id"]] - coords[event["retained_carbon_id"]]
    attack_direction = coords[event["attachment_id"]] - coords[event["attacked_carbon_id"]]
    original_oxygen_normal = np.array([0, -1, 0])
    assert np.dot(oxygen_direction, original_oxygen_normal) > 1
    assert np.dot(attack_direction, original_oxygen_normal) < -1
    assert np.dot(oxygen_direction, attack_direction) < 0


def test_real_peg_amine_count_opening_falls_back_to_another_available_site():
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    m = build_sheet(json.loads((root / "examples/finite_GO.json").read_text()))
    recipe = json.loads((root / "examples/open_epoxide_PEG2_amine.json").read_text())
    out = open_epoxides(m, recipe)
    assert sum(g["type"] == "opened_epoxide" for g in out["metadata"]["groups"]) == 1
    assert len(out["atoms"]) == len(m["atoms"]) + 18
    assert validate_model(out)["undercoordinated_atom_ids"] == []


def test_incoming_bond_points_backside_along_the_original_broken_c_o_axis():
    m = build_sheet({"nx": 4, "ny": 4, "groups": {"epoxide": 1}, "faces": "up"})
    out = open_epoxides(m, {"fragment": ammonia()})
    event = out["metadata"]["events"][-1]
    before = {a["id"]: np.array(a["position"]) for a in m["atoms"]}
    after = {a["id"]: np.array(a["position"]) for a in out["atoms"]}
    original = before[event["oxygen_id"]] - before[event["attacked_carbon_id"]]
    incoming = after[event["attachment_id"]] - after[event["attacked_carbon_id"]]
    assert np.dot(original, incoming) / (
        np.linalg.norm(original) * np.linalg.norm(incoming)
    ) == pytest.approx(-1)


def test_random_count_is_reproducible_and_derives_epoxide_fraction():
    model = build_sheet({"nx": 6, "ny": 6, "groups": {"epoxide": 3}, "seed": 17})
    recipe = {"count": 2, "fragment": ammonia(), "seed": 7}
    first = open_epoxides(model, recipe)
    assert first == open_epoxides(model, recipe)
    selected = first["metadata"]["functionalization_history"][-1]
    assert selected["mode"] == "random_count"
    assert selected["available_epoxides"] == 3
    assert selected["requested_count"] == selected["placed_count"] == 2
    assert "requested_fraction" not in selected
    assert selected["actual_fraction"] == pytest.approx(2 / 3)
    assert selected["seed"] == 7
    assert len(set(selected["selected_epoxide_ids"])) == 2
    assert len(first["atoms"]) == len(model["atoms"]) + 8
    assert sum(g["type"] == "epoxide" for g in first["metadata"]["groups"]) == 1


def test_random_seed_changes_sites_and_explicit_selection_stays_exact():
    model = build_sheet({"nx": 6, "ny": 6, "groups": {"epoxide": 6}, "seed": 17})
    chosen = []
    for seed in (1, 2):
        out = open_epoxides(model, {"count": 2, "fragment": ammonia(), "seed": seed})
        chosen.append(out["metadata"]["functionalization_history"][-1]["selected_epoxide_ids"])
    assert set(chosen[0]) != set(chosen[1])
    explicit = open_epoxides(model, {"id": chosen[0], "fragment": ammonia(), "seed": 3})
    history = explicit["metadata"]["functionalization_history"][-1]
    assert history["mode"] == "explicit_sites"
    assert history["selected_epoxide_ids"] == chosen[0]


@pytest.mark.parametrize("fraction", [0.5, 1.1, True, "0.5", float("nan"), float("inf")])
def test_fraction_is_output_only_and_cannot_be_supplied_as_input(fraction):
    model = build_sheet({"nx": 3, "ny": 3, "groups": {"epoxide": 1}})
    before = copy.deepcopy(model)
    with pytest.raises(ValueError, match="fraction"):
        open_epoxides(model, {"fraction_epoxides": fraction, "fragment": ammonia()})
    assert model == before


def test_explicit_id_count_must_match_the_requested_count():
    model = build_sheet({"nx": 3, "ny": 3, "groups": {"epoxide": 1}})
    with pytest.raises(ValueError, match="count must match the number of id entries"):
        open_epoxides(model, {"count": 2, "id": [1], "fragment": ammonia()})


@pytest.mark.parametrize("sites", [[[1]], [True], [1, 1], [1.0], [None], "1"])
def test_malformed_explicit_site_ids_raise_an_explanatory_error(sites):
    model = build_sheet({"nx": 3, "ny": 3, "groups": {"epoxide": 1}})
    with pytest.raises(ValueError, match="id"):
        open_epoxides(model, {"id": sites, "fragment": ammonia()})


@pytest.mark.parametrize("count", [0, 2])
def test_fraction_is_derived_for_zero_and_complete_loading(count):
    model = build_sheet({"nx": 5, "ny": 5, "groups": {"epoxide": 2}, "seed": 17})
    out = open_epoxides(model, {"count": count, "fragment": ammonia(), "seed": 1})
    assert len(out["atoms"]) == len(model["atoms"]) + 4 * count
    history = out["metadata"]["functionalization_history"][-1]
    assert history["placed_count"] == count
    assert history["actual_fraction"] == count / 2
