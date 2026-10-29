"""Graph identity, image conventions and lossless exchange regressions."""

import copy
import json

import pytest

from go_builder.model import load_model, save_model, validate_model, write_pdb


def pair():
    return {
        "schema_version": 1,
        "name": "pair",
        "atoms": [
            {"id": i, "element": "C", "position": p, "formal_charge": 0}
            for i, p in [(4, [0, 0, 0]), (19, [1.42, 0, 0])]
        ],
        "bonds": [{"a": 4, "b": 19, "order": 1.5, "image": [0, 0, 0]}],
        "cell": [[10, 0, 0], [0, 10, 0], [0, 0, 40]],
        "pbc": [False] * 3,
        "metadata": {},
    }


def test_json_preserves_sparse_ids_and_explicit_bonds(tmp_path):
    m = pair()
    save_model(m, tmp_path / "m.json")
    assert load_model(tmp_path / "m.json") == m
    with pytest.raises(FileExistsError):
        save_model(m, tmp_path / "m.json")
    write_pdb(m, tmp_path / "m.pdb")
    assert "CONECT    4   19" in (tmp_path / "m.pdb").read_text()


@pytest.mark.parametrize("mutation", ["duplicate", "dangling", "image", "nan", "charge"])
def test_invalid_graph_is_rejected(mutation):
    m = pair()
    if mutation == "duplicate":
        m["atoms"].append(copy.deepcopy(m["atoms"][0]))
    if mutation == "dangling":
        m["bonds"][0]["b"] = 42
    if mutation == "image":
        m["bonds"][0]["image"] = [1, 0, 0]
    if mutation == "nan":
        m["atoms"][0]["position"][0] = float("nan")
    if mutation == "charge":
        m["atoms"][0]["charge"] = float("inf")
    with pytest.raises(ValueError):
        validate_model(m)


def test_validation_does_not_change_or_renumber_graph():
    m = pair()
    before = json.dumps(m)
    report = validate_model(m)
    assert report["valid"] and report["atom_count"] == 2
    assert json.dumps(m) == before


def test_structurally_undercoordinated_import_is_flagged_for_review():
    # The exchange schema admits fragments, but never silently calls two
    # singly bonded aromatic carbons a saturated chemical structure.
    assert validate_model(pair())["undercoordinated_atom_ids"] == [4, 19]


def test_fused_aromatic_carbon_is_not_mistaken_for_hypervalent_carbon():
    from go_builder.build import build_sheet

    m = build_sheet({"nx": 2, "ny": 2})
    for atom in m["atoms"]:
        atom["metadata"] = {}
    report = validate_model(m)
    assert report["valid"]
    assert report["undercoordinated_atom_ids"] == []
