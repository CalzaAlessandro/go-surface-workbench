from copy import deepcopy
from pathlib import Path

import pytest

from go_builder.ff import assign_current_fragments, load_fragment_library, normalize_libraries

REFERENCE = Path(__file__).resolve().parents[1] / "resources/reference/current"


def fragment_model(name="HYD"):
    lib = load_fragment_library(REFERENCE / f"{name}_gaff2.lib", name)
    return {
        "schema_version": 1,
        "name": name,
        "pbc": [False] * 3,
        "cell": None,
        "atoms": [
            {
                "id": 10 + i,
                "name": a["name"],
                "element": {1: "H", 6: "C", 8: "O"}[a["atomic_number"]],
                "resname": name,
                "resid": 7,
                "position": [float(i), 0.0, 0.0],
                "charge": None,
                "atom_type": None,
                "formal_charge": 0,
                "metadata": {},
            }
            for i, a in enumerate(lib["atoms"])
        ],
        "bonds": [
            {"a": 9 + a, "b": 9 + b, "order": 1, "image": [0, 0, 0]} for a, b, _ in lib["bonds"]
        ],
        "metadata": {},
    }


def test_scoped_parser_ignores_extraneous_unit_and_audits_stale_perturbation():
    lib = load_fragment_library(REFERENCE / "EPH_gaff2.lib", "EPH")
    assert len(lib["atoms"]) == 4
    assert lib["atoms"][1]["atom_type"] == "cx"
    assert lib["atoms"][2]["name"] == "HK"
    assert any("unindexed" in x for x in lib["warnings"])
    assert any("perturbation" in x for x in lib["warnings"])


@pytest.mark.parametrize("name", ["GGG", "GRT", "HYD", "EPX", "EPH", "EPT", "CBX", "CBD", "CAB"])
def test_exact_current_fragment_charges_include_anchors(name):
    model = fragment_model(name)
    original = deepcopy(model)
    assigned = assign_current_fragments(model, REFERENCE)
    lib = load_fragment_library(REFERENCE / f"{name}_gaff2.lib", name)
    assert [a["charge"] for a in assigned["atoms"]] == [a["charge"] for a in lib["atoms"]]
    assert [a["atom_type"] for a in assigned["atoms"]] == [a["atom_type"] for a in lib["atoms"]]
    assert model == original
    if name == "CBX":
        assert assigned["atoms"][0]["charge"] == -0.115
    if name == "CBD":
        assert sum(a["charge"] for a in assigned["atoms"]) == pytest.approx(-1)


def test_missing_group_atom_and_internal_bond_are_rejected():
    model = fragment_model()
    model["atoms"].pop()
    model["bonds"].pop()
    with pytest.raises(ValueError, match="fragment|template"):
        assign_current_fragments(model, REFERENCE)
    model = fragment_model()
    model["bonds"].pop()
    with pytest.raises(ValueError, match="bond"):
        assign_current_fragments(model, REFERENCE)


def test_new_junction_requires_dedicated_template_even_with_existing_types():
    model = fragment_model()
    model["metadata"]["requires_parameterization"] = True
    for atom in model["atoms"]:
        atom["atom_type"] = "ca"
        atom["charge"] = 0
    with pytest.raises(ValueError, match="junction|parameterization"):
        assign_current_fragments(model, REFERENCE)


def test_normalization_copies_and_keeps_source_unchanged(tmp_path):
    source = REFERENCE / "EPH_gaff2.lib"
    before = source.read_bytes()
    report = normalize_libraries(REFERENCE, tmp_path / "normalized")
    normalized = (tmp_path / "normalized/EPH_gaff2.lib").read_text()
    assert "!entry.EPX." not in normalized
    assert '"HZ"' not in normalized
    assert source.read_bytes() == before
    assert report["EPH"]["source_sha256"] != report["EPH"]["normalized_sha256"]
    with pytest.raises(FileExistsError):
        normalize_libraries(REFERENCE, tmp_path / "normalized")


def test_all_normalized_libraries_roundtrip_with_parmed(tmp_path):
    import parmed as pmd

    output = tmp_path / "norm"
    report = normalize_libraries(REFERENCE, output)
    for name in report:
        parsed = pmd.amber.AmberOFFLibrary.parse(str(output / f"{name}_gaff2.lib"))
        assert list(parsed) == [name]
        source = load_fragment_library(REFERENCE / f"{name}_gaff2.lib", name)
        assert len(parsed[name].atoms) == len(source["atoms"])
        assert [a.charge for a in parsed[name].atoms] == pytest.approx(
            [a["charge"] for a in source["atoms"]]
        )
    assert any("boundbox" in x for x in report["CBD"]["warnings"])


def test_dedicated_template_does_not_accept_different_reaction_bond_orders(tmp_path):
    model = fragment_model("HYD")
    model["metadata"]["requires_parameterization"] = True
    frc = tmp_path / "junction.frcmod"
    frc.write_text("Dedicated externally supplied terms\n")
    template = {
        "atoms": [
            dict(
                a,
                charge=0.0,
                atom_type="c3" if a["element"] == "C" else "oh" if a["element"] == "O" else "ho",
            )
            for a in model["atoms"]
        ],
        "bonds": deepcopy(model["bonds"]),
        "frcmod": str(frc),
        "provenance": {"reference": "explicit test template"},
    }
    template["bonds"][0]["order"] = 2
    with pytest.raises(ValueError, match="bonds"):
        assign_current_fragments(model, REFERENCE, parameter_templates=template)


def test_source_epoxide_double_bond_is_audited():
    data = load_fragment_library(REFERENCE / "EPT_gaff2.lib", "EPT")
    assert any("double" in warning for warning in data["warnings"])
