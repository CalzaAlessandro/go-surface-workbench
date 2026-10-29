"""New candidates come from pinned bases, preserve type contracts and expose gaps."""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_FILES = (
    "gaff.dat",
    "gaff2.dat",
    "parm99.dat",
    "frcmod.ff99SB",
    "source_manifest.json",
)
MISSING_EXTERNAL = [
    name
    for name in EXTERNAL_FILES
    if not (ROOT / "references/system_forcefield_sources" / name).is_file()
]
requires_installed_amber_data = pytest.mark.skipif(
    bool(MISSING_EXTERNAL),
    reason="External installed Amber DAT sources are intentionally not vendored; missing: "
    + ", ".join(MISSING_EXTERNAL),
)


@requires_installed_amber_data
def test_exact_installed_bases_and_literature_version_distinction():
    from go_param.ff_candidates import parse_base_dat

    gaff = parse_base_dat(ROOT / "references/system_forcefield_sources/gaff.dat")
    gaff2 = parse_base_dat(ROOT / "references/system_forcefield_sources/gaff2.dat")
    assert "Version 1.81" in gaff["header"]
    assert "Version 2.2.20" in gaff2["header"]
    first = next(e for e in gaff["sections"]["BOND"] if e["types"] == ["ca", "ca"])
    second = next(e for e in gaff2["sections"]["BOND"] if e["types"] == ["ca", "ca"])
    assert first["values"] == [461.1, 1.3984]
    assert second["values"][0] > 0
    assert len(gaff2["sections"]["DIHE"]) > 1000


@requires_installed_amber_data
def test_candidates_have_native_and_isolated_ring_type_contracts(tmp_path):
    from go_param.ff_candidates import generate_candidates

    result = generate_candidates(ROOT, output_dir=tmp_path)
    assert set(result["candidates"]) == {
        "gaff_native",
        "gaff2_native",
        "amber99sb_go_analogy",
    }
    for name, manifest in result["candidates"].items():
        assert manifest["status"] == "unfitted_baseline_not_validated_for_production"
        assert manifest["parameter_counts"]["BOND"] > 0
        assert manifest["parameter_counts"]["ANGLE"] > 0
        assert manifest["parameter_counts"]["DIHE"] > 0
        assert manifest["parameter_counts"]["NONBON"] > 0
        assert (tmp_path / name / "candidate.frcmod").is_file()
        assert len(manifest["model_contracts"]) == 8
        assert manifest["gaps"] is not None
        assert manifest["requires_graph_specific_typing_validation"]
    assert result["candidates"]["gaff2_native"]["role_types"]["epoxide_O"] == "op"
    amber = result["candidates"]["amber99sb_go_analogy"]
    assert amber["role_types"]["epoxide_C"] == "GX"
    assert amber["role_types"]["epoxide_O"] == "GO"
    text = (tmp_path / "amber99sb_go_analogy/candidate.frcmod").read_text()
    assert "GX-GO-GX" in text
    assert "CY-CX-CX-CY" not in text
    assert (
        "addAtomTypes"
        in (tmp_path / "amber99sb_go_analogy/atom_types.leap").read_text()
    )


@requires_installed_amber_data
def test_amber_lj_equivalence_and_atom_identity_are_explicit(tmp_path):
    from go_param.ff_candidates import generate_candidates, parse_base_dat

    base = parse_base_dat(ROOT / "references/system_forcefield_sources/parm99.dat")
    assert base["lj_equivalences"]["CA"] == "C*"
    report = generate_candidates(ROOT, output_dir=tmp_path)
    contract = report["candidates"]["amber99sb_go_analogy"]["model_contracts"][
        "pyrene_COO"
    ]
    source = json.loads((ROOT / "models/pyrene_COO/model.json").read_text())
    assert [a["atom_id"] for a in contract["atoms"]] == [
        a["id"] for a in source["atoms"]
    ]
    assert (
        next(a for a in contract["atoms"] if a["role"] == "COO_C")["atom_type"] == "C"
    )
    assert {a["atom_type"] for a in contract["atoms"] if a["role"] == "COO_O"} == {"O2"}
    assert contract["partial_charges"] == "not_assigned_by_frcmod_generation"


@requires_installed_amber_data
def test_new_files_load_with_parameter_reader(tmp_path):
    AmberParameterSet = pytest.importorskip(
        "parmed.amber",
        reason="ParmEd is required for optional Amber parameter-reader integration",
    ).AmberParameterSet

    from go_param.ff_candidates import generate_candidates

    generate_candidates(ROOT, output_dir=tmp_path)
    for name in ("gaff_native", "gaff2_native", "amber99sb_go_analogy"):
        parameters = AmberParameterSet(str(tmp_path / name / "candidate.frcmod"))
        assert parameters.bond_types
        assert parameters.dihedral_types


def test_ring_context_collisions_are_reported_without_installed_data():
    from go_param.ff_candidates import _contract

    model = {
        "model_id": "synthetic_two_epoxides",
        "formal_charge": 0,
        "atoms": [
            {
                "id": i,
                "element": "C" if i in (1, 2, 3, 5) else "O",
                "role": "epoxide_C" if i in (1, 2, 3, 5) else "epoxide_O",
            }
            for i in range(1, 7)
        ],
        "bonds": [
            {"a": a, "b": b}
            for a, b in [(1, 2), (2, 3), (3, 5), (1, 4), (2, 4), (3, 6), (5, 6)]
        ],
    }
    contract = _contract(
        model,
        {"epoxide_C": "cx", "epoxide_O": "op"},
        {s: [] for s in ("BOND", "ANGLE", "DIHE", "IMPROPER")},
    )
    assert any(
        g.get("code") == "three_member_carbon_angle_used_outside_carbon_triangle"
        for g in contract["graph_term_gaps"]
    )
    assert any(
        g.get("code") == "cx_cx_bond_spans_distinct_epoxide_rings"
        for g in contract["graph_term_gaps"]
    )


def test_parse_base_dat_with_self_contained_synthetic_fixture(tmp_path):
    from go_param.ff_candidates import parse_base_dat

    path = tmp_path / "synthetic.dat"
    path.write_text("""SYNTHETIC PARSER FIXTURE - NOT FORCE FIELD DATA
ca 12.0 aromatic fixture
ha 1.0 hydrogen fixture

ha
ca-ca 100.0 1.4 fixture bond

ca-ca-ca 20.0 120.0 fixture angle

X -ca-ca-X 1 0.5 180.0 -2.0
X -ca-ca-X 1 0.2 0.0 1.0

X -X -ca-ha 1.0 180.0 2.0

C* ca

MOD4 RE
C* 1.9 0.08
ha 1.4 0.01

END
""")
    parsed = parse_base_dat(path)
    assert parsed["header"].startswith("SYNTHETIC")
    assert parsed["sections"]["BOND"][0]["values"] == [100.0, 1.4]
    assert parsed["sections"]["BOND"][0]["line_number"] == 6
    assert len(parsed["sections"]["DIHE"]) == 2
    assert parsed["sections"]["DIHE"][0]["values"][-1] == -2.0
    assert parsed["lj_equivalences"] == {"ca": "C*"}
