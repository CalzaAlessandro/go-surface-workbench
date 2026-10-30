import csv
import json
from pathlib import Path

import pytest

from go_param.context import analyze_dataset, write_outputs

ROOT = Path(__file__).resolve().parents[1]


def _input(tmp_path, model_name="pyrene_diol_epoxide", *, transform=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    model_path = ROOT / "models" / model_name / "model.json"
    model = json.loads(model_path.read_text())
    if transform:
        transform(model)
    model_copy = tmp_path / f"{model_name}.model.json"
    model_copy.write_text(json.dumps(model))
    charges = {
        "model_id": model["model_id"],
        "formal_charge": model["formal_charge"],
        "atoms": [
            {
                "atom_id": atom["id"],
                "element": atom["element"],
                "charge_e": (atom["id"] % 7 - 3) / 1000,
            }
            for atom in model["atoms"]
        ],
    }
    charges["atoms"][-1]["charge_e"] -= sum(a["charge_e"] for a in charges["atoms"])
    charges_path = tmp_path / f"{model_name}.charges.json"
    charges_path.write_text(json.dumps(charges))
    return {"model_path": str(model_copy), "charges_path": str(charges_path)}, model


def test_lists_all_three_functional_group_pairs_using_carbon_scaffold_paths(tmp_path):
    item, _ = _input(tmp_path)

    result = analyze_dataset([item])

    model = result["models"][0]
    assert len(model["functional_groups"]) == 3
    assert len(model["functional_group_pairs"]) == 3
    assert all(
        all(
            atom_id in model["scaffold_carbon_ids"]
            for atom_id in pair["graph_path_atom_ids"]
        )
        for pair in model["functional_group_pairs"]
    )
    assert all(
        pair["graph_distance_definition"]
        == "minimum number of C-C scaffold bonds between group scaffold-carbon anchors"
        for pair in model["functional_group_pairs"]
    )


def test_carboxyl_anchor_belongs_to_group_and_exocyclic_carbon_is_not_background(
    tmp_path,
):
    item, model = _input(tmp_path, "pyrene_COOH")

    result = analyze_dataset([item])

    output = result["models"][0]
    group = next(
        group for group in output["functional_groups"] if group["block_type"] == "COOH"
    )
    anchor = next(atom for atom in model["atoms"] if atom["role"] == "COOH_anchor")
    cooh_carbon = next(atom for atom in model["atoms"] if atom["role"] == "COOH_C")
    assert anchor["id"] in group["atom_ids"]
    assert anchor["id"] not in output["background_carbon_atom_ids"]
    assert cooh_carbon["id"] not in output["scaffold_carbon_ids"]
    assert cooh_carbon["id"] not in output["background_carbon_atom_ids"]


def test_background_charge_rows_keep_interior_and_hydrogen_edge_carbon_separate(
    tmp_path,
):
    item, _ = _input(tmp_path, "pyrene_pristine")

    result = analyze_dataset([item])["models"][0]

    rows = result["background_carbon_distance_charges"]
    by_kind = {row["background_kind"] for row in rows}
    assert by_kind == {"interior_scaffold", "hydrogen_terminated_edge"}
    assert result["coverage"]["background_kind_counts"] == {
        "interior_scaffold": 6,
        "hydrogen_terminated_edge": 10,
    }
    assert all(row["role"] for row in rows)


def test_declared_scaffold_ids_cannot_include_exocyclic_carboxyl_carbon(tmp_path):
    def corrupt_scaffold(model):
        model["metadata"]["scaffold_carbon_ids"] = [
            atom["id"] for atom in model["atoms"] if atom["element"] == "C"
        ]

    item, _ = _input(tmp_path, "pyrene_COOH", transform=corrupt_scaffold)

    with pytest.raises(ValueError, match="exocyclic"):
        analyze_dataset([item])


def test_exocyclic_carbonyl_groups_anchor_to_their_scaffold_neighbors(tmp_path):
    def mark_carbonyls_exocyclic(model):
        model["metadata"]["scaffold_carbon_ids"] = [
            atom["id"]
            for atom in model["atoms"]
            if atom["element"] == "C" and atom["role"] != "carbonyl_C"
        ]

    item, _ = _input(tmp_path, "pyrene_dione", transform=mark_carbonyls_exocyclic)

    result = analyze_dataset([item])["models"][0]

    assert len(result["functional_group_pairs"]) == 1
    pair = result["functional_group_pairs"][0]
    assert pair["graph_distance"] is not None
    assert all(
        atom_id in result["scaffold_carbon_ids"]
        for atom_id in pair["graph_path_atom_ids"]
    )
    assert all(
        group["scaffold_anchor_source"] == "directly bonded scaffold carbons"
        for group in result["functional_groups"]
    )


def test_graph_distance_does_not_change_when_only_coordinates_are_distorted(tmp_path):
    item_a, _ = _input(tmp_path / "a")

    # Make a separate output directory and stretch all coordinates nonuniformly.
    def distort(model):
        for atom in model["atoms"]:
            x, y, z = atom["position"]
            atom["position"] = [5 * x, 0.2 * y, 3 * z]

    item_b, _ = _input(tmp_path / "b", transform=distort)
    result_a = analyze_dataset([item_a])["models"][0]
    result_b = analyze_dataset([item_b])["models"][0]

    assert [p["graph_distance"] for p in result_a["functional_group_pairs"]] == [
        p["graph_distance"] for p in result_b["functional_group_pairs"]
    ]
    assert [
        p["carbon_centroid_distance_angstrom"]
        for p in result_a["functional_group_pairs"]
    ] != [
        p["carbon_centroid_distance_angstrom"]
        for p in result_b["functional_group_pairs"]
    ]


@pytest.mark.parametrize("corruption", ["atom_order", "element", "nonneutral"])
def test_charge_data_must_match_model_atoms_and_whole_charge(tmp_path, corruption):
    item, _model = _input(tmp_path)
    charges_path = Path(item["charges_path"])
    charges = json.loads(charges_path.read_text())
    if corruption == "atom_order":
        charges["atoms"][0]["atom_id"], charges["atoms"][1]["atom_id"] = 2, 1
    elif corruption == "element":
        charges["atoms"][0]["element"] = "O"
    else:
        charges["atoms"][0]["charge_e"] += 0.01
    charges_path.write_text(json.dumps(charges))

    with pytest.raises(ValueError):
        analyze_dataset([item])


def test_fit_evidence_is_preserved_as_provenance_not_acceptance(tmp_path):
    item, _ = _input(tmp_path)
    evidence_path = tmp_path / "fit.json"
    evidence = {
        "status": "diagnostic_only",
        "physical_conformers": 1,
        "orientations": 3,
        "protocol": {
            "method": "HF/6-31G* MK ESP; Amber hyperbolic RESP",
            "scope": "finite whole-model reference",
            "physical_conformers": 1,
            "orientations": 3,
            "parent_qm_optimization_frequencies_stability": "pending external gate",
        },
    }
    evidence_path.write_text(json.dumps(evidence))
    item["fit_evidence_path"] = str(evidence_path)

    result = analyze_dataset([item])["models"][0]

    assert result["provenance"]["fit_status"] == "diagnostic_only"
    assert result["provenance"]["qm_method"] == evidence["protocol"]["method"]
    assert result["provenance"]["scope"] == evidence["protocol"]["scope"]
    assert result["provenance"]["upstream_qm_gate"] == "pending external gate"
    assert result["provenance"]["fit_status_is_acceptance"] is False
    assert result["provenance"]["physical_conformers"] == 1
    assert result["provenance"]["orientations"] == 3


def test_diagnostic_model_scope_overrides_generic_protocol_scope_in_json_and_csv(
    tmp_path,
):
    def mark_model_diagnostic(model):
        model["metadata"].update(
            scope="rhf_instability_diagnostic_control_NOT_final_electronic_reference",
            wavefunction_stability_validated=False,
        )

    item, _ = _input(tmp_path / "input", transform=mark_model_diagnostic)
    evidence_path = tmp_path / "fit.json"
    evidence_path.write_text(
        json.dumps(
            {
                "status": "atomic_resp_fit_completed",
                "protocol": {
                    "scope": "finite whole-model atomic reference",
                    "method": "HF/6-31G* MK ESP; Amber hyperbolic RESP",
                },
            }
        )
    )
    item["fit_evidence_path"] = str(evidence_path)

    result = analyze_dataset([item])
    provenance = result["models"][0]["provenance"]
    output_dir = write_outputs(result, tmp_path / "outputs")
    with (output_dir / "per_atom_charges.csv").open(newline="") as handle:
        csv_provenance = next(csv.DictReader(handle))

    diagnostic_scope = (
        "rhf_instability_diagnostic_control_NOT_final_electronic_reference"
    )
    protocol_scope = "finite whole-model atomic reference"
    assert provenance["scope"] == diagnostic_scope
    assert provenance["protocol_scope"] == protocol_scope
    assert provenance["is_diagnostic"] is True
    assert provenance["reported_wavefunction_stability_validated"] is False
    assert provenance["fit_status_is_acceptance"] is False
    assert csv_provenance["scope"] == diagnostic_scope
    assert csv_provenance["protocol_scope"] == protocol_scope
    assert csv_provenance["is_diagnostic"] == "True"
    assert csv_provenance["reported_wavefunction_stability_validated"] == "False"


def test_separate_charge_fits_for_same_model_are_separate_contexts(tmp_path):
    item_a, _ = _input(tmp_path / "fit_a")
    item_b, _ = _input(tmp_path / "fit_b")
    result = analyze_dataset([item_a, item_b])

    assert result["context_count"] == 2
    assert result["model_count"] == 1
    assert result["models"][0]["context_id"] != result["models"][1]["context_id"]
