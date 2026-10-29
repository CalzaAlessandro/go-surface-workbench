"""An imported current basal epoxide can be annotated from explicit topology."""

from go_builder.build import build_sheet
from go_builder.functionalize import open_epoxides
from go_builder.topology import annotate_imported_epoxides


def test_explicit_import_annotation_can_be_functionalized_without_a_fake_cell():
    model = build_sheet({"nx": 5, "ny": 5, "groups": {"epoxide": 1}})
    model["cell"] = None
    model["metadata"] = {
        "representation": "amber_import",
        "bond_orders_unknown": True,
        "formal_charges_unknown": True,
    }
    for atom in model["atoms"]:
        atom["metadata"] = {}
    annotated = annotate_imported_epoxides(model)
    assert len(annotated["metadata"]["groups"]) == 1
    fragment = {
        "atoms": [
            {"id": 1, "element": "N", "position": [0, 0, 0], "formal_charge": 0},
            {"id": 2, "element": "H", "position": [0, 0, 1.01], "formal_charge": 0},
            {"id": 3, "element": "H", "position": [0.95, 0, -0.34], "formal_charge": 0},
            {"id": 4, "element": "H", "position": [-0.47, 0.82, -0.34], "formal_charge": 0},
        ],
        "bonds": [{"a": 1, "b": i, "order": 1, "image": [0, 0, 0]} for i in [2, 3, 4]],
        "attachment_atom": 1,
        "transfer_hydrogen": 2,
    }
    opened = open_epoxides(annotated, {"fragment": fragment, "count": 1})
    assert opened["cell"] is None
    assert len(opened["atoms"]) == len(model["atoms"]) + 4
    assert opened["metadata"]["requires_parameterization"] is True
