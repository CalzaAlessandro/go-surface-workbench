"""Compatibility gates for speed-only graph construction changes."""

import hashlib
import json

import numpy as np
import pytest

from go_builder.build import _collision_free_vectorized, build_sheet
from go_builder.errors import PlacementError
from go_builder.model import canonical_bond


def _model(elements, positions, *, pbc=(False, False, False), cell=None, bonds=()):
    atoms = [
        {
            "id": index,
            "element": element,
            "position": list(position),
        }
        for index, (element, position) in enumerate(zip(elements, positions, strict=True), 1)
    ]
    return {
        "atoms": atoms,
        "bonds": list(bonds),
        "pbc": list(pbc),
        "cell": cell or [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 40.0]],
    }


def _reference_collision_free(model, new_ids):
    """Original scalar collision loop retained as an independent test oracle."""
    from go_builder.build import _minimum_image

    bonded = {frozenset((b["a"], b["b"])) for b in model["bonds"]}
    atoms = model["atoms"]
    coords = np.array([a["position"] for a in atoms])
    for index, atom in enumerate(atoms):
        if atom["id"] not in new_ids:
            continue
        distances = np.linalg.norm(_minimum_image(coords - coords[index], model), axis=1)
        for j, other in enumerate(atoms):
            if j == index or frozenset((atom["id"], other["id"])) in bonded:
                continue
            cutoff = 0.8 if "H" in (atom["element"], other["element"]) else 1.25
            if distances[j] < cutoff:
                return False
    return True


def _collision_cases():
    return [
        (
            "heavy just inside cutoff",
            _model(["C", "C"], [[0, 0, 0], [1.249999, 0, 0]]),
            {1},
            False,
        ),
        (
            "heavy exactly at cutoff",
            _model(["C", "C"], [[0, 0, 0], [1.25, 0, 0]]),
            {1},
            True,
        ),
        (
            "hydrogen just inside cutoff",
            _model(["H", "C"], [[0, 0, 0], [0.799999, 0, 0]]),
            {1},
            False,
        ),
        (
            "hydrogen exactly at cutoff",
            _model(["H", "C"], [[0, 0, 0], [0.8, 0, 0]]),
            {1},
            True,
        ),
        (
            "bonded pair excluded",
            _model(["C", "O"], [[0, 0, 0], [0.2, 0, 0]], bonds=[canonical_bond(1, 2)]),
            {1},
            True,
        ),
        (
            "self excluded",
            _model(["O"], [[0, 0, 0]]),
            {1},
            True,
        ),
        (
            "periodic minimum image overlaps",
            _model(
                ["C", "C"],
                [[0.1, 0.1, 0], [9.6, 9.6, 0]],
                pbc=(True, True, False),
            ),
            {1},
            False,
        ),
        (
            "same coordinates without periodicity do not overlap",
            _model(["C", "C"], [[0.1, 0.1, 0], [9.6, 9.6, 0]]),
            {1},
            True,
        ),
        (
            "multiple new atoms preserve atom-order outcome",
            _model(["O", "H", "C"], [[0, 0, 0], [5, 0, 0], [0.7, 0, 0]]),
            {1, 2},
            False,
        ),
    ]


@pytest.mark.parametrize("label,model,new_ids,expected", _collision_cases())
def test_vectorized_collision_screen_matches_original_scalar_reference(
    label, model, new_ids, expected
):
    assert _reference_collision_free(model, new_ids) is expected
    assert _collision_free_vectorized(model, new_ids) is expected


def _model_hash(model):
    payload = json.dumps(model, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize(
    "recipe,expected_hash",
    [
        (
            {"nx": 8, "ny": 8, "seed": 17},
            "de2b9b8ef8af63edeed6dc52493753010aadcaebc05798503d18e0c985f2e5e2",
        ),
        (
            {"nx": 8, "ny": 8, "seed": 17, "holes": [{"size": 4, "seed_atom": 145}]},
            "61a8d986e9ade5b638b13b5afd0a3b412e840ba663cfdd6fbfd5d142cac18458",
        ),
        (
            {"nx": 8, "ny": 8, "seed": 17, "pbc": [True, True, False]},
            "dbd282b10e4a328faab18dd443afcaeb783318a1052dfc61f1f29f64749c7dbd",
        ),
        (
            {
                "nx": 6,
                "ny": 6,
                "seed": 17,
                "faces": "down",
                "groups": {
                    "OH": {"id": [58]},
                    "carbonyl": {"id": [26]},
                    "epoxide": 3,
                    "COOH": 1,
                    "COO": 1,
                },
            },
            "35d9c7759e8503c210e9c35062da7af18af4e7b510d98e052918dffe847f0a5f",
        ),
    ],
)
def test_same_seed_models_keep_baseline_serialized_hash(recipe, expected_hash):
    assert _model_hash(build_sheet(recipe)) == expected_hash


def test_placement_failure_report_remains_identical():
    with pytest.raises(PlacementError) as caught:
        build_sheet({"nx": 6, "ny": 6, "seed": 17, "groups": {"OH": 1000}})
    assert caught.value.as_report() == {
        "status": "INCOMPLETE",
        "message": "Requested functional group counts exceed eligible edge/basal sites",
        "requested": {"OH": 1000, "epoxide": 0, "carbonyl": 0, "COOH": 0, "COO": 0},
        "capacity": {"eligible_edge_carbons": 32, "eligible_basal_carbons": 98},
        "attempts": 0,
        "suggestions": ["Enlarge the sheet", "Reduce requested group counts", "Change the seed"],
    }
