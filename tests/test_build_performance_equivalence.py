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
            # Passo 1a raised the hard H steric cutoff from 0.80 A to 0.90 A.
            cutoff = 0.90 if "H" in (atom["element"], other["element"]) else 1.25
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
            _model(["H", "C"], [[0, 0, 0], [0.899999, 0, 0]]),
            {1},
            False,
        ),
        (
            "hydrogen exactly at cutoff",
            _model(["H", "C"], [[0, 0, 0], [0.90, 0, 0]]),
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
        # Updated for Passo 1a geometry upgrade: _cap_edges samples 20/-20 deg tilts when
        # armchair-bay H-H clearance is <= 1.80 A, and OH/epoxide/COOH/COO use physical
        # geometries with azimuth/torsion/pucker clearance optimization.
        (
            {"nx": 8, "ny": 8, "seed": 17},
            "a3a6d9eda2783669d516df2ea2845ad93248ccee84493c9f85a09a613a17d0fb",
        ),
        (
            {"nx": 8, "ny": 8, "seed": 17, "holes": [{"size": 4, "seed_atom": 145}]},
            "13b0bc8ff82914dda5da34200a8431ef5d0388dc5a4e0ccf7ee3a2aa4267b2a6",
        ),
        (
            {"nx": 8, "ny": 8, "seed": 17, "pbc": [True, True, False]},
            "dbd282b10e4a328faab18dd443afcaeb783318a1052dfc61f1f29f64749c7dbd",
        ),
        # Updated for Passo 1d-1: _orthonormal_edge_frame fixes non-orthogonal edge
        # vectors after sp3 puckering, and unified 1-4+ candidate scoring selects
        # conflict-free geometries (cases 0-2 remain bit-for-bit identical).
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
            "da1094d04a1c9abeaeaed83f7a90c008d74819850360f0760b5ee8532d01a325",
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
