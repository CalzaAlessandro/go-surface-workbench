"""Portable, explicit LayerOne and illustrative template creation recipes; no scientific engines imported."""

from __future__ import annotations

import copy
import hashlib
import json
from importlib.resources import files

LAYERONE_PRESET_NAMES = ("GOX8-P935", "REGOX8-P980", "REGOX8-P995", "REGOX8-P999")
TEMPLATE_PRESET_NAMES = ("template-go-ph7", "template-go-island", "template-rgo")
PRESET_NAMES = LAYERONE_PRESET_NAMES + TEMPLATE_PRESET_NAMES

_TEMPLATE_RECIPES: dict[str, dict] = {
    "template-go-ph7": {
        "name": "template-go-ph7",
        "lx_nm": 4.0,
        "ly_nm": 4.0,
        "shape": "rectangular",
        "pbc": [False, False, False],
        "faces": "both",
        "seed": 42,
        "composition": {
            "o_over_c": 0.35,
            "priority": ["OH", "epoxide", "COOH", "COO"],
            "ph": 7.4,
        },
        "placement": "random",
    },
    "template-go-island": {
        "name": "template-go-island",
        "lx_nm": 5.0,
        "ly_nm": 5.0,
        "shape": "rectangular",
        "pbc": [False, False, False],
        "faces": "both",
        "seed": 42,
        "holes": {"distribution": "mixed", "percent_vacancies": 8.0},
        "composition": {
            "o_over_c": 0.35,
            "priority": ["OH", "epoxide", "COOH", "COO"],
            "ph": 7.4,
        },
        "placement": {"mode": "island", "n_islands": 2},
    },
    "template-rgo": {
        "name": "template-rgo",
        "lx_nm": 4.0,
        "ly_nm": 4.0,
        "shape": "rectangular",
        "pbc": [False, False, False],
        "faces": "both",
        "seed": 42,
        "composition": {
            "o_over_c": 0.12,
            "priority": ["OH", "epoxide", "COOH", "COO"],
            "ph": 7.4,
        },
        "placement": "random",
    },
}

_TEMPLATE_DESCRIPTIONS: dict[str, str] = {
    "template-go-ph7": "Illustrative 4x4 nm oxidized GO sheet at pH 7.4 (O/C=0.35, random placement).",
    "template-go-island": "Illustrative 5x5 nm porous GO sheet with 8% mixed vacancies and 2 oxygen islands (O/C=0.35, pH 7.4).",
    "template-rgo": "Illustrative 4x4 nm reduced GO sheet at pH 7.4 (O/C=0.12, random placement).",
}


def _resource(name: str):
    if name not in LAYERONE_PRESET_NAMES:
        raise ValueError(f"Unknown LayerOne preset {name!r}; choose one of {', '.join(PRESET_NAMES)}")
    return files("go_builder").joinpath("data", "layerone", f"{name}.json")


def _source_info_resource():
    return files("go_builder").joinpath("data", "layerone", "source_info.json")


def load_source_info() -> dict:
    """Read public provenance and validation status for the preset collection."""
    value = json.loads(_source_info_resource().read_text(encoding="utf-8"))
    if not isinstance(value, dict) or not isinstance(value.get("presets"), dict):
        raise ValueError("LayerOne source_info.json must contain a presets object")
    return value


def load_preset(name: str) -> dict:
    """Read a fresh packaged recipe; preserve its explicit seed and composition."""
    if name in _TEMPLATE_RECIPES:
        return copy.deepcopy(_TEMPLATE_RECIPES[name])
    if name not in PRESET_NAMES:
        raise ValueError(f"Unknown preset {name!r}; choose one of {', '.join(PRESET_NAMES)}")
    recipe = json.loads(_resource(name).read_text(encoding="utf-8"))
    if not isinstance(recipe, dict):
        raise ValueError(f"Preset {name} must contain a recipe object")
    return recipe


def resolve_preset(name: str, *, seed: int | None = None) -> tuple[dict, dict]:
    """Return the construction recipe and provenance kept outside builder keys."""
    if name not in PRESET_NAMES:
        raise ValueError(f"Unknown preset {name!r}; choose one of {', '.join(PRESET_NAMES)}")
    if name in _TEMPLATE_RECIPES:
        source = copy.deepcopy(_TEMPLATE_RECIPES[name])
        recipe = copy.deepcopy(source)
        raw = json.dumps(source, sort_keys=True, separators=(",", ":")).encode("utf-8")
        if seed is not None:
            if type(seed) is not int or seed < 0:
                raise ValueError("seed must be a nonnegative integer")
            recipe["seed"] = seed
        comp = source["composition"]
        source_info = {
            "build_enabled": True,
            "readiness": "illustrative_starting_template_not_measured",
            "basis": "illustrative starting template, not a measured composition",
            "description": _TEMPLATE_DESCRIPTIONS[name],
            "target_o_over_c": comp["o_over_c"],
            "ph": comp["ph"],
        }
        provenance = {
            "name": name,
            "family": "Template",
            "resource": f"go_builder:templates/{name}",
            "recipe_sha256": hashlib.sha256(raw).hexdigest(),
            "source_recipe": source,
            "seed_override": seed,
            "source_info": source_info,
            "source_composition": {
                "basis": "illustrative_starting_template_not_measured",
                "target_O_C_atomic_ratio": comp["o_over_c"],
                "priority": list(comp["priority"]),
                "ph": comp["ph"],
            },
            "constructed_composition": {
                "basis": "illustrative starting template; group counts resolved by composition solver at build time",
                "lx_nm": source["lx_nm"],
                "ly_nm": source["ly_nm"],
                "holes": source.get("holes"),
                "placement": source.get("placement"),
            },
        }
        return recipe, provenance

    resource = _resource(name)
    raw = resource.read_bytes()
    source = json.loads(raw)
    recipe = json.loads(raw)
    source_info_raw = _source_info_resource().read_bytes()
    source_info = json.loads(source_info_raw)
    try:
        entry = source_info["presets"][name]
    except (KeyError, TypeError) as error:
        raise ValueError(f"No public source metadata recorded for preset {name}") from error
    if entry.get("recipe_file") != f"{name}.json":
        raise ValueError(f"Source metadata for {name} points to an unexpected recipe file")
    if seed is not None:
        if type(seed) is not int or seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        recipe["seed"] = seed
    provenance = {
        "name": name,
        "family": "LayerOne",
        "resource": f"go_builder:data/layerone/{name}.json",
        "recipe_sha256": hashlib.sha256(raw).hexdigest(),
        "source_recipe": source,
        "seed_override": seed,
        "source_info": entry,
        "source_info_sha256": hashlib.sha256(source_info_raw).hexdigest(),
        "source_composition": entry["source_composition"],
        "constructed_composition": {
            "basis": "construction recipe, not a measured composition",
            "groups": source.get("groups", {}),
            "nx": source.get("nx"),
            "ny": source.get("ny"),
        },
    }
    return recipe, provenance


def list_presets() -> list[dict]:
    """Describe packaged options without importing a graph/force-field engine."""
    result = []
    for name in PRESET_NAMES:
        recipe, provenance = resolve_preset(name)
        result.append(
            {
                "name": name,
                "family": provenance["family"],
                "resource": provenance["resource"],
                "default_seed": recipe.get("seed"),
                "source_composition": provenance["source_composition"],
                "constructed_composition": provenance["constructed_composition"],
                "build_enabled": provenance["source_info"]["build_enabled"],
                "readiness": provenance["source_info"]["readiness"],
                "source_info": provenance["source_info"],
            }
        )
    return result
