"""Portable, explicit LayerOne creation recipes; no scientific engines imported."""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files

PRESET_NAMES = ("GOX8-P935", "REGOX8-P980", "REGOX8-P995", "REGOX8-P999")


def _resource(name: str):
    if name not in PRESET_NAMES:
        raise ValueError(f"Unknown preset {name!r}; choose one of {', '.join(PRESET_NAMES)}")
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
    recipe = json.loads(_resource(name).read_text(encoding="utf-8"))
    if not isinstance(recipe, dict):
        raise ValueError(f"Preset {name} must contain a recipe object")
    return recipe


def resolve_preset(name: str, *, seed: int | None = None) -> tuple[dict, dict]:
    """Return the construction recipe and provenance kept outside builder keys."""
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
