#!/usr/bin/env python3
"""Prepare an exact-count builder recipe from C/O or O/C and group preferences."""

import argparse
import copy
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from go_builder.build import build_sheet
from go_builder.composition import resolve_composition
from go_builder.config import load_config
from go_builder.model import write_pdb

_ALLOWED_COMPOSITION_KEYS = {
    "C_over_O",
    "c_over_o",
    "O_over_C",
    "o_over_c",
    "preferences",
    "group_preferences",
    "priority",
    "mode",
    "enabled_groups",
    "ph",
    "pka",
    "deprotonation_fraction",
    "seed",
}

_SCAFFOLD_TOP_LEVEL_KEYS = {
    "name",
    "nx",
    "ny",
    "lx_nm",
    "ly_nm",
    "diameter_nm",
    "radius_nm",
    "shape",
    "edge_type",
    "seed",
    "pbc",
    "bond_length",
    "vacuum",
    "vacancies",
    "holes",
    "faces",
    "placement",
    "groups",
}


def _extract_recipe_and_spec(document: dict) -> tuple[dict, dict]:
    if "recipe" in document:
        recipe_raw = document["recipe"]
        if not isinstance(recipe_raw, dict):
            raise ValueError("Input config 'recipe' must be an object")
        conflicts = sorted(set(document) & _SCAFFOLD_TOP_LEVEL_KEYS)
        if conflicts:
            raise ValueError(
                f"Ambiguous config: 'recipe' cannot be combined with top-level keys {conflicts}"
            )
        if any(k in recipe_raw for k in ("composition", "c_over_o", "o_over_c")):
            raise ValueError(
                "Ambiguous config: nested 'recipe' must not contain composition keys"
            )
        if "composition" in document:
            if not isinstance(document["composition"], dict):
                raise ValueError("Input config 'composition' must be an object")
            extra_top = sorted(set(document) - {"recipe", "composition"})
            if extra_top:
                raise ValueError(
                    f"Ambiguous or unknown top-level keys alongside 'recipe' and 'composition': {extra_top}"
                )
            specification = dict(document["composition"])
        else:
            top_comp_keys = set(document) - {"recipe"}
            unknown = sorted(top_comp_keys - _ALLOWED_COMPOSITION_KEYS)
            if unknown:
                raise ValueError(f"Unknown top-level keys alongside 'recipe': {unknown}")
            specification = {k: document[k] for k in top_comp_keys}
            if not any(
                specification.get(k) is not None
                for k in ("c_over_o", "C_over_O", "o_over_c", "O_over_C")
            ):
                raise ValueError("Input config must contain a 'composition' object")
        base_recipe = copy.deepcopy(recipe_raw)
    else:
        if "composition" in document:
            if not isinstance(document["composition"], dict):
                raise ValueError("Input config 'composition' must be an object")
            if any(k in document for k in ("c_over_o", "C_over_O", "o_over_c", "O_over_C")):
                raise ValueError(
                    "Ambiguous config: specify either 'composition' or top-level 'c_over_o'/'o_over_c', not both"
                )
            specification = dict(document["composition"])
            base_recipe = {
                k: copy.deepcopy(v) for k, v in document.items() if k != "composition"
            }
        elif any(k in document for k in ("c_over_o", "C_over_O", "o_over_c", "O_over_C")):
            flat_comp_keys = {
                "c_over_o",
                "C_over_O",
                "o_over_c",
                "O_over_C",
                "preferences",
                "group_preferences",
                "priority",
                "mode",
                "enabled_groups",
                "ph",
                "pka",
                "deprotonation_fraction",
            }
            specification = {k: document[k] for k in flat_comp_keys if k in document}
            base_recipe = {
                k: copy.deepcopy(v) for k, v in document.items() if k not in flat_comp_keys
            }
        else:
            raise ValueError("Input config must contain a 'composition' object")

    extra = set(specification) - _ALLOWED_COMPOSITION_KEYS
    if extra:
        raise ValueError(f"Unknown composition keys: {sorted(extra)}")
    return base_recipe, specification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, type=Path, help="JSON or TOML with recipe and composition objects"
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Resolved exact-count recipe JSON"
    )
    parser.add_argument("--report", type=Path, help="Composition/provenance report JSON")
    parser.add_argument("--model", type=Path, help="Also build and validate the starting structure")
    parser.add_argument("--pdb", type=Path, help="Also export built structure as finite PDB")
    args = parser.parse_args()

    report_path = args.report or args.output.with_suffix(".composition.json")
    all_output_paths = [args.output, report_path]
    if args.model is not None:
        all_output_paths.append(args.model)
    if args.pdb is not None:
        all_output_paths.append(args.pdb)
    resolved_paths = [p.resolve() for p in all_output_paths]
    if len(set(resolved_paths)) != len(resolved_paths):
        raise ValueError("Output paths must be distinct")
    for path in all_output_paths:
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)

    document = load_config(args.input)
    base_recipe, specification = _extract_recipe_and_spec(document)

    recipe, report = resolve_composition(
        base_recipe,
        c_over_o=specification.get("c_over_o", specification.get("C_over_O")),
        preferences=specification.get("preferences", specification.get("group_preferences")),
        o_over_c=specification.get("o_over_c", specification.get("O_over_C")),
        priority=specification.get("priority"),
        mode=specification.get("mode", "balanced"),
        enabled_groups=specification.get("enabled_groups"),
        seed=specification.get("seed"),
        ph=specification.get("ph"),
        pka=specification.get("pka", 4.5),
        deprotonation_fraction=specification.get("deprotonation_fraction"),
    )
    outputs = {
        args.output: recipe,
        report_path: report,
    }
    model = None
    if args.model is not None or args.pdb is not None:
        model = build_sheet(recipe)
        if args.pdb is not None and any(any(b.get("image", [0, 0, 0])) for b in model["bonds"]):
            raise ValueError("PDB cannot preserve periodic image bonds; use JSON")
        counts = Counter(a["element"] for a in model["atoms"])
        assert (
            counts["C"] == report["predicted_carbon_atoms"]
            and counts["O"] == report["predicted_oxygen_atoms"]
        )
        report["status"] = "BUILT_AND_COMPOSITION_VERIFIED"
        report["actual_composition"] = dict(counts)
        if args.model is not None:
            outputs[args.model] = model
    for path, value in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as handle:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
    if args.pdb is not None and model is not None:
        write_pdb(model, args.pdb)
    print(
        json.dumps(
            {
                "C_over_O": report["predicted_C_over_O"],
                "O_over_C": report["predicted_O_over_C"],
                "groups": recipe["groups"],
                "status": report["status"],
            }
        )
    )


if __name__ == "__main__":
    main()
