"""Lightweight TOML and JSON configuration I/O and annotated templates (stdlib only)."""

from __future__ import annotations

import copy
import json
import math
import re
import tomllib
from pathlib import Path
from typing import Any

_BARE_KEY_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def _toml_escape_str(value: str) -> str:
    return json.dumps(value, ensure_ascii=False).replace("\x7f", "\\u007f")


def _format_toml_key(key: str) -> str:
    if _BARE_KEY_RE.match(key):
        return key
    return _toml_escape_str(key)


def _format_toml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"TOML does not support non-finite float: {value!r}")
        return repr(value)
    if isinstance(value, str):
        return _toml_escape_str(value)
    raise TypeError(f"Unsupported TOML scalar type: {type(value).__name__}")


def _is_dict_list(value: list[Any]) -> bool:
    return bool(value) and all(isinstance(item, dict) for item in value)


def _format_toml_list(values: list[Any]) -> str:
    parts: list[str] = []
    for item in values:
        if item is None:
            raise ValueError("TOML arrays cannot contain None/null values")
        if isinstance(item, list):
            parts.append(_format_toml_list(item))
        elif isinstance(item, dict):
            pairs = [
                f"{_format_toml_key(str(k))} = {_format_toml_value(v)}"
                for k, v in item.items()
                if v is not None
            ]
            parts.append("{ " + ", ".join(pairs) + " }")
        else:
            parts.append(_format_toml_scalar(item))
    return "[" + ", ".join(parts) + "]"


def _format_toml_value(value: Any) -> str:
    if isinstance(value, list):
        return _format_toml_list(value)
    return _format_toml_scalar(value)


def _emit_toml_table(lines: list[str], table: dict[str, Any], prefix: list[str]) -> None:
    scalars: list[tuple[str, Any]] = []
    subtables: list[tuple[str, dict[str, Any]]] = []
    array_tables: list[tuple[str, list[dict[str, Any]]]] = []

    for raw_key, value in table.items():
        if value is None:
            continue
        key = str(raw_key)
        if isinstance(value, dict):
            subtables.append((key, value))
        elif isinstance(value, list) and _is_dict_list(value):
            array_tables.append((key, value))
        else:
            scalars.append((key, value))

    for key, value in scalars:
        lines.append(f"{_format_toml_key(key)} = {_format_toml_value(value)}")

    for key, sub in subtables:
        path = [*prefix, key]
        if lines and lines[-1] != "":
            lines.append("")
        header = ".".join(_format_toml_key(p) for p in path)
        lines.append(f"[{header}]")
        _emit_toml_table(lines, sub, path)

    for key, items in array_tables:
        path = [*prefix, key]
        header = ".".join(_format_toml_key(p) for p in path)
        for item in items:
            if lines and lines[-1] != "":
                lines.append("")
            lines.append(f"[[{header}]]")
            _emit_toml_table(lines, item, path)


def dumps_toml(data: dict[str, Any]) -> str:
    """Serialize a recipe dictionary to deterministic TOML."""
    if not isinstance(data, dict):
        raise ValueError("TOML root must be a dictionary")
    lines: list[str] = []
    _emit_toml_table(lines, data, [])
    return "\n".join(lines) + "\n"


def _rename_alias(mapping: dict[str, Any], alias: str, canonical: str, context: str) -> None:
    if alias in mapping:
        if canonical in mapping:
            raise ValueError(f"Specify either {canonical!r} or {alias!r} in {context}, not both")
        mapping[canonical] = mapping.pop(alias)


def normalize_config(raw: dict[str, Any]) -> dict[str, Any]:
    """Strip documentation keys and normalize unambiguous config aliases."""
    if not isinstance(raw, dict):
        raise ValueError("Config root must be an object/dictionary")
    cfg = copy.deepcopy(raw)
    for key in list(cfg):
        if isinstance(key, str) and key.startswith("_"):
            cfg.pop(key)

    _rename_alias(cfg, "C_over_O", "c_over_o", "config")
    _rename_alias(cfg, "O_over_C", "o_over_c", "config")
    _rename_alias(cfg, "edge", "edge_type", "config")
    _rename_alias(cfg, "diameter", "diameter_nm", "config")
    _rename_alias(cfg, "radius", "radius_nm", "config")
    _rename_alias(cfg, "lx", "lx_nm", "config")
    _rename_alias(cfg, "ly", "ly_nm", "config")

    if isinstance(cfg.get("pbc"), bool):
        cfg["pbc"] = [True, True, False] if cfg["pbc"] else [False, False, False]

    comp = cfg.get("composition")
    if isinstance(comp, dict):
        _rename_alias(comp, "C_over_O", "c_over_o", "composition")
        _rename_alias(comp, "O_over_C", "o_over_c", "composition")
        _rename_alias(comp, "group_preferences", "preferences", "composition")
        _rename_alias(comp, "group_priority", "priority", "composition")
        _rename_alias(comp, "composition_mode", "mode", "composition")

    return cfg


def load_config(path: Path) -> dict[str, Any]:
    """Load a .toml or .json configuration file and normalize clear aliases."""
    path = Path(path)
    suffix = path.suffix.lower()
    text = path.read_text(encoding="utf-8")
    if suffix == ".toml":
        raw = tomllib.loads(text)
    elif suffix == ".json":
        raw = json.loads(text)
    else:
        try:
            raw = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Unsupported config extension {suffix!r}; use .toml or .json"
            ) from exc
    if not isinstance(raw, dict):
        raise ValueError(f"Expected a top-level object in config file: {path}")
    return normalize_config(raw)


def dump_config(recipe: dict[str, Any], path: Path) -> None:
    """Save a recipe dictionary as clean .toml or .json without overwriting."""
    if not isinstance(recipe, dict):
        raise ValueError("Recipe must be a dictionary")
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".toml":
        content = dumps_toml(recipe)
    elif suffix == ".json":
        content = json.dumps(recipe, indent=2, allow_nan=False) + "\n"
    else:
        raise ValueError(f"Unsupported config extension {suffix!r}; use .toml or .json")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        handle.write(content)


def render_annotated_toml_template() -> str:
    """Return a complete, annotated TOML configuration template."""
    return """# GO Builder Configuration Template (TOML)
# All units: dimensions in nm (lx_nm, ly_nm, diameter_nm, radius_nm);
# bond_length, vacuum, probe_radius_angstrom, and island radius/centers in Angstrom (A).

name = "go_sheet"                 # Structure label (default: "GO")
seed = 42                         # Reproducible RNG seed >= 0 (default: 0)

# --- Sheet Geometry & Dimensions ---
# Choose ONE sizing mode:
#   1) Physical rectangular / enclosing dimensions in nm: lx_nm & ly_nm
#   2) Nanoflake diameter/radius in nm (requires shape = "hexagonal" or "circular"): diameter_nm or radius_nm
#   3) Honeycomb tile counts (integers >= 2, default: nx = 6, ny = 6)
lx_nm = 4.0
ly_nm = 4.0
# nx = 6
# ny = 6
# diameter_nm = 4.0
# radius_nm = 2.0

shape = "rectangular"             # "rectangular" (default) | "hexagonal" | "circular"
# edge_type = "zigzag"            # Only for shape = "hexagonal": "zigzag" (default) | "armchair"
pbc = [false, false, false]       # Periodic flags [x, y, z]; z must be false (default: [false, false, false])
                                  # Note: hexagonal/circular shapes require [false, false, false]
faces = "both"                    # Basal functionalization faces: "both" (default) | "up" | "down"
# bond_length = 1.42              # C-C scaffold bond length in A, range [1.30, 1.55] (default: 1.42)
# vacuum = 40.0                   # Out-of-plane vacuum cell height in A, >= 15 (default: 40.0)

# --- Holes / Pores (optional) ---
# Either a single [holes] profile, or [[holes]] list with center/radius or size/seed_atom.
[holes]
distribution = "mixed"            # "small" (1-4 C) | "medium" (5-12 C) | "large" (13-30 C) | "mixed" | "random"
percent_vacancies = 5.0           # Target vacancy percentage in (0, 50]
# count = 2                       # Alternative or paired target hole count >= 1
# min_size = 2                    # Optional minimum carbons per hole
# max_size = 12                   # Optional maximum carbons per hole
avoid_outer_edge = true           # Prefer interior pores rather than outer-rim notches (default: true)
probe_radius_angstrom = 1.4       # Probe radius for accessible pore area diagnostics (default: 1.4 A)

# --- Composition & Protonation ---
# Choose EITHER [composition] (O/C or C/O target) OR explicit [groups] counts/percentages (never both).
[composition]
o_over_c = 0.35                   # Atomic O/C ratio (scalar or [low, high]); mutually exclusive with c_over_o
# c_over_o = 2.86                 # Alternative atomic C/O ratio (scalar or [low, high])
priority = ["OH", "epoxide", "COOH", "COO"]  # Ordered preference list; or use preferences = { OH = 2, epoxide = 1 }
mode = "balanced"                 # "balanced" (default) | "random"
ph = 7.4                          # Partition carboxylates via Henderson-Hasselbalch (or use deprotonation_fraction = 0.9)
pka = 4.5                         # Effective carboxyl pKa (default: 4.5)
# enabled_groups = ["OH", "epoxide", "carbonyl", "COOH", "COO"]

# --- Alternative: Explicit [groups] specification (uncomment and remove [composition]) ---
# Allowed groups: OH, epoxide, carbonyl, COOH, COO, epoxide_eph, epoxide_ept
# Each entry may be an integer count, or an object with count / percent / fraction / id:
# [groups]
# OH = 10
# epoxide = { percent = 10.0, basis = "eligible" }   # basis: "eligible" (default) | "carbons"
# COOH = { fraction = 0.1, basis = "eligible" }
# carbonyl = { id = [1, 5] }                         # Explicit scaffold carbon IDs (pairs [[a,b]] for epoxides)

# --- Spatial Placement ---
[placement]
mode = "random"                   # "random" (default) | "island" | "irregular_island" | "dispersed"
# Additional options when mode = "island" or "irregular_island":
# style = "compact"               # "compact" (default for "island") | "irregular" (branched Lerf-Klinowski)
# branching = 0.65                # Branching intensity in [0.0, 1.0] for irregular islands
# n_islands = 2                   # Number of island nuclei (default: 1)
# compactness = 0.85              # Clustering tightness in [0.0, 1.0] (default: 0.85)
# edge_mode = "random"            # "random" (default) | "near_islands"
# radius = 12.0                   # Optional max disk radius in A around each island center
# centers = [[10.0, 10.0], [25.0, 25.0]]  # Explicit [x, y] island centers in A (mutually exclusive with seed_atoms)
# seed_atoms = [25, 80]           # Explicit scaffold carbon IDs for island nuclei

# --- Closed-Shell Clar Electronic Filter ---
# Enforces bipartite sp2 sublattice balance (|N_A - N_B| = 0) and a perfect Kekule matching
# (closed-shell aromatic domains). Default ON in CLI; pass --no-clar or set enabled = false for radical tests.
[clar]
enabled = true
# min_huckel_gap = 0.0            # Minimum HOMO-LUMO pi gap in |beta| units (default: 0.0)
# max_attempts = 128
"""


def render_annotated_json_template() -> str:
    """Return a complete, valid JSON configuration template with inline option guide."""
    template = {
        "_documentation": {
            "units": "Sheet dimensions in nm (lx_nm, ly_nm, diameter_nm, radius_nm); bond_length, vacuum, probe_radius_angstrom, island radius/centers in Angstrom.",
            "sizing_modes": {
                "physical_nm": "lx_nm and ly_nm (positive floats)",
                "nanoflake_diameter": "diameter_nm or radius_nm (requires shape='hexagonal' or 'circular')",
                "tile_counts": "nx and ny integers >= 2 (default: nx=6, ny=6)",
            },
            "shape": {"allowed": ["rectangular", "hexagonal", "circular"], "default": "rectangular"},
            "edge_type": {"allowed": ["zigzag", "armchair"], "default": "zigzag (hexagonal only)"},
            "pbc": {"allowed": "[bool, bool, false]", "default": [False, False, False]},
            "faces": {"allowed": ["both", "up", "down"], "default": "both"},
            "seed": {"type": "integer >= 0", "default": 0},
            "bond_length": {"range_angstrom": [1.30, 1.55], "default": 1.42},
            "vacuum": {"min_angstrom": 15.0, "default": 40.0},
            "holes": {
                "distribution": ["small", "medium", "large", "mixed", "random"],
                "fields": ["percent_vacancies (0..50]", "count >= 1", "min_size", "max_size", "avoid_outer_edge (bool)", "probe_radius_angstrom (default 1.4)"],
            },
            "composition_vs_groups": "Use EITHER 'composition' (o_over_c or c_over_o, priority/preferences, mode='balanced'|'random', ph/deprotonation_fraction, pka=4.5, enabled_groups) OR explicit 'groups' ({group: count | {percent, basis} | {fraction, basis} | {id}}).",
            "placement": {
                "mode": ["random", "island", "irregular_island", "dispersed"],
                "island_fields": [
                    "style ('compact'|'irregular')",
                    "branching (0..1, irregular only, default 0.65)",
                    "n_islands",
                    "compactness (0..1, default 0.85)",
                    "edge_mode ('random'|'near_islands')",
                    "radius (A)",
                    "centers ([[x,y],...])",
                    "seed_atoms ([id,...])",
                ],
            },
            "clar": {
                "default": "Enabled by default in CLI and templates (disable via --no-clar or enabled=false)",
                "format": "boolean or {enabled, min_huckel_gap (>=0, default 0.0), max_attempts (>=1, default 128)}",
            },
        },
        "name": "go_sheet",
        "lx_nm": 4.0,
        "ly_nm": 4.0,
        "shape": "rectangular",
        "pbc": [False, False, False],
        "faces": "both",
        "seed": 42,
        "holes": {
            "distribution": "mixed",
            "percent_vacancies": 5.0,
            "avoid_outer_edge": True,
            "probe_radius_angstrom": 1.4,
        },
        "composition": {
            "o_over_c": 0.35,
            "priority": ["OH", "epoxide", "COOH", "COO"],
            "mode": "balanced",
            "ph": 7.4,
            "pka": 4.5,
        },
        "placement": {
            "mode": "random",
        },
        "clar": {
            "enabled": True,
        },
    }
    return json.dumps(template, indent=2, allow_nan=False) + "\n"
