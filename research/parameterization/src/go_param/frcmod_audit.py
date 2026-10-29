"""Read-only Amber override audit; text consistency is not scientific validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SECTIONS = ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON")
COUNTS = {
    "MASS": (1, 1),
    "BOND": (2, 2),
    "ANGLE": (3, 2),
    "DIHE": (4, 4),
    "IMPROPER": (4, 3),
    "NONBON": (1, 2),
}
UNITS = {
    "MASS": ["atomic mass unit", "optional polarizability field"],
    "BOND": ["kcal/mol/angstrom^2", "angstrom"],
    "ANGLE": ["kcal/mol/radian^2", "degree"],
    "DIHE": [
        "dimensionless divisor",
        "kcal/mol (before division)",
        "degree",
        "dimensionless signed continuation periodicity",
    ],
    "IMPROPER": ["kcal/mol", "degree", "dimensionless periodicity"],
    "NONBON": ["Rmin/2 angstrom (not sigma)", "epsilon kcal/mol"],
}


def _number(value):
    return float(value.replace("D", "E").replace("d", "e"))


def _key(section, types):
    key = tuple(types)
    return min(key, key[::-1]) if section in ("BOND", "ANGLE", "DIHE") else key


def parse_frcmod(path):
    """Preserve case, source line, comments and all Fourier components."""
    path = Path(path)
    data = path.read_bytes()
    sections = {name: [] for name in SECTIONS}
    issues = []
    section = None
    for line_number, raw in enumerate(data.decode(errors="replace").splitlines(), 1):
        stripped = raw.strip()
        if stripped in SECTIONS:
            section = stripped
            continue
        if stripped in ("END", "LJEDIT", "HBOND"):
            section = None
            continue
        if section is None or not stripped or stripped.startswith(("#", "!")):
            continue
        ntypes, nvalues = COUNTS[section]
        pattern = (
            r"^\s*([A-Za-z0-9?*+]+(?:\s*-\s*[A-Za-z0-9?*+]+){"
            + str(ntypes - 1)
            + r"})\s+(.+)$"
        )
        match = re.match(pattern, raw)
        if not match:
            issues.append(
                {
                    "section": section,
                    "line_number": line_number,
                    "code": "unparsed_parameter_line",
                    "severity": "error",
                    "raw": raw,
                }
            )
            continue
        types = [x.strip() for x in re.split(r"\s*-\s*", match[1])]
        tokens = match[2].split()
        try:
            values = [_number(x) for x in tokens[:nvalues]]
            if len(values) != nvalues:
                raise ValueError("missing numeric fields")
            consumed = nvalues
            if section == "MASS" and len(tokens) > 1:
                try:
                    values.append(_number(tokens[1]))
                    consumed = 2
                except ValueError:
                    pass
        except ValueError:
            issues.append(
                {
                    "section": section,
                    "line_number": line_number,
                    "code": "invalid_numeric_fields",
                    "severity": "error",
                    "raw": raw,
                }
            )
            continue
        sections[section].append(
            {
                "line_number": line_number,
                "types": types,
                "values": values,
                "comment": " ".join(tokens[consumed:]),
                "raw": raw,
            }
        )
    return {
        "file": str(path),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "sections": sections,
        "issues": issues,
        "units": UNITS,
    }


def audit_file(path):
    """Separate invalid numerics, ambiguity, and scientific-context review flags."""
    result = parse_frcmod(path)
    issues = result["issues"]
    grouped = defaultdict(list)
    series = []
    active = None

    def flag(entry, section, code, severity="review"):
        issues.append(
            {
                "section": section,
                "line_number": entry["line_number"],
                "types": entry["types"],
                "values": entry["values"],
                "comment": entry["comment"],
                "code": code,
                "severity": severity,
            }
        )

    for section, entries in result["sections"].items():
        for entry in entries:
            values, types = entry["values"], entry["types"]
            finite = all(math.isfinite(v) for v in values)
            if not finite:
                flag(entry, section, "nonfinite_numeric_value", "error")
            comment = entry["comment"]
            if re.search(
                r"ATTN|ATTENTION|need\w*\s+revision|placeholder|missing parameter",
                comment,
                re.IGNORECASE,
            ):
                flag(entry, section, "placeholder_comment")
            if re.search(
                r"same as|calculated|empirical|default value|general improper",
                comment,
                re.IGNORECASE,
            ):
                flag(entry, section, "analogy_or_default_comment")
            penalty = re.search(
                r"penalty score\s*=\s*([0-9.]+)", comment, re.IGNORECASE
            )
            if penalty and float(penalty[1]) >= 50:
                flag(entry, section, "high_analogy_penalty_not_accuracy_measure")
            if any(t in ("?", "DU", "du") for t in types):
                flag(entry, section, "placeholder_or_dummy_type_review")
            if section == "MASS" and finite and values[0] <= 0:
                flag(entry, section, "nonpositive_mass", "error")
            if section in ("BOND", "ANGLE") and finite:
                if values[0] < 0:
                    flag(entry, section, "negative_harmonic_force_constant", "error")
                elif values[0] == 0:
                    flag(
                        entry, section, "zero_harmonic_force_constant_context_required"
                    )
                if section == "BOND" and values[1] <= 0:
                    flag(entry, section, "nonpositive_bond_length", "error")
                if section == "ANGLE" and not 0 < values[1] <= 180:
                    flag(entry, section, "invalid_equilibrium_angle", "error")
            if section == "NONBON" and finite:
                if any(v < 0 for v in values):
                    flag(entry, section, "negative_lj_parameter", "error")
                elif 0 in values:
                    code = (
                        "contextual_zero_hydroxyl_hydrogen_lj"
                        if types[0].lower() == "ho" and values == [0.0, 0.0]
                        else "zero_lj_parameter_context_required"
                    )
                    flag(
                        entry,
                        section,
                        code,
                        "info" if code.startswith("contextual") else "review",
                    )
            if section in ("DIHE", "IMPROPER") and finite:
                if section == "DIHE" and values[0] <= 0:
                    flag(entry, section, "invalid_dihedral_divisor", "error")
                amplitude = values[1] if section == "DIHE" else values[0]
                if amplitude == 0:
                    flag(
                        entry,
                        section,
                        "zero_fourier_amplitude_context_required",
                        "info",
                    )
                if values[-1] == 0 or abs(values[-1]) != int(abs(values[-1])):
                    flag(entry, section, "invalid_or_noninteger_periodicity", "error")
            key = _key(section, types)
            if section != "DIHE":
                grouped[(section, key)].append(
                    {"line_numbers": [entry["line_number"]], "values": [values]}
                )
                continue
            if active is not None and active["key"] != key:
                flag(
                    active["entries"][-1],
                    "DIHE",
                    "broken_dihedral_continuation",
                    "error",
                )
                active = None
            if active is None:
                active = {"key": key, "entries": []}
                series.append(active)
            active["entries"].append(entry)
            if values[-1] >= 0 or not finite:
                active = None
    if active is not None:
        flag(
            active["entries"][-1], "DIHE", "unterminated_dihedral_continuation", "error"
        )
    for chain in series:
        # Sign is continuation bookkeeping; barrier coefficients remain unchanged.
        normalized = [
            e["values"][:-1] + [abs(e["values"][-1])] for e in chain["entries"]
        ]
        grouped[("DIHE", chain["key"])].append(
            {
                "line_numbers": [e["line_number"] for e in chain["entries"]],
                "values": normalized,
            }
        )
    duplicates = []
    for (section, key), definitions in grouped.items():
        if len(definitions) < 2:
            continue
        conflicting = any(
            d["values"] != definitions[0]["values"] for d in definitions[1:]
        )
        duplicates.append(
            {
                "section": section,
                "types": list(key),
                "conflicting": conflicting,
                "definitions": definitions,
            }
        )
        issues.append(
            {
                "section": section,
                "types": list(key),
                "line_number": definitions[-1]["line_numbers"][0],
                "code": "conflicting_definition"
                if conflicting
                else "identical_duplicate_definition",
                "severity": "review",
                "definitions": definitions,
            }
        )
    result["duplicate_groups"] = duplicates
    result["dihedral_series"] = [
        {
            "types": list(s["key"]),
            "component_count": len(s["entries"]),
            "line_numbers": [e["line_number"] for e in s["entries"]],
        }
        for s in series
    ]
    result["summary"] = {
        "section_counts": {s: len(e) for s, e in result["sections"].items()},
        "issue_counts": dict(Counter(i["severity"] for i in issues)),
        "zero_components_are_not_automatically_invalid": True,
    }
    return result


def audit_library_types(directory, override, *, suffix):
    """Read actual OFF atom tables; missing override is not a missing parameter."""
    libraries, issues, all_types = [], [], set()
    for path in sorted(Path(directory).glob("*" + suffix)):
        rows, current_unit = [], None
        raw = path.read_bytes()
        indexed_units = re.findall(
            r'^\s*"([^"]+)"\s*$',
            raw.decode().split("!entry.", 1)[0],
            flags=re.MULTILINE,
        )
        for line_number, line in enumerate(raw.decode().splitlines(), 1):
            match = re.match(r"!entry\.(\S+)\.unit\.atoms table", line)
            if match:
                current_unit = match[1]
                continue
            if line.startswith("!"):
                current_unit = None
            if current_unit is None:
                continue
            atom = re.match(r'\s*"([^" ]+)"\s+"([^" ]+)"', line)
            if atom:
                rows.append(
                    {
                        "unit": current_unit,
                        "atom_name": atom[1],
                        "atom_type": atom[2],
                        "line_number": line_number,
                    }
                )
                all_types.add(atom[2])
        libraries.append(
            {
                "file": str(path),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "atoms": rows,
                "indexed_units": indexed_units,
                "orphan_unit_blocks": sorted(
                    {a["unit"] for a in rows} - set(indexed_units)
                ),
                "runtime_load_behavior": "not_tested_do_not_assume_orphan_blocks_are_loaded",
            }
        )
        expected_lower = suffix == "_gaff2.lib"
        if any(
            (a["atom_type"].isupper() if expected_lower else a["atom_type"].islower())
            for a in rows
        ):
            issues.append(
                {
                    "file": str(path),
                    "code": "mixed_type_namespace",
                    "severity": "review",
                    "atoms": rows,
                }
            )
        if len({a["unit"] for a in rows}) > 1:
            issues.append(
                {
                    "file": str(path),
                    "code": "multiple_units_in_library_load_order_review",
                    "severity": "review",
                    "units": sorted({a["unit"] for a in rows}),
                }
            )
    mass = {e["types"][0] for e in override["sections"]["MASS"]}
    lj = {e["types"][0] for e in override["sections"]["NONBON"]}
    return {
        "libraries": libraries,
        "issues": issues,
        "atom_types_used": sorted(all_types),
        "atom_type_scope": "all raw unit tables including unindexed orphan blocks; runtime-loaded types not established",
        "override_coverage": {
            "types_without_mass_override": sorted(all_types - mass),
            "types_without_lj_override": sorted(all_types - lj),
            "base_force_field_completeness": "not_assessed",
            "interpretation": "Absent override can be inherited from the loaded base force field; full graph/path matching requires a built topology.",
        },
    }


def _safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: _safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe(v) for v in value]
    return value


def _compare_installed_bases(files, root):
    source = root / "references/system_forcefield_sources"
    if not (source / "gaff2.dat").is_file():
        return {
            "status": "not_available",
            "reason": "No pinned installed DAT sources supplied.",
        }
    from .ff_candidates import parse_base_dat

    bases = {
        name: parse_base_dat(source / name)
        for name in ("gaff.dat", "gaff2.dat", "parm99.dat")
    }
    comparisons = {}
    for override in files:
        filename = Path(override["file"]).name
        base_name = (
            "gaff2.dat"
            if filename == "GO_gaff2.frcmod"
            else "parm99.dat"
            if filename in ("GO_amber.frcmod", "GO.frcmod")
            else "gaff.dat"
        )
        base = bases[base_name]
        terms = []
        for section, rows in override["sections"].items():
            for row in rows:
                matches = [
                    e
                    for e in base["sections"][section]
                    if _key(section, e["types"]) == _key(section, row["types"])
                ]
                if section == "DIHE":
                    matches = [
                        e
                        for e in matches
                        if abs(e["values"][-1]) == abs(row["values"][-1])
                        and e["values"][2] % 360 == row["values"][2] % 360
                    ]
                equal = any(
                    all(abs(a - b) <= 1e-8 for a, b in zip(row["values"], e["values"]))
                    for e in matches
                )
                terms.append(
                    {
                        "section": section,
                        "types": row["types"],
                        "line_number": row["line_number"],
                        "override_values": row["values"],
                        "base_candidates": matches,
                        "status": "same_numeric_values"
                        if equal
                        else "different_numeric_values"
                        if matches
                        else "no_exact_key_match_not_a_missing_parameter_claim",
                    }
                )
        comparisons[filename] = {
            "base_file": base_name,
            "base_header": base["header"],
            "base_sha256": base["sha256"],
            "terms": terms,
            "counts": dict(Counter(t["status"] for t in terms)),
        }
    return {
        "status": "available_versioned_direct_type_comparison",
        "comparisons": comparisons,
        "limitations": [
            "Exact/reversed type-key comparisons only; wildcard fallback, aliases, load precedence and full topology coverage remain separate.",
            "OLD/current-old are compared to GAFF1.81 for provenance inspection, not assumed generated by that exact historical version.",
            "Base MASS comparison uses mass only; optional polarizability is not compared.",
        ],
    }


def generate_report(root=ROOT, *, output_dir=None):
    root = Path(root)
    directory = root / "references/local_snapshot/current"
    paths = [
        directory / name
        for name in ("GO_gaff2.frcmod", "GO_amber.frcmod", "GO_gaff2_old.frcmod")
    ]
    paths += sorted((root / "references/local_snapshot/old").glob("*.frcmod"))
    paths += [root / "references/gopy/GO_tutorial/FRCMOD/GO.frcmod"]
    files = [audit_file(p) for p in paths if p.is_file()]
    library_audits = {}
    for suffix, name in (
        ("_gaff2.lib", "GO_gaff2.frcmod"),
        ("_amber.lib", "GO_amber.frcmod"),
    ):
        override = next(f for f in files if Path(f["file"]).name == name)
        library_audits[suffix] = audit_library_types(directory, override, suffix=suffix)
    # Compare explicit override keys without treating different source contexts as competing validated models.
    current = next(f for f in files if Path(f["file"]).name == "GO_gaff2.frcmod")
    previous = next(f for f in files if Path(f["file"]).name == "GO_gaff2_old.frcmod")
    differences = []
    for section in ("MASS", "BOND", "ANGLE", "NONBON"):
        old = {_key(section, e["types"]): e for e in previous["sections"][section]}
        for entry in current["sections"][section]:
            prior = old.get(_key(section, entry["types"]))
            if prior and prior["values"] != entry["values"]:
                differences.append(
                    {
                        "section": section,
                        "types": entry["types"],
                        "current": entry,
                        "previous": prior,
                        "delta_current_minus_previous": [
                            a - b for a, b in zip(entry["values"], prior["values"])
                        ],
                    }
                )
    payload = {
        "scientific_status": "text_audit_not_force_field_validation",
        "files": files,
        "library_audits": library_audits,
        "current_vs_previous_override_differences": differences,
        "base_comparison": _compare_installed_bases(files, root),
        "historical_parameter_provenance": {
            "status": "user_reported_not_line_by_line_independently_verified",
            "scope": "current and OLD GAFF-related frcmod sets",
            "method": "parmchk generation followed by manual addition of missing terms using analogous parameters inspected in Amber DAT files",
            "interpretation": "Unmatched values are not presumed fabricated or QM-refitted; exact base version, analogy choice and graph context require validation.",
        },
        "limits": [
            "No original file is modified.",
            "No force-field winner follows from text inspection.",
            "Improper order is retained, not naively reverse-canonicalized.",
            "Conflicting definitions need actual LEaP load order and resolved topology to determine effective values.",
        ],
    }
    output = Path(output_dir) if output_dir is not None else root / "reports"
    output.mkdir(parents=True, exist_ok=True)
    (output / "frcmod_audit.json").write_text(
        json.dumps(_safe(payload), indent=2, allow_nan=False) + "\n"
    )
    with (output / "frcmod_audit_issues.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "file",
                "sha256",
                "line_number",
                "section",
                "types",
                "values",
                "severity",
                "code",
                "comment",
            ],
        )
        writer.writeheader()
        for file in files:
            for issue in file["issues"]:
                writer.writerow(
                    {
                        "file": file["file"],
                        "sha256": file["sha256"],
                        **{
                            k: issue.get(k, "")
                            for k in (
                                "line_number",
                                "section",
                                "severity",
                                "code",
                                "comment",
                            )
                        },
                        "types": "-".join(issue.get("types", [])),
                        "values": json.dumps(
                            _safe(issue.get("values", issue.get("definitions", [])))
                        ),
                    }
                )
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    generate_report(args.root, output_dir=args.output_dir)
