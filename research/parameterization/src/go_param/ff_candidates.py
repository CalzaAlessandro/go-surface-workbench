"""Generate pinned, unfitted native-base GO controls and explicit graph gaps."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from itertools import combinations, permutations
from pathlib import Path

from .frcmod_audit import COUNTS, SECTIONS, _key

ROOT = Path(__file__).resolve().parents[2]


def parse_base_dat(path):
    """Read the fixed Amber DAT section sequence in the exported pinned sources.

    Preserve original line numbers. Exclude prose after END; preserve LJ aliases.
    Reject unexpected block layout instead of guessing section assignments.
    """
    path = Path(path)
    raw = path.read_bytes()
    blocks, block = [], []
    for line_number, line in enumerate(raw.decode().splitlines(), 1):
        if line.strip():
            block.append((line_number, line))
        elif block:
            blocks.append(block)
            block = []
    if block:
        blocks.append(block)
    if len(blocks) < 7 or not blocks[6][0][1].strip().startswith("MOD4"):
        raise ValueError(f"Unexpected Amber DAT layout: {path}")
    sections = {name: [] for name in SECTIONS}
    for section, index, skip in (
        ("MASS", 0, 1),
        ("BOND", 1, 1),
        ("ANGLE", 2, 0),
        ("DIHE", 3, 0),
        ("IMPROPER", 4, 0),
        ("NONBON", 6, 1),
    ):
        ntypes, nvalues = COUNTS[section]
        pattern = (
            r"^\s*([A-Za-z0-9?*+]+(?:\s*-\s*[A-Za-z0-9?*+]+){"
            + str(ntypes - 1)
            + r"})\s+(.+)$"
        )
        for line_number, line in blocks[index][skip:]:
            match = re.match(pattern, line)
            if not match:
                raise ValueError(f"Cannot parse {section} {path}:{line_number}: {line}")
            tokens = match[2].split()
            values = [float(t.replace("D", "E")) for t in tokens[:nvalues]]
            if len(values) != nvalues:
                raise ValueError(f"Missing fields at {path}:{line_number}")
            sections[section].append(
                {
                    "types": [x.strip() for x in re.split(r"\s*-\s*", match[1])],
                    "values": values,
                    "line_number": line_number,
                    "comment": " ".join(tokens[nvalues:]),
                    "source_file": path.name,
                }
            )
    equivalences = {}
    for _, line in blocks[5]:
        if re.search(r"\d", line):
            continue  # Legacy 10-12 water record, not an LJ-equivalence list.
        names = line.split()
        if names:
            equivalences.update({name: names[0] for name in names[1:]})
    return {
        "file": str(path),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "header": blocks[0][0][1],
        "sections": sections,
        "lj_equivalences": equivalences,
    }


def _role_types(name):
    roles = {
        "scaffold_C": "ca",
        "edge_C": "ca",
        "edge_H": "ha",
        "hydroxyl_C": "c3",
        "hydroxyl_O": "oh",
        "hydroxyl_H": "ho",
        "epoxide_C": "cx",
        "epoxide_O": "op",
        "carbonyl_C": "c",
        "carbonyl_O": "o",
        "COOH_anchor": "ca",
        "COOH_C": "c",
        "COOH_O": "o",
        "COOH_OH": "oh",
        "COOH_H": "ho",
        "COO_anchor": "ca",
        "COO_C": "c",
        "COO_O": "o",
    }
    if name == "gaff2_native":
        roles["hydroxyl_C"] = "c6"
    if name == "amber99sb_go_analogy":
        translation = {
            "ca": "CA",
            "ha": "HA",
            "c3": "CT",
            "oh": "OH",
            "ho": "HO",
            "cx": "GX",
            "op": "GO",
            "c": "C",
            "o": "O",
        }
        roles = {role: translation[atom_type] for role, atom_type in roles.items()}
        roles["COO_O"] = "O2"
    return roles


def _subset(base, selected):
    sections = {
        section: [copy.deepcopy(e) for e in rows if set(e["types"]) <= selected | {"X"}]
        for section, rows in base["sections"].items()
    }
    # Parm99 LJ-equivalent aromatic types are explicitly materialized with provenance.
    present = {e["types"][0] for e in sections["NONBON"]}
    lookup = {e["types"][0]: e for e in base["sections"]["NONBON"]}
    for name in sorted(selected - present):
        parent = base["lj_equivalences"].get(name)
        if parent in lookup:
            entry = copy.deepcopy(lookup[parent])
            entry["types"] = [name]
            entry["equivalent_source_type"] = parent
            sections["NONBON"].append(entry)
    return sections


def _contract(model, roles, sections):
    atoms = [
        {
            "atom_id": a["id"],
            "element": a["element"],
            "role": a["role"],
            "atom_type": roles[a["role"]],
        }
        for a in model["atoms"]
    ]
    types = {a["atom_id"]: a["atom_type"] for a in atoms}
    ns = {i: set() for i in types}
    for b in model["bonds"]:
        ns[b["a"]].add(b["b"])
        ns[b["b"]].add(b["a"])
    queries = {"BOND": {}, "ANGLE": {}, "DIHE": {}, "IMPROPER": {}}
    for b in model["bonds"]:
        ids = (b["a"], b["b"])
        queries["BOND"][_key("BOND", [types[i] for i in ids])] = ids
        for first in ns[ids[0]] - {ids[1]}:
            for last in ns[ids[1]] - {ids[0]}:
                path = (first, *ids, last)
                if len(set(path)) == 4:
                    queries["DIHE"][_key("DIHE", [types[i] for i in path])] = path
    planar_roles = {
        "scaffold_C",
        "edge_C",
        "carbonyl_C",
        "COOH_C",
        "COO_C",
        "COOH_anchor",
        "COO_anchor",
    }
    for atom in atoms:
        center = atom["atom_id"]
        for pair in combinations(sorted(ns[center]), 2):
            ids = (pair[0], center, pair[1])
            queries["ANGLE"][_key("ANGLE", [types[i] for i in ids])] = ids
        if len(ns[center]) == 3 and atom["role"] in planar_roles:
            outer = sorted(ns[center])
            ids = (outer[0], outer[1], center, outer[2])
            queries["IMPROPER"][tuple(types[i] for i in ids)] = ids
    gaps = []
    for bond in model["bonds"]:
        a, b = bond["a"], bond["b"]
        if types[a] in ("cx", "GX") and types[b] == types[a]:
            shared_oxygen = any(
                next(x for x in atoms if x["atom_id"] == i)["element"] == "O"
                for i in ns[a] & ns[b]
            )
            if not shared_oxygen:
                gaps.append(
                    {
                        "model_id": model["model_id"],
                        "section": "BOND",
                        "types": [types[a], types[b]],
                        "example_atom_ids": [a, b],
                        "code": "cx_cx_bond_spans_distinct_epoxide_rings",
                        "status": "chemistry_context_unvalidated_even_if_numeric_term_exists",
                    }
                )
    for center, neighbors in ns.items():
        if types[center] not in ("cx", "GX"):
            continue
        for a, b in combinations(sorted(neighbors), 2):
            if types[a] == types[b] == types[center] and b not in ns[a]:
                gaps.append(
                    {
                        "model_id": model["model_id"],
                        "section": "ANGLE",
                        "types": [types[a], types[center], types[b]],
                        "example_atom_ids": [a, center, b],
                        "code": "three_member_carbon_angle_used_outside_carbon_triangle",
                        "status": "chemistry_context_unvalidated_even_if_numeric_term_exists",
                    }
                )
    for section, paths in queries.items():
        for signature, ids in paths.items():
            variants = [signature]
            if section == "IMPROPER":
                variants = [
                    (p[0], p[1], signature[2], p[2])
                    for p in permutations((signature[0], signature[1], signature[3]))
                ]
            elif section in ("BOND", "ANGLE", "DIHE"):
                variants.append(signature[::-1])
            covered = any(
                all(p == "X" or p == q for p, q in zip(entry["types"], v))
                for entry in sections[section]
                for v in variants
            )
            if not covered:
                gaps.append(
                    {
                        "model_id": model["model_id"],
                        "section": section,
                        "types": list(signature),
                        "example_atom_ids": list(ids),
                        "status": "unresolved_no_numeric_padding",
                    }
                )
    return {
        "atoms": atoms,
        "formal_charge": model["formal_charge"],
        "partial_charges": "not_assigned_by_frcmod_generation",
        "typing_status": "declared_role_contract_requires_native_graph_typing_review",
        "graph_term_gaps": gaps,
        "improper_matching": "coverage screen permits outer permutations with center third; actual LEaP matching still required",
    }


def _write_frcmod(sections, path, title):
    lines = [title]
    for section in SECTIONS:
        lines.append(section)
        for entry in sections[section]:
            types = "-".join(f"{t:<2s}" for t in entry["types"])
            if section == "DIHE":
                values = f"{int(entry['values'][0]):4d}" + "".join(
                    f"{v:15.8f}" for v in entry["values"][1:]
                )
            else:
                values = " ".join(f"{v:.8f}" for v in entry["values"])
            origin = f"{entry['source_file']}:{entry['line_number']}"
            if entry.get("equivalent_source_type"):
                origin += f" LJ-equivalent {entry['equivalent_source_type']}"
            lines.append(f"{types:<14s} {values}    source {origin}")
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


def generate_candidates(root=ROOT, *, output_dir=None):
    """Create three nonempty new controls without claiming a QM parameter refit."""
    root = Path(root)
    source = root / "references/system_forcefield_sources"
    bases = {
        name: parse_base_dat(source / file)
        for name, file in (
            ("gaff", "gaff.dat"),
            ("gaff2", "gaff2.dat"),
            ("amber99", "parm99.dat"),
        )
    }
    # Enforce the exported-source manifest, not filenames alone.
    source_manifest = json.loads((source / "source_manifest.json").read_text())
    expected = {e["copy"]: e["sha256"] for e in source_manifest["files"] if "copy" in e}
    for base in bases.values():
        if base["sha256"] != expected[Path(base["file"]).name]:
            raise ValueError(f"Pinned source hash mismatch: {base['file']}")
    ff99 = source / "frcmod.ff99SB"
    if hashlib.sha256(ff99.read_bytes()).hexdigest() != expected[ff99.name]:
        raise ValueError("Pinned ff99SB source hash mismatch")
    output = (
        Path(output_dir) if output_dir is not None else root / "parameters/candidates"
    )
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "scientific_status": "new_baseline_candidates_not_final_fitted_parameters",
        "candidates": {},
    }
    models = [
        json.loads(p.read_text())
        for p in sorted((root / "models").glob("*/model.json"))
    ]
    for name, base_name in (
        ("gaff_native", "gaff"),
        ("gaff2_native", "gaff2"),
        ("amber99sb_go_analogy", "amber99"),
    ):
        folder = output / name
        folder.mkdir(exist_ok=True)
        roles = _role_types(name)
        base = bases[base_name]
        sections = _subset(base, set(roles.values()) - {"GX", "GO"})
        custom_notes = []
        if name == "amber99sb_go_analogy":
            # New isolated ring namespace; no legacy CY/CX nitrile/sp confusion.
            for section in ("MASS", "NONBON"):
                for target, parent in (("GX", "CT"), ("GO", "OS")):
                    entry = copy.deepcopy(
                        next(
                            e
                            for e in base["sections"][section]
                            if e["types"] == [parent]
                        )
                    )
                    entry["types"] = [target]
                    entry["equivalent_source_type"] = parent
                    sections[section].append(entry)
            projection = {
                "cx": "GX",
                "op": "GO",
                "ca": "CA",
                "c3": "CT",
                "c": "C",
                "o": "O",
                "oh": "OH",
                "ho": "HO",
                "ha": "HA",
                "hc": "HC",
                "h1": "H1",
            }
            for section in ("BOND", "ANGLE", "DIHE", "IMPROPER"):
                for original in bases["gaff"]["sections"][section]:
                    if "cx" not in original["types"] and "op" not in original["types"]:
                        continue
                    if not set(original["types"]) <= set(projection) | {"X"}:
                        continue
                    entry = copy.deepcopy(original)
                    entry["source_types"] = entry["types"][:]
                    entry["types"] = [projection.get(t, t) for t in entry["types"]]
                    sections[section].append(entry)
            (folder / "atom_types.leap").write_text(
                'addAtomTypes {\n { "GX" "C" "sp3" }\n { "GO" "O" "sp3" }\n}\n'
            )
            custom_notes = [
                "GX/GO mass and LJ use Amber CT/OS analogies.",
                "Ring bonded terms are explicit GAFF1.81 cx/op projections; this is a hybrid analogy, not native AMBER99SB GO.",
                "ff99SB backbone correction contains nitrogen paths and does not alter these C/H/O-only selections.",
            ]
        contracts = {m["model_id"]: _contract(m, roles, sections) for m in models}
        gaps = [
            gap
            for contract in contracts.values()
            for gap in contract["graph_term_gaps"]
        ]
        path = folder / "candidate.frcmod"
        _write_frcmod(
            sections,
            path,
            f"{name}: PINNED UNFITTED BASELINE - NOT PRODUCTION VALIDATED",
        )
        dependencies = [
            {
                "file": Path(base["file"]).name,
                "sha256": base["sha256"],
                "header": base["header"],
            }
        ]
        if name == "amber99sb_go_analogy":
            dependencies += [
                {"file": ff99.name, "sha256": expected[ff99.name]},
                {
                    "file": "gaff.dat",
                    "sha256": bases["gaff"]["sha256"],
                    "scope": "ring bonded analogy only",
                },
            ]
        manifest = {
            "candidate": name,
            "status": "unfitted_baseline_not_validated_for_production",
            "base_force_field": base["header"],
            "dependencies": dependencies,
            "role_types": roles,
            "parameter_counts": {s: len(e) for s, e in sections.items()},
            "candidate_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "requires_graph_specific_typing_validation": True,
            "model_contracts": contracts,
            "gaps": gaps,
            "term_provenance": sections,
            "notes": custom_notes
            + [
                "Selected source terms are rematerialized controls, not newly optimized coefficients.",
                "No partial charges assigned. Per-graph parmchk2/LEaP and QM validation required.",
                "Role contracts do not replace native antechamber atom typing or bond-order/aromaticity review.",
                "Missing graph paths are not filled with zero/default coefficients.",
            ],
        }
        (folder / "manifest.json").write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n"
        )
        (folder / "gaps.json").write_text(json.dumps(gaps, indent=2) + "\n")
        report["candidates"][name] = manifest
    (output / "index.json").write_text(
        json.dumps(
            {
                "scientific_status": report["scientific_status"],
                "candidates": {
                    name: {
                        "path": f"{name}/candidate.frcmod",
                        "gaps": len(m["gaps"]),
                        "status": m["status"],
                    }
                    for name, m in report["candidates"].items()
                },
            },
            indent=2,
        )
        + "\n"
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    generate_candidates(args.root, output_dir=args.output_dir)
