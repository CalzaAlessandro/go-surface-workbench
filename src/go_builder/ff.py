"""Read-only current GAFF2 fragment assignment, with explicit provenance.

These charges reproduce the supplied fragments. They are not a validation of
whole-sheet electrostatics or of a new covalent junction.
"""

from __future__ import annotations

import hashlib
import math
import re
import shlex
from copy import deepcopy
from pathlib import Path

FRAGMENTS = ("GGG", "GRT", "HYD", "EPX", "EPH", "EPT", "CBX", "CBD", "CAB")
DEFAULT_REFERENCE = Path(__file__).resolve().parents[2] / "resources/reference/current"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_fragment_library(path, residue=None):
    """Parse one indexed OFF unit without accepting appended unindexed units."""
    path = Path(path)
    text = path.read_text()
    indexed = []
    for line in text.splitlines()[1:]:
        if line.startswith("!"):
            break
        indexed.extend(shlex.split(line))
    residue = residue or (indexed[0] if len(indexed) == 1 else None)
    if residue is None or residue not in indexed:
        raise ValueError(f"{path}: requested fragment {residue!r} is not uniquely indexed")
    sections = {}
    all_units = set()
    active = None
    for line in text.splitlines():
        match = re.match(r"!entry\.([^.]+)\.unit\.([^ ]+)", line)
        if match:
            unit, section = match.groups()
            all_units.add(unit)
            active = section if unit == residue else None
            if active:
                if active in sections:
                    raise ValueError(f"{path}: duplicate section {active}")
                sections[active] = {"header": line, "rows": []}
        elif line.startswith("!"):
            active = None
        elif active and line.strip():
            sections[active]["rows"].append(shlex.split(line))
    atoms = []
    for row in sections.get("atoms", {}).get("rows", []):
        if len(row) != 8:
            raise ValueError(f"{path}: malformed atom row")
        atoms.append(
            {
                "name": row[0],
                "atom_type": row[1],
                "atomic_number": int(row[6]),
                "charge": float(row[7]),
            }
        )
    if not atoms or len({a["name"] for a in atoms}) != len(atoms):
        raise ValueError(f"{path}: missing atoms or duplicate atom names in {residue}")
    if any(not math.isfinite(a["charge"]) for a in atoms):
        raise ValueError(f"{path}: nonfinite charge")
    bonds = [tuple(map(int, row)) for row in sections.get("connectivity", {}).get("rows", [])]
    if any(not (1 <= a <= len(atoms) and 1 <= b <= len(atoms) and a != b) for a, b, _ in bonds):
        raise ValueError(f"{path}: invalid connectivity indices")
    warnings = []
    if all_units - set(indexed):
        warnings.append(
            "unindexed appended unit(s) excluded: " + ",".join(sorted(all_units - set(indexed)))
        )
    perturb = sections.get("atomspertinfo", {}).get("rows", [])
    if len(perturb) != len(atoms) or any(
        r[0] != a["name"] or r[1] != a["atom_type"] or float(r[-1]) != a["charge"]
        for r, a in zip(perturb, atoms, strict=False)
    ):
        warnings.append(
            "stale perturbation atom name/type/charge table normalized from primary atoms"
        )
    child = sections.get("childsequence", {}).get("rows", [])
    if child != [[str(2)]]:
        warnings.append("unit childsequence normalized to 2 for one-residue unit")
    if residue in ("EPX", "EPH", "EPT") and any((flags & 15) == 2 for _, _, flags in bonds):
        warnings.append(
            "epoxide source contains a double-bond flag; explicit input graph bond orders take precedence"
        )
    for box_field in ("boundbox", "solventcap"):
        rows = sections.get(box_field, {}).get("rows", [])
        if len(rows) == 4 and float(rows[0][0]) < 0:
            warnings.append(
                f"{box_field}: missing unused fifth field padded with zero in normalized copy"
            )
        elif len(rows) != 5:
            raise ValueError(f"{path}: malformed {box_field} cannot be safely normalized")
    return {
        "residue": residue,
        "atoms": atoms,
        "bonds": bonds,
        "warnings": warnings,
        "source": str(path.resolve()),
        "source_sha256": _sha(path),
        "_sections": sections,
    }


def normalize_libraries(reference_dir, output_dir):
    """Write audited, scoped copies only. Never mutate supplied source files."""
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    libraries = {
        name: load_fragment_library(Path(reference_dir) / f"{name}_gaff2.lib", name)
        for name in FRAGMENTS
    }
    output_dir.mkdir(parents=True)
    report = {}
    for name, lib in libraries.items():
        lines = ["!!index array str", f' "{name}"']
        for section, data in lib["_sections"].items():
            lines.append(data["header"])
            rows = data["rows"]
            if section == "atomspertinfo":
                rows = [
                    [a["name"], a["atom_type"], "0", "-1", str(a["charge"])] for a in lib["atoms"]
                ]
            elif section in ("boundbox", "solventcap") and len(rows) == 4:
                rows = rows + [["0.0"]]
            elif section == "childsequence":
                rows = [["2"]]
            elif section == "name":
                rows = [[name]]
            elif section == "residues":
                rows = [[name, "1", str(len(lib["atoms"]) + 1), "1", "?", "0"]]
            # OFF columns named 'str' are quoted; every other column is numeric.
            types = re.findall(r"\b(str|int|dbl)\s+\w+", data["header"])
            if " single str" in data["header"]:
                types = ["str"]
            for row in rows:
                lines.append(
                    " "
                    + " ".join(
                        f'"{v}"' if i < len(types) and types[i] == "str" else str(v)
                        for i, v in enumerate(row)
                    )
                )
        target = output_dir / f"{name}_gaff2.lib"
        target.write_text("\n".join(lines) + "\n")
        report[name] = {
            k: v for k, v in lib.items() if not k.startswith("_") and k not in ("atoms", "bonds")
        }
        report[name].update(normalized_path=str(target.resolve()), normalized_sha256=_sha(target))
    return report


def assign_current_fragments(model, reference_dir=None, *, parameter_templates=None):
    """Return a copy with exact fragment charges/types, including native anchors.

    Each (resid,resname) must contain a complete fragment and its explicit
    internal bonds. A new junction requires a complete, reviewed model template
    (atom assignments, bond list, provenance, dedicated frcmod), never guessed
    charges. Dedicated templates are accepted by prepare_amber, not inferred here.
    """
    model = deepcopy(model)
    metadata = model.setdefault("metadata", {})
    metadata["parameterized"] = False
    if parameter_templates is not None:
        return _assign_dedicated(model, parameter_templates)
    if metadata.get("requires_parameterization") or any(
        a.get("metadata", {}).get("junction") for a in model["atoms"]
    ):
        raise ValueError(
            "New covalent junction requires dedicated charge/parameterization template and frcmod"
        )
    reference_dir = Path(reference_dir or DEFAULT_REFERENCE)
    grouped = {}
    for atom in model["atoms"]:
        grouped.setdefault((atom["resid"], atom["resname"]), []).append(atom)
    atom_residues = {a["id"]: (a["resid"], a["resname"]) for a in model["atoms"]}
    local_bonds = {key: set() for key in grouped}
    for bond in model["bonds"]:
        key = atom_residues[bond["a"]]
        if key == atom_residues[bond["b"]]:
            local_bonds[key].add(frozenset((bond["a"], bond["b"])))
    libraries = {}
    for (resid, name), atoms in grouped.items():
        if name not in FRAGMENTS:
            raise ValueError(
                f"No current fragment template for {name}; supply a dedicated parameter template"
            )
        if name not in libraries:
            libraries[name] = load_fragment_library(reference_dir / f"{name}_gaff2.lib", name)
        lib = libraries[name]
        by_name = {a["name"]: a for a in atoms}
        if len(by_name) != len(atoms) or set(by_name) != {a["name"] for a in lib["atoms"]}:
            raise ValueError(
                f"Incomplete or ambiguous fragment template {name}: {[a['name'] for a in atoms]}"
            )
        expected = {
            frozenset(
                (
                    by_name[lib["atoms"][a - 1]["name"]]["id"],
                    by_name[lib["atoms"][b - 1]["name"]]["id"],
                )
            )
            for a, b, _ in lib["bonds"]
        }
        present = local_bonds[(resid, name)]
        if expected != present:
            raise ValueError(f"Explicit internal bonds differ from fragment {name}")
        for record in lib["atoms"]:
            atom = by_name[record["name"]]
            symbol = {1: "H", 6: "C", 8: "O"}.get(record["atomic_number"])
            if atom["element"] != symbol:
                raise ValueError(f"Element differs from fragment template {name}.{record['name']}")
            atom.update(atom_type=record["atom_type"], charge=record["charge"])
            atom.setdefault("metadata", {})["charge_source"] = f"{name}_gaff2.lib:{record['name']}"
    metadata["force_field"] = {
        "name": "GAFF2 + supplied GO_gaff2.frcmod",
        "assignment": "current fragment transfer; whole-sheet charges unvalidated",
        "fragments": [
            {"fragment": name, "source_sha256": lib["source_sha256"]}
            for name, lib in libraries.items()
        ],
    }
    return model


def _assign_dedicated(model, template):
    required = ("atoms", "bonds", "frcmod", "provenance")
    if not isinstance(template, dict) or any(not template.get(k) for k in required):
        raise ValueError("Dedicated junction template requires atoms, bonds, frcmod and provenance")
    records = {a["id"]: a for a in template["atoms"]}
    if set(records) != {a["id"] for a in model["atoms"]}:
        raise ValueError("Dedicated junction template must cover the complete model atom IDs")

    def bond_signature(b):
        image = tuple(b.get("image", [0, 0, 0]))
        if b["a"] > b["b"]:
            image = tuple(-x for x in image)
        return (*sorted((b["a"], b["b"])), float(b.get("order", 1)), image)

    expected = {bond_signature(b) for b in template["bonds"]}
    if expected != {bond_signature(b) for b in model["bonds"]}:
        raise ValueError("Dedicated junction template bonds do not match model")
    if not Path(template["frcmod"]).is_file():
        raise ValueError("Dedicated junction frcmod missing")
    for atom in model["atoms"]:
        record = records[atom["id"]]
        if record.get("element") != atom["element"] or not record.get("atom_type"):
            raise ValueError("Dedicated junction template element/type mismatch")
        if not isinstance(record.get("charge"), (float, int)) or not math.isfinite(
            record["charge"]
        ):
            raise ValueError("Dedicated junction template requires finite explicit charges")
        atom.update(atom_type=record["atom_type"], charge=float(record["charge"]))
    model["metadata"]["junction_template"] = deepcopy(template)
    model["metadata"]["requires_parameterization"] = False
    return model
