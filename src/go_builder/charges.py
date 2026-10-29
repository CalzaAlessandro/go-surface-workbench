"""Traceable descriptive charge audits; source tables never establish accuracy.

OFF atom tables are read independently from unrelated box metadata. This is an
audit parser, not a replacement for LEaP's full library validation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shlex
import statistics
from collections import defaultdict
from pathlib import Path

REFERENCE = Path(__file__).resolve().parents[2] / "resources/reference"
ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I"}
MOTIFS = {
    "HYD": "hydroxyl",
    "EPX": "epoxide",
    "EPH": "epoxide",
    "EPT": "epoxide",
    "CAB": "carbonyl",
    "CBX": "carboxylic_acid",
    "CBD": "carboxylate",
    "CARB": "carboxylic_acid",
    "ETH": "ether",
    "GRA": "graphene",
    "GGG": "graphene",
    "GRAT": "H_terminated_edge",
    "GRT": "H_terminated_edge",
}


def _finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Nonfinite charge in source")
    return value


def classify_connectivity(atoms, bonds) -> dict:
    """Connectivity-only motifs, with ambiguity explicit when bond orders lack.

    Amber topologies do not establish bond multiplicity/protonation for a
    monovalent oxygen. Such environments retain an ambiguous label rather than
    being silently called carbonyls or carboxylates. Motifs can overlap.
    """
    by_id = {a["id"]: a for a in atoms}
    adjacent = {i: set() for i in by_id}
    for a, b in bonds:
        adjacent[a].add(b)
        adjacent[b].add(a)
    marks = {i: set() for i in by_id}
    for i, atom in by_id.items():
        if atom["element"] != "O":
            continue
        ns = adjacent[i]
        cs = {j for j in ns if by_id[j]["element"] == "C"}
        hs = {j for j in ns if by_id[j]["element"] == "H"}
        members = {i} | ns
        if len(cs) == 2 and len(ns) == 2:
            a, b = sorted(cs)
            motif = "epoxide" if b in adjacent[a] else "ether"
        elif len(cs) == 1:
            c = next(iter(cs))
            oxygens = {j for j in adjacent[c] if by_id[j]["element"] == "O"}
            if len(oxygens) == 2:
                has_proton = any(by_id[k]["element"] == "H" for j in oxygens for k in adjacent[j])
                motif = "carboxylic_acid" if has_proton else "carboxylate_or_unprotonated_CO2"
                members |= oxygens | {
                    k for j in oxygens for k in adjacent[j] if by_id[k]["element"] == "H"
                }
            elif len(hs) == 1 and len(ns) == 2:
                motif = "hydroxyl"
            elif len(ns) == 1:
                motif = "carbonyl_or_deprotonated_oxygen"
            else:
                motif = "other_oxygen_environment"
        else:
            motif = "other_oxygen_environment"
        for j in members:
            marks[j].add(motif)
    result = {}
    for i, atom in by_id.items():
        if marks[i]:
            result[i] = "+".join(sorted(marks[i]))
        elif atom["element"] == "C":
            result[i] = (
                "H_terminated_carbon"
                if any(by_id[j]["element"] == "H" for j in adjacent[i])
                else "carbon_scaffold"
            )
        elif atom["element"] == "H" and any(by_id[j]["element"] == "C" for j in adjacent[i]):
            result[i] = "H_terminated_carbon"
        else:
            result[i] = "unassigned"
    return result


def read_off_atoms(path) -> list[dict]:
    """Extract only actual atoms and intratemplate connectivity, with 1-based IDs.

    Atom identity is (unit, sequential atom index), never atom name: several
    historical templates contain duplicated names. Boundbox/perturbation/stray
    metadata are not used as atom definitions.
    """
    tables = defaultdict(list)
    section = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("!"):
            match = re.match(r"!entry\.([^.]+)\.unit\.(atoms|connectivity)\s+table\b", line)
            section = match.groups() if match else None
        elif section and line.strip():
            tables[section].append(shlex.split(line))
    rows = []
    for (unit, field), lines in tables.items():
        if field != "atoms":
            continue
        atoms = []
        for i, fields in enumerate(lines, 1):
            if len(fields) != 8:
                raise ValueError(f"{path}: invalid OFF atom row {i}")
            atomic_number = int(fields[6])
            atoms.append(
                {
                    "id": i,
                    "resname": unit,
                    "atom_name": fields[0],
                    "atom_type": fields[1],
                    "element": ELEMENTS.get(atomic_number, f"Z{atomic_number}"),
                    "charge": _finite(fields[7]),
                    "neighbors": [],
                    "motif": MOTIFS.get(unit, "unassigned"),
                    "context_scope": "intratemplate_neighbors_only",
                }
            )
        bonds = []
        for fields in tables.get((unit, "connectivity"), []):
            a, b = int(fields[0]) - 1, int(fields[1]) - 1
            if not 0 <= a < len(atoms) or not 0 <= b < len(atoms):
                raise ValueError(f"{path}: OFF bond outside atom table")
            atoms[a]["neighbors"].append(atoms[b]["element"])
            atoms[b]["neighbors"].append(atoms[a]["element"])
            bonds.append((a + 1, b + 1))
        motifs = classify_connectivity(atoms, bonds)
        for atom in atoms:
            atom["neighbors"].sort()
            atom["residue_motif_hint"] = atom["motif"]
            atom["motif"] = motifs[atom["id"]]
        rows.extend(atoms)
    if not rows:
        raise ValueError(f"{path}: no OFF atom table")
    return rows


def read_hiergo_charges(path) -> list[dict]:
    """Read HierGO Oplsatom2 preserving residue+atom-name charge context."""
    result = []
    for i, line in enumerate(Path(path).read_text().splitlines(), 1):
        body, _, description = line.partition("###")
        if not body.strip():
            continue
        fields = body.split()
        if len(fields) != 4 or not fields[0].startswith("opls_"):
            raise ValueError(f"{path}:{i}: invalid HierGO charge row")
        atom_type, resname, atom_name, charge = fields
        result.append(
            {
                "id": i,
                "atom_type": atom_type,
                "resname": resname,
                "atom_name": atom_name,
                "charge": _finite(charge),
                "element": atom_name[0],
                "motif": MOTIFS.get(resname, "unassigned"),
                "description": description.strip(),
                "neighbors": [],
                "context_scope": "HierGO_residue_and_atom_typing_rule",
            }
        )
    return result


def _topology_atoms(path, max_resid=None):
    import parmed as pmd

    topology = pmd.load_file(str(path))
    selected = [a for a in topology.atoms if max_resid is None or a.residue.idx < max_resid]
    ids = {a.idx for a in selected}
    result = []
    for atom in selected:
        result.append(
            {
                "id": atom.idx + 1,
                "resid": atom.residue.idx + 1,
                "resname": atom.residue.name,
                "atom_name": atom.name,
                "atom_type": atom.type,
                "element": ELEMENTS.get(atom.atomic_number, f"Z{atom.atomic_number}"),
                "charge": _finite(atom.charge),
                "motif": MOTIFS.get(atom.residue.name, "unassigned"),
                "neighbors": sorted(
                    ELEMENTS.get(a.atomic_number, f"Z{a.atomic_number}")
                    for a in atom.bond_partners
                    if a.idx in ids
                ),
                "context_scope": "compiled_topology_local_neighbors",
            }
        )
    motifs = classify_connectivity(
        result,
        [
            (b.atom1.idx + 1, b.atom2.idx + 1)
            for b in topology.bonds
            if b.atom1.idx in ids and b.atom2.idx in ids
        ],
    )
    for atom in result:
        atom["residue_motif_hint"] = atom["motif"]
        atom["motif"] = motifs[atom["id"]]
    return result


def _stats(values):
    return {
        "n": len(values),
        "sum_charge_e": sum(values),
        "mean_charge_e": statistics.mean(values),
        "std_charge_e": statistics.pstdev(values),
        "min_charge_e": min(values),
        "max_charge_e": max(values),
    }


def _csv(path, rows):
    if not rows:
        path.write_text("")
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(
            {k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()}
            for row in rows
        )


def audit_charge_sources(
    output_dir,
    *,
    current_dir=None,
    old_dir=None,
    opls_dir=None,
    current_topology=None,
    current_go_residue_count=494,
) -> dict:
    """Write atom/motif/source audits without ranking unmatched systems.

    Current and OLD directories default to bundled read-only source snapshots.
    current_topology is optional: select explicitly the first 494 GO residues of
    the supplied dry complex, excluding ligand residue 495. The selection and
    hash are recorded; no other topology is automatically cropped.
    """
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    current_dir = Path(current_dir) if current_dir else REFERENCE / "current"
    old_dir = Path(old_dir) if old_dir else REFERENCE / "old"
    opls_dir = Path(opls_dir) if opls_dir else REFERENCE / "opls"
    rows, provenance, datasets = [], [], {}

    def add(path, family, atoms, scope, selection="all"):
        dataset = f"{family}:{path.name}"
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        provenance.append(
            {"dataset": dataset, "path": str(path.resolve()), "sha256": sha, "selection": selection}
        )
        charges = [a["charge"] for a in atoms]
        if not charges:
            raise ValueError(f"Empty charge dataset: {path}")
        counts = defaultdict(int)
        for atom in atoms:
            counts[atom["element"]] += 1
            atom.update(source=family, dataset=dataset, source_sha256=sha)
            rows.append(atom)
        datasets[dataset] = {
            "scope": scope,
            "n_records": len(atoms),
            "element_counts": dict(counts),
            "total_charge_e": sum(charges) if scope != "contextual_parameter_table" else None,
            "statistics": _stats(charges),
            "selection": selection,
        }

    paths = sorted(current_dir.glob("*.lib"))
    if not paths:
        raise ValueError(f"No OFF libraries in {current_dir}")
    for path in paths:
        if path.name.endswith("_gaff2.lib"):
            family = "CURRENT"
        elif path.name.endswith("_amber.lib"):
            family = "BUNDLED_LEGACY_AMBER"
        else:
            family = "BUNDLED_AUXILIARY"
        add(path, family, read_off_atoms(path), "fragment_atom_table")
    old_paths = sorted(old_dir.glob("*.parm7"))
    if not old_paths:
        raise ValueError(f"No OLD compiled topologies in {old_dir}")
    for path in old_paths:
        add(path, "OLD", _topology_atoms(path), "compiled_system")
    if current_topology:
        if not isinstance(current_go_residue_count, int) or current_go_residue_count < 1:
            raise ValueError("current_go_residue_count must be a positive integer")
        path = Path(current_topology)
        add(
            path,
            "CURRENT_COMPILED",
            _topology_atoms(path, current_go_residue_count),
            "compiled_GO_selection",
            f"1-based residues 1-{current_go_residue_count}; all later residues excluded",
        )
    add(
        opls_dir / "Oplsatom2",
        "HIERGO",
        read_hiergo_charges(opls_dir / "Oplsatom2"),
        "contextual_parameter_table",
    )
    hiergo_provenance = json.loads((opls_dir / "PROVENANCE.json").read_text())
    for record in hiergo_provenance["files"]:
        relative = Path(record["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Invalid HierGO provenance relative path")
        copied = opls_dir / relative
        if hashlib.sha256(copied.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"HierGO provenance hash mismatch: {relative}")
    grouped = defaultdict(list)
    for row in rows:
        key = (
            row["dataset"],
            row["motif"],
            row["element"],
            row["atom_type"],
            row["atom_name"],
            row["context_scope"],
            ",".join(row["neighbors"]),
        )
        grouped[key].append(row["charge"])
    motifs = [
        dict(
            zip(
                (
                    "dataset",
                    "motif",
                    "element",
                    "atom_type",
                    "atom_name",
                    "context_scope",
                    "neighbors",
                ),
                key,
                strict=True,
            ),
            **_stats(values),
        )
        for key, values in sorted(grouped.items())
    ]
    result = {
        "schema_version": 1,
        "ranking": None,
        "datasets": datasets,
        "current_charge_method": {
            "origin": "User: standard Amber parameterization of separate GO fragments",
            "qm_method": None,
            "basis": None,
            "charge_fit": None,
            "status": "Recover original Gaussian/ESP/RESP logs before attributing a method",
        },
        "provenance": provenance,
        "hiergo_provenance": hiergo_provenance,
        "comparison_status": "Descriptive only: systems and templates are unmatched in atom count, connectivity, protonation, and net charge.",
        "limits": [
            "CURRENT denotes the nine *_gaff2.lib fragment tables. Historical *_amber.lib and auxiliary EHH.lib are separate source families.",
            "Atom names and force-field type labels do not establish chemical equivalence.",
            "Fragment sums and contextual parameter-table sums are not whole-surface charge constraints.",
            "Amber/GAFF2 variants are separate datasets; duplicate templates are not independent samples.",
            "OFF audit ignores unrelated metadata errors; full topology preparation must still validate library semantics.",
            "Residue motif labels are hints; local connectivity is reported separately, including unassigned OLD chemistry.",
            "Charge magnitudes do not measure accuracy. No winner without common QM reference and conformation checks.",
        ],
        "validation_plan": [
            "Select matched finite H-capped graphene motifs: pristine edge, hydroxyl, epoxide, carbonyl, COOH, COO− and opened epoxide adduct.",
            "Use identical atom IDs, coordinates, protonation and integer total charge for every candidate on each fragment.",
            "Recover current fragment RESP provenance; prepare a declared HF/6-31G* MK reference and test additional conformations separately.",
            "Record any inability to assign a source to the matched chemistry; never pad missing charges or silently neutralize.",
            "Compare ESP RMSE and dipoles about one origin using compare_esp; charged-fragment dipoles depend on origin.",
            "Then validate Amber stability, geometry and solvation/binding sensitivity with matched GB radii and sampling.",
        ],
    }
    output_dir.mkdir(parents=True)
    _csv(output_dir / "atoms.csv", rows)
    _csv(output_dir / "motifs.csv", motifs)
    _csv(output_dir / "sources.csv", [dict(dataset=k, **v) for k, v in datasets.items()])
    (output_dir / "audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
