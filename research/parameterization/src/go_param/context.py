"""Analyze atomic RESP charge context without averaging functional groups.

The unit of analysis is one model/charge fit. Multiple rigid ESP rotations in a
fit are preserved as provenance and do not become independent samples.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import deque
from itertools import combinations
from pathlib import Path

from .modular import partition_model
from .resp import _validate_model

GRAPH_DISTANCE_DEFINITION = (
    "minimum number of C-C scaffold bonds between group scaffold-carbon anchors"
)
CHARGE_ROUNDING_DIGITS = 6


def _read_json(value):
    if isinstance(value, (str, Path)):
        return json.loads(Path(value).read_text())
    return value


def _charge_data(model, charge_data):
    atoms = model["atoms"]
    rows = charge_data.get("atoms")
    if not isinstance(rows, list) or len(rows) != len(atoms):
        raise ValueError("Charge atom count differs from model")
    if charge_data.get("model_id", model["model_id"]) != model["model_id"]:
        raise ValueError("Charge model_id differs from model")
    if (
        charge_data.get("formal_charge", model["formal_charge"])
        != model["formal_charge"]
    ):
        raise ValueError("Charge-file formal charge differs from model")

    charges = {}
    for index, (atom, row) in enumerate(zip(atoms, rows), 1):
        if row.get("atom_id") != atom["id"] or atom["id"] != index:
            raise ValueError("Charge atom IDs/order differ from model")
        if row.get("element") != atom["element"]:
            raise ValueError(f"Charge element differs from model at atom {index}")
        if "role" in row and row["role"] != atom["role"]:
            raise ValueError(f"Charge role differs from model at atom {index}")
        try:
            q = float(row["charge_e"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid/missing charge at atom {index}") from exc
        if not math.isfinite(q):
            raise ValueError(f"Nonfinite charge at atom {index}")
        charges[atom["id"]] = q

    total = sum(charges.values())
    # Allow worst-case accumulated rounding of six-decimal atom charges, with a
    # one-microelectron floor for files that round the total separately.
    tolerance = max(1e-6, len(atoms) * 0.5 * 10 ** (-CHARGE_ROUNDING_DIGITS) + 1e-12)
    error = total - model["formal_charge"]
    if abs(error) > tolerance:
        raise ValueError(
            f"Atomic charge sum {total:.9f} differs from formal charge "
            f"{model['formal_charge']} by {error:.9f} e (tolerance {tolerance:.9f} e)"
        )
    return charges, {
        "expected_formal_charge_e": model["formal_charge"],
        "observed_charge_sum_e": total,
        "error_e": error,
        "rounding_tolerance_e": tolerance,
        "passed": True,
    }


def _scaffold_carbon_ids(model):
    atoms = {atom["id"]: atom for atom in model["atoms"]}
    metadata = model.get("metadata", {})
    declared = metadata.get("scaffold_carbon_ids")
    if declared is not None:
        ids = list(declared)
        if (
            not ids
            or any(type(i) is not int for i in ids)
            or len(ids) != len(set(ids))
            or any(i not in atoms or atoms[i]["element"] != "C" for i in ids)
        ):
            raise ValueError(
                "metadata.scaffold_carbon_ids must name unique model carbons"
            )
        exocyclic_roles = {"COOH_C", "COO_C"}
        if any(atoms[i]["role"] in exocyclic_roles for i in ids):
            raise ValueError(
                "metadata.scaffold_carbon_ids includes an exocyclic functional-group carbon"
            )
        return sorted(ids), "metadata.scaffold_carbon_ids"

    # Current pyrene models explicitly tag every original scaffold carbon in
    # the roles array, including edge and functionalized anchor carbons.
    tagged = []
    ambiguous = []
    for atom in model["atoms"]:
        if atom["element"] != "C":
            continue
        roles = atom.get("roles")
        if roles is None and atom.get("role") != "scaffold_C":
            ambiguous.append(atom["id"])
        elif "scaffold_C" in (roles if roles is not None else [atom["role"]]):
            tagged.append(atom["id"])
    if ambiguous:
        raise ValueError(
            "Cannot identify scaffold carbons unambiguously from roles for atoms "
            f"{ambiguous}; provide metadata.scaffold_carbon_ids"
        )
    if not tagged:
        raise ValueError(
            "Cannot identify scaffold carbons unambiguously: provide "
            "metadata.scaffold_carbon_ids or explicit scaffold_C roles"
        )
    return sorted(tagged), "explicit scaffold_C role fallback"


def _scaffold_graph(model, scaffold_ids):
    scaffold = set(scaffold_ids)
    graph = {atom_id: set() for atom_id in scaffold_ids}
    for bond in model["bonds"]:
        a, b = bond["a"], bond["b"]
        if a in scaffold and b in scaffold:
            graph[a].add(b)
            graph[b].add(a)
    return graph


def _bfs(graph, start):
    distances = {start: 0}
    previous = {start: None}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for neighbor in sorted(graph[node]):
            if neighbor not in distances:
                distances[neighbor] = distances[node] + 1
                previous[neighbor] = node
                queue.append(neighbor)
    return distances, previous


def _shortest_anchor_path(graph, starts, ends):
    """Shortest scaffold-only path; lengths count C-C edges."""
    ends = set(ends)
    best = None
    for start in sorted(starts):
        distances, previous = _bfs(graph, start)
        for end in sorted(ends.intersection(distances)):
            candidate = (distances[end], start, end, previous)
            if best is None or candidate[:3] < best[:3]:
                best = candidate
    if best is None:
        return None, None, None, None
    distance, start, end, previous = best
    path = []
    node = end
    while node is not None:
        path.append(node)
        if node == start:
            break
        node = previous[node]
    path.reverse()
    return distance, start, end, path


def _group_record(group, group_id, atoms, charges, scaffold_ids, bonds):
    ids = group["atom_ids"]
    carbon_ids = [i for i in ids if atoms[i]["element"] == "C"]
    anchors = sorted(set(ids).intersection(scaffold_ids))
    anchor_source = "group scaffold carbons"
    if not anchors:
        group_atom_ids = set(ids)
        anchors = sorted(
            {
                bond["a"] if bond["b"] in group_atom_ids else bond["b"]
                for bond in bonds
                if (bond["a"] in group_atom_ids and bond["b"] in scaffold_ids)
                or (bond["b"] in group_atom_ids and bond["a"] in scaffold_ids)
            }
        )
        anchor_source = "directly bonded scaffold carbons"
    return {
        "group_id": group_id,
        "block_type": group["block_type"],
        "atom_ids": ids,
        "roles": sorted(set(group["roles"])),
        "net_charge_e": sum(charges[i] for i in ids),
        "target_formal_charge_e": group["target_charge_e"],
        "scaffold_carbon_anchor_ids": anchors,
        "scaffold_anchor_source": anchor_source if anchors else None,
        "scaffold_anchor_charge_e_sum": sum(charges[i] for i in anchors),
        "carbon_atom_ids_for_centroid": carbon_ids,
        "carbon_centroid_angstrom": (
            [
                sum(atoms[i]["position"][axis] for i in carbon_ids) / len(carbon_ids)
                for axis in range(3)
            ]
            if carbon_ids
            else None
        ),
    }


def _named_values(value, key_name):
    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key.lower() == key_name.lower():
                found.append(item)
            found.extend(_named_values(item, key_name))
    elif isinstance(value, list):
        for item in value:
            found.extend(_named_values(item, key_name))
    return found


def _diagnostic_label(value):
    return isinstance(value, str) and any(
        token in value.lower()
        for token in ("diagnostic", "instability_control", "not_final")
    )


def _provenance(item, model, fit_evidence):
    metadata = model.get("metadata", {})
    protocol = fit_evidence.get("protocol", {}) if fit_evidence else {}
    evidence_provenance = fit_evidence.get("provenance", {}) if fit_evidence else {}
    metadata_scopes = [
        value for value in _named_values(metadata, "scope") if isinstance(value, str)
    ]
    evidence_scopes = [
        value
        for value in _named_values(fit_evidence, "scope")
        if isinstance(value, str)
    ]
    all_scope_labels = metadata_scopes + evidence_scopes
    diagnostic_scope = next(
        (value for value in all_scope_labels if _diagnostic_label(value)), None
    )
    status_labels = _named_values(metadata, "geometry_status")
    status_labels += _named_values(metadata, "status")
    status_labels += _named_values(fit_evidence, "status")
    stability_values = _named_values(metadata, "wavefunction_stability_validated")
    if not stability_values:
        stability_values = _named_values(
            fit_evidence, "wavefunction_stability_validated"
        )
    reported_stability = next(
        (value for value in stability_values if type(value) is bool), None
    )
    model_scope = metadata.get("scope")
    protocol_scope = protocol.get("scope")
    fallback_scope = (
        protocol_scope
        or (fit_evidence.get("scope") if fit_evidence else None)
        or model_scope
    )
    effective_scope = diagnostic_scope or fallback_scope
    face_annotation = {
        key: value for key, value in metadata.items() if "face" in key.lower()
    }
    return {
        "model_path": str(item["model_path"]),
        "charges_path": str(item["charges_path"]),
        "fit_evidence_path": (
            str(item["fit_evidence_path"]) if item.get("fit_evidence_path") else None
        ),
        "model_id": model["model_id"],
        "origin": evidence_provenance.get(
            "model_source_path",
            evidence_provenance.get("source_path", metadata.get("origin")),
        ),
        "model_metadata": metadata,
        "face_metadata_annotation_only": face_annotation,
        "fit_status": fit_evidence.get("status") if fit_evidence else None,
        "qm_method": protocol.get(
            "method",
            fit_evidence.get("qm_method", metadata.get("qm_method"))
            if fit_evidence
            else metadata.get("qm_method"),
        ),
        "scope": effective_scope,
        "protocol_scope": protocol_scope,
        "model_scope": model_scope,
        "is_diagnostic": any(_diagnostic_label(label) for label in all_scope_labels)
        or any(_diagnostic_label(label) for label in status_labels),
        "reported_wavefunction_stability_validated": reported_stability,
        "upstream_qm_gate": protocol.get(
            "parent_qm_optimization_frequencies_stability",
            fit_evidence.get("upstream_qm_gate") if fit_evidence else None,
        ),
        "physical_conformers": (
            protocol.get(
                "physical_conformers", fit_evidence.get("physical_conformers", 1)
            )
            if fit_evidence
            else 1
        ),
        "orientations": (
            protocol.get("orientations", fit_evidence.get("orientations"))
            if fit_evidence
            else None
        ),
        "analysis_context_count": 1,
        "context_unit": "one atomic charge fit; rigid rotations are provenance, not independent contexts",
        "fit_status_is_acceptance": False,
        "fit_evidence_metadata_only": fit_evidence,
    }


def _analyze_one(item):
    if not isinstance(item, dict) or not {"model_path", "charges_path"}.issubset(item):
        raise ValueError("Each dataset entry requires model_path and charges_path")
    model = _read_json(item["model_path"])
    _validate_model(model)
    partition = partition_model(model)
    charges_data = _read_json(item["charges_path"])
    charges, charge_check = _charge_data(model, charges_data)
    fit_evidence = (
        _read_json(item["fit_evidence_path"]) if item.get("fit_evidence_path") else None
    )

    atoms = {atom["id"]: atom for atom in model["atoms"]}
    scaffold_ids, scaffold_source = _scaffold_carbon_ids(model)
    scaffold_set = set(scaffold_ids)
    graph = _scaffold_graph(model, scaffold_ids)

    group_rows = []
    for group in partition["groups"]:
        if group["block_type"] == "GRT":
            continue
        group_id = f"{group['block_type']}:{sum(g['block_type'] == group['block_type'] for g in group_rows) + 1}"
        group_rows.append(
            _group_record(group, group_id, atoms, charges, scaffold_set, model["bonds"])
        )

    pairs = []
    for left, right in combinations(group_rows, 2):
        distance, start, end, path = _shortest_anchor_path(
            graph,
            left["scaffold_carbon_anchor_ids"],
            right["scaffold_carbon_anchor_ids"],
        )
        p1, p2 = left["carbon_centroid_angstrom"], right["carbon_centroid_angstrom"]
        centroid_distance = (
            math.dist(p1, p2) if p1 is not None and p2 is not None else None
        )
        pairs.append(
            {
                "group_a_id": left["group_id"],
                "group_a_type": left["block_type"],
                "group_a_net_charge_e": left["net_charge_e"],
                "group_a_scaffold_anchor_charge_e_sum": left[
                    "scaffold_anchor_charge_e_sum"
                ],
                "group_b_id": right["group_id"],
                "group_b_type": right["block_type"],
                "group_b_net_charge_e": right["net_charge_e"],
                "group_b_scaffold_anchor_charge_e_sum": right[
                    "scaffold_anchor_charge_e_sum"
                ],
                "combined_group_net_charge_e": left["net_charge_e"]
                + right["net_charge_e"],
                "graph_distance": distance,
                "graph_distance_definition": GRAPH_DISTANCE_DEFINITION,
                "graph_anchor_atom_ids": [start, end] if path is not None else None,
                "graph_path_atom_ids": path,
                "carbon_centroid_distance_angstrom": centroid_distance,
                "centroid_definition": "Euclidean distance between centroids of all carbon atoms in each partitioned functional-group block; no periodic boundary correction",
            }
        )

    fg_atom_ids = {atom_id for group in group_rows for atom_id in group["atom_ids"]}
    fg_anchors = {
        group["group_id"]: group["scaffold_carbon_anchor_ids"] for group in group_rows
    }
    background_ids = [
        atom_id
        for atom_id in scaffold_ids
        if atom_id not in fg_atom_ids and atoms[atom_id]["element"] == "C"
    ]
    background = []
    for atom_id in background_ids:
        atom = atoms[atom_id]
        roles = atom.get("roles", [atom["role"]])
        background_kind = (
            "hydrogen_terminated_edge"
            if atom["role"] == "edge_C" or "edge_C" in roles
            else "interior_scaffold"
        )
        distances, _ = _bfs(graph, atom_id)
        available = [
            (distances[anchor], group_id)
            for group_id, anchor_ids in fg_anchors.items()
            for anchor in anchor_ids
            if anchor in distances
        ]
        min_distance = min((d for d, _ in available), default=None)
        nearest_groups = sorted(gid for d, gid in available if d == min_distance)
        background.append(
            {
                "atom_id": atom_id,
                "element": atoms[atom_id]["element"],
                "role": atom["role"],
                "roles": roles,
                "background_kind": background_kind,
                "charge_e": charges[atom_id],
                "nearest_functional_group_graph_distance": min_distance,
                "nearest_functional_group_ids": nearest_groups
                if min_distance is not None
                else None,
                "graph_distance_definition": GRAPH_DISTANCE_DEFINITION,
            }
        )

    atom_rows = []
    group_by_atom = {}
    for group in group_rows:
        for atom_id in group["atom_ids"]:
            group_by_atom.setdefault(atom_id, []).append(group["group_id"])
    for atom in model["atoms"]:
        atom_rows.append(
            {
                "atom_id": atom["id"],
                "element": atom["element"],
                "name": atom.get("name"),
                "role": atom["role"],
                "roles": atom.get("roles", [atom["role"]]),
                "position_angstrom": atom["position"],
                "charge_e": charges[atom["id"]],
                "functional_group_ids": group_by_atom.get(atom["id"], []),
                "is_scaffold_carbon": atom["id"] in scaffold_set,
            }
        )

    return {
        "model_id": model["model_id"],
        "provenance": _provenance(item, model, fit_evidence),
        "charge_check": charge_check,
        "coverage": {
            "model_atom_count": len(atoms),
            "charge_atom_count": len(charges),
            "atom_identity_order_elements_match": True,
            "all_atom_charges_finite": True,
            "whole_model_charge_passed": True,
            "functional_group_count_excluding_GRT": len(group_rows),
            "functional_group_pair_count": len(pairs),
            "scaffold_carbon_count": len(scaffold_ids),
            "background_scaffold_carbon_count": len(background),
            "background_kind_counts": {
                kind: sum(row["background_kind"] == kind for row in background)
                for kind in ("interior_scaffold", "hydrogen_terminated_edge")
            },
            "scaffold_carbon_source": scaffold_source,
            "all_pair_distances_computed": all(
                p["graph_distance"] is not None for p in pairs
            ),
        },
        "scaffold_carbon_ids": scaffold_ids,
        "background_carbon_atom_ids": background_ids,
        "per_atom_charges": atom_rows,
        "functional_groups": group_rows,
        "functional_group_pairs": pairs,
        "background_carbon_distance_charges": background,
    }


def analyze_dataset(spec):
    """Analyze a JSON-list-compatible specification of model/charge fits."""
    entries = _read_json(spec)
    if not isinstance(entries, list) or not entries:
        raise ValueError("Dataset specification must be a nonempty JSON list")
    models = []
    for index, item in enumerate(entries, 1):
        model = _analyze_one(item)
        model["context_index"] = index
        model["context_id"] = f"{model['model_id']}::fit_{index:03d}"
        model["provenance"]["context_id"] = model["context_id"]
        models.append(model)
    ids = [model["model_id"] for model in models]
    return {
        "schema_version": 1,
        "analysis": "atomic partial-charge chemical context; no group averaging",
        "graph_distance_definition": GRAPH_DISTANCE_DEFINITION,
        "model_count": len(set(ids)),
        "context_count": len(models),
        "coverage": {
            "all_models_passed_charge_identity_and_total_checks": all(
                model["charge_check"]["passed"] for model in models
            ),
            "all_requested_functional_group_pairs_computed": all(
                model["coverage"]["all_pair_distances_computed"] for model in models
            ),
            "functional_group_pair_count": sum(
                len(model["functional_group_pairs"]) for model in models
            ),
        },
        "models": models,
    }


def _write_csv(path, rows, fieldnames=None):
    rows = list(rows)
    fields = fieldnames or sorted({key for row in rows for key in row})
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, allow_nan=False)
                    if isinstance(value, (list, dict))
                    else ("" if value is None else value)
                    for key, value in row.items()
                }
            )


def write_outputs(result, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "context_analysis.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    for key, filename, fields in (
        (
            "per_atom_charges",
            "per_atom_charges.csv",
            [
                "context_id",
                "model_id",
                "atom_id",
                "element",
                "name",
                "role",
                "roles",
                "position_angstrom",
                "charge_e",
                "functional_group_ids",
                "is_scaffold_carbon",
            ],
        ),
        (
            "functional_groups",
            "functional_group_charges.csv",
            [
                "context_id",
                "model_id",
                "group_id",
                "block_type",
                "atom_ids",
                "roles",
                "net_charge_e",
                "target_formal_charge_e",
                "scaffold_carbon_anchor_ids",
                "scaffold_anchor_source",
                "scaffold_anchor_charge_e_sum",
                "carbon_atom_ids_for_centroid",
                "carbon_centroid_angstrom",
            ],
        ),
        (
            "functional_group_pairs",
            "functional_group_pairs.csv",
            [
                "context_id",
                "model_id",
                "group_a_id",
                "group_a_type",
                "group_a_net_charge_e",
                "group_a_scaffold_anchor_charge_e_sum",
                "group_b_id",
                "group_b_type",
                "group_b_net_charge_e",
                "group_b_scaffold_anchor_charge_e_sum",
                "combined_group_net_charge_e",
                "graph_distance",
                "graph_distance_definition",
                "graph_anchor_atom_ids",
                "graph_path_atom_ids",
                "carbon_centroid_distance_angstrom",
                "centroid_definition",
            ],
        ),
        (
            "background_carbon_distance_charges",
            "background_carbon_distance_charges.csv",
            [
                "context_id",
                "model_id",
                "atom_id",
                "element",
                "role",
                "roles",
                "charge_e",
                "background_kind",
                "nearest_functional_group_graph_distance",
                "nearest_functional_group_ids",
                "graph_distance_definition",
            ],
        ),
    ):
        provenance_fields = [
            "scope",
            "protocol_scope",
            "model_scope",
            "is_diagnostic",
            "reported_wavefunction_stability_validated",
            "fit_status_is_acceptance",
        ]
        fields = [
            *fields[:2],
            *provenance_fields,
            *(field for field in fields[2:] if field not in provenance_fields),
        ]
        _write_csv(
            output_dir / filename,
            (
                {
                    "context_id": model["context_id"],
                    "model_id": model["model_id"],
                    **{
                        field: model["provenance"][field] for field in provenance_fields
                    },
                    **row,
                }
                for model in result["models"]
                for row in model[key]
            ),
            fields,
        )
    return output_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "spec_json", help="JSON list of model_path/charges_path/evidence entries"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Directory for JSON and CSV outputs"
    )
    args = parser.parse_args(argv)
    result = analyze_dataset(args.spec_json)
    write_outputs(result, args.output_dir)
    print(
        json.dumps(
            {"model_count": result["model_count"], "output_dir": args.output_dir}
        )
    )


if __name__ == "__main__":
    main()
