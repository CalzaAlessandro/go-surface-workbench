"""Explicit, atom-conserving anti epoxide opening for starting geometries.

A neutral amine/alcohol/thiol donor supplies its own explicit transfer hydrogen.
This operation changes the graph and places an anti initial geometry; it is
not a reaction-path, transition-state, kinetic, or force-field calculation.
"""

from __future__ import annotations

import copy
import math
import random

import numpy as np

from .build import _add_atom, _collision_free, _integer
from .errors import PlacementError
from .model import adjacency, canonical_bond, electronic_summary, validate_model


def _validated_fragment(fragment):
    if (
        not isinstance(fragment, dict)
        or not {"atoms", "bonds", "attachment_atom", "transfer_hydrogen"} <= fragment.keys()
    ):
        raise ValueError(
            "Provide explicit fragment atoms/bonds, attachment_atom and transfer_hydrogen"
        )
    donor = fragment["attachment_atom"]
    proton = fragment["transfer_hydrogen"]
    graph = {
        "schema_version": 1,
        "name": "nucleophile",
        "atoms": copy.deepcopy(fragment["atoms"]),
        "bonds": copy.deepcopy(fragment["bonds"]),
        "cell": [[100, 0, 0], [0, 100, 0], [0, 0, 100]],
        "pbc": [False] * 3,
        "metadata": copy.deepcopy(fragment.get("metadata", {})),
    }
    for atom in graph["atoms"]:
        if "formal_charge" not in atom:
            raise ValueError("Fragment atoms require explicit formal_charge")
    report = validate_model(graph)
    if report["connected_components"] != 1:
        raise ValueError("Nucleophile fragment must be connected")
    atoms = {a["id"]: a for a in graph["atoms"]}
    ns = adjacency(graph)
    if donor not in atoms or proton not in atoms or donor == proton:
        raise ValueError("Attachment and transfer-H IDs must exist and be distinct")
    if (
        atoms[proton]["element"] != "H"
        or atoms[proton]["formal_charge"] != 0
        or len(ns[proton]) != 1
        or ns[proton][0][0] != donor
    ):
        raise ValueError("Transfer hydrogen must be neutral H bonded to attachment atom")
    element = atoms[donor]["element"]
    required = {"N": 3, "O": 2, "S": 2}.get(element)
    if (
        required is None
        or atoms[donor]["formal_charge"] != 0
        or len(ns[donor]) != required
        or any(n[2] != 1 for n in ns[donor])
    ):
        raise ValueError(
            "Supported donors are explicit neutral amines, alcohols or thiols with a transferable attached H"
        )
    # Reject incomplete fragments instead of silently inventing missing H or
    # protonation/charge states. Charged remote groups are allowed if explicit.
    for identifier, atom in atoms.items():
        order = sum(n[2] for n in ns[identifier])
        charge = atom["formal_charge"]
        element = atom["element"]
        allowed = {
            "H": {0: {1}},
            "C": {0: {4, 4.5}},
            "N": {0: {3}, 1: {4}},
            "O": {0: {2}, -1: {1}},
            "S": {0: {2}},
            "F": {0: {1}},
            "Cl": {0: {1}},
            "Br": {0: {1}},
            "I": {0: {1}},
        }
        if order not in allowed.get(element, {}).get(charge, set()):
            raise ValueError(
                f"Unsupported/incomplete fragment valence or charge at atom {identifier}"
            )
    return graph, donor, proton


def _rotation(source, target):
    source = np.array(source, dtype=float)
    target = np.array(target, dtype=float)
    source /= np.linalg.norm(source)
    target /= np.linalg.norm(target)
    cross = np.cross(source, target)
    cosine = float(np.dot(source, target))
    if cosine > 1 - 1e-10:
        return np.eye(3)
    if cosine < -1 + 1e-10:
        axis = np.cross(source, [1, 0, 0])
        if np.linalg.norm(axis) < 1e-8:
            axis = np.cross(source, [0, 1, 0])
        axis /= np.linalg.norm(axis)
        return 2 * np.outer(axis, axis) - np.eye(3)
    matrix = np.array(
        [[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]]
    )
    return np.eye(3) + matrix + matrix @ matrix / (1 + cosine)


def _open_one(model, group, fragment, donor, proton, attacked, angle):
    atoms = {a["id"]: a for a in model["atoms"]}
    oxygen = group["oxygen_id"]
    carbons = group["site_ids"]
    retained = next(i for i in carbons if i != attacked)
    ns = adjacency(model)
    if {n[0] for n in ns[oxygen]} != set(carbons) or atoms[oxygen]["element"] != "O":
        raise ValueError("Epoxide metadata does not match intact explicit C-O-C graph")
    if any(len(ns[c]) != 4 for c in carbons):
        raise ValueError("Epoxide carbon must have four explicit neighbors before opening")
    cell = np.zeros((3, 3)) if model["cell"] is None else np.array(model["cell"])
    oxygen_image = next(image for j, image, _ in ns[retained] if j == oxygen)
    retained_position = np.array(atoms[retained]["position"])
    original_vector = np.array(atoms[oxygen]["position"]) + oxygen_image @ cell - retained_position
    attack_oxygen_image = next(image for j, image, _ in ns[attacked] if j == oxygen)
    original_attacked_co = (
        np.array(atoms[oxygen]["position"])
        + attack_oxygen_image @ cell
        - np.array(atoms[attacked]["position"])
    )
    incoming_axis = original_attacked_co / np.linalg.norm(original_attacked_co)
    carbon_vectors = [
        np.array(atoms[j]["position"]) + image @ cell - retained_position
        for j, image, _ in ns[retained]
        if atoms[j]["element"] == "C"
    ]
    if len(carbon_vectors) != 3:
        raise ValueError("Epoxide opening needs three sheet-carbon neighbors per basal carbon")
    normal = np.cross(carbon_vectors[1] - carbon_vectors[0], carbon_vectors[2] - carbon_vectors[0])
    norm = np.linalg.norm(normal)
    if norm < 1e-8:
        raise ValueError("Epoxide local sheet normal is undefined")
    normal /= norm
    projection = float(np.dot(original_vector, normal))
    if abs(projection) < 0.1:
        raise ValueError("Epoxide face is undefined")
    if projection < 0:
        normal *= -1
    face = (1 if normal[2] > 0 else -1) if abs(normal[2]) > 1e-8 else group.get("face", 1)
    # Keep original oxygen on its original face and retain the existing
    # periodic bond image to the carbon that keeps the C-O bond.
    atoms[oxygen]["position"] = (retained_position + normal * 1.43 - oxygen_image @ cell).tolist()
    model["bonds"] = [b for b in model["bonds"] if {b["a"], b["b"]} != {attacked, oxygen}]
    attachment_length = {"N": 1.47, "O": 1.43, "S": 1.81}[
        next(a["element"] for a in fragment["atoms"] if a["id"] == donor)
    ]
    attachment_position = np.array(atoms[attacked]["position"]) - incoming_axis * attachment_length
    fragment_atoms = {a["id"]: a for a in fragment["atoms"]}
    origin = np.array(fragment_atoms[donor]["position"])
    hvector = np.array(fragment_atoms[proton]["position"]) - origin
    # The consumed donor-H direction points toward the carbon. Rotation about
    # that axis searches fragment torsions while fixing the anti face.
    rotate = _rotation(hvector, incoming_axis)
    skew = np.array(
        [
            [0, -incoming_axis[2], incoming_axis[1]],
            [incoming_axis[2], 0, -incoming_axis[0]],
            [-incoming_axis[1], incoming_axis[0], 0],
        ]
    )
    torsion = np.eye(3) + math.sin(angle) * skew + (1 - math.cos(angle)) * (skew @ skew)
    mapping = {}
    used_names = set()
    resid = max(a.get("resid", 1) for a in model["atoms"]) + 1
    for atom in fragment["atoms"]:
        position = attachment_position + torsion @ rotate @ (np.array(atom["position"]) - origin)
        name = atom.get("name") or f"{atom['element']}{atom['id']}"
        if name in used_names:
            name = f"{atom['element']}{atom['id']}"
        if name in used_names:
            raise ValueError("Fragment atom names must be unique within its residue")
        used_names.add(name)
        identifier = _add_atom(
            model,
            atom["element"],
            position,
            name,
            "NUC",
            resid,
            role="nucleophile",
            junction=True,
            fragment_atom_id=atom["id"],
            epoxide_group_id=group["id"],
        )
        mapping[atom["id"]] = identifier
        model["atoms"][-1]["formal_charge"] = atom["formal_charge"]
        model["atoms"][-1]["metadata"]["fragment_atom_metadata"] = copy.deepcopy(
            atom.get("metadata", {})
        )
    for bond in fragment["bonds"]:
        if proton in (bond["a"], bond["b"]):
            continue
        new_bond = canonical_bond(mapping[bond["a"]], mapping[bond["b"]], bond["order"])
        new_bond["metadata"] = {"fragment_bond_metadata": copy.deepcopy(bond.get("metadata", {}))}
        model["bonds"].append(new_bond)
    attachment = mapping[donor]
    transferred = mapping[proton]
    model["bonds"].append(canonical_bond(attacked, attachment))
    # All atoms including transfer H are retained. The H changes its explicit
    # covalent partner and is positioned in the oxygen coordinate frame.
    new_atoms = {a["id"]: a for a in model["atoms"]}
    carbon_image = next(image for j, image, _ in ns[retained] if j == attacked)
    tangent = np.array(atoms[attacked]["position"]) + carbon_image @ cell - retained_position
    tangent -= normal * np.dot(tangent, normal)
    tangent /= np.linalg.norm(tangent)
    new_atoms[transferred]["position"] = (
        np.array(new_atoms[oxygen]["position"]) + 0.84 * tangent + 0.485 * normal
    ).tolist()
    new_atoms[transferred]["metadata"]["role"] = "transferred_hydrogen"
    model["bonds"].append(canonical_bond(oxygen, transferred))
    for identifier in carbons + [oxygen]:
        new_atoms[identifier]["atom_type"] = None
        new_atoms[identifier]["charge"] = None
        new_atoms[identifier]["metadata"].update(junction=True, group="opened_epoxide")
        new_atoms[identifier]["resname"] = "OPN"
    group.update(
        type="opened_epoxide",
        attachment_id=attachment,
        attacked_carbon_id=attacked,
        retained_carbon_id=retained,
        transferred_hydrogen_id=transferred,
        atom_ids=group["atom_ids"] + list(mapping.values()),
    )
    event = {
        "operation": "open_epoxide",
        "group_id": group["id"],
        "attacked_carbon_id": attacked,
        "retained_carbon_id": retained,
        "oxygen_id": oxygen,
        "attachment_id": attachment,
        "transferred_hydrogen_id": transferred,
        "oxygen_face": face,
        "attack_face": -face,
        "original_oxygen_side_normal": normal.tolist(),
        "fragment_atom_mapping": {str(k): v for k, v in mapping.items()},
        "geometry_status": "constructed_unrelaxed",
        "fragment_preparation": copy.deepcopy(fragment.get("metadata", {})),
    }
    model["metadata"]["events"].append(event)
    return _collision_free(model, set(mapping.values()) | {oxygen})


def open_epoxides(model: dict, recipe: dict) -> dict:
    """Return a copy with explicit anti opening and transfer of donor-H to O.

    Set `count` and optionally `id` (explicit group IDs); provide `fragment` containing
    atoms, bonds, attachment_atom and transfer_hydrogen. Fragment atoms need
    explicit positions and formal_charge. Complete neutral N/O/S donor
    chemistry is supported; no reaction is inferred from a molecule name.
    New junctions require parameterization before Amber/GROMACS preparation.
    """
    validate_model(model)
    if not isinstance(recipe, dict):
        raise ValueError("Opening recipe must contain count, optional id, fragment and seed")
    unknown = set(recipe) - {"count", "id", "fragment", "seed"}
    if unknown:
        raise ValueError(
            f"Unknown opening recipe keys: {sorted(unknown)}; use count and optional id"
        )
    fragment, donor, proton = _validated_fragment(recipe.get("fragment"))
    seed = _integer(recipe.get("seed", model.get("metadata", {}).get("seed", 0)), "seed")
    rng = random.Random(seed)
    groups = [g for g in model.get("metadata", {}).get("groups", []) if g["type"] == "epoxide"]
    explicit_sites = recipe.get("id") is not None
    if explicit_sites:
        selected = recipe["id"]
        ids = {g["id"] for g in groups}
        if (
            not isinstance(selected, list)
            or any(type(i) is not int or i not in ids for i in selected)
            or len(set(selected)) != len(selected)
        ):
            raise ValueError("id must list distinct available epoxide group IDs")
        requested_count = len(selected)
        if "count" in recipe and _integer(recipe["count"], "count") != requested_count:
            raise ValueError("count must match the number of id entries")
    else:
        count = _integer(recipe.get("count", 1), "count")
        if count > len(groups):
            raise PlacementError(
                "Requested opening count exceeds available epoxides",
                requested={"grafts": count},
                capacity={"available_epoxides": len(groups)},
                attempts=0,
                suggestions=[
                    "Reduce count or build a larger sheet with more epoxide groups.",
                    "Check that previously opened rings are not counted as available sites.",
                ],
            )
        requested_count = count
        # A count selects any feasible site. Keep searching other intact
        # epoxides when the first site's supplied fragment geometry clashes.
        selected = rng.sample([g["id"] for g in groups], len(groups)) if count else []
    opened = []
    attempted_orientations = 0
    out = copy.deepcopy(model)
    out["metadata"].setdefault("events", [])
    out["metadata"].setdefault("next_atom_id", max(a["id"] for a in out["atoms"]) + 1)
    for group_id in selected:
        success = False
        original_group = next(g for g in out["metadata"]["groups"] if g["id"] == group_id)
        attack_sites = list(original_group["site_ids"])
        rng.shuffle(attack_sites)
        offset = rng.random() * math.pi / 60
        # Complete 3-degree torsion coverage of both regioisomeric attacks.
        # A finite grid can still miss a narrower window; failures remain
        # explicit and do not relax the collision criterion.
        trials = [
            (attacked, offset + k * math.pi / 60) for k in range(120) for attacked in attack_sites
        ]
        for attacked, angle in trials:
            attempted_orientations += 1
            trial = copy.deepcopy(out)
            group = next(g for g in trial["metadata"]["groups"] if g["id"] == group_id)
            if _open_one(trial, group, fragment, donor, proton, attacked, angle):
                validate_model(trial)
                out = trial
                success = True
                break
        if success:
            opened.append(group_id)
            if len(opened) == requested_count:
                break
        elif explicit_sites:
            raise PlacementError(
                f"No collision-free anti opening found for epoxide group {group_id}",
                requested={"grafts": requested_count},
                achieved={"grafts": len(opened)},
                capacity={"available_epoxides": len(groups)},
                attempts=attempted_orientations,
                suggestions=[
                    "Choose different id sites, a lower loading or a larger sheet.",
                    "Try another fragment conformation; rigid attachment orientations were exhausted.",
                ],
            )
    if len(opened) != requested_count:
        raise PlacementError(
            f"No collision-free exact opening found: requested {requested_count}, "
            f"found {len(opened)} within site/torsion search",
            requested={"grafts": requested_count},
            achieved={"grafts": len(opened)},
            capacity={"available_epoxides": len(groups)},
            attempts=attempted_orientations,
            suggestions=[
                "Increase the sheet size or reduce the requested loading.",
                "Try another seed or fragment conformation; sites and rigid torsions were searched.",
            ],
        )
    selected = opened
    if selected:
        out["metadata"]["requires_parameterization"] = True
        out["metadata"]["parameterization_status"] = "new_covalent_junction_unparameterized"
        out["metadata"]["geometry_status"] = "constructed_unrelaxed"
        # A topology/charge assignment on the starting graph is invalidated.
        out["metadata"].pop("parameterized", None)
    out["metadata"].setdefault("functionalization_history", []).append(
        {
            "reaction": "anti_epoxide_opening",
            "mode": "explicit_sites" if explicit_sites else "random_count",
            "seed": seed,
            "available_epoxides": len(groups),
            "requested_count": requested_count,
            "placed_count": len(selected),
            "actual_fraction": len(selected) / len(groups) if groups else 0.0,
            "selected_epoxide_ids": selected,
            "placement": "shuffled_sites_with_steric_rejection_and_rigid_fragment_torsion_search",
        }
    )
    out["metadata"].update(electronic_summary(out))
    validate_model(out)
    expected_atoms = len(model["atoms"]) + len(selected) * len(fragment["atoms"])
    expected_charge = sum(a.get("formal_charge", 0) for a in model["atoms"]) + len(selected) * sum(
        a["formal_charge"] for a in fragment["atoms"]
    )
    if (
        len(out["atoms"]) != expected_atoms
        or sum(a.get("formal_charge", 0) for a in out["atoms"]) != expected_charge
    ):
        raise ValueError("Opening failed atom/charge conservation")
    return out
