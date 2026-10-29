"""Reproducible explicit honeycomb graphs and functional-group construction.

Original implementation of analytical four-carbon rectangular honeycomb tiles.
The rectangular-cell/edge chemistry design is informed by HierGO/GOPY; no
source code is copied. This constructs starting geometry, not an equilibrated
GO ensemble, a charge model, or a validated force field.
"""

from __future__ import annotations

import copy
import math
import random

import numpy as np

from .errors import PlacementError
from .model import adjacency, canonical_bond, electronic_summary, validate_model

GROUPS = ("OH", "epoxide", "carbonyl", "COOH", "COO")


def _integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _atom(identifier, element, position, name, resname="GGG", resid=1, **metadata):
    return {
        "id": identifier,
        "element": element,
        "name": name,
        "resname": resname,
        "resid": resid,
        "position": [float(v) for v in position],
        "atom_type": None,
        "charge": None,
        "formal_charge": 0,
        "metadata": metadata,
    }


def _add_atom(model, element, position, name, resname, resid, **metadata):
    identifier = model["metadata"]["next_atom_id"]
    model["metadata"]["next_atom_id"] += 1
    model["atoms"].append(_atom(identifier, element, position, name, resname, resid, **metadata))
    return identifier


def _coordinates(model):
    return {a["id"]: np.array(a["position"], dtype=float) for a in model["atoms"]}


def _displacement(model, a, b, image):
    coords = _coordinates(model)
    return coords[b] + np.asarray(image) @ np.asarray(model["cell"]) - coords[a]


def _minimum_image(delta, model):
    result = np.array(delta, dtype=float)
    if not any(model["pbc"]):
        return result
    cell = np.array(model["cell"])
    fractional = result @ np.linalg.inv(cell)
    for k, periodic in enumerate(model["pbc"]):
        if periodic:
            fractional[..., k] -= np.rint(fractional[..., k])
    return fractional @ cell


def _collision_free_vectorized(model, new_ids):
    """Reject severe nonbonded overlaps with O(N) temporary arrays."""
    atoms = model["atoms"]
    atom_ids = np.fromiter((atom["id"] for atom in atoms), dtype=np.int64, count=len(atoms))
    elements = np.asarray([atom["element"] for atom in atoms])
    coords = np.asarray([atom["position"] for atom in atoms], dtype=float)
    bonded_neighbors = {}
    for bond in model["bonds"]:
        bonded_neighbors.setdefault(bond["a"], set()).add(bond["b"])
        bonded_neighbors.setdefault(bond["b"], set()).add(bond["a"])

    for index, atom in enumerate(atoms):
        if atom["id"] not in new_ids:
            continue
        distances = np.linalg.norm(_minimum_image(coords - coords[index], model), axis=1)
        eligible = np.ones(len(atoms), dtype=bool)
        eligible[index] = False
        bonded_ids = bonded_neighbors.get(atom["id"])
        if bonded_ids:
            eligible &= ~np.isin(atom_ids, tuple(bonded_ids))
        cutoffs = np.where((elements == "H") | (atom["element"] == "H"), 0.8, 1.25)
        if np.any(distances[eligible] < cutoffs[eligible]):
            return False
    return True


def _collision_free(model, new_ids):
    """Reject severe nonbonded overlaps, using periodic minimum images.

    This is a steric screen only. Passing it does not establish optimized
    bond angles, torsions, energetics, or suitability for dynamics.
    """
    return _collision_free_vectorized(model, new_ids)


def _connected(model):
    ns = adjacency(model)
    if not ns:
        return False
    visited = set()
    stack = [next(iter(ns))]
    while stack:
        identifier = stack.pop()
        if identifier in visited:
            continue
        visited.add(identifier)
        stack.extend(n[0] for n in ns[identifier])
    return len(visited) == len(ns)


def _remove_atoms(model, removed):
    model["atoms"] = [a for a in model["atoms"] if a["id"] not in removed]
    model["bonds"] = [b for b in model["bonds"] if not {b["a"], b["b"]} & removed]


def _prune_tips(model):
    removed = set()
    while True:
        ns = adjacency(model)
        tips = {i for i, neighbors in ns.items() if len(neighbors) < 2}
        if not tips:
            return removed
        removed.update(tips)
        _remove_atoms(model, tips)
        if not model["atoms"]:
            return removed


def _tile(nx, ny, length, pbc, vacuum, seed, name):
    h = math.sqrt(3) * length / 2
    template = [(0, 0), (length, 0), (1.5 * length, h), (2.5 * length, h)]
    cell = [[3 * length * nx, 0, 0], [0, 2 * h * ny, 0], [0, 0, vacuum]]
    model = {
        "schema_version": 1,
        "name": name,
        "atoms": [],
        "bonds": [],
        "cell": cell,
        "pbc": pbc,
        "metadata": {
            "seed": seed,
            "events": [],
            "groups": [],
            "geometry_status": "constructed_unrelaxed",
            "next_atom_id": 4 * nx * ny + 1,
            "construction": "analytical_rectangular_honeycomb",
        },
    }

    def identifier(i, j, k):
        return 4 * (i * ny + j) + k + 1

    for i in range(nx):
        for j in range(ny):
            for k, (x, y) in enumerate(template):
                model["atoms"].append(
                    _atom(
                        identifier(i, j, k),
                        "C",
                        [x + 3 * length * i, y + 2 * h * j, 0],
                        "CX",
                        role="sheet_carbon",
                        hybridization="sp2",
                        sublattice=(0 if k in (0, 2) else 1),
                    )
                )
    # Only these analytic lattice edges exist. No distance-based bond inference.
    edges = [(0, 1, 0, 0), (1, 2, 0, 0), (1, 2, 0, -1), (2, 3, 0, 0), (3, 0, 1, 0), (3, 0, 1, 1)]
    for i in range(nx):
        for j in range(ny):
            for ka, kb, di, dj in edges:
                ii, jj = i + di, j + dj
                if not pbc[0] and not 0 <= ii < nx:
                    continue
                if not pbc[1] and not 0 <= jj < ny:
                    continue
                image = [ii // nx, jj // ny, 0]
                model["bonds"].append(
                    canonical_bond(
                        identifier(i, j, ka), identifier(ii % nx, jj % ny, kb), 1.5, image
                    )
                )
    return model


def _apply_holes(model, recipe, rng):
    removed = set()
    ids = {a["id"] for a in model["atoms"]}
    vacancies = recipe.get("vacancies", 0)
    if type(vacancies) is int:
        _integer(vacancies, "vacancies")
        if vacancies >= len(ids):
            raise ValueError("Too many vacancies")
        removed.update(rng.sample(sorted(ids), vacancies))
    elif isinstance(vacancies, list):
        if any(type(i) is not int or i not in ids for i in vacancies) or len(set(vacancies)) != len(
            vacancies
        ):
            raise ValueError("Vacancies must list distinct existing carbon IDs")
        removed.update(vacancies)
    else:
        raise ValueError("vacancies must be a count or a list of original carbon IDs")
    holes = recipe.get("holes", [])
    if not isinstance(holes, list) or any(not isinstance(h, dict) for h in holes):
        raise ValueError("holes must be a list of hole specification objects")
    for hole in holes:
        available = ids - removed
        if "center" in hole and "radius" in hole:
            center = np.array(hole["center"], dtype=float)
            radius = float(hole["radius"])
            if (
                center.shape != (2,)
                or not np.isfinite(center).all()
                or not math.isfinite(radius)
                or radius <= 0
            ):
                raise ValueError("Hole requires finite xy center and positive radius")
            for atom in model["atoms"]:
                delta = np.array(atom["position"]) - np.r_[center, 0]
                if np.linalg.norm(_minimum_image(delta, model)[:2]) <= radius:
                    removed.add(atom["id"])
        elif "size" in hole:
            size = _integer(hole["size"], "hole size", 1)
            if size >= len(available):
                raise ValueError("Cluster hole consumes the sheet")
            start = hole.get("seed_atom", rng.choice(sorted(available)))
            if start not in available:
                raise ValueError("Hole seed is not an available carbon")
            ns = adjacency(model)
            cluster = {start}
            while len(cluster) < size:
                frontier = sorted({n[0] for i in cluster for n in ns[i]} & available - cluster)
                if not frontier:
                    raise ValueError("Cannot grow connected hole of requested size")
                cluster.add(rng.choice(frontier))
            removed.update(cluster)
        else:
            raise ValueError("Hole requires center/radius or size/seed_atom")
    _remove_atoms(model, removed)
    removed.update(_prune_tips(model))
    if not model["atoms"] or not _connected(model):
        raise ValueError("Vacancies/holes leave an empty or disconnected sheet")
    model["metadata"]["removed_carbon_ids"] = sorted(removed)
    model["metadata"]["requested_vacancies"] = copy.deepcopy(vacancies)


def _max_matching(model, candidates, rng):
    """Exact cardinality matching on the bipartite honeycomb, randomized ties."""
    by_id = {a["id"]: a for a in model["atoms"]}
    ns = adjacency(model)
    left = [i for i in candidates if by_id[i]["metadata"]["sublattice"] == 0]
    rng.shuffle(left)
    match = {}

    def augment(i, seen):
        options = [n[0] for n in ns[i] if n[0] in candidates]
        rng.shuffle(options)
        for j in options:
            if j in seen:
                continue
            seen.add(j)
            if j not in match or augment(match[j], seen):
                match[j] = i
                return True
        return False

    for i in left:
        augment(i, set())
    return [(i, j) for j, i in match.items()]


def _set_site(model, identifier, name, resname, resid, group_id, group, sp3=False):
    atom = next(a for a in model["atoms"] if a["id"] == identifier)
    atom.update(name=name, resname=resname, resid=resid)
    atom["metadata"].update(group_id=group_id, group=group)
    if sp3:
        atom["metadata"]["hybridization"] = "sp3"
    # Saturated/edge carbonyl sites cannot retain aromatic bond order.
    if sp3 or group == "carbonyl":
        for bond in model["bonds"]:
            if identifier in (bond["a"], bond["b"]):
                bond["order"] = 1.0


def _outward(model, identifier, *, neighbors=None, coordinates=None, cell=None):
    if neighbors is None:
        neighbors = adjacency(model)
    if coordinates is None:
        coordinates = _coordinates(model)
    if cell is None:
        cell = np.array(model["cell"])
    ns = neighbors[identifier]
    coords = coordinates
    direction = -sum(
        (coords[j] + image @ cell - coords[identifier] for j, image, _ in ns),
        start=np.zeros(3),
    )
    norm = np.linalg.norm(direction)
    if norm < 1e-8:
        raise ValueError("Undefined edge outward direction")
    return direction / norm


def _place_group(model, group, sites, face, group_id):
    resid = group_id + 1
    coords = _coordinates(model)
    origin = coords[sites[0]]
    initial_ids = {a["id"] for a in model["atoms"]}
    resname = {"OH": "HYD", "epoxide": "EPX", "carbonyl": "CAB", "COOH": "CBX", "COO": "CBD"}[group]
    atom_ids = list(sites)

    def add(element, position, name, charge=0):
        identifier = _add_atom(
            model,
            element,
            position,
            name,
            resname,
            resid,
            role="group",
            group_id=group_id,
            group=group,
        )
        model["atoms"][-1]["formal_charge"] = charge
        atom_ids.append(identifier)
        return identifier

    def bond(a, b, order=1, image=(0, 0, 0)):
        model["bonds"].append(canonical_bond(a, b, order, image))

    if group == "OH":
        _set_site(model, sites[0], "CL", resname, resid, group_id, group, True)
        oxygen = add("O", origin + [0, 0, face * 1.43], "OL")
        hydrogen = add("H", origin + [0.83, 0, face * 1.93], "HL")
        bond(sites[0], oxygen)
        bond(oxygen, hydrogen)
    elif group == "epoxide":
        a, b = sites
        image = next(image for j, image, _ in adjacency(model)[a] if j == b)
        midpoint = (origin + coords[b] + image @ np.array(model["cell"])) / 2
        half = np.linalg.norm(coords[b] + image @ np.array(model["cell"]) - origin) / 2
        height = math.sqrt(1.46**2 - half**2)
        oxygen = add("O", midpoint + [0, 0, face * height], "OE")
        _set_site(model, a, "CY", resname, resid, group_id, group, True)
        _set_site(model, b, "CZ", resname, resid, group_id, group, True)
        bond(a, oxygen)
        bond(b, oxygen, image=-image)
    else:
        direction = _outward(model, sites[0])
        perpendicular = np.cross([0, 0, 1], direction)
        _set_site(
            model, sites[0], "CN" if group == "carbonyl" else "CY", resname, resid, group_id, group
        )
        if group == "carbonyl":
            oxygen = add("O", origin + 1.23 * direction, "ON")
            bond(sites[0], oxygen, 2)
        else:
            center = origin + 1.52 * direction
            carbon = add("C", center, "C4")
            bond(sites[0], carbon)
            oxygen = add(
                "O", center + 1.21 * (0.5 * direction + math.sqrt(3) / 2 * perpendicular), "OJ"
            )
            other = add(
                "O",
                center + 1.32 * (0.5 * direction - math.sqrt(3) / 2 * perpendicular),
                "OK",
                -1 if group == "COO" else 0,
            )
            bond(carbon, oxygen, 2)
            bond(carbon, other)
            if group == "COOH":
                hydrogen = add(
                    "H",
                    center
                    + 1.32 * (0.5 * direction - math.sqrt(3) / 2 * perpendicular)
                    + 0.97 * direction,
                    "HK",
                )
                bond(other, hydrogen)
    model["metadata"]["groups"].append(
        {
            "id": group_id,
            "type": group,
            "site_ids": list(sites),
            "atom_ids": atom_ids,
            "oxygen_id": oxygen,
            "face": face,
        }
    )
    return _collision_free(model, {a["id"] for a in model["atoms"]} - initial_ids)


def _group_requests(requested):
    counts, fixed = {}, {}
    for group in GROUPS:
        value = requested.get(group, 0)
        label = f"groups.{group}"
        if not isinstance(value, dict):
            counts[group] = _integer(value, label)
            continue
        if not value or set(value) - {"count", "id"}:
            raise ValueError(f"{label} requires count and/or id only")
        if "id" in value:
            sites = value["id"]
            if not isinstance(sites, list):
                raise ValueError(f"{label}.id must be a list of original carbon sites")
            if group == "epoxide":
                if any(not isinstance(pair, list) or len(pair) != 2 for pair in sites):
                    raise ValueError(f"{label}.id must list C-C pairs [[a,b],...]")
                fixed[group] = [list(pair) for pair in sites]
            else:
                fixed[group] = [[site] for site in sites]
            count = _integer(value.get("count", len(sites)), f"{label}.count")
            if count != len(sites):
                raise ValueError(f"{label}.count must equal the number of sites in id")
        else:
            count = _integer(value["count"], f"{label}.count")
        counts[group] = count
    return counts, fixed


def _decorate(base, counts, fixed, faces, rng):
    ns = adjacency(base)
    edge = [i for i, n in ns.items() if len(n) == 2]
    basal = {i for i, n in ns.items() if len(n) == 3}
    reserved = set()
    for group, sites in fixed.items():
        eligible = basal if group in ("OH", "epoxide") else set(edge)
        for site in sites:
            for identifier in site:
                if type(identifier) is not int or identifier not in ns:
                    raise ValueError(f"groups.{group}.id requires existing original carbon IDs")
                if identifier in reserved:
                    raise ValueError("Functional group site IDs must be distinct and disjoint")
                if identifier not in eligible:
                    degree = 3 if group in ("OH", "epoxide") else 2
                    raise ValueError(f"groups.{group}.id requires degree-{degree} carbon sites")
                reserved.add(identifier)
            if group == "epoxide" and not any(j == site[1] for j, _, _ in ns[site[0]]):
                raise ValueError("groups.epoxide.id requires actual basal C-C bonds")
    # Reserve every fixed site before choosing any random group, including
    # groups placed earlier in the construction sequence.
    edge = [i for i in edge if i not in reserved]
    basal -= reserved
    random_counts = {group: 0 if group in fixed else count for group, count in counts.items()}
    edge_count = sum(random_counts[g] for g in ("carbonyl", "COOH", "COO"))
    capacity = {"eligible_edge_carbons": len(edge), "eligible_basal_carbons": len(basal)}
    suggestions = ["Enlarge the sheet", "Reduce requested group counts", "Change the seed"]
    if edge_count > len(edge) or random_counts["OH"] + 2 * random_counts["epoxide"] > len(basal):
        raise PlacementError(
            "Requested functional group counts exceed eligible edge/basal sites",
            requested=counts,
            capacity=capacity,
            attempts=0,
            suggestions=suggestions,
        )
    matching = _max_matching(base, basal, rng)
    capacity["max_disjoint_epoxide_sites"] = len(matching)
    if random_counts["epoxide"] > len(matching):
        raise PlacementError(
            "Requested epoxide count exceeds disjoint eligible C-C bonds",
            requested=counts,
            capacity=capacity,
            attempts=0,
            suggestions=suggestions,
        )
    # Retry stereochemical placements and randomized tie choices. Count is exact
    # on success; search failure is explicit and never silently undershoots.
    best_achieved = dict.fromkeys(GROUPS, 0)
    for _ in range(128):
        model = copy.deepcopy(base)
        rng.shuffle(edge)
        matching = _max_matching(base, basal, rng)
        rng.shuffle(matching)
        epoxy_sites = fixed.get("epoxide", matching[: random_counts["epoxide"]])
        used = {i for pair in epoxy_sites for i in pair}
        oh_sites = sorted(basal - used)
        rng.shuffle(oh_sites)
        assignments = []
        cursor = 0
        for group in ("COOH", "COO", "carbonyl"):
            sites = fixed.get(
                group, [[site] for site in edge[cursor : cursor + random_counts[group]]]
            )
            assignments += [(group, site) for site in sites]
            cursor += random_counts[group]
        assignments += [("epoxide", list(pair)) for pair in epoxy_sites]
        assignments += [
            ("OH", site) for site in fixed.get("OH", [[i] for i in oh_sites[: random_counts["OH"]]])
        ]
        success = True
        achieved = dict.fromkeys(GROUPS, 0)
        for group_id, (group, sites) in enumerate(assignments, 1):
            face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
            if not _place_group(model, group, sites, face, group_id):
                success = False
                break
            achieved[group] += 1
        if success:
            return model
        if sum(achieved.values()) > sum(best_achieved.values()):
            best_achieved = achieved
    raise PlacementError(
        "No collision-free exact placement found in 128 attempts; "
        f"best tried placement counts: {best_achieved}; adjust composition, seed or size",
        requested=counts,
        achieved=best_achieved,
        capacity=capacity,
        attempts=128,
        suggestions=suggestions,
    )


def _cap_edges(model):
    ns = adjacency(model)
    coordinates = _coordinates(model)
    cell = np.array(model["cell"])
    atoms_by_id = {atom["id"]: atom for atom in model["atoms"]}
    next_resid = max(a.get("resid", 1) for a in model["atoms"]) + 1
    occupied = {i for g in model["metadata"]["groups"] for i in g["site_ids"]}
    for identifier, neighbors in ns.items():
        atom = atoms_by_id[identifier]
        if (
            atom.get("metadata", {}).get("role") != "sheet_carbon"
            or len(neighbors) != 2
            or identifier in occupied
        ):
            continue
        resid = next_resid
        next_resid += 1
        atom.update(name="CH", resname="GRT", resid=resid)
        outward = _outward(
            model, identifier, neighbors=ns, coordinates=coordinates, cell=cell
        )
        hydrogen = _add_atom(
            model,
            "H",
            np.array(atom["position"]) + 1.09 * outward,
            "HC",
            "GRT",
            resid,
            role="edge_hydrogen",
            anchor_id=identifier,
        )
        model["bonds"].append(canonical_bond(identifier, hydrogen))
        # Small holes have colliding planar inward H caps. Explicitly tilt
        # starting caps when necessary; preserve the C-H bond and report it.
        for angle in (0, 35, -35, 50, -50):
            radians = math.radians(angle)
            direction = math.cos(radians) * outward + np.array([0, 0, math.sin(radians)])
            model["atoms"][-1]["position"] = (
                np.array(atom["position"]) + 1.09 * direction
            ).tolist()
            if _collision_free(model, {hydrogen}):
                model["atoms"][-1]["metadata"]["cap_tilt_degrees"] = angle
                break
        else:
            raise ValueError("No nonoverlapping edge H cap; enlarge or modify the hole")


def build_sheet(recipe: dict) -> dict:
    """Build exact-count GO from rectangular honeycomb cells.

    Recipe: nx/ny >=2 (4 C per untrimmed cell), seed, pbc=[x,y,false],
    bond_length=1.42 Å, vacuum=40 Å, vacancies=count or original ID list,
    holes=[{center:[x,y],radius:r} or {size:n,seed_atom:id}],
    groups={OH,epoxide,carbonyl,COOH,COO: integer counts}, faces=both/up/down.
    Each group may instead be {count:N,id:[...]}: OH/carbonyl/COOH/COO id
    lists original sheet-carbon IDs; epoxide id lists bonded C-C pairs
    [[a,b],...]. Count may be omitted with id and then equals its length;
    when both are given they must agree. Counts must be nonnegative integers.
    Integer counts and {count:N} choose random sites. Explicit IDs must exist
    after holes/tip removal, be distinct across all groups, and are reserved
    before random selection. Global retries retain these sites and enforce
    the same faces and steric checks; failure never substitutes other sites.
    Basal OH/epoxide require degree-three carbon; other groups require an
    edge/defect degree-two site. Tips are removed, remaining edges H-capped.
    """
    if not isinstance(recipe, dict):
        raise ValueError("Recipe must be a dictionary")
    allowed = {
        "nx",
        "ny",
        "seed",
        "pbc",
        "name",
        "bond_length",
        "vacuum",
        "vacancies",
        "holes",
        "groups",
        "faces",
    }
    if set(recipe) - allowed:
        raise ValueError(f"Unknown recipe keys: {sorted(set(recipe) - allowed)}")
    nx, ny = _integer(recipe.get("nx", 6), "nx", 2), _integer(recipe.get("ny", 6), "ny", 2)
    seed = _integer(recipe.get("seed", 0), "seed")
    pbc = recipe.get("pbc", [False, False, False])
    if (
        not isinstance(pbc, list)
        or len(pbc) != 3
        or any(type(x) is not bool for x in pbc)
        or pbc[2]
    ):
        raise ValueError("pbc must be [bool,bool,false]")
    length, vacuum = float(recipe.get("bond_length", 1.42)), float(recipe.get("vacuum", 40))
    if (
        not math.isfinite(length)
        or not 1.2 <= length <= 1.6
        or not math.isfinite(vacuum)
        or vacuum < 15
    ):
        raise ValueError("bond_length must be 1.2–1.6 Å and vacuum >=15 Å")
    requested = recipe.get("groups", {})
    if not isinstance(requested, dict) or set(requested) - set(GROUPS):
        raise ValueError("Unknown functional group")
    counts, fixed = _group_requests(requested)
    faces = recipe.get("faces", "both")
    if faces not in ("both", "up", "down"):
        raise ValueError("faces must be both/up/down")
    rng = random.Random(seed)
    model = _tile(nx, ny, length, list(pbc), vacuum, seed, str(recipe.get("name", "GO")))
    _apply_holes(model, recipe, rng)
    model = _decorate(model, counts, fixed, faces, rng)
    _cap_edges(model)
    # GGG is a one-atom library fragment, so each pristine carbon needs a
    # distinct residue while retaining its original atom ID and graph order.
    next_resid = max(a.get("resid", 1) for a in model["atoms"]) + 1
    for atom in model["atoms"]:
        if atom["resname"] == "GGG":
            atom["resid"] = next_resid
            next_resid += 1
    model["metadata"]["recipe"] = copy.deepcopy(recipe)
    model["metadata"]["events"].append(
        {
            "operation": "build_sheet",
            "group_counts": counts,
            "removed_carbon_ids": model["metadata"]["removed_carbon_ids"],
        }
    )
    model["metadata"].update(electronic_summary(model))
    validate_model(model)
    if not _collision_free(model, {a["id"] for a in model["atoms"]}):
        raise ValueError("Hydrogen capping creates severe overlap; adjust defects/size")
    return model
