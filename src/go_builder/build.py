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
from collections import deque
from numbers import Real

import numpy as np

from .errors import PlacementError
from .model import (
    _build_spatial_cells,
    _compute_scaffold_diameter_nm,
    _minimum_image_delta,
    adjacency,
    assert_export_ready,
    canonical_bond,
    check_geometry_integrity,
    electronic_summary,
)

LEGACY_GROUPS = ("OH", "epoxide", "carbonyl", "COOH", "COO")
PRIORITY2_GROUPS = ("epoxide_eph", "epoxide_ept")
GROUPS = LEGACY_GROUPS + PRIORITY2_GROUPS


def _integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _finite_real(value, label):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(float(value)):
        raise ValueError(f"{label} must be a finite number")
    return float(value)


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
    cell = np.asarray(model["cell"], dtype=float)
    inv_cell = np.linalg.inv(cell)
    return _minimum_image_delta(result, model["pbc"], cell, inv_cell)


class _SpatialGrid:
    """Incremental 3D PBC-aware voxel cell index for O(1) local neighbor queries."""

    def __init__(self, coords_list, ids_list, elems_list, pbc, cell, inv_cell, cell_size=6.5):
        self.coords = [list(map(float, c)) for c in coords_list]
        self.ids = list(ids_list)
        self.elems = list(elems_list)
        self.pbc = list(pbc)
        self.cell = np.asarray(cell, dtype=float)
        self.inv_cell = inv_cell
        self.has_pbc = any(self.pbc)
        self.cell_size = float(cell_size)
        self.id_to_slot = {aid: idx for idx, aid in enumerate(self.ids)}
        self.slot_key: list[tuple[int, int, int]] = []
        self.cells: dict[tuple[int, int, int], list[int]] = {}
        self.n_bins: list[int | None] = [None, None, None]
        self.h_perp: list[float] = [0.0, 0.0, 0.0]
        if self.has_pbc:
            for k in range(3):
                h_k = 1.0 / float(np.linalg.norm(self.inv_cell[:, k]))
                self.h_perp[k] = h_k
                if self.pbc[k]:
                    self.n_bins[k] = max(1, int(math.floor(h_k / self.cell_size)))

        for slot, pos in enumerate(self.coords):
            key = self._key(pos)
            self.slot_key.append(key)
            self.cells.setdefault(key, []).append(slot)

    def _key(self, pos) -> tuple[int, int, int]:
        if not self.has_pbc:
            cs = self.cell_size
            return (
                int(math.floor(pos[0] / cs)),
                int(math.floor(pos[1] / cs)),
                int(math.floor(pos[2] / cs)),
            )
        frac = np.asarray(pos, dtype=float) @ self.inv_cell
        b = [0, 0, 0]
        for k in range(3):
            m_k = self.n_bins[k]
            if m_k is not None:
                fk = float(frac[k]) - math.floor(float(frac[k]))
                b[k] = int(math.floor(fk * m_k)) % m_k
            else:
                b[k] = int(math.floor(float(frac[k]) * (self.h_perp[k] / self.cell_size)))
        return (b[0], b[1], b[2])

    def add(self, aid: int, pos, elem: str) -> None:
        slot = len(self.coords)
        pos_list = [float(pos[0]), float(pos[1]), float(pos[2])]
        self.coords.append(pos_list)
        self.ids.append(int(aid))
        self.elems.append(str(elem))
        self.id_to_slot[int(aid)] = slot
        key = self._key(pos_list)
        self.slot_key.append(key)
        self.cells.setdefault(key, []).append(slot)

    def update_pos(self, aid: int, new_pos) -> None:
        slot = self.id_to_slot[int(aid)]
        pos_list = [float(new_pos[0]), float(new_pos[1]), float(new_pos[2])]
        self.coords[slot] = pos_list
        new_key = self._key(pos_list)
        old_key = self.slot_key[slot]
        if new_key != old_key:
            self.cells[old_key].remove(slot)
            self.cells.setdefault(new_key, []).append(slot)
            self.slot_key[slot] = new_key

    def query_radius(self, origin, radius: float):
        if not self.coords:
            return (
                np.empty((0, 3), dtype=float),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype="<U2"),
            )
        ck = self._key(origin)
        offs = []
        for k in range(3):
            b_k = ck[k]
            m_k = self.n_bins[k]
            if m_k is not None:
                offs.append(sorted({(b_k + d) % m_k for d in (-1, 0, 1)}))
            else:
                offs.append((b_k - 1, b_k, b_k + 1))

        cand_slots = []
        for kx in offs[0]:
            for ky in offs[1]:
                for kz in offs[2]:
                    bucket = self.cells.get((kx, ky, kz))
                    if bucket:
                        cand_slots.extend(bucket)
        if not cand_slots:
            return (
                np.empty((0, 3), dtype=float),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype="<U2"),
            )
        sub_coords = np.asarray([self.coords[s] for s in cand_slots], dtype=float)
        sub_ids = np.asarray([self.ids[s] for s in cand_slots], dtype=np.int64)
        sub_elems = np.asarray([self.elems[s] for s in cand_slots])
        delta = sub_coords - np.asarray(origin, dtype=float)
        if self.has_pbc:
            delta = _minimum_image_delta(delta, self.pbc, self.cell, self.inv_cell)
        dists = np.linalg.norm(delta, axis=1)
        mask = dists <= radius
        return sub_coords[mask], sub_ids[mask], sub_elems[mask]


def _make_spatial_grid(model, cell_size=6.5):
    atoms = model["atoms"]
    pbc = model["pbc"]
    cell = np.asarray(model["cell"], dtype=float)
    inv_cell = np.linalg.inv(cell) if any(pbc) else None
    return _SpatialGrid(
        [a["position"] for a in atoms],
        [a["id"] for a in atoms],
        [a["element"] for a in atoms],
        pbc,
        cell,
        inv_cell,
        cell_size=cell_size,
    )


def _collision_free_vectorized(
    model, new_ids, *, min_h_distance: float = 0.90, min_heavy_distance: float = 1.25
):
    """Reject severe nonbonded overlaps with O(N) spatial cell-list screening."""
    atoms = model["atoms"]
    if not atoms or not new_ids:
        return True
    n_atoms = len(atoms)
    atom_ids = np.fromiter((atom["id"] for atom in atoms), dtype=np.int64, count=n_atoms)
    elements = np.asarray([atom["element"] for atom in atoms])
    coords = np.asarray([atom["position"] for atom in atoms], dtype=float)
    bonded_neighbors = {}
    for bond in model["bonds"]:
        bonded_neighbors.setdefault(bond["a"], set()).add(bond["b"])
        bonded_neighbors.setdefault(bond["b"], set()).add(bond["a"])

    if n_atoms <= 256:
        for index, atom in enumerate(atoms):
            if atom["id"] not in new_ids:
                continue
            distances = np.linalg.norm(_minimum_image(coords - coords[index], model), axis=1)
            eligible = np.ones(n_atoms, dtype=bool)
            eligible[index] = False
            bonded_ids = bonded_neighbors.get(atom["id"])
            if bonded_ids:
                eligible &= ~np.isin(atom_ids, tuple(bonded_ids))
            cutoffs = np.where(
                (elements == "H") | (atom["element"] == "H"), min_h_distance, min_heavy_distance
            )
            if np.any(distances[eligible] < cutoffs[eligible]):
                return False
        return True

    pbc = model["pbc"]
    cell = np.zeros((3, 3), dtype=float) if model["cell"] is None else np.asarray(model["cell"], dtype=float)
    inv_cell = np.linalg.inv(cell) if any(pbc) else None
    cutoff_max = max(2.0, min_h_distance, min_heavy_distance)
    cells, n_bins = _build_spatial_cells(coords, pbc, cell, inv_cell, cutoff_max)
    atom_cell = {}
    for key, idxs in cells.items():
        for idx in idxs.tolist():
            atom_cell[idx] = key

    for index, atom in enumerate(atoms):
        aid = atom["id"]
        if aid not in new_ids:
            continue
        ck = atom_cell[index]
        offs = []
        for k in range(3):
            b_k = ck[k]
            m_k = n_bins[k]
            if m_k is not None:
                offs.append(sorted({(b_k + d) % m_k for d in (-1, 0, 1)}))
            else:
                offs.append((b_k - 1, b_k, b_k + 1))
        cand_chunks = []
        for kx in offs[0]:
            for ky in offs[1]:
                for kz in offs[2]:
                    bucket = cells.get((kx, ky, kz))
                    if bucket is not None:
                        cand_chunks.append(bucket)
        if not cand_chunks:
            continue
        cand_idx = np.concatenate(cand_chunks)
        cand_idx = cand_idx[cand_idx != index]
        if not len(cand_idx):
            continue
        bonded_ids = bonded_neighbors.get(aid)
        if bonded_ids:
            cand_idx = cand_idx[~np.isin(atom_ids[cand_idx], tuple(bonded_ids))]
        if not len(cand_idx):
            continue
        delta = _minimum_image_delta(coords[cand_idx] - coords[index], pbc, cell, inv_cell)
        distances = np.linalg.norm(delta, axis=1)
        cutoffs = np.where(
            (elements[cand_idx] == "H") | (atom["element"] == "H"),
            min_h_distance,
            min_heavy_distance,
        )
        if np.any(distances < cutoffs):
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


_HOLE_PROFILE_KEYS = {
    "distribution",
    "percent_vacancies",
    "count",
    "min_size",
    "max_size",
    "avoid_outer_edge",
    "probe_radius_angstrom",
}


def _sample_hole_size(distribution, min_size, max_size, available_count, rng):
    cap = max(2, available_count // 5)
    if min_size is not None and max_size is not None:
        lo = min_size
        hi = min(max_size, max(1, available_count - 4))
        return rng.randint(lo, max(lo, hi))
    if distribution == "small":
        lo = min_size if min_size is not None else 1
        hi = max_size if max_size is not None else 4
        hi = min(hi, max(1, available_count // 3))
        lo = min(lo, hi)
        return rng.randint(lo, max(lo, hi))
    if distribution == "medium":
        lo = min_size if min_size is not None else 5
        hi = max_size if max_size is not None else 12
        hi = min(hi, max(2, available_count // 4))
        lo = min(lo, hi)
        return rng.randint(lo, max(lo, hi))
    if distribution == "large":
        lo = min_size if min_size is not None else min(13, cap)
        hi = max_size if max_size is not None else min(30, cap)
        lo = min(lo, hi)
        return rng.randint(lo, max(lo, hi))
    if distribution == "mixed":
        if min_size is not None or max_size is not None:
            lo = min_size if min_size is not None else 1
            hi = max_size if max_size is not None else min(30, cap)
            lo = min(lo, hi)
            return rng.randint(lo, max(lo, hi))
        tier = rng.random()
        if tier < 0.55:
            return rng.randint(1, min(4, max(1, available_count // 3)))
        if tier < 0.85:
            lo_m = min(5, max(2, available_count // 4))
            hi_m = min(12, max(2, available_count // 4))
            return rng.randint(lo_m, max(lo_m, hi_m))
        lo_l = min(13, cap)
        hi_l = min(30, cap)
        return rng.randint(lo_l, max(lo_l, hi_l))
    lo = min_size if min_size is not None else 1
    hi = max_size if max_size is not None else min(30, cap)
    lo = min(lo, hi)
    return rng.randint(lo, max(lo, hi))


def _estimate_pore_accessible_area(
    model,
    orig_coords: dict[int, np.ndarray],
    surviving_ids: set[int],
    removed_ids: set[int],
    *,
    bond_length: float = 1.42,
    probe_radius: float = 1.4,
    grid_step: float = 0.25,
) -> float:
    """Deterministic 2D grid estimate (0.25 A) of probe-accessible pore area in A^2."""
    if not removed_ids or not surviving_ids:
        return 0.0
    pbc = model["pbc"]
    cell = np.asarray(model["cell"], dtype=float)
    inv_cell = np.linalg.inv(cell) if any(pbc) else None
    has_pbc = any(pbc)
    r_excl = 1.70 + float(probe_radius)
    r_pad = float(bond_length) + 0.05

    lx = float(cell[0, 0])
    ly = float(cell[1, 1])
    if pbc[0] and lx > 0.0:
        nx_grid = max(1, int(round(lx / grid_step)))
        step_x = lx / nx_grid
    else:
        nx_grid = None
        step_x = grid_step
    if pbc[1] and ly > 0.0:
        ny_grid = max(1, int(round(ly / grid_step)))
        step_y = ly / ny_grid
    else:
        ny_grid = None
        step_y = grid_step

    r_pad_sq = r_pad * r_pad
    grid_points: set[tuple[int, int]] = set()
    for cid in sorted(removed_ids):
        pos = orig_coords[cid]
        rx, ry = float(pos[0]), float(pos[1])
        ix_lo = math.floor((rx - r_pad) / step_x)
        ix_hi = math.ceil((rx + r_pad) / step_x)
        iy_lo = math.floor((ry - r_pad) / step_y)
        iy_hi = math.ceil((ry + r_pad) / step_y)
        for ix in range(ix_lo, ix_hi + 1):
            px = ix * step_x
            dx = px - rx
            for iy in range(iy_lo, iy_hi + 1):
                py = iy * step_y
                dy = py - ry
                if dx * dx + dy * dy <= r_pad_sq:
                    gx = ix % nx_grid if nx_grid is not None else ix
                    gy = iy % ny_grid if ny_grid is not None else iy
                    grid_points.add((gx, gy))

    if not grid_points:
        return 0.0

    pts = np.asarray([[gx * step_x, gy * step_y, 0.0] for gx, gy in sorted(grid_points)], dtype=float)
    surv_coords = np.asarray([orig_coords[cid] for cid in sorted(surviving_ids)], dtype=float)

    if len(surv_coords) <= 256:
        # Vectorized chunked distance check
        accessible = 0
        chunk = 1024
        for start in range(0, len(pts), chunk):
            sub = pts[start : start + chunk]
            delta = surv_coords[None, :, :] - sub[:, None, :]
            if has_pbc:
                delta = _minimum_image_delta(delta, pbc, cell, inv_cell)
            dists = np.linalg.norm(delta[..., :2], axis=2)
            accessible += int(np.sum(np.min(dists, axis=1) >= r_excl))
        return float(accessible * step_x * step_y)

    grid = _SpatialGrid(
        surv_coords.tolist(),
        list(range(len(surv_coords))),
        ["C"] * len(surv_coords),
        pbc,
        cell,
        inv_cell,
        cell_size=max(4.0, r_excl + 0.5),
    )
    accessible = 0
    for pt in pts:
        near_coords, _, _ = grid.query_radius(pt, r_excl)
        if len(near_coords) == 0:
            accessible += 1
        else:
            delta = near_coords - pt
            if has_pbc:
                delta = _minimum_image_delta(delta, pbc, cell, inv_cell)
            if float(np.min(np.linalg.norm(delta[..., :2], axis=1))) >= r_excl:
                accessible += 1
    return float(accessible * step_x * step_y)


def _grow_and_validate_profile_hole(
    model,
    ns,
    coords,
    available,
    target_size,
    min_required_size,
    avoid_outer_edge,
    rng,
    *,
    max_delta_allowed=None,
    outer_edge_ids=None,
):
    forbidden = (outer_edge_ids or set()) if avoid_outer_edge else set()
    eligible_pool = available - forbidden
    if len(eligible_pool) < min_required_size:
        return None, 0, 0, 0

    cand_deg = {
        cid: sum(1 for nbr, _, _ in ns[cid] if nbr in available) for cid in eligible_pool
    }
    if avoid_outer_edge:
        deg3 = [cid for cid in sorted(eligible_pool) if cand_deg[cid] == 3]
        deep = [
            cid
            for cid in deg3
            if all(cand_deg.get(nbr, 0) == 3 for nbr, _, _ in ns[cid])
        ]
        deep_set = set(deep)
        deg3_edge = [cid for cid in deg3 if cid not in deep_set]
        deg3_set = set(deg3)
        rim = [cid for cid in sorted(eligible_pool) if cid not in deg3_set]
        rng.shuffle(deep)
        rng.shuffle(deg3_edge)
        rng.shuffle(rim)
        seed_candidates = deep + deg3_edge + rim
    else:
        seed_candidates = sorted(eligible_pool)
        rng.shuffle(seed_candidates)

    pbc = model["pbc"]
    max_tries = min(len(seed_candidates), 80)
    fallback_best = None
    for attempt_idx, start in enumerate(seed_candidates[:max_tries]):
        cluster = {start}
        start_pos = coords[start]
        while len(cluster) < target_size:
            frontier = sorted(
                {n[0] for i in cluster for n in ns[i]} & eligible_pool - cluster
            )
            if not frontier:
                break
            best_cand = None
            best_score = None
            for cand in frontier:
                in_cluster = sum(1 for n, _, _ in ns[cand] if n in cluster)
                delta = coords[cand] - start_pos
                if any(pbc):
                    delta = _minimum_image(delta, model)
                dist = float(np.linalg.norm(delta[:2]))
                edge_pen = 1.5 if (avoid_outer_edge and cand_deg[cand] < 3) else 0.0
                score = (
                    round(dist - 0.85 * in_cluster + edge_pen, 6),
                    rng.random(),
                )
                if best_score is None or score < best_score:
                    best_score = score
                    best_cand = cand
            cluster.add(best_cand)

        if len(cluster) < min_required_size:
            continue

        rem = available - cluster
        pruned = set()
        while rem:
            tips = {
                cid
                for cid in rem
                if sum(1 for nbr, _, _ in ns[cid] if nbr in rem) < 2
            }
            if not tips:
                break
            pruned.update(tips)
            rem -= tips

        if forbidden and (pruned & forbidden):
            continue

        if len(rem) < 4:
            continue

        visited = set()
        stack = [next(iter(rem))]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            for nbr, _, _ in ns[curr]:
                if nbr in rem and nbr not in visited:
                    stack.append(nbr)
        if len(visited) != len(rem):
            continue

        total_hole_removed = len(cluster) + len(pruned)
        if max_delta_allowed is not None and total_hole_removed > max_delta_allowed:
            continue
        if not pruned:
            return cluster, total_hole_removed, len(cluster), 0, cluster
        if fallback_best is None or len(pruned) < fallback_best[0]:
            fallback_best = (
                len(pruned),
                cluster | pruned,
                total_hole_removed,
                len(cluster),
                len(pruned),
                cluster,
            )
        if attempt_idx >= 15:
            return (
                fallback_best[1],
                fallback_best[2],
                fallback_best[3],
                fallback_best[4],
                fallback_best[5],
            )
    if fallback_best is not None:
        return (
            fallback_best[1],
            fallback_best[2],
            fallback_best[3],
            fallback_best[4],
            fallback_best[5],
        )
    return None, 0, 0, 0, set()


def _analyze_physical_pores(
    deleted_ids: set[int],
    trimmed_initial_ids: set[int],
    orig_ns: dict[int, list[tuple[int, np.ndarray, float]]],
    orig_outer_edge_ids: set[int],
) -> dict:
    """Distinguish deleted-atom clusters from physical openings in the surviving sheet."""
    if not deleted_ids:
        return {
            "deleted_atom_clusters": 0,
            "final_pore_count": 0,
            "internal_pores": 0,
            "outer_edge_notches": 0,
            "final_pore_sizes": [],
            "opening_ring_sizes": [],
        }

    init_adj: dict[int, set[int]] = {
        cid: {nbr for nbr, _, _ in orig_ns[cid] if nbr in trimmed_initial_ids}
        for cid in trimmed_initial_ids
    }
    init_adj_img: dict[int, list[tuple[int, tuple[int, int, int]]]] = {
        cid: [
            (nbr, (int(img[0]), int(img[1]), int(img[2])))
            for nbr, img, _ in orig_ns[cid]
            if nbr in trimmed_initial_ids
        ]
        for cid in trimmed_initial_ids
    }

    # 1. Connected components of deleted carbons on the trimmed initial scaffold graph
    deleted_clusters: list[set[int]] = []
    unseen = set(deleted_ids)
    while unseen:
        seed_id = min(unseen)
        unseen.remove(seed_id)
        comp = {seed_id}
        stack = [seed_id]
        while stack:
            curr = stack.pop()
            for nbr in init_adj[curr]:
                if nbr in unseen:
                    unseen.remove(nbr)
                    comp.add(nbr)
                    stack.append(nbr)
        deleted_clusters.append(comp)

    # 2. Physical openings: merge deleted atoms that are adjacent OR share a true 6-ring
    parent: dict[int, int] = {cid: cid for cid in deleted_ids}

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x: int, y: int) -> None:
        rx, ry = _find(x), _find(y)
        if rx != ry:
            if rx < ry:
                parent[ry] = rx
            else:
                parent[rx] = ry

    for u in sorted(deleted_ids):
        for v in init_adj[u]:
            if v in deleted_ids:
                _union(u, v)

    incident_rings: list[tuple[int, tuple[int, ...]]] = []
    for u in sorted(deleted_ids):
        nbrs_u = init_adj_img[u]
        for idx_a in range(len(nbrs_u)):
            a, img_ua = nbrs_u[idx_a]
            for idx_b in range(idx_a + 1, len(nbrs_u)):
                b, img_ub = nbrs_u[idx_b]
                for x, img_ax in init_adj_img[a]:
                    if x == u or x == b:
                        continue
                    for y, img_by in init_adj_img[b]:
                        if y in (u, a, x):
                            continue
                        for w, img_xw in init_adj_img[x]:
                            if w in (u, a, b, y):
                                continue
                            for y_nbr, img_yw in init_adj_img[y]:
                                if y_nbr != w:
                                    continue
                                if (
                                    img_ua[0] + img_ax[0] + img_xw[0] - img_yw[0] - img_by[0] - img_ub[0] == 0
                                    and img_ua[1] + img_ax[1] + img_xw[1] - img_yw[1] - img_by[1] - img_ub[1] == 0
                                    and img_ua[2] + img_ax[2] + img_xw[2] - img_yw[2] - img_by[2] - img_ub[2] == 0
                                ):
                                    ring = (u, a, x, w, y, b)
                                    for v in ring:
                                        if v in deleted_ids:
                                            _union(u, v)
                                    incident_rings.append((u, ring))

    openings_by_root: dict[int, set[int]] = {}
    for u in deleted_ids:
        openings_by_root.setdefault(_find(u), set()).add(u)
    ring_verts_by_root: dict[int, set[int]] = {root: set() for root in openings_by_root}
    for u, ring in incident_rings:
        ring_verts_by_root[_find(u)].update(ring)

    physical_openings = list(openings_by_root.values())
    final_pore_sizes = sorted((len(c) for c in physical_openings), reverse=True)
    opening_ring_sizes = sorted(
        (len(ring_verts_by_root[root] - comp) for root, comp in openings_by_root.items()),
        reverse=True,
    )
    outer_notch_count = sum(1 for c in physical_openings if c & orig_outer_edge_ids)
    internal_pore_count = len(physical_openings) - outer_notch_count
    return {
        "deleted_atom_clusters": len(deleted_clusters),
        "final_pore_count": len(physical_openings),
        "internal_pores": internal_pore_count,
        "outer_edge_notches": outer_notch_count,
        "final_pore_sizes": final_pore_sizes,
        "opening_ring_sizes": opening_ring_sizes,
    }


def _apply_holes(model, recipe, rng, *, bond_length: float = 1.42):
    removed = set()
    ids = {a["id"] for a in model["atoms"]}

    # Snapshot initial scaffold topology, coordinates, and outer edge ONCE
    # before applying any vacancies, legacy holes, or profile holes.
    orig_ns = adjacency(model)
    orig_coords = _coordinates(model)
    trimmed_initial_ids = set(ids)
    while trimmed_initial_ids:
        init_tips = {
            cid
            for cid in trimmed_initial_ids
            if sum(1 for nbr, _, _ in orig_ns[cid] if nbr in trimmed_initial_ids) < 2
        }
        if not init_tips:
            break
        trimmed_initial_ids -= init_tips
    initial_carbons = len(trimmed_initial_ids)
    pbc_flags = model["pbc"]
    if not (pbc_flags[0] and pbc_flags[1]):
        orig_outer_edge_ids = {
            cid
            for cid in trimmed_initial_ids
            if sum(1 for nbr, _, _ in orig_ns[cid] if nbr in trimmed_initial_ids) < 3
        }
    else:
        orig_outer_edge_ids = set()

    direct_removed_ids: set[int] = set()
    vacancies = recipe.get("vacancies", 0)
    if type(vacancies) is int:
        _integer(vacancies, "vacancies")
        if vacancies >= len(ids):
            raise ValueError("Too many vacancies")
        picked = rng.sample(sorted(ids), vacancies)
        removed.update(picked)
        direct_removed_ids.update(picked)
    elif isinstance(vacancies, list):
        if any(type(i) is not int or i not in ids for i in vacancies) or len(set(vacancies)) != len(
            vacancies
        ):
            raise ValueError("Vacancies must list distinct existing carbon IDs")
        removed.update(vacancies)
        direct_removed_ids.update(vacancies)
    else:
        raise ValueError("vacancies must be a count or a list of original carbon IDs")

    holes_raw = recipe.get("holes", [])
    if isinstance(holes_raw, dict):
        if not holes_raw:
            raise ValueError("Hole profile object must not be empty")
        holes = [holes_raw]
    elif isinstance(holes_raw, list) and all(isinstance(h, dict) for h in holes_raw):
        holes = holes_raw
    else:
        raise ValueError("holes must be a hole profile object or a list of hole specification objects")

    used_profile = False
    pre_tips_done = False
    legacy_hole_events = 0
    profile_distributions = []
    profile_req_pvs = []
    profile_req_counts = []
    profile_hole_sizes = []
    profile_cluster_sizes = []
    profile_pruned_counts = []
    probe_radius_ang = 1.4

    # Canonicalize mixed hole lists so deterministic explicit holes (center/radius or
    # explicit seed_atom) always run before random-seed legacy holes and stochastic
    # profile holes. This preserves relative order within each class while guaranteeing
    # order-invariance when combining explicit and profile hole specs.
    def _hole_priority(h: dict) -> int:
        if set(h) & _HOLE_PROFILE_KEYS:
            return 2
        if ("center" in h and "radius" in h) or ("size" in h and "seed_atom" in h):
            return 0
        return 1

    holes = sorted(holes, key=_hole_priority)

    for hole in holes:
        available = ids - removed
        is_profile = bool(set(hole) & _HOLE_PROFILE_KEYS)
        if is_profile:
            extra = set(hole) - _HOLE_PROFILE_KEYS
            if extra:
                raise ValueError(f"Unknown hole profile keys: {sorted(extra)}")
            distribution = hole.get("distribution", "mixed")
            if distribution not in ("small", "medium", "large", "mixed", "random"):
                raise ValueError(
                    "holes.distribution must be 'small', 'medium', 'large', 'mixed', or 'random'"
                )
            pv = hole.get("percent_vacancies")
            if pv is not None:
                pv = _finite_real(pv, "holes.percent_vacancies")
                if not 0.0 < pv <= 50.0:
                    raise ValueError("holes.percent_vacancies must be in (0, 50]")
            req_count = hole.get("count")
            if req_count is not None:
                req_count = _integer(req_count, "holes.count", 1)
            if pv is None and req_count is None:
                raise ValueError("Hole profile requires percent_vacancies and/or count")
            min_size = hole.get("min_size")
            max_size = hole.get("max_size")
            if min_size is not None:
                min_size = _integer(min_size, "holes.min_size", 1)
            if max_size is not None:
                max_size = _integer(max_size, "holes.max_size", 1)
            if min_size is not None and max_size is not None and min_size > max_size:
                raise ValueError("holes.min_size must be <= holes.max_size")
            avoid_outer_edge = hole.get("avoid_outer_edge", True)
            if type(avoid_outer_edge) is not bool:
                raise ValueError("holes.avoid_outer_edge must be a boolean")
            if "probe_radius_angstrom" in hole and hole["probe_radius_angstrom"] is not None:
                probe_radius_ang = _finite_real(
                    hole["probe_radius_angstrom"], "holes.probe_radius_angstrom"
                )
                if probe_radius_ang <= 0.0:
                    raise ValueError("holes.probe_radius_angstrom must be positive")

            if not pre_tips_done:
                # Trim pre-existing dangling edge tips from _tile (and any prior
                # vacancies/legacy holes) using the initial snapshot adjacency.
                pre_rem = ids - removed
                while pre_rem:
                    pre_tips = {
                        cid
                        for cid in pre_rem
                        if sum(1 for nbr, _, _ in orig_ns[cid] if nbr in pre_rem) < 2
                    }
                    if not pre_tips:
                        break
                    removed.update(pre_tips)
                    pre_rem -= pre_tips
                pre_tips_done = True

            used_profile = True
            profile_distributions.append(distribution)
            profile_req_pvs.append(pv)
            profile_req_counts.append(req_count)

            target_removals = (
                max(1, int(round(initial_carbons * pv / 100.0))) if pv is not None else None
            )

            available_now = ids - removed
            max_removable = len(available_now) - 4
            if avoid_outer_edge:
                max_removable = min(max_removable, len(available_now - orig_outer_edge_ids))
            if target_removals is not None:
                if min_size is not None and target_removals < min_size:
                    raise ValueError(
                        "Contradictory hole constraints: holes.min_size exceeds vacancy budget implied by holes.percent_vacancies"
                    )
                if req_count is not None:
                    if min_size is not None and req_count * min_size > target_removals:
                        raise ValueError(
                            "Contradictory hole constraints: holes.count * holes.min_size exceeds vacancy budget implied by holes.percent_vacancies"
                        )
                    if max_size is not None and req_count * max_size < target_removals:
                        raise ValueError(
                            "Contradictory hole constraints: holes.count * holes.max_size is below vacancy budget implied by holes.percent_vacancies"
                        )
            if req_count is not None:
                min_each = min_size if min_size is not None else 1
                if req_count * min_each > max_removable:
                    raise ValueError(
                        "Contradictory hole constraints: requested count * min_size exceeds removable interior carbons"
                    )

            removed_here = 0
            holes_here = 0

            while True:
                if req_count is not None and holes_here >= req_count:
                    break
                if req_count is None and target_removals is not None and removed_here >= target_removals:
                    break

                available = ids - removed
                if len(available) < 8:
                    if req_count is not None or holes_here == 0:
                        raise ValueError(
                            "Cannot create requested hole profile while keeping sheet connected and non-empty"
                        )
                    break

                sampled = _sample_hole_size(
                    distribution, min_size, max_size, len(available), rng
                )
                max_delta_allowed = None
                min_req = min_size if min_size is not None else 1

                if target_removals is not None and req_count is None:
                    remaining_budget = target_removals - removed_here
                    if holes_here > 0 and remaining_budget <= max(1, int(round(0.10 * target_removals))):
                        break
                    if holes_here > 0 and min_size is not None and remaining_budget < min_size:
                        break
                    if holes_here == 0 and target_removals <= 6:
                        budget_cap = max(target_removals + 2, int(round(1.5 * target_removals)))
                    else:
                        budget_cap = remaining_budget + max(3, int(round(0.25 * target_removals)))
                    if sampled > budget_cap:
                        if holes_here > 0 and min_size is not None and min_size > budget_cap:
                            break
                        target_size = max(min_req, remaining_budget)
                    else:
                        target_size = sampled
                    if min_size is not None and target_size < min_size:
                        if min_size <= budget_cap:
                            target_size = min_size
                        elif holes_here > 0:
                            break
                        else:
                            raise ValueError(
                                "Contradictory hole constraints: holes.min_size cannot fit within vacancy budget"
                            )
                    if max_size is not None and target_size > max_size:
                        target_size = max_size
                    max_delta_allowed = budget_cap
                elif target_removals is not None and req_count is not None:
                    remaining_budget = max(1, target_removals - removed_here)
                    holes_left = req_count - holes_here
                    ideal_each = max(1, int(round(remaining_budget / holes_left)))
                    lo_bound = min_size if min_size is not None else 1
                    hi_bound = max_size if max_size is not None else max(lo_bound, sampled)
                    target_size = max(lo_bound, min(hi_bound, ideal_each))
                else:
                    target_size = sampled

                hole_atoms, hole_sz, cluster_sz, pruned_sz, cluster_atoms = _grow_and_validate_profile_hole(
                    model,
                    orig_ns,
                    orig_coords,
                    available,
                    target_size,
                    min_req,
                    avoid_outer_edge,
                    rng,
                    max_delta_allowed=max_delta_allowed,
                    outer_edge_ids=orig_outer_edge_ids,
                )
                if hole_atoms is None and max_delta_allowed is not None and target_size > min_req:
                    # Retry with tighter cluster size to absorb tip-pruning overhead
                    hole_atoms, hole_sz, cluster_sz, pruned_sz, cluster_atoms = _grow_and_validate_profile_hole(
                        model,
                        orig_ns,
                        orig_coords,
                        available,
                        max(min_req, target_size - 2),
                        min_req,
                        avoid_outer_edge,
                        rng,
                        max_delta_allowed=max_delta_allowed,
                        outer_edge_ids=orig_outer_edge_ids,
                    )
                if hole_atoms is None:
                    raise ValueError(
                        "Cannot create requested hole profile while keeping sheet connected and non-empty"
                    )
                removed.update(hole_atoms)
                direct_removed_ids.update(cluster_atoms)
                removed_here += hole_sz
                holes_here += 1
                profile_hole_sizes.append(hole_sz)
                profile_cluster_sizes.append(cluster_sz)
                profile_pruned_counts.append(pruned_sz)

        elif "center" in hole and "radius" in hole:
            center = np.array(hole["center"], dtype=float)
            radius = float(hole["radius"])
            if (
                center.shape != (2,)
                or not np.isfinite(center).all()
                or not math.isfinite(radius)
                or radius <= 0
            ):
                raise ValueError("Hole requires finite xy center and positive radius")
            legacy_hits = set()
            for atom in model["atoms"]:
                aid = atom["id"]
                delta = np.array(atom["position"]) - np.r_[center, 0]
                if np.linalg.norm(_minimum_image(delta, model)[:2]) <= radius:
                    removed.add(aid)
                    if aid in available:
                        legacy_hits.add(aid)
            direct_removed_ids.update(legacy_hits)
            legacy_hole_events += 1
        elif "size" in hole:
            size = _integer(hole["size"], "hole size", 1)
            if size >= len(available):
                raise ValueError("Cluster hole consumes the sheet")
            start = hole.get("seed_atom", rng.choice(sorted(available)))
            if start not in available:
                raise ValueError("Hole seed is not an available carbon")
            cluster = {start}
            while len(cluster) < size:
                frontier = sorted({n[0] for i in cluster for n in orig_ns[i]} & available - cluster)
                if not frontier:
                    raise ValueError("Cannot grow connected hole of requested size")
                cluster.add(rng.choice(frontier))
            removed.update(cluster)
            direct_removed_ids.update(cluster)
            legacy_hole_events += 1
        else:
            raise ValueError("Hole requires center/radius, size/seed_atom, or a hole profile")
    _remove_atoms(model, removed)
    extra_pruned = _prune_tips(model)
    removed.update(extra_pruned)
    if not model["atoms"] or not _connected(model):
        raise ValueError("Vacancies/holes leave an empty or disconnected sheet")
    model["metadata"]["removed_carbon_ids"] = sorted(removed)
    model["metadata"]["requested_vacancies"] = copy.deepcopy(vacancies)
    if used_profile:
        req_pv = profile_req_pvs[0] if len(profile_req_pvs) == 1 else (
            sum(x for x in profile_req_pvs if x is not None) if any(x is not None for x in profile_req_pvs) else None
        )
        req_cnt = profile_req_counts[0] if len(profile_req_counts) == 1 else (
            sum(x for x in profile_req_counts if x is not None) if any(x is not None for x in profile_req_counts) else None
        )
        surviving_ids = {a["id"] for a in model["atoms"]}
        deleted_from_trimmed = trimmed_initial_ids - surviving_ids
        total_removed = len(deleted_from_trimmed)
        direct_in_trimmed = direct_removed_ids & deleted_from_trimmed
        total_pruned_tips = len(deleted_from_trimmed - direct_in_trimmed)

        pore_stats = _analyze_physical_pores(
            deleted_from_trimmed,
            trimmed_initial_ids,
            orig_ns,
            orig_outer_edge_ids,
        )
        pore_area_a2 = _estimate_pore_accessible_area(
            model,
            orig_coords,
            surviving_ids,
            deleted_from_trimmed,
            bond_length=bond_length,
            probe_radius=probe_radius_ang,
            grid_step=0.25,
        )

        excavation_events = len(profile_hole_sizes) + legacy_hole_events
        model["metadata"]["hole_summary"] = {
            "distribution": (
                profile_distributions[0]
                if len(profile_distributions) == 1
                else profile_distributions
            ),
            "requested_percent_vacancies": req_pv,
            "requested_hole_count": req_cnt,
            "holes_created": len(profile_hole_sizes),
            "hole_sizes": profile_hole_sizes,
            "excavation_events": excavation_events,
            "excavation_sizes": profile_cluster_sizes,
            "pruned_tip_carbons": total_pruned_tips,
            "deleted_atom_clusters": pore_stats["deleted_atom_clusters"],
            "final_pore_count": pore_stats["final_pore_count"],
            "internal_pores": pore_stats["internal_pores"],
            "outer_edge_notches": pore_stats["outer_edge_notches"],
            "final_pore_components": pore_stats["final_pore_count"],
            "final_pore_sizes": pore_stats["final_pore_sizes"],
            "opening_ring_sizes": pore_stats["opening_ring_sizes"],
            "internal_pore_count": pore_stats["internal_pores"],
            "outer_edge_notch_count": pore_stats["outer_edge_notches"],
            "probe_radius_angstrom": float(probe_radius_ang),
            "pore_accessible_area_angstrom2": round(pore_area_a2, 4),
            "pore_accessible_area_nm2": round(pore_area_a2 / 100.0, 6),
            "initial_carbons": initial_carbons,
            "removed_carbons": total_removed,
            "actual_vacancy_percent": 100.0 * total_removed / initial_carbons,
        }


def _hopcroft_karp_matching(left_nodes: list[int], adj_map: dict[int, list[int]]) -> list[tuple[int, int]]:
    """Iterative O(E sqrt(V)) Hopcroft-Karp maximum bipartite matching (zero recursion)."""
    if not left_nodes:
        return []
    inf = len(left_nodes) + 10
    pair_u: dict[int, int | None] = dict.fromkeys(left_nodes)
    pair_v: dict[int, int] = {}
    dist: dict[int, int] = {}

    while True:
        queue: deque[int] = deque()
        for u in left_nodes:
            if pair_u[u] is None:
                dist[u] = 0
                queue.append(u)
            else:
                dist[u] = inf
        dist_nil = inf
        while queue:
            u = queue.popleft()
            if dist[u] < dist_nil:
                d_next = dist[u] + 1
                for v in adj_map[u]:
                    u2 = pair_v.get(v)
                    if u2 is None:
                        if d_next < dist_nil:
                            dist_nil = d_next
                    elif dist[u2] == inf:
                        dist[u2] = d_next
                        queue.append(u2)
        if dist_nil == inf:
            break

        ptr = dict.fromkeys(left_nodes, 0)
        for u0 in left_nodes:
            if pair_u[u0] is not None:
                continue
            stack_u = [u0]
            stack_v: list[int | None] = [None]
            while stack_u:
                u = stack_u[-1]
                nbrs = adj_map[u]
                idx = ptr[u]
                if idx >= len(nbrs):
                    dist[u] = inf
                    stack_u.pop()
                    stack_v.pop()
                    continue
                v = nbrs[idx]
                ptr[u] = idx + 1
                u2 = pair_v.get(v)
                if u2 is None:
                    if dist[u] + 1 == dist_nil:
                        pair_v[v] = u
                        pair_u[u] = v
                        for k in range(len(stack_u) - 1, 0, -1):
                            v_edge = stack_v[k]
                            u_prev = stack_u[k - 1]
                            pair_v[v_edge] = u_prev
                            pair_u[u_prev] = v_edge
                        for u_path in stack_u:
                            dist[u_path] = inf
                        break
                elif dist[u2] == dist[u] + 1:
                    stack_u.append(u2)
                    stack_v.append(v)

    return [(u, v) for u, v in pair_u.items() if v is not None]


def _iterative_kuhn_matching(
    ns, left_nodes: list[int], allowed_right: set[int], rng
) -> list[tuple[int, int]]:
    """Zero-recursion iterative Kuhn augmenting-path matching preserving exact RNG sequence."""
    match: dict[int, int] = {}
    for start in left_nodes:
        seen: set[int] = set()
        opts0 = [n[0] for n in ns[start] if n[0] in allowed_right]
        rng.shuffle(opts0)
        stack: list[list] = [[start, opts0, 0]]
        while stack:
            u, opts, idx = stack[-1]
            if idx >= len(opts):
                stack.pop()
                continue
            j = opts[idx]
            stack[-1][2] = idx + 1
            if j in seen:
                continue
            seen.add(j)
            if j not in match:
                match[j] = u
                for k in range(len(stack) - 2, -1, -1):
                    p_u, p_opts, p_next = stack[k]
                    match[p_opts[p_next - 1]] = p_u
                break
            next_u = match[j]
            opts_n = [n[0] for n in ns[next_u] if n[0] in allowed_right]
            rng.shuffle(opts_n)
            stack.append([next_u, opts_n, 0])
    return [(i, j) for j, i in match.items()]


def _max_matching(model, candidates, rng):
    """Exact cardinality matching on the bipartite honeycomb, randomized ties (iterative)."""
    cand_set = candidates if isinstance(candidates, set) else set(candidates)
    if len(cand_set) < 2:
        return []
    by_id = {a["id"]: a for a in model["atoms"]}
    ns = adjacency(model)
    left = [i for i in candidates if by_id[i]["metadata"]["sublattice"] == 0]
    rng.shuffle(left)
    if len(left) <= 1000:
        return _iterative_kuhn_matching(ns, left, cand_set, rng)
    adj_map: dict[int, list[int]] = {}
    for u in left:
        options = [n[0] for n in ns[u] if n[0] in cand_set]
        rng.shuffle(options)
        adj_map[u] = options
    return _hopcroft_karp_matching(left, adj_map)


def _bipartite_matching(ns, left_nodes, right_set, rng):
    """Exact cardinality matching between disjoint left_nodes and right_set (iterative)."""
    left = list(left_nodes)
    if not left or not right_set:
        return []
    rng.shuffle(left)
    if len(left) <= 1000:
        return _iterative_kuhn_matching(ns, left, right_set, rng)
    adj_map: dict[int, list[int]] = {}
    for u in left:
        options = [n[0] for n in ns[u] if n[0] in right_set]
        rng.shuffle(options)
        adj_map[u] = options
    return _hopcroft_karp_matching(left, adj_map)


def _set_site(
    model,
    identifier,
    name,
    resname,
    resid,
    group_id,
    group,
    sp3=False,
    *,
    by_id=None,
    bond_by_atom=None,
):
    atom = by_id[identifier] if by_id is not None else next(
        a for a in model["atoms"] if a["id"] == identifier
    )
    atom.update(name=name, resname=resname, resid=resid)
    atom["metadata"].update(group_id=group_id, group=group)
    if sp3:
        atom["metadata"]["hybridization"] = "sp3"
    # Saturated/edge carbonyl sites cannot retain aromatic bond order.
    if sp3 or group == "carbonyl":
        bonds = bond_by_atom[identifier] if bond_by_atom is not None else model["bonds"]
        for bond in bonds:
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


def _orthonormal_edge_frame(outward):
    """Construct a strictly right-handed orthonormal triad {u_out, u_perp, u_n}."""
    u_out = np.asarray(outward, dtype=float)
    out_norm = float(np.linalg.norm(u_out))
    if out_norm < 1e-8:
        raise ValueError("Undefined edge outward direction")
    if abs(u_out[2]) < 1e-15:
        if abs(out_norm - 1.0) > 1e-12:
            u_out = u_out / out_norm
        u_n = np.array([0.0, 0.0, 1.0], dtype=float)
        u_perp = np.cross(u_n, u_out)
        u_perp = u_perp / float(np.linalg.norm(u_perp))
        return u_out, u_perp, u_n
    u_out = u_out / out_norm
    z_ref = np.array([0.0, 0.0, 1.0], dtype=float)
    perp_raw = np.cross(z_ref, u_out)
    perp_norm = float(np.linalg.norm(perp_raw))
    if perp_norm < 1e-8:
        perp_raw = np.cross(np.array([0.0, 1.0, 0.0], dtype=float), u_out)
        perp_norm = float(np.linalg.norm(perp_raw))
    u_perp = perp_raw / perp_norm
    u_n = np.cross(u_out, u_perp)
    u_n = u_n / float(np.linalg.norm(u_n))
    u_perp = np.cross(u_n, u_out)
    u_perp = u_perp / float(np.linalg.norm(u_perp))
    return u_out, u_perp, u_n


def _local_candidates_check(
    model,
    origin,
    existing_coords,
    existing_ids,
    existing_elements,
    *,
    pbc,
    cell,
    inv_cell,
    radius=6.5,
):
    if not any(pbc):
        delta = existing_coords - origin
    else:
        fractional = (existing_coords - origin) @ inv_cell
        for k, periodic in enumerate(pbc):
            if periodic:
                fractional[..., k] -= np.rint(fractional[..., k])
        delta = fractional @ cell
    dists = np.linalg.norm(delta, axis=1)
    mask = dists <= radius
    return (
        existing_coords[mask],
        existing_ids[mask],
        existing_elements[mask],
    )


def _score_candidate_atoms(
    cand_atoms,
    local_coords,
    local_ids,
    local_elements,
    *,
    pbc,
    cell,
    inv_cell,
    virtual_h_coords=None,
):
    """Return (valid, min_14_plus, min_non12) for candidate atoms against local + virtual H.

    Each entry in cand_atoms is (pos, element, ids_12_set, ids_13_set).
    Enforces hard non-1-2 cutoffs (0.90 A for H, 1.25 A heavy-heavy) and
    1-4+ cutoffs (1.20 A for X-H, 1.60 A heavy-heavy).
    """
    min_14 = float("inf")
    min_non12 = float("inf")

    for pos, elem, ids_12, ids_13 in cand_atoms:
        if len(local_coords):
            delta = local_coords - pos
            if any(pbc):
                fractional = delta @ inv_cell
                for k, periodic in enumerate(pbc):
                    if periodic:
                        fractional[..., k] -= np.rint(fractional[..., k])
                delta = fractional @ cell
            dists = np.linalg.norm(delta, axis=1)

            non12 = np.ones(len(local_ids), dtype=bool)
            if ids_12:
                non12 &= ~np.isin(local_ids, tuple(ids_12))
            if np.any(non12):
                is_h_pair = (local_elements == "H") | (elem == "H")
                hard_cutoffs = np.where(is_h_pair, 0.90, 1.25)
                if np.any(dists[non12] < hard_cutoffs[non12]):
                    return False, 0.0, 0.0

                is_14 = non12.copy()
                if ids_13:
                    is_14 &= ~np.isin(local_ids, tuple(ids_13))
                if np.any(is_14):
                    cutoffs_14 = np.where(is_h_pair, 1.20, 1.60)
                    if np.any(dists[is_14] < cutoffs_14[is_14]):
                        return False, 0.0, 0.0
                    local_min_14 = float(np.min(dists[is_14]))
                    if local_min_14 < min_14:
                        min_14 = local_min_14

                local_min_non12 = float(np.min(dists[non12]))
                if local_min_non12 < min_non12:
                    min_non12 = local_min_non12

        if virtual_h_coords is not None and len(virtual_h_coords):
            v_delta = virtual_h_coords - pos
            if any(pbc):
                v_frac = v_delta @ inv_cell
                for k, periodic in enumerate(pbc):
                    if periodic:
                        v_frac[..., k] -= np.rint(v_frac[..., k])
                v_delta = v_frac @ cell
            v_dists = np.linalg.norm(v_delta, axis=1)
            v_min = float(np.min(v_dists))
            if v_min < 1.20:
                return False, 0.0, 0.0
            if v_min < min_14:
                min_14 = v_min
            if v_min < min_non12:
                min_non12 = v_min

    return True, min_14, min_non12


def _virtual_edge_hydrogens(
    atoms,
    adj,
    coords_by_id,
    cell,
    occupied_edges,
    exclude_site,
    origin,
    model,
    *,
    spatial_grid=None,
    by_id=None,
):
    """Estimate GRT edge H positions on adjacent undecorated edge carbons."""
    excluded = {exclude_site} if isinstance(exclude_site, int) else set(exclude_site)
    positions = []
    if spatial_grid is not None and by_id is not None:
        _, cand_ids, _ = spatial_grid.query_radius(origin, 5.5)
        for aid_np in cand_ids:
            aid = int(aid_np)
            if aid in excluded or aid in occupied_edges or len(adj[aid]) != 2:
                continue
            atom = by_id[aid]
            if atom.get("metadata", {}).get("role") != "sheet_carbon":
                continue
            pos = coords_by_id[aid]
            outward = _outward(model, aid, neighbors=adj, coordinates=coords_by_id, cell=cell)
            positions.append(pos + 1.09 * outward)
    else:
        for atom in atoms:
            aid = atom["id"]
            if (
                aid in excluded
                or aid in occupied_edges
                or atom.get("metadata", {}).get("role") != "sheet_carbon"
                or len(adj[aid]) != 2
            ):
                continue
            pos = coords_by_id[aid]
            if np.linalg.norm(_minimum_image(pos - origin, model)) > 5.5:
                continue
            outward = _outward(model, aid, neighbors=adj, coordinates=coords_by_id, cell=cell)
            positions.append(pos + 1.09 * outward)
    if not positions:
        return None
    return np.asarray(positions, dtype=float)


def _place_group(
    model,
    group,
    sites,
    face,
    group_id,
    *,
    faces="both",
    adj=None,
    coords_by_id=None,
    occupied_edges=None,
    spatial_grid=None,
    by_id=None,
    bond_by_atom=None,
):
    if adj is None:
        adj = adjacency(model)
    if coords_by_id is None:
        coords_by_id = _coordinates(model)
    if occupied_edges is None:
        occupied_edges = {i for g in model["metadata"]["groups"] for i in g["site_ids"]}

    atoms = model["atoms"]
    pbc = model["pbc"]
    cell = np.asarray(model["cell"], dtype=float)
    inv_cell = np.linalg.inv(cell) if any(pbc) else None
    if by_id is None:
        by_id = {a["id"]: a for a in atoms}

    resid = group_id + 1
    origin = coords_by_id[sites[0]]
    resname = {
        "OH": "HYD",
        "epoxide": "EPX",
        "epoxide_eph": "EPH",
        "epoxide_ept": "EPT",
        "carbonyl": "CAB",
        "COOH": "CBX",
        "COO": "CBD",
    }[group]

    if spatial_grid is not None:
        local_coords, local_ids, local_elements = spatial_grid.query_radius(origin, 6.5)
    else:
        existing_coords = np.asarray([a["position"] for a in atoms], dtype=float)
        existing_ids = np.fromiter((a["id"] for a in atoms), dtype=np.int64, count=len(atoms))
        existing_elements = np.asarray([a["element"] for a in atoms])
        local_coords, local_ids, local_elements = _local_candidates_check(
            model,
            origin,
            existing_coords,
            existing_ids,
            existing_elements,
            pbc=pbc,
            cell=cell,
            inv_cell=inv_cell,
        )

    face_options = [face, -face] if faces == "both" else [face]
    chosen_face = face
    site_updates: dict[int, np.ndarray] = {}
    new_atoms_spec: list[tuple[str, np.ndarray, str, int]] = []
    bond_spec: list[tuple[str | int, str | int, float, tuple[int, int, int]]] = []
    oxygen_tag: str = ""

    if group == "OH":
        site_id = sites[0]
        site_neighbors = [nbr for nbr, _, _ in adj[site_id]]
        site_12 = set(site_neighbors)
        site_13 = {
            nbr2 for nbr in site_neighbors for nbr2, _, _ in adj[nbr] if nbr2 != site_id
        }
        # Exclude the site itself from local arrays when testing shifted site + O/H
        keep = local_ids != site_id
        l_coords = local_coords[keep]
        l_ids = local_ids[keep]
        l_elems = local_elements[keep]

        oh_theta = math.radians(108.5)
        sin_t = math.sin(oh_theta)
        neg_cos_t = -math.cos(oh_theta)

        found = False
        for face_try in face_options:
            c_pos = origin + np.array([0.0, 0.0, face_try * 0.25])
            # Verify sp3 puckering preserves physical C-C bond lengths (1.40-1.56 A)
            cc_ok = True
            for nbr, image, _ in adj[site_id]:
                d_cc = np.linalg.norm(coords_by_id[nbr] + image @ cell - c_pos)
                if not 1.40 <= d_cc <= 1.56:
                    cc_ok = False
                    break
            if not cc_ok:
                continue

            o_pos = c_pos + np.array([0.0, 0.0, face_try * 1.43])
            c_and_o_ok, _, _ = _score_candidate_atoms(
                [
                    (c_pos, "C", site_12, site_13),
                    (o_pos, "O", {site_id}, site_12),
                ],
                l_coords,
                l_ids,
                l_elems,
                pbc=pbc,
                cell=cell,
                inv_cell=inv_cell,
            )
            if not c_and_o_ok:
                continue

            best_az = None
            for az_deg in range(0, 360, 30):
                phi = math.radians(az_deg)
                h_pos = o_pos + 0.97 * np.array(
                    [sin_t * math.cos(phi), sin_t * math.sin(phi), face_try * neg_cos_t]
                )
                ok, min_14, min_non12 = _score_candidate_atoms(
                    [(h_pos, "H", set(), {site_id})],
                    l_coords,
                    l_ids,
                    l_elems,
                    pbc=pbc,
                    cell=cell,
                    inv_cell=inv_cell,
                )
                if ok:
                    score = (round(min_14, 9), round(min_non12, 9))
                    if best_az is None or score > best_az[0]:
                        best_az = (score, h_pos)
            if best_az is not None:
                chosen_face = face_try
                site_updates[site_id] = c_pos
                new_atoms_spec = [
                    ("O", o_pos, "OL", 0),
                    ("H", best_az[1], "HL", 0),
                ]
                bond_spec = [
                    (site_id, "OL", 1.0, (0, 0, 0)),
                    ("OL", "HL", 1.0, (0, 0, 0)),
                ]
                oxygen_tag = "OL"
                found = True
                break
        if not found:
            return False

    elif group == "epoxide":
        a, b = sites
        image = next(img for j, img, _ in adj[a] if j == b)
        a_neighbors = [nbr for nbr, _, _ in adj[a]]
        b_neighbors = [nbr for nbr, _, _ in adj[b]]
        a_12 = set(a_neighbors)
        b_12 = set(b_neighbors)
        a_13 = {n2 for n1 in a_neighbors for n2, _, _ in adj[n1] if n2 != a} - b_12
        b_13 = {n2 for n1 in b_neighbors for n2, _, _ in adj[n1] if n2 != b} - a_12
        o_13 = (a_12 | b_12) - {a, b}

        keep = (local_ids != a) & (local_ids != b)
        l_coords = local_coords[keep]
        l_ids = local_ids[keep]
        l_elems = local_elements[keep]

        found = False
        for face_try in face_options:
            shift = np.array([0.0, 0.0, face_try * 0.18])
            c_a = coords_by_id[a] + shift
            c_b = coords_by_id[b] + shift
            cc_ok = True
            for site_i, c_new in ((a, c_a), (b, c_b)):
                for nbr, img_n, _ in adj[site_i]:
                    other_pos = (
                        c_b if nbr == b else c_a if nbr == a else coords_by_id[nbr]
                    )
                    d_cc = np.linalg.norm(other_pos + img_n @ cell - c_new)
                    if not 1.40 <= d_cc <= 1.56:
                        cc_ok = False
                        break
                if not cc_ok:
                    break
            if not cc_ok:
                continue

            midpoint = (c_a + c_b + image @ cell) / 2.0
            half = np.linalg.norm(c_b + image @ cell - c_a) / 2.0
            if half >= 1.45:
                continue
            height = math.sqrt(1.45**2 - half**2)
            o_pos = midpoint + np.array([0.0, 0.0, face_try * height])

            ok, _, _ = _score_candidate_atoms(
                [
                    (c_a, "C", a_12, a_13),
                    (c_b, "C", b_12, b_13),
                    (o_pos, "O", {a, b}, o_13),
                ],
                l_coords,
                l_ids,
                l_elems,
                pbc=pbc,
                cell=cell,
                inv_cell=inv_cell,
            )
            if ok:
                chosen_face = face_try
                site_updates[a] = c_a
                site_updates[b] = c_b
                new_atoms_spec = [("O", o_pos, "OE", 0)]
                bond_spec = [
                    (a, "OE", 1.0, (0, 0, 0)),
                    (b, "OE", 1.0, tuple(-image)),
                ]
                oxygen_tag = "OE"
                found = True
                break
        if not found:
            return False

    elif group == "epoxide_eph":
        a, b = sites  # a = basal (CZ), b = edge (CK)
        image = next(img for j, img, _ in adj[a] if j == b)
        a_neighbors = [nbr for nbr, _, _ in adj[a]]
        b_neighbors = [nbr for nbr, _, _ in adj[b]]
        a_12 = set(a_neighbors)
        b_12 = set(b_neighbors)
        a_13 = {n2 for n1 in a_neighbors for n2, _, _ in adj[n1] if n2 != a} - b_12
        b_13 = {n2 for n1 in b_neighbors for n2, _, _ in adj[n1] if n2 != b} - a_12
        o_13 = (a_12 | b_12) - {a, b}

        keep = (local_ids != a) & (local_ids != b)
        l_coords = local_coords[keep]
        l_ids = local_ids[keep]
        l_elems = local_elements[keep]
        outward_b = _outward(model, b, neighbors=adj, coordinates=coords_by_id, cell=cell)
        u_out_b, _, u_n_b = _orthonormal_edge_frame(outward_b)
        virtual_h = _virtual_edge_hydrogens(
            atoms,
            adj,
            coords_by_id,
            cell,
            occupied_edges,
            {a, b},
            coords_by_id[b],
            model,
            spatial_grid=spatial_grid,
            by_id=by_id,
        )

        found = False
        for face_try in face_options:
            shift = np.array([0.0, 0.0, face_try * 0.18])
            c_a = coords_by_id[a] + shift
            c_b = coords_by_id[b] + shift
            cc_ok = True
            for site_i, c_new in ((a, c_a), (b, c_b)):
                for nbr, img_n, _ in adj[site_i]:
                    other_pos = (
                        c_b if nbr == b else c_a if nbr == a else coords_by_id[nbr]
                    )
                    d_cc = np.linalg.norm(other_pos + img_n @ cell - c_new)
                    if not 1.40 <= d_cc <= 1.56:
                        cc_ok = False
                        break
                if not cc_ok:
                    break
            if not cc_ok:
                continue

            midpoint = (c_a + c_b + image @ cell) / 2.0
            half = np.linalg.norm(c_b + image @ cell - c_a) / 2.0
            if half >= 1.45:
                continue
            height = math.sqrt(1.45**2 - half**2)
            o_pos = midpoint + np.array([0.0, 0.0, face_try * height])

            co_ok, _, _ = _score_candidate_atoms(
                [
                    (c_a, "C", a_12, a_13),
                    (c_b, "C", b_12, b_13),
                    (o_pos, "O", {a, b}, o_13),
                ],
                l_coords,
                l_ids,
                l_elems,
                pbc=pbc,
                cell=cell,
                inv_cell=inv_cell,
            )
            if not co_ok:
                continue

            best_hk = None
            for tilt_deg in (35, 25, 45, 20, 50):
                rad = math.radians(tilt_deg)
                h_dir = math.cos(rad) * u_out_b + (-face_try * math.sin(rad)) * u_n_b
                hk_pos = c_b + 1.09 * h_dir
                if np.linalg.norm(_minimum_image(hk_pos - o_pos, model)) < 1.20:
                    continue
                if np.linalg.norm(_minimum_image(hk_pos - (c_a - image @ cell), model)) < 1.20:
                    continue
                ok, min_14, min_non12 = _score_candidate_atoms(
                    [(hk_pos, "H", set(), b_12)],
                    l_coords,
                    l_ids,
                    l_elems,
                    pbc=pbc,
                    cell=cell,
                    inv_cell=inv_cell,
                    virtual_h_coords=virtual_h,
                )
                if not ok:
                    continue
                score = (round(min_14, 9), round(min_non12, 9))
                if tilt_deg == 35 and min_14 >= 1.80:
                    best_hk = (score, hk_pos)
                    break
                if best_hk is None or score > best_hk[0]:
                    best_hk = (score, hk_pos)
            if best_hk is not None:
                chosen_face = face_try
                site_updates[a] = c_a
                site_updates[b] = c_b
                new_atoms_spec = [("H", best_hk[1], "HK", 0), ("O", o_pos, "OE", 0)]
                bond_spec = [
                    (a, "OE", 1.0, (0, 0, 0)),
                    (b, "OE", 1.0, tuple(-image)),
                    (b, "HK", 1.0, (0, 0, 0)),
                ]
                oxygen_tag = "OE"
                found = True
                break
        if not found:
            return False

    elif group == "epoxide_ept":
        a, b = sites  # a = CP (edge), b = CK (edge)
        image = next(img for j, img, _ in adj[a] if j == b)
        a_neighbors = [nbr for nbr, _, _ in adj[a]]
        b_neighbors = [nbr for nbr, _, _ in adj[b]]
        a_12 = set(a_neighbors)
        b_12 = set(b_neighbors)
        a_13 = {n2 for n1 in a_neighbors for n2, _, _ in adj[n1] if n2 != a} - b_12
        b_13 = {n2 for n1 in b_neighbors for n2, _, _ in adj[n1] if n2 != b} - a_12
        o_13 = (a_12 | b_12) - {a, b}

        keep = (local_ids != a) & (local_ids != b)
        l_coords = local_coords[keep]
        l_ids = local_ids[keep]
        l_elems = local_elements[keep]
        outward_a = _outward(model, a, neighbors=adj, coordinates=coords_by_id, cell=cell)
        outward_b = _outward(model, b, neighbors=adj, coordinates=coords_by_id, cell=cell)
        u_out_a, _, u_n_a = _orthonormal_edge_frame(outward_a)
        u_out_b, _, u_n_b = _orthonormal_edge_frame(outward_b)
        edge_mid = (coords_by_id[a] + coords_by_id[b] + image @ cell) / 2.0
        virtual_h = _virtual_edge_hydrogens(
            atoms,
            adj,
            coords_by_id,
            cell,
            occupied_edges,
            {a, b},
            edge_mid,
            model,
            spatial_grid=spatial_grid,
            by_id=by_id,
        )

        found = False
        for face_try in face_options:
            shift = np.array([0.0, 0.0, face_try * 0.18])
            c_a = coords_by_id[a] + shift
            c_b = coords_by_id[b] + shift
            cc_ok = True
            for site_i, c_new in ((a, c_a), (b, c_b)):
                for nbr, img_n, _ in adj[site_i]:
                    other_pos = (
                        c_b if nbr == b else c_a if nbr == a else coords_by_id[nbr]
                    )
                    d_cc = np.linalg.norm(other_pos + img_n @ cell - c_new)
                    if not 1.40 <= d_cc <= 1.56:
                        cc_ok = False
                        break
                if not cc_ok:
                    break
            if not cc_ok:
                continue

            midpoint = (c_a + c_b + image @ cell) / 2.0
            half = np.linalg.norm(c_b + image @ cell - c_a) / 2.0
            if half >= 1.45:
                continue
            height = math.sqrt(1.45**2 - half**2)
            o_pos = midpoint + np.array([0.0, 0.0, face_try * height])

            co_ok, _, _ = _score_candidate_atoms(
                [
                    (c_a, "C", a_12, a_13),
                    (c_b, "C", b_12, b_13),
                    (o_pos, "O", {a, b}, o_13),
                ],
                l_coords,
                l_ids,
                l_elems,
                pbc=pbc,
                cell=cell,
                inv_cell=inv_cell,
            )
            if not co_ok:
                continue

            hp_candidates = []
            hk_candidates = []
            for tilt_deg in (35, 25, 45, 20, 50):
                rad = math.radians(tilt_deg)
                dir_a = math.cos(rad) * u_out_a + (-face_try * math.sin(rad)) * u_n_a
                hp_pos = c_a + 1.09 * dir_a
                if (
                    np.linalg.norm(_minimum_image(hp_pos - o_pos, model)) >= 1.20
                    and np.linalg.norm(_minimum_image(hp_pos - (c_b + image @ cell), model))
                    >= 1.20
                ):
                    ok_a, min14_a, min12_a = _score_candidate_atoms(
                        [(hp_pos, "H", set(), a_12)],
                        l_coords,
                        l_ids,
                        l_elems,
                        pbc=pbc,
                        cell=cell,
                        inv_cell=inv_cell,
                        virtual_h_coords=virtual_h,
                    )
                    if ok_a:
                        hp_candidates.append((tilt_deg, hp_pos, min14_a, min12_a))

                dir_b = math.cos(rad) * u_out_b + (-face_try * math.sin(rad)) * u_n_b
                hk_pos = c_b + 1.09 * dir_b
                if (
                    np.linalg.norm(_minimum_image(hk_pos - o_pos, model)) >= 1.20
                    and np.linalg.norm(_minimum_image(hk_pos - (c_a - image @ cell), model))
                    >= 1.20
                ):
                    ok_b, min14_b, min12_b = _score_candidate_atoms(
                        [(hk_pos, "H", set(), b_12)],
                        l_coords,
                        l_ids,
                        l_elems,
                        pbc=pbc,
                        cell=cell,
                        inv_cell=inv_cell,
                        virtual_h_coords=virtual_h,
                    )
                    if ok_b:
                        hk_candidates.append((tilt_deg, hk_pos, min14_b, min12_b))

            best_pair = None
            for ta, hp_pos, min14_a, min12_a in hp_candidates:
                for tb, hk_pos, min14_b, min12_b in hk_candidates:
                    d_hh = float(np.linalg.norm(_minimum_image(hp_pos - hk_pos, model)))
                    if d_hh < 1.20:
                        continue
                    score = (
                        round(min(min14_a, min14_b, d_hh), 9),
                        round(min(min12_a, min12_b, d_hh), 9),
                    )
                    if ta == 35 and tb == 35 and min(min14_a, min14_b, d_hh) >= 1.80:
                        best_pair = (score, hp_pos, hk_pos)
                        break
                    if best_pair is None or score > best_pair[0]:
                        best_pair = (score, hp_pos, hk_pos)
                if best_pair is not None and ta == 35 and best_pair[0][0] >= 1.80:
                    break
            if best_pair is not None:
                chosen_face = face_try
                site_updates[a] = c_a
                site_updates[b] = c_b
                new_atoms_spec = [
                    ("H", best_pair[1], "HP", 0),
                    ("H", best_pair[2], "HK", 0),
                    ("O", o_pos, "OE", 0),
                ]
                bond_spec = [
                    (a, "OE", 1.0, (0, 0, 0)),
                    (b, "OE", 1.0, tuple(-image)),
                    (a, "HP", 1.0, (0, 0, 0)),
                    (b, "HK", 1.0, (0, 0, 0)),
                ]
                oxygen_tag = "OE"
                found = True
                break
        if not found:
            return False

    else:
        site_id = sites[0]
        raw_direction = _outward(model, site_id, neighbors=adj, coordinates=coords_by_id, cell=cell)
        direction, in_plane_perp, z_axis = _orthonormal_edge_frame(raw_direction)
        site_13 = {nbr for nbr, _, _ in adj[site_id]}

        if group == "carbonyl":
            o_pos = origin + 1.23 * direction
            ok, _, _ = _score_candidate_atoms(
                [(o_pos, "O", {site_id}, site_13)],
                local_coords,
                local_ids,
                local_elements,
                pbc=pbc,
                cell=cell,
                inv_cell=inv_cell,
            )
            if not ok:
                return False
            new_atoms_spec = [("O", o_pos, "ON", 0)]
            bond_spec = [(site_id, "ON", 2.0, (0, 0, 0))]
            oxygen_tag = "ON"
        else:
            center = origin + 1.52 * direction
            virtual_h = _virtual_edge_hydrogens(
                atoms,
                adj,
                coords_by_id,
                cell,
                occupied_edges,
                site_id,
                origin,
                model,
                spatial_grid=spatial_grid,
                by_id=by_id,
            )
            # Check C4 once (independent of dihedral tau)
            c4_ok, c4_min14, c4_min_non12 = _score_candidate_atoms(
                [(center, "C", {site_id}, site_13)],
                local_coords,
                local_ids,
                local_elements,
                pbc=pbc,
                cell=cell,
                inv_cell=inv_cell,
                virtual_h_coords=virtual_h,
            )
            if not c4_ok:
                return False

            torsions = (0, 30, -30, 45, -45, 60, -60, 75, -75, 90)
            best_choice = None
            clean_threshold = 1.80

            if group == "COO":
                half_angle = math.radians(63.0)  # O-C-O = 126 deg
                cos_h = math.cos(half_angle)
                sin_h = math.sin(half_angle)
                for tau_deg in torsions:
                    tau = math.radians(tau_deg)
                    u_perp = math.cos(tau) * in_plane_perp + math.sin(tau) * z_axis
                    oj_pos = center + 1.26 * (cos_h * direction + sin_h * u_perp)
                    ok_pos = center + 1.26 * (cos_h * direction - sin_h * u_perp)
                    ok, min_14, min_non12 = _score_candidate_atoms(
                        [
                            (oj_pos, "O", set(), {site_id}),
                            (ok_pos, "O", set(), {site_id}),
                        ],
                        local_coords,
                        local_ids,
                        local_elements,
                        pbc=pbc,
                        cell=cell,
                        inv_cell=inv_cell,
                        virtual_h_coords=virtual_h,
                    )
                    if not ok:
                        continue
                    clearance = min(c4_min14, min_14)
                    if clearance >= clean_threshold:
                        best_choice = (clearance, oj_pos, ok_pos, None)
                        break
                    if best_choice is None or clearance > best_choice[0] + 1e-9:
                        best_choice = (clearance, oj_pos, ok_pos, None)
            else:
                # COOH: d(C=O)=1.21, d(C-OH)=1.34, d(O-H)=0.97, C-O-H=108 deg
                a_double = math.radians(58.0)
                a_single = math.radians(66.0)
                coh_rad = math.radians(108.0)
                cos_d, sin_d = math.cos(a_double), math.sin(a_double)
                cos_s, sin_s = math.cos(a_single), math.sin(a_single)
                cos_oh = math.cos(math.pi - coh_rad)
                sin_oh = math.sin(math.pi - coh_rad)

                for conformer in ("syn", "anti"):
                    sign_h = 1.0 if conformer == "syn" else -1.0
                    for tau_deg in torsions:
                        tau = math.radians(tau_deg)
                        u_perp = math.cos(tau) * in_plane_perp + math.sin(tau) * z_axis
                        oj_pos = center + 1.21 * (cos_d * direction + sin_d * u_perp)
                        e_co = cos_s * direction - sin_s * u_perp
                        ok_pos = center + 1.34 * e_co
                        n_syn = sin_s * direction + cos_s * u_perp
                        e_oh = cos_oh * e_co + sign_h * sin_oh * n_syn
                        hk_pos = ok_pos + 0.97 * e_oh
                        # Avoid severe intramolecular O=C...H clash in anti/syn
                        if np.linalg.norm(oj_pos - hk_pos) < 1.20:
                            continue
                        ok, min_14, min_non12 = _score_candidate_atoms(
                            [
                                (oj_pos, "O", set(), {site_id}),
                                (ok_pos, "O", set(), {site_id}),
                                (hk_pos, "H", set(), set()),
                            ],
                            local_coords,
                            local_ids,
                            local_elements,
                            pbc=pbc,
                            cell=cell,
                            inv_cell=inv_cell,
                            virtual_h_coords=virtual_h,
                        )
                        if not ok:
                            continue
                        clearance = min(c4_min14, min_14)
                        if conformer == "syn" and clearance >= clean_threshold:
                            best_choice = (clearance, oj_pos, ok_pos, hk_pos)
                            break
                        if best_choice is None or clearance > best_choice[0] + 0.02:
                            best_choice = (clearance, oj_pos, ok_pos, hk_pos)
                    if (
                        conformer == "syn"
                        and best_choice is not None
                        and best_choice[0] >= clean_threshold
                    ):
                        break

            if best_choice is None:
                return False
            _, oj_pos, ok_pos, hk_pos = best_choice
            new_atoms_spec = [
                ("C", center, "C4", 0),
                ("O", oj_pos, "OJ", 0),
                ("O", ok_pos, "OK", -1 if group == "COO" else 0),
            ]
            bond_spec = [
                (site_id, "C4", 1.0, (0, 0, 0)),
                ("C4", "OJ", 2.0, (0, 0, 0)),
                ("C4", "OK", 1.0, (0, 0, 0)),
            ]
            if group == "COOH":
                new_atoms_spec.append(("H", hk_pos, "HK", 0))
                bond_spec.append(("OK", "HK", 1.0, (0, 0, 0)))
            oxygen_tag = "OJ"

    # Commit validated geometry and topology to model
    if group == "OH":
        _set_site(
            model, sites[0], "CL", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
    elif group == "epoxide":
        _set_site(
            model, sites[0], "CY", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
        _set_site(
            model, sites[1], "CZ", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
    elif group == "epoxide_eph":
        _set_site(
            model, sites[0], "CZ", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
        _set_site(
            model, sites[1], "CK", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
    elif group == "epoxide_ept":
        _set_site(
            model, sites[0], "CP", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
        _set_site(
            model, sites[1], "CK", resname, resid, group_id, group, True,
            by_id=by_id, bond_by_atom=bond_by_atom,
        )
    else:
        _set_site(
            model,
            sites[0],
            "CN" if group == "carbonyl" else "CY",
            resname,
            resid,
            group_id,
            group,
            by_id=by_id,
            bond_by_atom=bond_by_atom,
        )

    for sid, new_pos in site_updates.items():
        pos_list = [float(v) for v in new_pos]
        by_id[sid]["position"] = pos_list
        coords_by_id[sid] = np.asarray(pos_list, dtype=float)
        if spatial_grid is not None:
            spatial_grid.update_pos(sid, pos_list)

    atom_ids = list(sites)
    tag_to_id: dict[str, int] = {}
    for elem, pos, name, charge in new_atoms_spec:
        identifier = _add_atom(
            model,
            elem,
            pos,
            name,
            resname,
            resid,
            role="group",
            group_id=group_id,
            group=group,
        )
        new_atom = model["atoms"][-1]
        new_atom["formal_charge"] = charge
        by_id[identifier] = new_atom
        coords_by_id[identifier] = np.asarray(new_atom["position"], dtype=float)
        adj[identifier] = []
        if bond_by_atom is not None:
            bond_by_atom[identifier] = []
        if spatial_grid is not None:
            spatial_grid.add(identifier, new_atom["position"], elem)
        atom_ids.append(identifier)
        tag_to_id[name] = identifier

    for a_ref, b_ref, order, image in bond_spec:
        a_id = tag_to_id[a_ref] if isinstance(a_ref, str) else a_ref
        b_id = tag_to_id[b_ref] if isinstance(b_ref, str) else b_ref
        new_bond = canonical_bond(a_id, b_id, order, image)
        model["bonds"].append(new_bond)
        if bond_by_atom is not None:
            bond_by_atom[a_id].append(new_bond)
            bond_by_atom[b_id].append(new_bond)
        img_arr = np.asarray(image, dtype=int)
        adj[a_id].append((b_id, img_arr, float(order)))
        adj[b_id].append((a_id, -img_arr, float(order)))

    if occupied_edges is not None:
        occupied_edges.update(sites)

    model["metadata"]["groups"].append(
        {
            "id": group_id,
            "type": group,
            "site_ids": list(sites),
            "atom_ids": atom_ids,
            "oxygen_id": tag_to_id[oxygen_tag],
            "face": chosen_face,
        }
    )
    return True


def _parse_group_requests(requested):
    raw_specs, fixed = {}, {}
    has_relative = False
    for group in GROUPS:
        value = requested.get(group, 0)
        label = f"groups.{group}"
        if not isinstance(value, dict):
            raw_specs[group] = _integer(value, label)
            continue
        allowed_keys = {"count", "id", "percent", "fraction", "basis"}
        if not value or set(value) - allowed_keys:
            raise ValueError(
                f"{label} accepts count/id or percent/fraction with optional basis only"
            )
        if "id" in value:
            if any(k in value for k in ("percent", "fraction", "basis")):
                raise ValueError(f"{label}.id cannot be combined with percent, fraction, or basis")
            sites = value["id"]
            if not isinstance(sites, list):
                raise ValueError(f"{label}.id must be a list of original carbon sites")
            if group in ("epoxide", "epoxide_eph", "epoxide_ept"):
                if any(not isinstance(pair, list) or len(pair) != 2 for pair in sites):
                    raise ValueError(f"{label}.id must list C-C pairs [[a,b],...]")
                fixed[group] = [list(pair) for pair in sites]
            else:
                fixed[group] = [[site] for site in sites]
            count = _integer(value.get("count", len(sites)), f"{label}.count")
            if count != len(sites):
                raise ValueError(f"{label}.count must equal the number of sites in id")
            raw_specs[group] = count
        elif "percent" in value or "fraction" in value:
            if "count" in value:
                raise ValueError(f"{label}.count cannot be combined with percent or fraction")
            if "percent" in value and "fraction" in value:
                raise ValueError(f"{label} cannot specify both percent and fraction")
            if "percent" in value:
                pct = _finite_real(value["percent"], f"{label}.percent")
                if not 0.0 <= pct <= 100.0:
                    raise ValueError(f"{label}.percent must be between 0 and 100")
                fraction = pct / 100.0
            else:
                fraction = _finite_real(value["fraction"], f"{label}.fraction")
                if not 0.0 <= fraction <= 1.0:
                    raise ValueError(f"{label}.fraction must be between 0 and 1")
            basis = value.get("basis", "eligible")
            if basis not in ("eligible", "carbons"):
                raise ValueError(f"{label}.basis must be 'eligible' or 'carbons'")
            raw_specs[group] = {"fraction": fraction, "basis": basis}
            has_relative = True
        else:
            if "basis" in value:
                raise ValueError(f"{label}.basis requires percent or fraction")
            if "count" not in value:
                raise ValueError(f"{label} requires count, id, percent, or fraction")
            raw_specs[group] = _integer(value["count"], f"{label}.count")
    return raw_specs, fixed, has_relative


def _resolve_group_counts(base, raw_specs, has_relative):
    if not has_relative:
        return dict(raw_specs)
    ns = adjacency(base)
    n_carbons = len(base["atoms"])
    edge = {i for i, n in ns.items() if len(n) == 2}
    basal = {i for i, n in ns.items() if len(n) == 3}
    det_rng = random.Random(0)
    epx_cap = None
    eph_cap = None
    ept_cap = None
    counts = {}
    for group in GROUPS:
        spec = raw_specs[group]
        if not isinstance(spec, dict):
            counts[group] = spec
            continue
        fraction = spec["fraction"]
        basis = spec["basis"]
        if basis == "carbons":
            ref_capacity = n_carbons
        else:
            if group == "OH":
                ref_capacity = len(basal)
            elif group == "epoxide":
                if epx_cap is None:
                    epx_cap = len(_max_matching(base, basal, det_rng))
                ref_capacity = epx_cap
            elif group in ("carbonyl", "COOH", "COO"):
                ref_capacity = len(edge)
            elif group == "epoxide_eph":
                if eph_cap is None:
                    eph_cap = len(_bipartite_matching(ns, sorted(basal), edge, det_rng))
                ref_capacity = eph_cap
            else:  # epoxide_ept
                if ept_cap is None:
                    ept_cap = len(_max_matching(base, edge, det_rng))
                ref_capacity = ept_cap
        counts[group] = int(round(fraction * ref_capacity))
    return counts


def _parse_placement(raw_placement):
    if raw_placement is None:
        return None
    if isinstance(raw_placement, str):
        raw = {"mode": raw_placement}
    elif isinstance(raw_placement, dict):
        raw = dict(raw_placement)
    else:
        raise ValueError("placement must be a string or an object")
    mode = raw.get("mode")
    if mode not in ("random", "island", "dispersed"):
        raise ValueError("placement.mode must be 'random', 'island', or 'dispersed'")
    if mode in ("random", "dispersed"):
        extra = set(raw) - {"mode"}
        if extra:
            raise ValueError(f"Unknown placement keys for mode {mode!r}: {sorted(extra)}")
        return {"mode": mode}

    allowed_island = {
        "mode",
        "n_islands",
        "centers",
        "seed_atoms",
        "radius",
        "compactness",
        "edge_mode",
    }
    extra = set(raw) - allowed_island
    if extra:
        raise ValueError(f"Unknown placement keys for mode 'island': {sorted(extra)}")
    centers = raw.get("centers")
    seed_atoms = raw.get("seed_atoms")
    if centers is not None and seed_atoms is not None:
        raise ValueError("placement.centers and placement.seed_atoms are mutually exclusive")
    parsed_centers = None
    parsed_seeds = None
    if centers is not None:
        if not isinstance(centers, list) or not centers:
            raise ValueError("placement.centers must be a non-empty list of [x, y] coordinates")
        parsed_centers = []
        for idx, pt in enumerate(centers):
            if not isinstance(pt, (list, tuple)) or len(pt) != 2:
                raise ValueError(f"placement.centers[{idx}] must be [x, y]")
            parsed_centers.append(
                [
                    _finite_real(pt[0], f"placement.centers[{idx}][0]"),
                    _finite_real(pt[1], f"placement.centers[{idx}][1]"),
                ]
            )
        if "n_islands" in raw:
            n_isl = _integer(raw["n_islands"], "placement.n_islands", 1)
            if n_isl != len(parsed_centers):
                raise ValueError("placement.n_islands must equal len(placement.centers)")
        n_islands = len(parsed_centers)
    elif seed_atoms is not None:
        if not isinstance(seed_atoms, list) or not seed_atoms:
            raise ValueError("placement.seed_atoms must be a non-empty list of carbon IDs")
        if any(type(s) is not int or s < 1 for s in seed_atoms) or len(set(seed_atoms)) != len(
            seed_atoms
        ):
            raise ValueError("placement.seed_atoms must contain distinct positive integer IDs")
        parsed_seeds = list(seed_atoms)
        if "n_islands" in raw:
            n_isl = _integer(raw["n_islands"], "placement.n_islands", 1)
            if n_isl != len(parsed_seeds):
                raise ValueError("placement.n_islands must equal len(placement.seed_atoms)")
        n_islands = len(parsed_seeds)
    else:
        n_islands = _integer(raw.get("n_islands", 1), "placement.n_islands", 1)

    radius = raw.get("radius")
    if radius is not None:
        radius = _finite_real(radius, "placement.radius")
        if radius <= 0.0:
            raise ValueError("placement.radius must be positive")
    compactness = _finite_real(raw.get("compactness", 0.85), "placement.compactness")
    if not 0.0 <= compactness <= 1.0:
        raise ValueError("placement.compactness must be between 0 and 1")
    edge_mode = raw.get("edge_mode", "random")
    if edge_mode not in ("random", "near_islands"):
        raise ValueError("placement.edge_mode must be 'random' or 'near_islands'")
    return {
        "mode": "island",
        "n_islands": n_islands,
        "centers": parsed_centers,
        "seed_atoms": parsed_seeds,
        "radius": radius,
        "compactness": compactness,
        "edge_mode": edge_mode,
    }


def _min_img_xy_dist(pos_a, pos_b, model):
    delta = np.array([float(pos_a[0]) - float(pos_b[0]), float(pos_a[1]) - float(pos_b[1]), 0.0])
    return float(np.linalg.norm(_minimum_image(delta, model)[:2]))


def _pair_midpoint(coords_by_id, a, b, model):
    pos_a = coords_by_id[a]
    pos_b = coords_by_id[b]
    return pos_a + 0.5 * _minimum_image(pos_b - pos_a, model)


def _min_dist_to_centers(pos, centers, model):
    return min(_min_img_xy_dist(pos, center, model) for center in centers)


def _choose_island_nuclei(pool, n_islands, coords_by_id, model, rng):
    if n_islands > len(pool):
        raise ValueError(
            f"placement.n_islands ({n_islands}) exceeds available sheet sites ({len(pool)})"
        )
    chosen = [rng.choice(pool)]
    if n_islands == 1:
        return chosen
    min_dists = {cid: _min_img_xy_dist(coords_by_id[cid], coords_by_id[chosen[0]], model) for cid in pool if cid != chosen[0]}
    while len(chosen) < n_islands:
        max_d = max(min_dists.values())
        ties = [cid for cid, d in min_dists.items() if d >= max_d - 1e-6]
        pick = rng.choice(ties)
        chosen.append(pick)
        del min_dists[pick]
        p_pick = coords_by_id[pick]
        for cid in list(min_dists):
            d = _min_img_xy_dist(coords_by_id[cid], p_pick, model)
            if d < min_dists[cid]:
                min_dists[cid] = d
    return chosen


def _order_by_island(items, pos_fn, island_centers, compactness, model, rng):
    if not items:
        return []
    dists = [_min_dist_to_centers(pos_fn(item), island_centers, model) for item in items]
    d_min = min(dists)
    span = max(max(dists) - d_min, 1e-9)
    scored = []
    for item, d in zip(items, dists, strict=True):
        jitter = rng.random()
        if compactness >= 1.0:
            key = (round(d, 9), jitter)
        else:
            score = compactness * ((d - d_min) / span) + (1.0 - compactness) * jitter
            key = (score, jitter)
        scored.append((key, item))
    scored.sort(key=lambda pair: pair[0])
    return [item for _, item in scored]


def _compute_mean_basal_group_nn_distance(model, coords_by_id):
    basal_groups = [g for g in model["metadata"]["groups"] if g["type"] in ("OH", "epoxide")]
    if len(basal_groups) < 2:
        return None
    centers = []
    for g in basal_groups:
        sids = g["site_ids"]
        if len(sids) == 1:
            centers.append(coords_by_id[sids[0]][:2])
        else:
            centers.append(_pair_midpoint(coords_by_id, sids[0], sids[1], model)[:2])
    nn_dists = []
    for i, pos_i in enumerate(centers):
        best_d = min(
            _min_img_xy_dist(pos_i, centers[j], model)
            for j in range(len(centers))
            if j != i
        )
        nn_dists.append(best_d)
    return float(np.mean(nn_dists))


def _clone_model_fast(base: dict) -> dict:
    """Fast structural copy of a sheet model avoiding copy.deepcopy overhead."""
    return {
        "schema_version": base["schema_version"],
        "name": base["name"],
        "pbc": list(base["pbc"]),
        "cell": [list(row) for row in base["cell"]] if base["cell"] is not None else None,
        "atoms": [
            {
                "id": a["id"],
                "element": a["element"],
                "name": a["name"],
                "resname": a["resname"],
                "resid": a["resid"],
                "position": list(a["position"]),
                "atom_type": a["atom_type"],
                "charge": a["charge"],
                "formal_charge": a["formal_charge"],
                "metadata": dict(a["metadata"]),
            }
            for a in base["atoms"]
        ],
        "bonds": [
            {
                "a": b["a"],
                "b": b["b"],
                "order": b["order"],
                "image": list(b["image"]),
            }
            for b in base["bonds"]
        ],
        "metadata": copy.deepcopy(base["metadata"]),
    }


def _decorate(base, counts, fixed, faces, rng, *, placement_spec=None, has_p2=False, max_attempts=128):
    ns = adjacency(base)
    edge = [i for i, n in ns.items() if len(n) == 2]
    basal = {i for i, n in ns.items() if len(n) == 3}
    base_coords = _coordinates(base)

    if placement_spec and placement_spec["mode"] == "island" and placement_spec["seed_atoms"] is not None:
        for sid in placement_spec["seed_atoms"]:
            if sid not in ns:
                raise ValueError(
                    f"placement.seed_atoms ID {sid} is not present in the sheet after holes/vacancies"
                )

    reserved = set()
    normalized_fixed = {}
    for group, sites in fixed.items():
        normalized_sites = []
        for site in sites:
            for identifier in site:
                if type(identifier) is not int or identifier not in ns:
                    raise ValueError(f"groups.{group}.id requires existing original carbon IDs")
                if identifier in reserved:
                    raise ValueError("Functional group site IDs must be distinct and disjoint")
            if group in ("OH", "epoxide"):
                for identifier in site:
                    if identifier not in basal:
                        raise ValueError(f"groups.{group}.id requires degree-3 carbon sites")
                if group == "epoxide" and not any(j == site[1] for j, _, _ in ns[site[0]]):
                    raise ValueError("groups.epoxide.id requires actual basal C-C bonds")
            elif group == "epoxide_eph":
                if not any(j == site[1] for j, _, _ in ns[site[0]]):
                    raise ValueError("groups.epoxide_eph.id requires actual C-C bonds")
                deg0, deg1 = len(ns[site[0]]), len(ns[site[1]])
                if deg0 == 2 and deg1 == 3:
                    site = [site[1], site[0]]
                elif not (deg0 == 3 and deg1 == 2):
                    raise ValueError(
                        "groups.epoxide_eph.id requires one degree-3 basal and one degree-2 edge carbon"
                    )
            elif group == "epoxide_ept":
                for identifier in site:
                    if identifier not in edge:
                        raise ValueError("groups.epoxide_ept.id requires degree-2 carbon sites")
                if not any(j == site[1] for j, _, _ in ns[site[0]]):
                    raise ValueError("groups.epoxide_ept.id requires actual edge C-C bonds")
            else:
                for identifier in site:
                    if identifier not in edge:
                        raise ValueError(f"groups.{group}.id requires degree-2 carbon sites")
            reserved.update(site)
            normalized_sites.append(list(site))
        normalized_fixed[group] = normalized_sites
    fixed = normalized_fixed

    # Reserve every fixed site before choosing any random group, including
    # groups placed earlier in the construction sequence.
    edge = [i for i in edge if i not in reserved]
    basal -= reserved
    random_counts = {group: 0 if group in fixed else count for group, count in counts.items()}
    edge_count = (
        sum(random_counts[g] for g in ("carbonyl", "COOH", "COO"))
        + random_counts["epoxide_eph"]
        + 2 * random_counts["epoxide_ept"]
    )
    basal_count = random_counts["OH"] + 2 * random_counts["epoxide"] + random_counts["epoxide_eph"]

    report_keys = GROUPS if has_p2 else LEGACY_GROUPS
    requested_report = {k: counts[k] for k in report_keys}
    capacity = {"eligible_edge_carbons": len(edge), "eligible_basal_carbons": len(basal)}
    suggestions = ["Enlarge the sheet", "Reduce requested group counts", "Change the seed"]
    if edge_count > len(edge) or basal_count > len(basal):
        raise PlacementError(
            "Requested functional group counts exceed eligible edge/basal sites",
            requested=requested_report,
            capacity=capacity,
            attempts=0,
            suggestions=suggestions,
        )
    has_epx_work = random_counts["epoxide"] > 0
    if has_epx_work:
        matching = _max_matching(base, basal, rng)
        capacity["max_disjoint_epoxide_sites"] = len(matching)
        if random_counts["epoxide"] > len(matching):
            raise PlacementError(
                "Requested epoxide count exceeds disjoint eligible C-C bonds",
                requested=requested_report,
                capacity=capacity,
                attempts=0,
                suggestions=suggestions,
            )

    if has_p2:
        if random_counts["epoxide_eph"] > 0:
            eph_cap_match = _bipartite_matching(ns, sorted(basal), set(edge), rng)
            capacity["max_disjoint_eph_sites"] = len(eph_cap_match)
            if random_counts["epoxide_eph"] > len(eph_cap_match):
                raise PlacementError(
                    "Requested epoxide_eph count exceeds disjoint (basal, edge) C-C bonds",
                    requested=requested_report,
                    capacity=capacity,
                    attempts=0,
                    suggestions=suggestions,
                )
        if random_counts["epoxide_ept"] > 0:
            ept_cap_match = _max_matching(base, set(edge), rng)
            capacity["max_disjoint_ept_sites"] = len(ept_cap_match)
            if random_counts["epoxide_ept"] > len(ept_cap_match):
                raise PlacementError(
                    "Requested epoxide_ept count exceeds disjoint edge C-C bonds",
                    requested=requested_report,
                    capacity=capacity,
                    attempts=0,
                    suggestions=suggestions,
                )

    mode = placement_spec["mode"] if placement_spec is not None else "random"
    island_radius = placement_spec["radius"] if mode == "island" else None
    fixed_island_centers = None
    if mode == "island":
        if placement_spec["centers"] is not None:
            fixed_island_centers = placement_spec["centers"]
        elif placement_spec["seed_atoms"] is not None:
            fixed_island_centers = [
                [float(base_coords[sid][0]), float(base_coords[sid][1])]
                for sid in placement_spec["seed_atoms"]
            ]

    if mode == "island" and island_radius is not None and fixed_island_centers is not None:
        disk_basal = {
            s
            for s in basal
            if _min_dist_to_centers(base_coords[s], fixed_island_centers, base) <= island_radius
        }
        capacity["island_eligible_basal_carbons"] = len(disk_basal)
        if has_epx_work:
            disk_match = _max_matching(base, disk_basal, rng)
            capacity["island_max_disjoint_epoxide_sites"] = len(disk_match)
        else:
            disk_match = []
        if (
            random_counts["OH"] + 2 * random_counts["epoxide"] > len(disk_basal)
            or random_counts["epoxide"] > len(disk_match)
        ):
            raise PlacementError(
                "Island radius is too small for requested functional group counts",
                requested=requested_report,
                capacity=capacity,
                attempts=0,
                suggestions=["Increase placement.radius", "Reduce basal group counts", "Add island centers"],
            )

    # Retry stereochemical placements and randomized tie choices. Count is exact
    # on success; search failure is explicit and never silently undershoots.
    best_achieved = dict.fromkeys(report_keys, 0)
    basal_sorted = sorted(basal)
    has_eph_work = random_counts["epoxide_eph"] > 0 or bool(fixed.get("epoxide_eph"))
    has_ept_work = random_counts["epoxide_ept"] > 0 or bool(fixed.get("epoxide_ept"))
    has_any_groups = any(counts.values()) or bool(fixed)
    radius_too_small_seen = False

    for _ in range(max_attempts):
        model = _clone_model_fast(base)
        adj = adjacency(model)
        coords_by_id = _coordinates(model)
        by_id = {a["id"]: a for a in model["atoms"]}
        bond_by_atom: dict[int, list[dict]] = {aid: [] for aid in by_id}
        for bond in model["bonds"]:
            bond_by_atom[bond["a"]].append(bond)
            bond_by_atom[bond["b"]].append(bond)
        spatial_grid = _make_spatial_grid(model)
        occupied_edges = set(reserved)

        # Keep exact baseline RNG sequence when placement == "random" and P2 counts == 0
        edge_pool = list(edge)
        rng.shuffle(edge_pool)

        island_centers = fixed_island_centers
        chosen_seed_atoms = placement_spec["seed_atoms"] if mode == "island" else None
        if mode == "island":
            if island_centers is None:
                nuclei_pool = basal_sorted if basal_sorted else sorted(ns)
                chosen_seed_atoms = _choose_island_nuclei(
                    nuclei_pool, placement_spec["n_islands"], base_coords, base, rng
                )
                island_centers = [
                    [float(base_coords[sid][0]), float(base_coords[sid][1])]
                    for sid in chosen_seed_atoms
                ]
            if island_radius is not None:
                cand_basal = {
                    s
                    for s in basal
                    if _min_dist_to_centers(base_coords[s], island_centers, base) <= island_radius
                }
                matching = _max_matching(base, cand_basal, rng) if has_epx_work else []
                if (
                    random_counts["OH"] + 2 * random_counts["epoxide"] > len(cand_basal)
                    or random_counts["epoxide"] > len(matching)
                ):
                    radius_too_small_seen = True
                    continue
                cand_basal_sorted = sorted(cand_basal)
            else:
                cand_basal_sorted = basal_sorted
                matching = _max_matching(base, basal, rng) if has_epx_work else []

            compactness = placement_spec["compactness"]
            if has_epx_work:
                matching = _order_by_island(
                    matching,
                    lambda p: _pair_midpoint(base_coords, p[0], p[1], base),
                    island_centers,
                    compactness,
                    base,
                    rng,
                )
            oh_pool = _order_by_island(
                cand_basal_sorted,
                lambda s: base_coords[s],
                island_centers,
                compactness,
                base,
                rng,
            )
            if placement_spec["edge_mode"] == "near_islands":
                edge_pool = _order_by_island(
                    edge_pool,
                    lambda s: base_coords[s],
                    island_centers,
                    compactness,
                    base,
                    rng,
                )
        else:
            if has_epx_work:
                matching = _max_matching(base, basal, rng)
                rng.shuffle(matching)
            else:
                matching = []
            oh_pool = list(basal_sorted)
            rng.shuffle(oh_pool)

        ept_matching = []
        if has_ept_work and random_counts["epoxide_ept"] > 0:
            ept_matching = _max_matching(base, set(edge), rng)
            rng.shuffle(ept_matching)
            if mode == "island" and placement_spec["edge_mode"] == "near_islands":
                ept_matching = _order_by_island(
                    ept_matching,
                    lambda p: _pair_midpoint(base_coords, p[0], p[1], base),
                    island_centers,
                    placement_spec["compactness"],
                    base,
                    rng,
                )

        used_edge = set()
        used_basal = set()
        group_id = 1
        achieved = dict.fromkeys(report_keys, 0)
        success = True

        placed_basal_pts = []
        placed_edge_pts = []
        if mode == "dispersed":
            for gname in ("OH", "epoxide"):
                for sites in fixed.get(gname, ()):
                    if len(sites) == 1:
                        placed_basal_pts.append(base_coords[sites[0]][:2])
                    else:
                        placed_basal_pts.append(_pair_midpoint(base_coords, sites[0], sites[1], base)[:2])
            for gname in ("COOH", "COO", "carbonyl", "epoxide_eph", "epoxide_ept"):
                for sites in fixed.get(gname, ()):
                    if len(sites) == 1:
                        placed_edge_pts.append(base_coords[sites[0]][:2])
                    else:
                        placed_edge_pts.append(_pair_midpoint(base_coords, sites[0], sites[1], base)[:2])

        def _place_dispersed(candidates, pos_fn, placed_pts, needed, group_name, is_pair, update_sets):
            nonlocal group_id
            remaining = []
            for cand in candidates:
                if is_pair:
                    if cand[0] in used_basal or cand[1] in used_basal or cand[0] in used_edge or cand[1] in used_edge:
                        continue
                else:
                    if cand in (used_basal if group_name == "OH" else used_edge):
                        continue
                remaining.append((cand, pos_fn(cand)[:2], rng.random()))
            min_d = [
                min((_min_img_xy_dist(pt, p, base) for p in placed_pts), default=float("inf"))
                for _, pt, _ in remaining
            ]
            placed_here = 0
            while placed_here < needed and remaining:
                if not placed_pts:
                    order = list(range(len(remaining)))
                else:
                    order = sorted(
                        range(len(remaining)),
                        key=lambda idx: (-round(min_d[idx], 9), remaining[idx][2]),
                    )
                chosen_idx = None
                for idx in order:
                    cand = remaining[idx][0]
                    sites_arg = list(cand) if is_pair else [cand]
                    face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                    if _place_group(
                        model,
                        group_name,
                        sites_arg,
                        face,
                        group_id,
                        faces=faces,
                        adj=adj,
                        coords_by_id=coords_by_id,
                        occupied_edges=occupied_edges,
                        spatial_grid=spatial_grid,
                        by_id=by_id,
                        bond_by_atom=bond_by_atom,
                    ):
                        chosen_idx = idx
                        break
                if chosen_idx is None:
                    return False
                cand, pt, _ = remaining[chosen_idx]
                update_sets(cand)
                group_id += 1
                achieved[group_name] += 1
                placed_here += 1
                placed_pts.append(pt)
                # Remove chosen and update incremental min_d
                new_remaining = []
                new_min_d = []
                for idx, item in enumerate(remaining):
                    if idx == chosen_idx:
                        continue
                    c_item = item[0]
                    if is_pair:
                        if c_item[0] in used_basal or c_item[1] in used_basal or c_item[0] in used_edge or c_item[1] in used_edge:
                            continue
                    d_new = _min_img_xy_dist(item[1], pt, base)
                    new_remaining.append(item)
                    new_min_d.append(min(min_d[idx], d_new))
                remaining = new_remaining
                min_d = new_min_d
            return placed_here == needed

        # 1. Place Priority-2 epoxide_ept (only if requested/fixed)
        if has_ept_work:
            if "epoxide_ept" in fixed:
                for sites in fixed["epoxide_ept"]:
                    face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                    if not _place_group(
                        model,
                        "epoxide_ept",
                        list(sites),
                        face,
                        group_id,
                        faces=faces,
                        adj=adj,
                        coords_by_id=coords_by_id,
                        occupied_edges=occupied_edges,
                        spatial_grid=spatial_grid,
                        by_id=by_id,
                        bond_by_atom=bond_by_atom,
                    ):
                        success = False
                        break
                    used_edge.update(sites)
                    group_id += 1
                    achieved["epoxide_ept"] += 1
            elif random_counts["epoxide_ept"] > 0:
                needed = random_counts["epoxide_ept"]
                if mode == "dispersed":
                    if not _place_dispersed(
                        ept_matching,
                        lambda p: _pair_midpoint(base_coords, p[0], p[1], base),
                        placed_edge_pts,
                        needed,
                        "epoxide_ept",
                        True,
                        lambda p: used_edge.update(p),
                    ):
                        success = False
                else:
                    placed_here = 0
                    for pair in ept_matching:
                        if pair[0] in used_edge or pair[1] in used_edge:
                            continue
                        face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                        if _place_group(
                            model,
                            "epoxide_ept",
                            list(pair),
                            face,
                            group_id,
                            faces=faces,
                            adj=adj,
                            coords_by_id=coords_by_id,
                            occupied_edges=occupied_edges,
                            spatial_grid=spatial_grid,
                            by_id=by_id,
                            bond_by_atom=bond_by_atom,
                        ):
                            used_edge.update(pair)
                            group_id += 1
                            achieved["epoxide_ept"] += 1
                            placed_here += 1
                            if placed_here == needed:
                                break
                    if placed_here < needed:
                        success = False

        # 2. Place Priority-2 epoxide_eph (only if requested/fixed)
        if success and has_eph_work:
            if "epoxide_eph" in fixed:
                for sites in fixed["epoxide_eph"]:
                    face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                    if not _place_group(
                        model,
                        "epoxide_eph",
                        list(sites),
                        face,
                        group_id,
                        faces=faces,
                        adj=adj,
                        coords_by_id=coords_by_id,
                        occupied_edges=occupied_edges,
                        spatial_grid=spatial_grid,
                        by_id=by_id,
                        bond_by_atom=bond_by_atom,
                    ):
                        success = False
                        break
                    used_basal.add(sites[0])
                    used_edge.add(sites[1])
                    group_id += 1
                    achieved["epoxide_eph"] += 1
            elif random_counts["epoxide_eph"] > 0:
                needed = random_counts["epoxide_eph"]
                avail_basal_eph = sorted(basal - used_basal)
                avail_edge_eph = set(edge) - used_edge
                eph_matching = _bipartite_matching(ns, avail_basal_eph, avail_edge_eph, rng)
                rng.shuffle(eph_matching)
                if mode == "island" and (
                    placement_spec["edge_mode"] == "near_islands" or island_radius is None
                ):
                    eph_matching = _order_by_island(
                        eph_matching,
                        lambda p: _pair_midpoint(base_coords, p[0], p[1], base),
                        island_centers,
                        placement_spec["compactness"],
                        base,
                        rng,
                    )
                if mode == "dispersed":
                    if not _place_dispersed(
                        eph_matching,
                        lambda p: _pair_midpoint(base_coords, p[0], p[1], base),
                        placed_edge_pts,
                        needed,
                        "epoxide_eph",
                        True,
                        lambda p: (used_basal.add(p[0]), used_edge.add(p[1])),
                    ):
                        success = False
                else:
                    placed_here = 0
                    for pair in eph_matching:
                        if pair[0] in used_basal or pair[1] in used_edge:
                            continue
                        face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                        if _place_group(
                            model,
                            "epoxide_eph",
                            list(pair),
                            face,
                            group_id,
                            faces=faces,
                            adj=adj,
                            coords_by_id=coords_by_id,
                            occupied_edges=occupied_edges,
                            spatial_grid=spatial_grid,
                            by_id=by_id,
                            bond_by_atom=bond_by_atom,
                        ):
                            used_basal.add(pair[0])
                            used_edge.add(pair[1])
                            group_id += 1
                            achieved["epoxide_eph"] += 1
                            placed_here += 1
                            if placed_here == needed:
                                break
                    if placed_here < needed:
                        success = False

        # 3. Place edge groups (COOH, COO, carbonyl)
        if success:
            for group in ("COOH", "COO", "carbonyl"):
                if group in fixed:
                    for sites in fixed[group]:
                        face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                        if not _place_group(
                            model,
                            group,
                            sites,
                            face,
                            group_id,
                            faces=faces,
                            adj=adj,
                            coords_by_id=coords_by_id,
                            occupied_edges=occupied_edges,
                            spatial_grid=spatial_grid,
                            by_id=by_id,
                            bond_by_atom=bond_by_atom,
                        ):
                            success = False
                            break
                        group_id += 1
                        achieved[group] += 1
                    if not success:
                        break
                else:
                    needed = random_counts[group]
                    if needed:
                        if mode == "dispersed":
                            if not _place_dispersed(
                                edge_pool,
                                lambda s: base_coords[s],
                                placed_edge_pts,
                                needed,
                                group,
                                False,
                                lambda s: used_edge.add(s),
                            ):
                                success = False
                                break
                        else:
                            placed_here = 0
                            for site in edge_pool:
                                if site in used_edge:
                                    continue
                                face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                                if _place_group(
                                    model,
                                    group,
                                    [site],
                                    face,
                                    group_id,
                                    faces=faces,
                                    adj=adj,
                                    coords_by_id=coords_by_id,
                                    occupied_edges=occupied_edges,
                                    spatial_grid=spatial_grid,
                                    by_id=by_id,
                                    bond_by_atom=bond_by_atom,
                                ):
                                    used_edge.add(site)
                                    group_id += 1
                                    achieved[group] += 1
                                    placed_here += 1
                                    if placed_here == needed:
                                        break
                            if placed_here < needed:
                                success = False
                                break

        # 4. Place epoxide groups
        if success:
            if "epoxide" in fixed:
                for sites in fixed["epoxide"]:
                    face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                    if not _place_group(
                        model,
                        "epoxide",
                        list(sites),
                        face,
                        group_id,
                        faces=faces,
                        adj=adj,
                        coords_by_id=coords_by_id,
                        occupied_edges=occupied_edges,
                        spatial_grid=spatial_grid,
                        by_id=by_id,
                        bond_by_atom=bond_by_atom,
                    ):
                        success = False
                        break
                    used_basal.update(sites)
                    group_id += 1
                    achieved["epoxide"] += 1
            else:
                needed = random_counts["epoxide"]
                if needed:
                    if mode == "dispersed":
                        if not _place_dispersed(
                            matching,
                            lambda p: _pair_midpoint(base_coords, p[0], p[1], base),
                            placed_basal_pts,
                            needed,
                            "epoxide",
                            True,
                            lambda p: used_basal.update(p),
                        ):
                            success = False
                    else:
                        placed_here = 0
                        for pair in matching:
                            if pair[0] in used_basal or pair[1] in used_basal:
                                continue
                            face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                            if _place_group(
                                model,
                                "epoxide",
                                list(pair),
                                face,
                                group_id,
                                faces=faces,
                                adj=adj,
                                coords_by_id=coords_by_id,
                                occupied_edges=occupied_edges,
                                spatial_grid=spatial_grid,
                                by_id=by_id,
                                bond_by_atom=bond_by_atom,
                            ):
                                used_basal.update(pair)
                                group_id += 1
                                achieved["epoxide"] += 1
                                placed_here += 1
                                if placed_here == needed:
                                    break
                        if placed_here < needed:
                            success = False

        # 5. Place OH groups
        if success:
            if "OH" in fixed:
                for sites in fixed["OH"]:
                    face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                    if not _place_group(
                        model,
                        "OH",
                        list(sites),
                        face,
                        group_id,
                        faces=faces,
                        adj=adj,
                        coords_by_id=coords_by_id,
                        occupied_edges=occupied_edges,
                        spatial_grid=spatial_grid,
                        by_id=by_id,
                        bond_by_atom=bond_by_atom,
                    ):
                        success = False
                        break
                    used_basal.update(sites)
                    group_id += 1
                    achieved["OH"] += 1
            else:
                needed = random_counts["OH"]
                if needed:
                    if mode == "dispersed":
                        if not _place_dispersed(
                            oh_pool,
                            lambda s: base_coords[s],
                            placed_basal_pts,
                            needed,
                            "OH",
                            False,
                            lambda s: used_basal.add(s),
                        ):
                            success = False
                    else:
                        placed_here = 0
                        for site in oh_pool:
                            if site in used_basal:
                                continue
                            face = 1 if faces == "up" else -1 if faces == "down" else rng.choice([-1, 1])
                            if _place_group(
                                model,
                                "OH",
                                [site],
                                face,
                                group_id,
                                faces=faces,
                                adj=adj,
                                coords_by_id=coords_by_id,
                                occupied_edges=occupied_edges,
                                spatial_grid=spatial_grid,
                                by_id=by_id,
                                bond_by_atom=bond_by_atom,
                            ):
                                used_basal.add(site)
                                group_id += 1
                                achieved["OH"] += 1
                                placed_here += 1
                                if placed_here == needed:
                                    break
                        if placed_here < needed:
                            success = False

        if success:
            try:
                _cap_edges(
                    model,
                    spatial_grid=spatial_grid,
                    by_id=by_id,
                    adj=adj,
                    coords_by_id=coords_by_id,
                )
            except ValueError:
                if not has_any_groups:
                    raise
                success = False
            if success and has_any_groups:
                integrity = check_geometry_integrity(model)
                if not integrity["passed"]:
                    success = False

        if success:
            if placement_spec is not None:
                pmeta = {"mode": mode}
                if mode == "island":
                    pmeta.update(
                        {
                            "n_islands": len(island_centers),
                            "island_centers": [[float(x), float(y)] for x, y in island_centers],
                            "seed_atoms": (
                                list(chosen_seed_atoms) if chosen_seed_atoms is not None else None
                            ),
                            "radius": placement_spec["radius"],
                            "compactness": placement_spec["compactness"],
                            "edge_mode": placement_spec["edge_mode"],
                        }
                    )
                mean_nn = _compute_mean_basal_group_nn_distance(model, base_coords)
                if mean_nn is not None:
                    pmeta["mean_basal_group_nn_distance"] = mean_nn
                model["metadata"]["placement"] = pmeta
            return model
        if sum(achieved.values()) > sum(best_achieved.values()):
            best_achieved = achieved

    if radius_too_small_seen and sum(best_achieved.values()) == 0:
        raise PlacementError(
            "Island radius is too small for requested functional group counts",
            requested=requested_report,
            capacity=capacity,
            attempts=max_attempts,
            suggestions=["Increase placement.radius", "Reduce basal group counts", "Add island centers"],
        )
    raise PlacementError(
        f"No collision-free exact placement found in {max_attempts} attempts; "
        f"best tried placement counts: {best_achieved}; adjust composition, seed or size",
        requested=requested_report,
        achieved=best_achieved,
        capacity=capacity,
        attempts=max_attempts,
        suggestions=suggestions,
    )


def _cap_edges(model, *, spatial_grid=None, by_id=None, adj=None, coords_by_id=None):
    if adj is None:
        adj = adjacency(model)
    if coords_by_id is None:
        coords_by_id = _coordinates(model)
    cell = np.array(model["cell"])
    pbc = model["pbc"]
    inv_cell = np.linalg.inv(cell) if any(pbc) else None
    if by_id is None:
        by_id = {atom["id"]: atom for atom in model["atoms"]}
    next_resid = max(a.get("resid", 1) for a in model["atoms"]) + 1
    occupied = {i for g in model["metadata"]["groups"] for i in g["site_ids"]}
    if spatial_grid is None:
        spatial_grid = _make_spatial_grid(model)

    edge_items = [
        (identifier, neighbors)
        for identifier, neighbors in adj.items()
        if by_id[identifier].get("metadata", {}).get("role") == "sheet_carbon"
        and len(neighbors) == 2
        and identifier not in occupied
    ]

    for identifier, neighbors in edge_items:
        atom = by_id[identifier]
        resid = next_resid
        next_resid += 1
        atom.update(name="CH", resname="GRT", resid=resid)
        outward = _outward(
            model, identifier, neighbors=adj, coordinates=coords_by_id, cell=cell
        )
        u_out, _, u_n = _orthonormal_edge_frame(outward)
        anchor_pos = coords_by_id[identifier]
        l_coords, l_ids, _ = spatial_grid.query_radius(anchor_pos, 5.5)
        ids_13 = {nbr for nbr, _, _ in neighbors}
        non12_mask = l_ids != identifier
        mask_14 = non12_mask & ~np.isin(l_ids, tuple(ids_13))

        # Small holes or nearby carboxyls can collide with planar H caps.
        # Sample tilts and enforce 1-4+ clearance (>= 1.20 A), preferring 0 deg if > 1.8 A.
        best_choice = None
        for idx_angle, angle in enumerate((0, 20, -20, 35, -35, 50, -50, 65, -65)):
            radians = math.radians(angle)
            direction = math.cos(radians) * u_out + math.sin(radians) * u_n
            pos = anchor_pos + 1.09 * direction
            if len(l_coords):
                delta = l_coords - pos
                if any(pbc):
                    frac = delta @ inv_cell
                    for k, periodic in enumerate(pbc):
                        if periodic:
                            frac[..., k] -= np.rint(frac[..., k])
                    delta = frac @ cell
                dists = np.linalg.norm(delta, axis=1)
                if np.any(dists[non12_mask] < 0.90):
                    continue
                min_non12 = float(np.min(dists[non12_mask])) if np.any(non12_mask) else float("inf")
                min_14 = float(np.min(dists[mask_14])) if np.any(mask_14) else float("inf")
            else:
                min_non12 = float("inf")
                min_14 = float("inf")
            if min_14 < 1.20:
                continue
            if angle == 0 and min_non12 > 1.80:
                best_choice = (min_non12, -idx_angle, angle, pos)
                break
            candidate = (round(min_non12, 9), -idx_angle, angle, pos)
            if best_choice is None or candidate > best_choice:
                best_choice = candidate

        if best_choice is None:
            raise ValueError("No nonoverlapping edge H cap; enlarge or modify the hole")

        _, _, chosen_angle, chosen_pos = best_choice
        hydrogen = _add_atom(
            model,
            "H",
            chosen_pos,
            "HC",
            "GRT",
            resid,
            role="edge_hydrogen",
            anchor_id=identifier,
            cap_tilt_degrees=chosen_angle,
        )
        new_h = model["atoms"][-1]
        by_id[hydrogen] = new_h
        coords_by_id[hydrogen] = np.asarray(new_h["position"], dtype=float)
        spatial_grid.add(hydrogen, new_h["position"], "H")
        model["bonds"].append(canonical_bond(identifier, hydrogen))


def _parse_composition_spec(recipe):
    has_comp = "composition" in recipe and recipe["composition"] is not None
    has_c_over_o = "c_over_o" in recipe and recipe["c_over_o"] is not None
    has_o_over_c = "o_over_c" in recipe and recipe["o_over_c"] is not None
    if not (has_comp or has_c_over_o or has_o_over_c):
        return None
    if recipe.get("groups"):
        raise ValueError(
            "Cannot specify both non-empty groups and composition/c_over_o/o_over_c"
        )
    if has_comp and (has_c_over_o or has_o_over_c):
        raise ValueError(
            "Specify either 'composition' or top-level 'c_over_o'/'o_over_c', not both"
        )
    if has_c_over_o and has_o_over_c:
        raise ValueError("Specify either top-level c_over_o or o_over_c, not both")

    if has_c_over_o:
        return {
            "c_over_o": recipe["c_over_o"],
            "o_over_c": None,
            "preferences": None,
            "priority": None,
            "mode": "balanced",
            "enabled_groups": None,
            "ph": None,
            "pka": 4.5,
            "deprotonation_fraction": None,
            "seed": None,
        }
    if has_o_over_c:
        return {
            "c_over_o": None,
            "o_over_c": recipe["o_over_c"],
            "preferences": None,
            "priority": None,
            "mode": "balanced",
            "enabled_groups": None,
            "ph": None,
            "pka": 4.5,
            "deprotonation_fraction": None,
            "seed": None,
        }

    comp = recipe["composition"]
    if not isinstance(comp, dict):
        raise ValueError("composition must be a dictionary")
    allowed_comp = {
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
    extra = set(comp) - allowed_comp
    if extra:
        raise ValueError(f"Unknown composition keys: {sorted(extra)}")
    if "C_over_O" in comp and "c_over_o" in comp:
        raise ValueError("Specify either C_over_O or c_over_o, not both")
    if "O_over_C" in comp and "o_over_c" in comp:
        raise ValueError("Specify either O_over_C or o_over_c, not both")
    c_over_o_val = comp.get("C_over_O", comp.get("c_over_o"))
    o_over_c_val = comp.get("O_over_C", comp.get("o_over_c"))
    if c_over_o_val is not None and o_over_c_val is not None:
        raise ValueError("Specify either C_over_O/c_over_o or O_over_C/o_over_c, not both")
    if c_over_o_val is None and o_over_c_val is None:
        raise ValueError("composition requires C_over_O/c_over_o or O_over_C/o_over_c")
    if "preferences" in comp and "group_preferences" in comp:
        raise ValueError("Specify either preferences or group_preferences, not both")
    prefs = comp.get("preferences", comp.get("group_preferences"))
    return {
        "c_over_o": c_over_o_val,
        "o_over_c": o_over_c_val,
        "preferences": prefs,
        "priority": comp.get("priority"),
        "mode": comp.get("mode", "balanced"),
        "enabled_groups": comp.get("enabled_groups"),
        "ph": comp.get("ph"),
        "pka": comp.get("pka", 4.5),
        "deprotonation_fraction": comp.get("deprotonation_fraction"),
        "seed": comp.get("seed"),
    }


def _apply_shape(
    model,
    *,
    nx: int,
    ny: int,
    length: float,
    shape: str,
    edge_type: str | None,
    cut_radius: float | None,
):
    if shape == "rectangular":
        return
    h = math.sqrt(3.0) * length / 2.0
    span_x = (3.0 * nx - 0.5) * length
    span_y = (2.0 * ny - 1.0) * h
    target_cx = 0.5 * span_x
    target_cy = 0.5 * span_y
    ring_centers = []
    for i in range(nx):
        for j in range(ny):
            ring_centers.append((3.0 * length * i + 0.5 * length, 2.0 * h * j + h))
            if j > 0:
                ring_centers.append((3.0 * length * i + 2.0 * length, 2.0 * h * j))
    cx, cy = min(
        ring_centers,
        key=lambda pt: (pt[0] - target_cx) ** 2 + (pt[1] - target_cy) ** 2,
    )
    if cut_radius is None:
        cut_radius = min(cx, span_x - cx, cy, span_y - cy) + 1e-6

    if shape == "circular":
        keep_ids = {
            a["id"]
            for a in model["atoms"]
            if math.hypot(a["position"][0] - cx, a["position"][1] - cy) <= cut_radius + 1e-6
        }
    elif shape == "hexagonal":
        offset_deg = 0.0 if (edge_type or "zigzag") == "zigzag" else 30.0
        apothem = cut_radius * math.sqrt(3.0) / 2.0
        normals = [
            (
                math.cos(math.radians(offset_deg + 60.0 * k)),
                math.sin(math.radians(offset_deg + 60.0 * k)),
            )
            for k in range(6)
        ]
        keep_ids = {
            a["id"]
            for a in model["atoms"]
            if all(
                (a["position"][0] - cx) * ux + (a["position"][1] - cy) * uy <= apothem + 1e-6
                for ux, uy in normals
            )
        }
    else:
        raise ValueError(f"Unknown shape {shape!r}")

    removed = {a["id"] for a in model["atoms"]} - keep_ids
    _remove_atoms(model, removed)
    _prune_tips(model)
    if len(model["atoms"]) < 4 or not _connected(model):
        raise ValueError(
            "Requested nanoflake shape/size leaves fewer than 4 carbons or a disconnected sheet"
        )
    model["metadata"]["construction"] = f"analytical_{shape}_honeycomb"


def _scaffold_morphology_metrics(model, length: float) -> dict:
    ns = adjacency(model)
    n_c = len(model["atoms"])
    edge_set = {cid for cid, nbrs in ns.items() if len(nbrs) == 2}
    perimeter_carbons = len(edge_set)
    armchair_edge_pairs = 0
    zigzag_edge_carbons = 0
    seen_edge = set()
    for cid in sorted(edge_set):
        if cid in seen_edge:
            continue
        comp = []
        stack = [cid]
        seen_edge.add(cid)
        while stack:
            curr = stack.pop()
            comp.append(curr)
            for nbr, _, _ in ns[curr]:
                if nbr in edge_set and nbr not in seen_edge:
                    seen_edge.add(nbr)
                    stack.append(nbr)
        armchair_edge_pairs += len(comp) // 2
        zigzag_edge_carbons += len(comp) % 2

    area_nm2 = (n_c * (3.0 * math.sqrt(3.0) / 4.0) * (length ** 2)) / 100.0
    perimeter_nm = (
        3.0 * length * armchair_edge_pairs + math.sqrt(3.0) * length * zigzag_edge_carbons
    ) / 10.0
    edge_to_area = (
        (perimeter_nm / area_nm2) if (area_nm2 > 0.0 and perimeter_carbons > 0) else 0.0
    )
    return {
        "area_nm2": round(area_nm2, 6),
        "perimeter_carbons": perimeter_carbons,
        "armchair_edge_pairs": armchair_edge_pairs,
        "zigzag_edge_carbons": zigzag_edge_carbons,
        "edge_to_area_ratio_nm_inv": round(edge_to_area, 6),
    }


def _model_dimensions_nm(model) -> tuple[list[float], list[float]]:
    cell = np.asarray(model["cell"], dtype=float)
    cell_nm = [round(float(np.linalg.norm(cell[k])) / 10.0, 6) for k in range(3)]
    coords = np.asarray([a["position"] for a in model["atoms"]], dtype=float)
    mins = np.min(coords, axis=0)
    maxs = np.max(coords, axis=0)
    extent_nm = [round(float(maxs[k] - mins[k]) / 10.0, 6) for k in range(3)]
    return cell_nm, extent_nm


def build_sheet(recipe: dict) -> dict:
    """Build exact-count GO from honeycomb cells (rectangular, hexagonal, or circular).

    Recipe: nx/ny >=2 (or lx_nm/ly_nm or diameter_nm/radius_nm), shape="rectangular"|"hexagonal"|"circular",
    edge_type="zigzag"|"armchair" (for hexagonal), seed, pbc=[x,y,false],
    bond_length=1.42 Å, vacuum=40 Å, vacancies=count or original ID list,
    holes=[{center:[x,y],radius:r} or {size:n,seed_atom:id}] or hole profile,
    groups={OH,epoxide,carbonyl,COOH,COO,epoxide_eph,epoxide_ept} or
    composition/c_over_o/o_over_c,
    faces=both/up/down, placement="random"|"island"|"dispersed"|{mode:...}.
    """
    if not isinstance(recipe, dict):
        raise ValueError("Recipe must be a dictionary")
    allowed = {
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
        "name",
        "bond_length",
        "vacuum",
        "vacancies",
        "holes",
        "groups",
        "faces",
        "placement",
        "composition",
        "c_over_o",
        "o_over_c",
    }
    if set(recipe) - allowed:
        raise ValueError(f"Unknown recipe keys: {sorted(set(recipe) - allowed)}")

    shape = recipe.get("shape", "rectangular")
    if shape not in ("rectangular", "hexagonal", "circular"):
        raise ValueError("shape must be 'rectangular', 'hexagonal', or 'circular'")
    edge_type = recipe.get("edge_type")
    if edge_type is not None:
        if edge_type not in ("zigzag", "armchair"):
            raise ValueError("edge_type must be 'zigzag' or 'armchair'")
        if shape != "hexagonal":
            raise ValueError("edge_type is only supported with shape='hexagonal'")
    effective_edge_type = (edge_type or "zigzag") if shape == "hexagonal" else "mixed"

    seed = _integer(recipe.get("seed", 0), "seed")
    pbc = recipe.get("pbc", [False, False, False])
    if (
        not isinstance(pbc, list)
        or len(pbc) != 3
        or any(type(x) is not bool for x in pbc)
        or pbc[2]
    ):
        raise ValueError("pbc must be [bool,bool,false]")
    if shape in ("hexagonal", "circular") and any(pbc):
        raise ValueError(f"shape='{shape}' requires non-periodic pbc=[False, False, False]")

    length, vacuum = float(recipe.get("bond_length", 1.42)), float(recipe.get("vacuum", 40))
    if (
        not math.isfinite(length)
        or not 1.30 <= length <= 1.55
        or not math.isfinite(vacuum)
        or vacuum < 15
    ):
        raise ValueError("bond_length must be 1.30–1.55 Å and vacuum >=15 Å")

    has_nx_ny = "nx" in recipe or "ny" in recipe
    lx_nm_raw = recipe.get("lx_nm")
    ly_nm_raw = recipe.get("ly_nm")
    has_lx_ly = lx_nm_raw is not None or ly_nm_raw is not None
    diam_raw = recipe.get("diameter_nm")
    rad_raw = recipe.get("radius_nm")
    has_diam = diam_raw is not None
    has_rad = rad_raw is not None

    if has_diam and has_rad:
        raise ValueError("Specify either diameter_nm or radius_nm, not both")
    if has_nx_ny and (has_lx_ly or has_diam or has_rad):
        raise ValueError("nx/ny cannot be combined with lx_nm/ly_nm/diameter_nm/radius_nm")
    if has_lx_ly and (has_diam or has_rad):
        raise ValueError("lx_nm/ly_nm cannot be combined with diameter_nm/radius_nm")
    if shape == "rectangular" and (has_diam or has_rad):
        raise ValueError("diameter_nm/radius_nm require shape='circular' or 'hexagonal'")

    h_lat = math.sqrt(3.0) * length / 2.0
    dx_cell_nm = (3.0 * length) / 10.0
    dy_cell_nm = (2.0 * h_lat) / 10.0
    use_nm_sizing = has_lx_ly or has_diam or has_rad
    use_morphology_meta = (
        shape != "rectangular" or "shape" in recipe or "edge_type" in recipe or use_nm_sizing
    )
    dim_meta = None
    cut_radius_a = None

    if has_diam or has_rad:
        if has_diam:
            diam_nm = _finite_real(diam_raw, "diameter_nm")
            if diam_nm <= 0.0:
                raise ValueError("diameter_nm must be positive")
            rad_nm = 0.5 * diam_nm
        else:
            rad_nm = _finite_real(rad_raw, "radius_nm")
            if rad_nm <= 0.0:
                raise ValueError("radius_nm must be positive")
            diam_nm = 2.0 * rad_nm
        cut_radius_a = rad_nm * 10.0
        nx = max(2, math.ceil((2.0 * cut_radius_a + 2.0 * length) / (3.0 * length)))
        ny = max(2, math.ceil((2.0 * cut_radius_a + 2.0 * length) / (2.0 * h_lat)))
        dim_meta = {
            "requested_diameter_nm": diam_nm,
            "requested_radius_nm": rad_nm,
            "rounding_policy": "ceil_enclosing_tile_min2",
            "nx": nx,
            "ny": ny,
        }
    elif has_lx_ly:
        if lx_nm_raw is None or ly_nm_raw is None:
            raise ValueError("Both lx_nm and ly_nm must be specified together")
        lx_nm = _finite_real(lx_nm_raw, "lx_nm")
        ly_nm = _finite_real(ly_nm_raw, "ly_nm")
        if lx_nm <= 0.0 or ly_nm <= 0.0:
            raise ValueError("lx_nm and ly_nm must be positive")
        if shape == "rectangular" and any(pbc):
            rounding_policy = "ceil_cell_multiple_min2"
            nx = max(2, math.ceil(lx_nm / dx_cell_nm))
            ny = max(2, math.ceil(ly_nm / dy_cell_nm))
        elif shape == "rectangular":
            rounding_policy = "round_cell_multiple_min2"
            nx = max(2, round(lx_nm / dx_cell_nm))
            ny = max(2, round(ly_nm / dy_cell_nm))
        else:
            rounding_policy = "ceil_enclosing_tile_min2"
            cut_radius_a = 0.5 * min(lx_nm, ly_nm) * 10.0
            nx = max(2, math.ceil((lx_nm * 10.0 + 2.0 * length) / (3.0 * length)))
            ny = max(2, math.ceil((ly_nm * 10.0 + 2.0 * length) / (2.0 * h_lat)))
        dim_meta = {
            "requested_lx_nm": lx_nm,
            "requested_ly_nm": ly_nm,
            "rounding_policy": rounding_policy,
            "nx": nx,
            "ny": ny,
        }
    else:
        nx, ny = _integer(recipe.get("nx", 6), "nx", 2), _integer(recipe.get("ny", 6), "ny", 2)
        if shape != "rectangular":
            dim_meta = {"rounding_policy": "inscribed_in_nx_ny", "nx": nx, "ny": ny}

    requested = recipe.get("groups", {})
    if not isinstance(requested, dict) or set(requested) - set(GROUPS):
        raise ValueError("Unknown functional group")
    comp_spec = _parse_composition_spec(recipe)
    raw_specs, fixed, has_relative = _parse_group_requests(requested)
    faces = recipe.get("faces", "both")
    if faces not in ("both", "up", "down"):
        raise ValueError("faces must be both/up/down")
    placement_spec = _parse_placement(recipe.get("placement"))
    rng = random.Random(seed)
    model = _tile(nx, ny, length, list(pbc), vacuum, seed, str(recipe.get("name", "GO")))
    _apply_shape(
        model,
        nx=nx,
        ny=ny,
        length=length,
        shape=shape,
        edge_type=effective_edge_type,
        cut_radius=cut_radius_a,
    )
    _apply_holes(model, recipe, rng, bond_length=length)

    morph_metrics = (
        _scaffold_morphology_metrics(model, length) if use_morphology_meta else None
    )

    comp_report = None
    if comp_spec is not None:
        from .composition import resolve_composition

        clean_recipe = copy.deepcopy(recipe)
        clean_recipe.pop("composition", None)
        clean_recipe.pop("c_over_o", None)
        clean_recipe.pop("o_over_c", None)
        clean_recipe["groups"] = {}
        base_comp_seed = seed if comp_spec["seed"] is None else comp_spec["seed"]
        tried_counts = set()
        last_error = None
        decorated_model = None
        final_counts = None
        weight_multipliers: dict[str, float] = {}
        has_coupled_carboxyl = (
            comp_spec["ph"] is not None or comp_spec["deprotonation_fraction"] is not None
        )

        for comp_attempt in range(12):
            resolved_counts = None
            attempt_report = None
            for sub_try in range(4):
                w_mult = weight_multipliers if weight_multipliers else None
                if comp_attempt == 0 and sub_try == 0:
                    _, attempt_report = resolve_composition(
                        clean_recipe,
                        c_over_o=comp_spec["c_over_o"],
                        preferences=comp_spec["preferences"],
                        o_over_c=comp_spec["o_over_c"],
                        priority=comp_spec["priority"],
                        mode=comp_spec["mode"],
                        enabled_groups=comp_spec["enabled_groups"],
                        seed=comp_spec["seed"],
                        ph=comp_spec["ph"],
                        pka=comp_spec["pka"],
                        deprotonation_fraction=comp_spec["deprotonation_fraction"],
                        _scaffold=model,
                    )
                elif comp_spec["mode"] == "random":
                    attempt_seed = base_comp_seed + (comp_attempt * 13 + sub_try) * 7919
                    _, attempt_report = resolve_composition(
                        clean_recipe,
                        c_over_o=comp_spec["c_over_o"],
                        preferences=comp_spec["preferences"],
                        o_over_c=comp_spec["o_over_c"],
                        priority=comp_spec["priority"],
                        mode=comp_spec["mode"],
                        enabled_groups=comp_spec["enabled_groups"],
                        seed=attempt_seed,
                        ph=comp_spec["ph"],
                        pka=comp_spec["pka"],
                        deprotonation_fraction=comp_spec["deprotonation_fraction"],
                        _scaffold=model,
                        _weight_multipliers=w_mult,
                    )
                else:
                    jitter_seed = (
                        None
                        if (sub_try == 0 and w_mult is not None)
                        else (seed + 1) * 10007 + comp_attempt * 131 + sub_try
                    )
                    _, attempt_report = resolve_composition(
                        clean_recipe,
                        c_over_o=comp_spec["c_over_o"],
                        preferences=comp_spec["preferences"],
                        o_over_c=comp_spec["o_over_c"],
                        priority=comp_spec["priority"],
                        mode=comp_spec["mode"],
                        enabled_groups=comp_spec["enabled_groups"],
                        seed=comp_spec["seed"],
                        ph=comp_spec["ph"],
                        pka=comp_spec["pka"],
                        deprotonation_fraction=comp_spec["deprotonation_fraction"],
                        _scaffold=model,
                        _jitter_seed=jitter_seed,
                        _weight_multipliers=w_mult,
                    )
                cand_counts = attempt_report["group_counts"]
                sig = tuple(cand_counts[g] for g in GROUPS)
                if sig not in tried_counts or comp_attempt == 0:
                    resolved_counts = cand_counts
                    break
            if resolved_counts is None:
                resolved_counts = attempt_report["group_counts"]
            tried_counts.add(tuple(resolved_counts[g] for g in GROUPS))
            has_p2 = (
                resolved_counts["epoxide_eph"] > 0 or resolved_counts["epoxide_ept"] > 0
            )
            attempt_rng = rng if comp_attempt == 0 else random.Random(seed + comp_attempt * 1009)
            try:
                decorated_model = _decorate(
                    model,
                    resolved_counts,
                    {},
                    faces,
                    attempt_rng,
                    placement_spec=placement_spec,
                    has_p2=has_p2,
                    max_attempts=16 if comp_attempt < 11 else 128,
                )
                final_counts = resolved_counts
                comp_report = attempt_report
                comp_report["composition_attempts"] = comp_attempt + 1
                break
            except PlacementError as exc:
                last_error = exc
                msg = str(exc)
                bottleneck_groups: list[str] = []
                if exc.attempts == 0:
                    if "epoxide count exceeds" in msg:
                        bottleneck_groups = ["epoxide"]
                    elif "epoxide_eph count exceeds" in msg:
                        bottleneck_groups = ["epoxide_eph"]
                    elif "epoxide_ept count exceeds" in msg:
                        bottleneck_groups = ["epoxide_ept"]
                    elif "Island radius is too small" in msg:
                        bottleneck_groups = [g for g in ("OH", "epoxide") if resolved_counts.get(g, 0) > 0]
                    elif "edge/basal sites" in msg:
                        cap = exc.capacity or {}
                        edge_req = (
                            resolved_counts.get("carbonyl", 0)
                            + resolved_counts.get("COOH", 0)
                            + resolved_counts.get("COO", 0)
                            + resolved_counts.get("epoxide_eph", 0)
                            + 2 * resolved_counts.get("epoxide_ept", 0)
                        )
                        if edge_req > cap.get("eligible_edge_carbons", 0):
                            bottleneck_groups = [
                                g
                                for g in ("COOH", "COO", "carbonyl", "epoxide_eph", "epoxide_ept")
                                if resolved_counts.get(g, 0) > 0
                            ]
                        else:
                            bottleneck_groups = [
                                g for g in ("OH", "epoxide", "epoxide_eph") if resolved_counts.get(g, 0) > 0
                            ]
                else:
                    achieved = exc.achieved or {}
                    for g in (
                        "epoxide_ept",
                        "epoxide_eph",
                        "COOH",
                        "COO",
                        "carbonyl",
                        "epoxide",
                        "OH",
                    ):
                        req_g = resolved_counts.get(g, 0)
                        if req_g > 0 and achieved.get(g, 0) < req_g:
                            bottleneck_groups = [g]
                            break
                    if not bottleneck_groups:
                        # All groups were placed, so edge capping or final 1-4+ integrity failed.
                        edge_nonzero = [
                            g
                            for g in ("COO", "COOH", "carbonyl", "epoxide_eph", "epoxide_ept")
                            if resolved_counts.get(g, 0) > 0
                        ]
                        if edge_nonzero:
                            bottleneck_groups = edge_nonzero
                        else:
                            nonzero = [(resolved_counts[g], g) for g in GROUPS if resolved_counts.get(g, 0) > 0]
                            if nonzero:
                                bottleneck_groups = [max(nonzero)[1]]

                scaled = set()
                for bg in bottleneck_groups:
                    if has_coupled_carboxyl and bg in ("COOH", "COO"):
                        scaled.update(("COOH", "COO"))
                    else:
                        scaled.add(bg)
                for sg in scaled:
                    weight_multipliers[sg] = weight_multipliers.get(sg, 1.0) * 0.45

        if decorated_model is None:
            raise last_error
        model = decorated_model
        counts = final_counts
    else:
        counts = _resolve_group_counts(model, raw_specs, has_relative)
        has_p2 = (
            any(g in requested for g in PRIORITY2_GROUPS)
            or counts["epoxide_eph"] > 0
            or counts["epoxide_ept"] > 0
        )
        model = _decorate(
            model,
            counts,
            fixed,
            faces,
            rng,
            placement_spec=placement_spec,
            has_p2=has_p2,
        )

    # GGG is a one-atom library fragment, so each pristine carbon needs a
    # distinct residue while retaining its original atom ID and graph order.
    next_resid = max(a.get("resid", 1) for a in model["atoms"]) + 1
    for atom in model["atoms"]:
        if atom["resname"] == "GGG":
            atom["resid"] = next_resid
            next_resid += 1
    model["metadata"]["recipe"] = copy.deepcopy(recipe)
    event_counts = {g: counts[g] for g in (GROUPS if has_p2 else LEGACY_GROUPS)}
    if has_relative:
        model["metadata"]["resolved_group_counts"] = dict(event_counts)
    if comp_report is not None:
        actual_c = sum(1 for a in model["atoms"] if a["element"] == "C")
        actual_o = sum(1 for a in model["atoms"] if a["element"] == "O")
        comp_report["actual_C"] = actual_c
        comp_report["actual_O"] = actual_o
        comp_report["actual_C_over_O"] = actual_c / actual_o
        comp_report["actual_O_over_C"] = actual_o / actual_c
        model["metadata"]["composition"] = comp_report
    if dim_meta is not None or use_morphology_meta:
        cell_nm, extent_nm = _model_dimensions_nm(model)
        scaffold_ids = {
            a["id"]
            for a in model["atoms"]
            if a["element"] == "C" and a.get("metadata", {}).get("role") == "sheet_carbon"
        }
        scaffold_diam_nm = _compute_scaffold_diameter_nm(model["atoms"], scaffold_ids)
        if dim_meta is not None:
            dim_meta["cell_nm"] = cell_nm
            dim_meta["extent_nm"] = extent_nm
            dim_meta["scaffold_diameter_nm"] = scaffold_diam_nm
            model["metadata"]["dimensions"] = dim_meta
        if use_morphology_meta:
            model["metadata"]["morphology"] = {
                "shape": shape,
                "edge_type": effective_edge_type,
                "extent_nm": extent_nm,
                "scaffold_diameter_nm": scaffold_diam_nm,
                **morph_metrics,
            }
    model["metadata"]["events"].append(
        {
            "operation": "build_sheet",
            "group_counts": event_counts,
            "removed_carbon_ids": model["metadata"]["removed_carbon_ids"],
        }
    )
    model["metadata"].update(electronic_summary(model))
    assert_export_ready(model)
    return model
