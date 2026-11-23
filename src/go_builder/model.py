"""Explicit molecular graph validation and lossless JSON exchange.

Coordinates are Å, atom IDs are stable integers. PDB is a finite visualization
export only; JSON preserves image bonds and is the authoritative graph.
"""

from __future__ import annotations

import itertools
import json
import math
import re
from collections import deque
from pathlib import Path

import numpy as np

ATOMIC_NUMBERS = dict(
    zip(
        "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split(),
        range(1, 119),
        strict=True,
    )
)
ATOMIC_NUMBERS["EP"] = 0


def canonical_bond(a, b, order=1.0, image=(0, 0, 0)):
    """Store a->b physical image displacement with ascending atom IDs."""
    image = [int(x) for x in image]
    if a > b:
        a, b, image = b, a, [-x for x in image]
    return {"a": a, "b": b, "order": float(order), "image": image}


def adjacency(model):
    neighbors = {a["id"]: [] for a in model["atoms"]}
    for bond in model["bonds"]:
        a, b, image = bond["a"], bond["b"], np.array(bond.get("image", [0, 0, 0]))
        neighbors[a].append((b, image, bond["order"]))
        neighbors[b].append((a, -image, bond["order"]))
    return neighbors


def electronic_summary(model):
    """Count electrons from explicit formal charge, without assigning spin."""
    if model.get("metadata", {}).get("formal_charges_unknown"):
        return {
            "electron_count": None,
            "electron_parity": "unknown",
            "electronic_state_not_validated": True,
            "pi_bonding_not_validated": True,
        }
    electrons = sum(
        ATOMIC_NUMBERS[a["element"]] - a.get("formal_charge", 0) for a in model["atoms"]
    )
    return {
        "electron_count": electrons,
        "electron_parity": "odd" if electrons % 2 else "even",
        "electronic_state_not_validated": True,
        "pi_bonding_not_validated": True,
    }


def validate_model(model: dict) -> dict:
    """Validate graph structure without mutating it; raise ValueError on errors.

    Native graphene has delocalized pi bonds. Its carbon chemistry is checked
    by explicit coordination and hybridization, not naive aromatic-order sums.
    This is structural validation, not force-field or electronic validation.
    """
    try:
        if model["schema_version"] != 1:
            raise ValueError("Unsupported schema_version")
        atoms, bonds = model["atoms"], model["bonds"]
        pbc = model["pbc"]
        if len(pbc) != 3 or any(type(x) is not bool for x in pbc):
            raise ValueError("pbc must contain three explicit boolean flags")
        if model["cell"] is None:
            if any(pbc):
                raise ValueError("Periodic graphs require a cell")
            cell = np.zeros((3, 3))  # Algebra only: the serialized cell remains None.
        else:
            cell = np.asarray(model["cell"], dtype=float)
            if (
                cell.shape != (3, 3)
                or not np.isfinite(cell).all()
                or abs(np.linalg.det(cell)) < 1e-8
            ):
                raise ValueError("Cell must be finite and nonsingular")
        imported = model.get("metadata", {}).get("representation") == "amber_import"
        orders_known = not model.get("metadata", {}).get("bond_orders_unknown", False)
        if not atoms:
            raise ValueError("Empty atom graph")
        ids = [a["id"] for a in atoms]
        if any(type(i) is not int or i < 1 for i in ids) or len(ids) != len(set(ids)):
            raise ValueError("Atom IDs must be unique positive integers")
        by_id = {a["id"]: a for a in atoms}
        for atom in atoms:
            if atom["element"] not in ATOMIC_NUMBERS:
                raise ValueError("Unsupported or missing element")
            position = np.asarray(atom["position"], dtype=float)
            if position.shape != (3,) or not np.isfinite(position).all():
                raise ValueError("Atom position must contain three finite coordinates")
            if atom.get("charge") is not None and not math.isfinite(float(atom["charge"])):
                raise ValueError("Nonfinite charge")
            if type(atom.get("formal_charge", 0)) is not int:
                raise ValueError("Formal charge must be integer")
        seen = set()
        for bond in bonds:
            a, b = bond["a"], bond["b"]
            if a not in by_id or b not in by_id or a >= b:
                raise ValueError("Bonds require existing, distinct, ascending atom IDs")
            image = bond.get("image", [0, 0, 0])
            if len(image) != 3 or any(type(x) is not int for x in image):
                raise ValueError("Bond images must be three integers")
            if any(image[k] and not pbc[k] for k in range(3)):
                raise ValueError("Bond image crosses nonperiodic cell direction")
            key = a, b, tuple(image)
            if key in seen:
                raise ValueError("Duplicate bond")
            seen.add(key)
            if bond["order"] not in (1, 1.5, 2, 3):
                raise ValueError("Unsupported bond order")
            d = (
                np.array(by_id[b]["position"])
                + np.array(image) @ cell
                - np.array(by_id[a]["position"])
            )
            if not imported and not 0.3 < np.linalg.norm(d) < 4.0:
                raise ValueError("Implausible bonded distance")
        neighbors = adjacency(model)
        undercoordinated = []
        for atom in atoms:
            ns = neighbors[atom["id"]]
            degree = len(ns)
            element = atom["element"]
            role = atom.get("metadata", {}).get("role")
            hybrid = atom.get("metadata", {}).get("hybridization")
            if role == "sheet_carbon" and degree != (4 if hybrid == "sp3" else 3):
                raise ValueError(f"Invalid graphene carbon coordination at {atom['id']}")
            max_degree = {"H": 1, "O": 2, "N": 4, "C": 4, "F": 1, "Cl": 1, "Br": 1, "I": 1}.get(
                element, 6
            )
            if degree > max_degree:
                raise ValueError(f"Excess coordination at {atom['id']}")
            if element == "H" and degree != 1:
                raise ValueError("Hydrogen must have one explicit bond")
            valence = sum(n[2] for n in ns)
            if element == "H" and orders_known and valence != 1:
                raise ValueError("Hydrogen requires one single bond")
            target = {"H": 1, "C": 4, "N": 3, "O": 2, "S": 2}.get(element)
            if atom.get("formal_charge", 0) == -1 and element == "O":
                target = 1
            if atom.get("formal_charge", 0) == 1 and element == "N":
                target = 4
            if (
                orders_known
                and role != "sheet_carbon"
                and target is not None
                and valence < target - 0.01
            ):
                undercoordinated.append(atom["id"])
            if orders_known and element == "O" and valence > 2.01:
                raise ValueError("Excess oxygen valence")
            fused_aromatic = degree == 3 and all(n[2] == 1.5 for n in ns)
            if (
                orders_known
                and element == "C"
                and role != "sheet_carbon"
                and valence > 4.01
                and not fused_aromatic
            ):
                raise ValueError("Excess carbon valence")
        remaining = set(ids)
        components = 0
        while remaining:
            components += 1
            stack = [remaining.pop()]
            while stack:
                for other, _, _ in neighbors[stack.pop()]:
                    if other in remaining:
                        remaining.remove(other)
                        stack.append(other)
        return {
            "valid": True,
            "atom_count": len(atoms),
            "bond_count": len(bonds),
            "connected_components": components,
            "formal_charge": sum(a.get("formal_charge", 0) for a in atoms),
            "undercoordinated_atom_ids": sorted(undercoordinated),
            **electronic_summary(model),
            "warnings": [
                "Constructed graph validation does not establish equilibration or parameter validity."
            ],
        }
    except (KeyError, TypeError, IndexError) as exc:
        raise ValueError(f"Malformed graph: {exc}") from exc


def _reduce_periodic_basis(basis: np.ndarray) -> np.ndarray:
    """Gauss/greedy integer reduction of 1..3 periodic lattice vectors."""
    red = np.array(basis, dtype=float, copy=True)
    m = len(red)
    if m <= 1:
        return red
    for _ in range(32):
        changed = False
        for i in range(m):
            for j in range(m):
                if i == j:
                    continue
                norm2 = float(np.dot(red[j], red[j]))
                if norm2 <= 1e-16:
                    continue
                mu = round(float(np.dot(red[i], red[j])) / norm2)
                if mu != 0:
                    cand = red[i] - mu * red[j]
                    if float(np.dot(cand, cand)) < float(np.dot(red[i], red[i])) - 1e-12:
                        red[i] = cand
                        changed = True
        if not changed:
            break
    return red


def _minimum_image_delta(
    delta: np.ndarray, pbc: list[bool], cell: np.ndarray, inv_cell: np.ndarray | None
) -> np.ndarray:
    if not any(pbc):
        return delta
    # Fast exact path for orthogonal cells
    if np.max(np.abs(cell - np.diag(np.diag(cell)))) <= 1e-12:
        fractional = delta @ inv_cell
        for k, periodic in enumerate(pbc):
            if periodic:
                fractional[..., k] -= np.rint(fractional[..., k])
        return fractional @ cell

    p_axes = [k for k, flag in enumerate(pbc) if flag]
    b_red = _reduce_periodic_basis(cell[p_axes])
    d_p = len(p_axes)
    gram_inv = np.linalg.inv(b_red @ b_red.T)
    orig_shape = delta.shape
    flat = delta.reshape(-1, 3)
    coeff = (flat @ b_red.T) @ gram_inv
    r0 = flat - np.rint(coeff) @ b_red
    shifts = np.asarray(list(itertools.product((-1, 0, 1), repeat=d_p)), dtype=float)
    shift_vecs = shifts @ b_red
    cands = r0[None, :, :] + shift_vecs[:, None, :]
    norms2 = np.sum(cands * cands, axis=2)
    best_idx = np.argmin(norms2, axis=0)
    best = cands[best_idx, np.arange(len(flat))]
    return best.reshape(orig_shape)


_MOTIF_BOND_RANGES: dict[tuple[str, str], tuple[float, float]] = {
    ("H", "O"): (0.90, 1.06),
    ("C", "H"): (1.00, 1.18),
    ("C", "O"): (1.12, 1.56),
    ("C", "C"): (1.30, 1.65),
}

_CO_MOTIF_RANGES: dict[str, tuple[float, float]] = {
    "C=O": (1.12, 1.32),
    "COO-": (1.18, 1.35),
    "COOH_C-OH": (1.25, 1.45),
    "C-O_sp3": (1.32, 1.56),
}


def _classify_co_bond(
    ia: int,
    ib: int,
    bond: dict,
    atoms: list[dict],
    bonded_sets: list[set[int]],
) -> tuple[str, tuple[float, float]]:
    ic, io = (ia, ib) if atoms[ia]["element"] == "C" else (ib, ia)
    ac = atoms[ic]
    ao = atoms[io]
    order = float(bond.get("order", 1.0))
    c_meta = ac.get("metadata", {})
    o_meta = ao.get("metadata", {})
    c_group = c_meta.get("group")
    o_group = o_meta.get("group")
    c_res = ac.get("resname")
    o_res = ao.get("resname")

    if o_group == "COO" or o_res == "CBD" or c_group == "COO" or c_res == "CBD":
        return "COO-", _CO_MOTIF_RANGES["COO-"]
    if o_group == "COOH" or o_res == "CBX" or c_group == "COOH" or c_res == "CBX":
        if any(atoms[nbr]["element"] == "H" for nbr in bonded_sets[io]) or ao.get("name") == "OK":
            return "COOH_C-OH", _CO_MOTIF_RANGES["COOH_C-OH"]
        return "C=O", _CO_MOTIF_RANGES["C=O"]
    if o_group == "carbonyl" or o_res == "CAB" or c_group == "carbonyl" or c_res == "CAB" or order == 2.0:
        return "C=O", _CO_MOTIF_RANGES["C=O"]
    if o_group in ("OH", "epoxide", "epoxide_eph", "epoxide_ept") or o_res in (
        "HYD",
        "EPX",
        "EPH",
        "EPT",
    ):
        return "C-O_sp3", _CO_MOTIF_RANGES["C-O_sp3"]

    c_o_nbrs = [nbr for nbr in bonded_sets[ic] if atoms[nbr]["element"] == "O"]
    if len(bonded_sets[ic]) == 3 and len(c_o_nbrs) >= 2:
        if any(atoms[nbr].get("formal_charge", 0) == -1 for nbr in c_o_nbrs) or ac.get(
            "formal_charge", 0
        ) == -1:
            return "COO-", _CO_MOTIF_RANGES["COO-"]
        if any(atoms[nbr]["element"] == "H" for nbr in bonded_sets[io]):
            return "COOH_C-OH", _CO_MOTIF_RANGES["COOH_C-OH"]
        if len(bonded_sets[io]) >= 2:
            return "C-O_sp3", _CO_MOTIF_RANGES["C-O_sp3"]
        return "C=O", _CO_MOTIF_RANGES["C=O"]
    return "C-O_sp3", _CO_MOTIF_RANGES["C-O_sp3"]


def _build_spatial_cells(
    coords: np.ndarray,
    pbc: list[bool],
    cell: np.ndarray,
    inv_cell: np.ndarray | None,
    cutoff: float,
) -> tuple[dict[tuple[int, int, int], np.ndarray], list[int | None]]:
    """Assign each atom to a 3D cell of width >= cutoff (PBC-aware, triclinic-safe)."""
    n_atoms = len(coords)
    if not any(pbc):
        cell_idx = np.floor(coords / cutoff).astype(np.int64)
        n_bins: list[int | None] = [None, None, None]
    else:
        frac = coords @ inv_cell
        cell_idx = np.empty((n_atoms, 3), dtype=np.int64)
        n_bins = [None, None, None]
        for k in range(3):
            h_k = 1.0 / float(np.linalg.norm(inv_cell[:, k]))
            if pbc[k]:
                f_k = frac[:, k] - np.floor(frac[:, k])
                m_k = max(1, int(math.floor(h_k / cutoff)))
                cell_idx[:, k] = np.floor(f_k * m_k).astype(np.int64) % m_k
                n_bins[k] = m_k
            else:
                cell_idx[:, k] = np.floor(frac[:, k] * (h_k / cutoff)).astype(np.int64)

    buckets: dict[tuple[int, int, int], list[int]] = {}
    for idx, row in enumerate(cell_idx.tolist()):
        key = (row[0], row[1], row[2])
        bucket = buckets.get(key)
        if bucket is None:
            buckets[key] = [idx]
        else:
            bucket.append(idx)
    cells = {key: np.asarray(idxs, dtype=np.int64) for key, idxs in buckets.items()}
    return cells, n_bins


def _record_pair_violation(
    i: int,
    j: int,
    d: float,
    is_h_pair: bool,
    is_13: bool,
    min_h_distance: float,
    min_heavy_distance: float,
    min_14_h_distance: float,
    min_14_heavy_distance: float,
    atom_ids: np.ndarray,
    elements: np.ndarray,
    clashes: list[dict],
    max_clashes: int,
) -> bool:
    hard_thresh = min_h_distance if is_h_pair else min_heavy_distance
    if d < hard_thresh:
        if len(clashes) < max_clashes:
            kind = "non12_h" if is_h_pair else "non12_heavy"
            clashes.append(
                {
                    "a": int(atom_ids[i]),
                    "b": int(atom_ids[j]),
                    "elements": [str(elements[i]), str(elements[j])],
                    "distance": float(d),
                    "threshold": float(hard_thresh),
                    "cutoff": float(hard_thresh),
                    "kind": kind,
                }
            )
        return True
    if not is_13:
        thresh_14 = min_14_h_distance if is_h_pair else min_14_heavy_distance
        if d < thresh_14:
            if len(clashes) < max_clashes:
                kind = "14plus_h" if is_h_pair else "14plus_heavy"
                clashes.append(
                    {
                        "a": int(atom_ids[i]),
                        "b": int(atom_ids[j]),
                        "elements": [str(elements[i]), str(elements[j])],
                        "distance": float(d),
                        "threshold": float(thresh_14),
                        "cutoff": float(thresh_14),
                        "kind": kind,
                    }
                )
            return True
    return False


def check_geometry_integrity(
    model: dict,
    *,
    min_h_distance: float = 0.90,
    min_heavy_distance: float = 1.25,
    min_14_h_distance: float = 1.20,
    min_14_heavy_distance: float = 1.60,
    check_bond_lengths: bool = True,
    max_clashes: int = 20,
) -> dict:
    """Verify nonbonded steric clearances (1-2/1-3/1-4+), bond lengths, coordination, and charge."""
    if (
        min_h_distance <= 0
        or min_heavy_distance <= 0
        or min_14_h_distance <= 0
        or min_14_heavy_distance <= 0
    ):
        raise ValueError("Distance thresholds must be positive")
    validation = validate_model(model)
    atoms = model["atoms"]
    bonds = model["bonds"]
    pbc = model["pbc"]
    imported = model.get("metadata", {}).get("representation") == "amber_import"
    cell = (
        np.zeros((3, 3), dtype=float)
        if model["cell"] is None
        else np.asarray(model["cell"], dtype=float)
    )
    inv_cell = np.linalg.inv(cell) if any(pbc) else None

    coords = np.asarray([atom["position"] for atom in atoms], dtype=float)
    atom_ids = np.fromiter((atom["id"] for atom in atoms), dtype=np.int64, count=len(atoms))
    elements = np.asarray([atom["element"] for atom in atoms])
    is_h = elements == "H"
    is_ep = elements == "EP"
    id_to_index = {int(aid): idx for idx, aid in enumerate(atom_ids)}

    if bonds:
        idx_a = np.fromiter(
            (id_to_index[b["a"]] for b in bonds), dtype=np.int64, count=len(bonds)
        )
        idx_b = np.fromiter(
            (id_to_index[b["b"]] for b in bonds), dtype=np.int64, count=len(bonds)
        )
        images = np.asarray([b.get("image", [0, 0, 0]) for b in bonds], dtype=float)
        bond_vecs = coords[idx_b] + images @ cell - coords[idx_a]
        bond_lengths = np.linalg.norm(bond_vecs, axis=1)
        min_bond_length = float(np.min(bond_lengths))
        max_bond_length = float(np.max(bond_lengths))
    else:
        idx_a = np.empty(0, dtype=np.int64)
        idx_b = np.empty(0, dtype=np.int64)
        bond_lengths = np.empty(0, dtype=float)
        min_bond_length = None
        max_bond_length = None

    n_atoms = len(atoms)
    bonded_sets: list[set[int]] = [set() for _ in range(n_atoms)]
    for ia, ib in zip(idx_a.tolist(), idx_b.tolist(), strict=True):
        bonded_sets[ia].add(ib)
        bonded_sets[ib].add(ia)

    pairs_13_sets: list[set[int]] = [set() for _ in range(n_atoms)]
    for center in range(n_atoms):
        nbrs = list(bonded_sets[center])
        n_nbr = len(nbrs)
        for u_idx in range(n_nbr):
            u = nbrs[u_idx]
            u_bonded = bonded_sets[u]
            for v_idx in range(u_idx + 1, n_nbr):
                v = nbrs[v_idx]
                if v not in u_bonded:
                    pairs_13_sets[u].add(v)
                    pairs_13_sets[v].add(u)

    min_nonbonded_h: float | None = None
    min_nonbonded_heavy: float | None = None
    min_14_h: float | None = None
    min_14_heavy: float | None = None
    nonbonded_minima_cutoff: float | None = None
    clash_count = 0
    clashes: list[dict] = []

    if n_atoms <= 256:
        for i in range(n_atoms):
            if i + 1 >= n_atoms:
                continue
            if imported and is_ep[i]:
                continue
            delta = coords[i + 1 :] - coords[i]
            delta = _minimum_image_delta(delta, pbc, cell, inv_cell)
            dists = np.linalg.norm(delta, axis=1)
            b_set = bonded_sets[i]
            s13 = pairs_13_sets[i]
            for offset, d_val in enumerate(dists.tolist()):
                j = i + 1 + offset
                if j in b_set:
                    continue
                if imported and is_ep[j]:
                    continue
                is_h_pair = bool(is_h[i] or is_h[j])
                is_13 = j in s13
                if is_h_pair:
                    if min_nonbonded_h is None or d_val < min_nonbonded_h:
                        min_nonbonded_h = d_val
                    if not is_13 and (min_14_h is None or d_val < min_14_h):
                        min_14_h = d_val
                else:
                    if min_nonbonded_heavy is None or d_val < min_nonbonded_heavy:
                        min_nonbonded_heavy = d_val
                    if not is_13 and (min_14_heavy is None or d_val < min_14_heavy):
                        min_14_heavy = d_val
                if _record_pair_violation(
                    i,
                    j,
                    d_val,
                    is_h_pair,
                    is_13,
                    min_h_distance,
                    min_heavy_distance,
                    min_14_h_distance,
                    min_14_heavy_distance,
                    atom_ids,
                    elements,
                    clashes,
                    max_clashes,
                ):
                    clash_count += 1
    else:
        cell_cutoff = max(
            3.2, min_h_distance, min_heavy_distance, min_14_h_distance, min_14_heavy_distance
        )
        nonbonded_minima_cutoff = float(cell_cutoff)
        cells, n_bins = _build_spatial_cells(coords, pbc, cell, inv_cell, cell_cutoff)
        for c_key, atoms_c in cells.items():
            offs = []
            for k in range(3):
                b_k = c_key[k]
                m_k = n_bins[k]
                if m_k is not None:
                    offs.append(sorted({(b_k + d) % m_k for d in (-1, 0, 1)}))
                else:
                    offs.append((b_k - 1, b_k, b_k + 1))

            for nx_b in offs[0]:
                for ny_b in offs[1]:
                    for nz_b in offs[2]:
                        n_key = (nx_b, ny_b, nz_b)
                        if n_key < c_key:
                            continue
                        atoms_n = cells.get(n_key)
                        if atoms_n is None:
                            continue
                        if n_key == c_key:
                            if len(atoms_c) < 2:
                                continue
                            delta = coords[atoms_c][None, :, :] - coords[atoms_c][:, None, :]
                            delta = _minimum_image_delta(delta, pbc, cell, inv_cell)
                            dists = np.linalg.norm(delta, axis=2)
                            u_idx, v_idx = np.triu_indices(len(atoms_c), k=1)
                            pair_dists = dists[u_idx, v_idx]
                            close_mask = pair_dists <= cell_cutoff
                            for k_pos in np.flatnonzero(close_mask).tolist():
                                i = int(atoms_c[u_idx[k_pos]])
                                j = int(atoms_c[v_idx[k_pos]])
                                if i > j:
                                    i, j = j, i
                                if j in bonded_sets[i]:
                                    continue
                                if imported and (is_ep[i] or is_ep[j]):
                                    continue
                                d_val = float(pair_dists[k_pos])
                                is_h_pair = bool(is_h[i] or is_h[j])
                                is_13 = j in pairs_13_sets[i]
                                if is_h_pair:
                                    if min_nonbonded_h is None or d_val < min_nonbonded_h:
                                        min_nonbonded_h = d_val
                                    if not is_13 and (min_14_h is None or d_val < min_14_h):
                                        min_14_h = d_val
                                else:
                                    if min_nonbonded_heavy is None or d_val < min_nonbonded_heavy:
                                        min_nonbonded_heavy = d_val
                                    if not is_13 and (
                                        min_14_heavy is None or d_val < min_14_heavy
                                    ):
                                        min_14_heavy = d_val
                                if _record_pair_violation(
                                    i,
                                    j,
                                    d_val,
                                    is_h_pair,
                                    is_13,
                                    min_h_distance,
                                    min_heavy_distance,
                                    min_14_h_distance,
                                    min_14_heavy_distance,
                                    atom_ids,
                                    elements,
                                    clashes,
                                    max_clashes,
                                ):
                                    clash_count += 1
                        else:
                            delta = coords[atoms_n][None, :, :] - coords[atoms_c][:, None, :]
                            delta = _minimum_image_delta(delta, pbc, cell, inv_cell)
                            dists = np.linalg.norm(delta, axis=2)
                            u_idx, v_idx = np.nonzero(dists <= cell_cutoff)
                            for k_pos in range(len(u_idx)):
                                r_u = int(u_idx[k_pos])
                                r_v = int(v_idx[k_pos])
                                i = int(atoms_c[r_u])
                                j = int(atoms_n[r_v])
                                if i > j:
                                    i, j = j, i
                                if j in bonded_sets[i]:
                                    continue
                                if imported and (is_ep[i] or is_ep[j]):
                                    continue
                                d_val = float(dists[r_u, r_v])
                                is_h_pair = bool(is_h[i] or is_h[j])
                                is_13 = j in pairs_13_sets[i]
                                if is_h_pair:
                                    if min_nonbonded_h is None or d_val < min_nonbonded_h:
                                        min_nonbonded_h = d_val
                                    if not is_13 and (min_14_h is None or d_val < min_14_h):
                                        min_14_h = d_val
                                else:
                                    if min_nonbonded_heavy is None or d_val < min_nonbonded_heavy:
                                        min_nonbonded_heavy = d_val
                                    if not is_13 and (
                                        min_14_heavy is None or d_val < min_14_heavy
                                    ):
                                        min_14_heavy = d_val
                                if _record_pair_violation(
                                    i,
                                    j,
                                    d_val,
                                    is_h_pair,
                                    is_13,
                                    min_h_distance,
                                    min_heavy_distance,
                                    min_14_h_distance,
                                    min_14_heavy_distance,
                                    atom_ids,
                                    elements,
                                    clashes,
                                    max_clashes,
                                ):
                                    clash_count += 1

    enforce_bonds = bool(check_bond_lengths and not imported and bonds)
    bond_lengths_valid = True
    bond_length_violations: list[dict] = []
    bond_tol = 1e-5
    if enforce_bonds:
        bond_lengths_valid = (0.85 - bond_tol) <= min_bond_length and max_bond_length <= (1.85 + bond_tol)
        for k in range(len(bonds)):
            d_val = float(bond_lengths[k])
            ia = int(idx_a[k])
            ib = int(idx_b[k])
            ea = str(elements[ia])
            eb = str(elements[ib])
            pair_key = (ea, eb) if ea <= eb else (eb, ea)
            if pair_key == ("C", "O"):
                motif_name, (lo, hi) = _classify_co_bond(ia, ib, bonds[k], atoms, bonded_sets)
            else:
                motif_name = f"{pair_key[0]}-{pair_key[1]}"
                lo, hi = _MOTIF_BOND_RANGES.get(pair_key, (0.85, 1.85))
            lo = max(0.85, lo)
            hi = min(1.85, hi)
            if not ((lo - bond_tol) <= d_val <= (hi + bond_tol)):
                bond_lengths_valid = False
                if len(bond_length_violations) < max_clashes:
                    bond_length_violations.append(
                        {
                            "a": int(atom_ids[ia]),
                            "b": int(atom_ids[ib]),
                            "elements": [ea, eb],
                            "motif": motif_name,
                            "distance": d_val,
                            "expected_range": [lo, hi],
                        }
                    )

    undercoordinated = validation["undercoordinated_atom_ids"]
    formal_charges_unknown = bool(model.get("metadata", {}).get("formal_charges_unknown"))
    formal_charge = None if formal_charges_unknown else validation["formal_charge"]
    charge_valid = formal_charges_unknown or (type(formal_charge) is int)
    elec = electronic_summary(model)
    warnings: list[str] = []
    if elec["electron_parity"] == "odd":
        warnings.append(
            f"Odd total electron count ({elec['electron_count']}): open-shell system "
            "(valid for MM and open-shell/doublet QM, incompatible with closed-shell singlet QM)."
        )

    passed = (
        clash_count == 0
        and len(undercoordinated) == 0
        and bond_lengths_valid
        and charge_valid
    )

    return {
        "passed": bool(passed),
        "atom_count": n_atoms,
        "bond_count": len(bonds),
        "min_bond_length": min_bond_length,
        "max_bond_length": max_bond_length,
        "bond_lengths_valid": bool(bond_lengths_valid),
        "bond_length_violations": bond_length_violations,
        "min_nonbonded_h_distance": min_nonbonded_h,
        "min_nonbonded_heavy_distance": min_nonbonded_heavy,
        "min_14_h_distance": min_14_h,
        "min_14_heavy_distance": min_14_heavy,
        "nonbonded_minima_cutoff": nonbonded_minima_cutoff,
        "nonbonded_minima_cutoff_limited": bool(nonbonded_minima_cutoff is not None),
        "clash_count": clash_count,
        "clashes": clashes,
        "undercoordinated_atom_ids": undercoordinated,
        "formal_charge": formal_charge,
        "electron_count": elec["electron_count"],
        "electron_parity": elec["electron_parity"],
        "warnings": warnings,
        "thresholds": {
            "min_h_distance": float(min_h_distance),
            "min_heavy_distance": float(min_heavy_distance),
            "min_14_h_distance": float(min_14_h_distance),
            "min_14_heavy_distance": float(min_14_heavy_distance),
            "min_bond_length": 0.85 if enforce_bonds else None,
            "max_bond_length": 1.85 if enforce_bonds else None,
            "motif_bond_ranges": (
                {
                    **{f"{a}-{b}": [lo, hi] for (a, b), (lo, hi) in _MOTIF_BOND_RANGES.items()},
                    **{motif: [lo, hi] for motif, (lo, hi) in _CO_MOTIF_RANGES.items()},
                }
                if enforce_bonds
                else None
            ),
        },
    }


def _rotating_calipers_diameter_sq(poly: np.ndarray) -> float:
    """Compute exact squared diameter of a CCW-ordered 2D convex polygon in O(h) time."""
    m = len(poly)
    if m < 2:
        return 0.0
    if m == 2:
        dx = float(poly[1, 0] - poly[0, 0])
        dy = float(poly[1, 1] - poly[0, 1])
        return dx * dx + dy * dy

    def _dist2(a: int, b: int) -> float:
        dx = float(poly[a, 0] - poly[b, 0])
        dy = float(poly[a, 1] - poly[b, 1])
        return dx * dx + dy * dy

    def _area2(i_idx: int, ni_idx: int, k_idx: int) -> float:
        ex = float(poly[ni_idx, 0] - poly[i_idx, 0])
        ey = float(poly[ni_idx, 1] - poly[i_idx, 1])
        return abs(
            ex * float(poly[k_idx, 1] - poly[i_idx, 1])
            - ey * float(poly[k_idx, 0] - poly[i_idx, 0])
        )

    max_d2 = 0.0
    j = 1
    for i in range(m):
        ni = (i + 1) % m
        cur = _area2(i, ni, j)
        steps = 0
        while steps < m:
            nj = (j + 1) % m
            nxt = _area2(i, ni, nj)
            if nxt > cur + 1e-12:
                j = nj
                cur = nxt
                steps += 1
            else:
                break
        nj = (j + 1) % m
        for a, b in ((i, j), (ni, j), (i, nj), (ni, nj)):
            d2 = _dist2(a, b)
            if d2 > max_d2:
                max_d2 = d2
    return max_d2


def _compute_scaffold_diameter_nm(atoms: list[dict], scaffold_ids: set[int]) -> float:
    target_ids = scaffold_ids or {a["id"] for a in atoms if a["element"] == "C"}
    if len(target_ids) < 2:
        return 0.0
    sc = np.asarray([a["position"] for a in atoms if a["id"] in target_ids], dtype=float)
    # Center coordinates before any distance/hull math to avoid catastrophic cancellation
    # under large rigid translations (e.g., 10^8 A shifts).
    sc = sc - np.mean(sc, axis=0)
    if len(sc) == 2:
        diff = sc[1] - sc[0]
        return round(math.sqrt(float(np.dot(diff, diff))) / 10.0, 6)

    from scipy.spatial import ConvexHull, QhullError

    z_span = float(np.ptp(sc[:, 2]))
    if z_span < 1e-9:
        try:
            hull = ConvexHull(sc[:, :2])
            poly = np.asarray(sc[hull.vertices, :2], dtype=float)
            max_d2 = _rotating_calipers_diameter_sq(poly)
            return round(math.sqrt(max(0.0, max_d2)) / 10.0, 6)
        except QhullError:
            pass
        pts = sc
    else:
        if len(sc) >= 4:
            try:
                hull = ConvexHull(sc)
                pts = sc[hull.vertices]
            except QhullError:
                pts = sc
        else:
            pts = sc

    max_d2 = 0.0
    chunk = 256
    for start in range(0, len(pts), chunk):
        sub = pts[start : start + chunk]
        diff = sub[:, None, :] - pts[None, :, :]
        d2 = np.einsum("ijk,ijk->ij", diff, diff)
        cur_max = float(np.max(d2))
        if cur_max > max_d2:
            max_d2 = cur_max
    return round(math.sqrt(max(0.0, max_d2)) / 10.0, 6)


def compute_chemistry_summary(model: dict) -> dict:
    """Compute structural chemistry, sp2 domains, and QM diagnostics on any valid model."""
    atoms = model["atoms"]
    metadata = model.get("metadata", {})
    n_atoms = len(atoms)

    atom_counts: dict[str, int] = {}
    for atom in atoms:
        el = str(atom["element"])
        atom_counts[el] = atom_counts.get(el, 0) + 1
    atom_counts = dict(sorted(atom_counts.items()))

    # Hill-style formula: C first, H second if C present, then alphabetical
    formula_parts: list[str] = []
    if "C" in atom_counts:
        for el in ["C"] + (["H"] if "H" in atom_counts else []) + [
            k for k in atom_counts if k not in ("C", "H")
        ]:
            cnt = atom_counts[el]
            formula_parts.append(f"{el}{cnt if cnt > 1 else ''}")
    else:
        for el, cnt in atom_counts.items():
            formula_parts.append(f"{el}{cnt if cnt > 1 else ''}")
    formula = "".join(formula_parts)

    atomic_fractions = {
        el: (cnt / n_atoms if n_atoms > 0 else 0.0) for el, cnt in atom_counts.items()
    }
    n_c = atom_counts.get("C", 0)
    n_o = atom_counts.get("O", 0)
    c_over_o = (n_c / n_o) if n_o > 0 else None
    o_over_c = (n_o / n_c) if n_c > 0 else None

    formal_charges_unknown = bool(metadata.get("formal_charges_unknown"))
    formal_charge = (
        None
        if formal_charges_unknown
        else sum(int(a.get("formal_charge", 0)) for a in atoms)
    )

    coords = (
        np.asarray([a["position"] for a in atoms], dtype=float)
        if atoms
        else np.zeros((0, 3), dtype=float)
    )
    if len(coords):
        spans = (np.max(coords, axis=0) - np.min(coords, axis=0)) / 10.0
        extent_nm = [float(spans[0]), float(spans[1]), float(spans[2])]
    else:
        extent_nm = [0.0, 0.0, 0.0]

    cell = model.get("cell")
    if cell is not None:
        cell_arr = np.asarray(cell, dtype=float)
        cell_nm = [float(np.linalg.norm(cell_arr[k]) / 10.0) for k in range(3)]
    else:
        cell_nm = None

    by_id = {a["id"]: a for a in atoms}
    adj = adjacency(model)

    # Identify scaffold carbons and sp2/sp3 counts + chemical environments
    scaffold_ids: set[int] = set()
    sp2_scaffold_ids: set[int] = set()
    sp3_scaffold_ids: set[int] = set()

    env_counts = {
        "C-C_sp2": 0,
        "C-O_sp3": 0,
        "C=O_carbonyl": 0,
        "O-C=O_carboxyl": 0,
    }
    carboxyl_c4 = 0
    carboxyl_cy = 0
    other_c = 0

    for atom in atoms:
        if atom["element"] != "C":
            continue
        aid = atom["id"]
        ameta = atom.get("metadata", {})
        role = ameta.get("role")
        group = ameta.get("group")
        resname = atom.get("resname")
        ns = adj[aid]
        deg = len(ns)

        is_sheet = role == "sheet_carbon"
        if role is None:
            # Heuristic fallback for imported/external models
            c_nbrs = sum(1 for nbr, _, _ in ns if by_id[nbr]["element"] == "C")
            o_nbrs = sum(1 for nbr, _, _ in ns if by_id[nbr]["element"] == "O")
            is_sheet = c_nbrs >= 2 and not (o_nbrs >= 2 and deg == 3)

        if is_sheet:
            scaffold_ids.add(aid)
            hybrid = ameta.get("hybridization")
            is_sp3 = (hybrid == "sp3") or (hybrid is None and deg == 4)
            if is_sp3:
                sp3_scaffold_ids.add(aid)
            else:
                sp2_scaffold_ids.add(aid)

        # Carbon chemical environment classification
        if group in ("OH", "epoxide", "epoxide_eph", "epoxide_ept") or resname in (
            "HYD",
            "EPX",
            "EPH",
            "EPT",
        ):
            env_counts["C-O_sp3"] += 1
        elif group == "carbonyl" or resname == "CAB":
            env_counts["C=O_carbonyl"] += 1
        elif (group in ("COOH", "COO") or resname in ("CBX", "CBD")) and role == "group":
            env_counts["O-C=O_carboxyl"] += 1
            carboxyl_c4 += 1
        elif (group in ("COOH", "COO") or resname in ("CBX", "CBD")) and is_sheet:
            env_counts["C-C_sp2"] += 1
            carboxyl_cy += 1
        elif is_sheet:
            # Check fallback for imported models with O neighbors
            if role is None:
                o_nbrs = [nbr for nbr, _, _ in ns if by_id[nbr]["element"] == "O"]
                if len(o_nbrs) >= 2 and deg == 3:
                    env_counts["O-C=O_carboxyl"] += 1
                    carboxyl_c4 += 1
                elif any(order == 2.0 for nbr, _, order in ns if by_id[nbr]["element"] == "O"):
                    env_counts["C=O_carbonyl"] += 1
                elif deg == 4 and o_nbrs:
                    env_counts["C-O_sp3"] += 1
                elif any(
                    len([n2 for n2, _, _ in adj[nbr] if by_id[n2]["element"] == "O"]) >= 2
                    for nbr, _, _ in ns
                    if by_id[nbr]["element"] == "C"
                ):
                    env_counts["C-C_sp2"] += 1
                    carboxyl_cy += 1
                else:
                    env_counts["C-C_sp2"] += 1
            else:
                env_counts["C-C_sp2"] += 1
        else:
            o_nbrs = sum(1 for nbr, _, _ in ns if by_id[nbr]["element"] == "O")
            if o_nbrs >= 2 and deg == 3:
                env_counts["O-C=O_carboxyl"] += 1
                carboxyl_c4 += 1
            else:
                other_c += 1

    if other_c > 0:
        env_counts["other"] = other_c

    scaffold_diameter_nm = _compute_scaffold_diameter_nm(atoms, scaffold_ids)

    n_scaffold = len(scaffold_ids)
    n_sp2 = len(sp2_scaffold_ids)
    n_sp3 = len(sp3_scaffold_ids)
    sp2_fraction = (n_sp2 / n_scaffold) if n_scaffold > 0 else 0.0
    sp3_fraction = (n_sp3 / n_scaffold) if n_scaffold > 0 else 0.0

    # Connected sp2 domains on the explicit PBC-aware scaffold graph
    remaining_sp2 = set(sp2_scaffold_ids)
    domain_sizes: list[int] = []
    while remaining_sp2:
        start = remaining_sp2.pop()
        stack = [start]
        size = 0
        while stack:
            curr = stack.pop()
            size += 1
            for nbr, _, _ in adj[curr]:
                if nbr in remaining_sp2:
                    remaining_sp2.remove(nbr)
                    stack.append(nbr)
        domain_sizes.append(size)
    domain_sizes.sort(reverse=True)

    sp2_domains = {
        "count": len(domain_sizes),
        "largest_domain_carbons": domain_sizes[0] if domain_sizes else 0,
        "mean_domain_carbons": (sum(domain_sizes) / len(domain_sizes)) if domain_sizes else 0.0,
        "size_distribution": domain_sizes,
    }

    env_fractions = {
        k: (v / n_c if n_c > 0 else 0.0) for k, v in env_counts.items()
    }

    elec = electronic_summary(model)
    has_sublattice = bool(scaffold_ids) and all(
        by_id[aid].get("metadata", {}).get("sublattice") in (0, 1) for aid in scaffold_ids
    )
    if has_sublattice:
        n_a = sum(
            1 for aid in sp2_scaffold_ids if by_id[aid]["metadata"]["sublattice"] == 0
        )
        n_b = sum(
            1 for aid in sp2_scaffold_ids if by_id[aid]["metadata"]["sublattice"] == 1
        )
        imbalance = abs(n_a - n_b)
        sublattice_counts = {"A": n_a, "B": n_b}
    else:
        sublattice_counts = None
        imbalance = None

    return {
        "formula": formula,
        "atom_counts": atom_counts,
        "atomic_fractions": atomic_fractions,
        "c_over_o": c_over_o,
        "o_over_c": o_over_c,
        "formal_charge": formal_charge,
        "dimensions_nm": {
            "cell_nm": cell_nm,
            "extent_nm": extent_nm,
            "scaffold_diameter_nm": scaffold_diameter_nm,
            "periodic": [bool(x) for x in model.get("pbc", [False, False, False])],
        },
        "scaffold_carbons": n_scaffold,
        "sp2_carbons": n_sp2,
        "sp3_carbons": n_sp3,
        "sp2_fraction": sp2_fraction,
        "sp3_fraction": sp3_fraction,
        "sp2_domains": sp2_domains,
        "chemical_environments": {
            "basis": "structural_graph_carbon_environments_not_xps_intensities",
            "total_carbons": n_c,
            "counts": env_counts,
            "fractions": env_fractions,
            "carboxyl_breakdown": {
                "carboxyl_C4": carboxyl_c4,
                "anchor_CY": carboxyl_cy,
            },
        },
        "qm_diagnostics": {
            "electron_count": elec["electron_count"],
            "electron_parity": elec["electron_parity"],
            "sp2_sublattice_counts": sublattice_counts,
            "sp2_sublattice_imbalance": imbalance,
            "caveat": (
                "sp2_sublattice_imbalance (|N_A - N_B|) is purely a topological scaffold count "
                "(Lieb theorem applies only to neutral bipartite pi-systems and does not generally predict "
                "spin multiplicity for oxidized, charged, carbonyl, or periodic systems); "
                "it does not replace electronic-structure QM."
            ),
        },
    }


def assert_export_ready(
    model: dict,
    *,
    allow_undercoordinated: bool = False,
    check_bond_lengths: bool = True,
) -> dict:
    """Assert geometric and valence integrity prior to exporting coordinates."""
    report = check_geometry_integrity(model, check_bond_lengths=check_bond_lengths)
    reasons = []
    if report["clash_count"] > 0:
        first = report["clashes"][0]
        thresh = first.get("threshold", first.get("cutoff"))
        kind = first.get("kind", "nonbonded")
        reasons.append(
            f"{report['clash_count']} nonbonded clash(es) "
            f"(e.g. {kind} atoms {first['a']}-{first['b']} at {first['distance']:.3f} A < {thresh:.2f} A)"
        )
    if not allow_undercoordinated and report["undercoordinated_atom_ids"]:
        sample = report["undercoordinated_atom_ids"][:10]
        reasons.append(f"undercoordinated atom IDs: {sample}")
    if not report["bond_lengths_valid"]:
        violations = report.get("bond_length_violations", [])
        if violations:
            v0 = violations[0]
            reasons.append(
                f"implausible bond length(s) (e.g. {v0['elements'][0]}-{v0['elements'][1]} "
                f"atoms {v0['a']}-{v0['b']} at {v0['distance']:.3f} A outside {v0['expected_range']})"
            )
        else:
            reasons.append(
                f"implausible bond lengths outside [0.85, 1.85] A "
                f"(range [{report['min_bond_length']:.3f}, {report['max_bond_length']:.3f}] A)"
            )
    if reasons:
        raise ValueError(f"Model failed pre-export integrity gate: {'; '.join(reasons)}")
    return report


def load_model(path) -> dict:
    with Path(path).open() as handle:
        model = json.load(handle)
    validate_model(model)
    return model


def save_model(model: dict, path) -> None:
    validate_model(model)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        json.dump(model, handle, indent=2, allow_nan=False)
        handle.write("\n")


def write_pdb(model: dict, path, *, allow_undercoordinated: bool = False) -> None:
    """Write explicit finite CONECT records; never flatten periodic seam bonds."""
    assert_export_ready(model, allow_undercoordinated=allow_undercoordinated)
    if any(any(b.get("image", [0, 0, 0])) for b in model["bonds"]):
        raise ValueError("PDB cannot preserve periodic image bonds; use JSON")
    if any(a["id"] > 99999 or a.get("resid", 1) > 9999 for a in model["atoms"]):
        raise ValueError("PDB identifier range exceeded; use JSON")
    lines = ["REMARK Unrelaxed constructed graph; authoritative connectivity is in JSON.\n"]
    for atom in model["atoms"]:
        x, y, z = atom["position"]
        if any(not (-999.999 <= float(v) <= 9999.999) for v in (x, y, z)):
            raise ValueError("PDB coordinate range exceeded")
        q = int(atom.get("formal_charge", 0))
        if q == 0:
            q_str = "  "
        elif -9 <= q <= 9:
            q_str = f"{abs(q)}{'+' if q > 0 else '-'}"
        else:
            raise ValueError(f"PDB formal charge range exceeded at atom {atom['id']}: {q}")
        raw_name = str(atom.get("name") or atom["element"]).strip() or str(atom["element"])
        elem = str(atom["element"]).strip()
        if len(elem) == 1 and len(raw_name) < 4:
            pdb_name = f" {raw_name:<3.3s}"
        else:
            pdb_name = f"{raw_name:<4.4s}"
        lines.append(
            f"HETATM{atom['id']:5d} {pdb_name} "
            f"{atom.get('resname', 'GO'):>3.3s} A{atom.get('resid', 1):4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom['element']:>2s}{q_str}\n"
        )
    conect_adj: dict[int, list[int]] = {}
    for bond in model["bonds"]:
        a_id, b_id = int(bond["a"]), int(bond["b"])
        conect_adj.setdefault(a_id, []).append(b_id)
        conect_adj.setdefault(b_id, []).append(a_id)
    for src_id in sorted(conect_adj):
        nbrs = sorted(conect_adj[src_id])
        for i in range(0, len(nbrs), 4):
            chunk = nbrs[i : i + 4]
            lines.append("CONECT" + f"{src_id:5d}" + "".join(f"{n:5d}" for n in chunk) + "\n")
    lines.append("END\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)


def write_xyz(model: dict, path) -> None:
    """Write standard XYZ coordinates (Angstrom); reject periodic image bonds."""
    report = assert_export_ready(model)
    if any(any(b.get("image", [0, 0, 0])) for b in model["bonds"]):
        raise ValueError("XYZ cannot preserve periodic image bonds; use write_extxyz or JSON")
    metadata = model.get("metadata", {})
    summary = electronic_summary(model)
    name = model.get("name", "GO")
    formal_charge = report["formal_charge"]
    electron_count = summary.get("electron_count")
    geometry_status = metadata.get("geometry_status", "unknown")
    comment = (
        f"name={name} formal_charge={formal_charge} "
        f"electron_count={electron_count} geometry_status={geometry_status}"
    )
    lines = [f"{len(model['atoms'])}\n", f"{comment}\n"]
    for atom in model["atoms"]:
        x, y, z = atom["position"]
        lines.append(f"{atom['element']:<2s} {x:14.6f} {y:14.6f} {z:14.6f}\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)


def _escape_extxyz_str(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def write_extxyz(model: dict, path) -> None:
    """Write Extended XYZ (ASE / xTB / CP2K / ORCA compatible) for periodic or finite sheets."""
    report = assert_export_ready(model)
    metadata = model.get("metadata", {})
    summary = electronic_summary(model)
    raw_name = "GO" if model.get("name") is None else str(model.get("name"))
    raw_geom = str(metadata.get("geometry_status", "unknown"))
    if "\n" in raw_name or "\r" in raw_name or "\n" in raw_geom or "\r" in raw_geom:
        raise ValueError("Extended XYZ single-line fields cannot contain newlines")
    name = _escape_extxyz_str(raw_name)
    formal_charge = report["formal_charge"]
    electron_count = summary.get("electron_count")
    geometry_status = _escape_extxyz_str(raw_geom)
    pbc_str = " ".join("T" if flag else "F" for flag in model["pbc"])
    header_parts = []
    if model.get("cell") is not None:
        cell_arr = np.asarray(model["cell"], dtype=float)
        lattice_vals = " ".join(f"{float(v):.8f}" for v in cell_arr.reshape(-1))
        header_parts.append(f'Lattice="{lattice_vals}"')
    header_parts.extend(
        [
            "Properties=species:S:1:pos:R:3",
            f'pbc="{pbc_str}"',
            f'name="{name}"',
            f"formal_charge={formal_charge}",
            f"electron_count={electron_count}",
            f'geometry_status="{geometry_status}"',
        ]
    )
    lines = [f"{len(model['atoms'])}\n", " ".join(header_parts) + "\n"]
    for atom in model["atoms"]:
        x, y, z = atom["position"]
        lines.append(f"{atom['element']:<2s} {x:14.6f} {y:14.6f} {z:14.6f}\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)


def _hopcroft_karp_matching_dict(
    left_nodes: list[int], adj_map: dict[int, list[int]]
) -> dict[int, int]:
    """Iterative O(E sqrt(V)) Hopcroft-Karp bipartite matching returning symmetric mate dict."""
    if not left_nodes:
        return {}
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

    mate: dict[int, int] = {}
    for u, v in pair_u.items():
        if v is not None:
            mate[u] = v
            mate[v] = u
    return mate


def _blossom_max_matching(nodes: list[int], adj_map: dict[int, list[int]]) -> dict[int, int]:
    """Iterative Edmonds Blossom maximum cardinality matching on general graphs."""
    mate: dict[int, int] = {}
    for root in nodes:
        if root in mate:
            continue
        parent: dict[int, int] = {}
        base: dict[int, int] = {v: v for v in nodes}
        used: set[int] = {root}
        blossom_mark: set[int] = set()
        queue: deque[int] = deque([root])

        def lca(a: int, b: int) -> int:
            seen: set[int] = set()
            while True:
                a = base[a]
                seen.add(a)
                if a not in mate:
                    break
                a = parent[mate[a]]
            while True:
                b = base[b]
                if b in seen:
                    return b
                b = parent[mate[b]]

        def mark_path(v: int, b: int, children: int) -> None:
            while base[v] != b:
                blossom_mark.add(base[v])
                blossom_mark.add(base[mate[v]])
                parent[v] = children
                children = mate[v]
                v = parent[mate[v]]

        found_end: int | None = None
        while queue and found_end is None:
            v = queue.popleft()
            for to in adj_map.get(v, ()):
                if base[v] == base[to] or mate.get(v) == to:
                    continue
                if to in used:
                    curbase = lca(v, to)
                    blossom_mark.clear()
                    mark_path(v, curbase, to)
                    mark_path(to, curbase, v)
                    for i in nodes:
                        if base[i] in blossom_mark:
                            base[i] = curbase
                            if i not in used:
                                used.add(i)
                                queue.append(i)
                elif to not in parent:
                    parent[to] = v
                    if to not in mate:
                        found_end = to
                        break
                    nxt = mate[to]
                    used.add(nxt)
                    queue.append(nxt)
        if found_end is not None:
            v = found_end
            while v is not None:
                pv = parent[v]
                ppv = mate.get(pv)
                mate[v] = pv
                mate[pv] = v
                v = ppv
    return mate


def _kekule_matching_for_sdf(
    atoms: list[dict], bonds: list[dict]
) -> tuple[dict[int, int], set[int]]:
    """Kekulize sp2 aromatic-order (1.5) carbon subgraph via maximum matching."""
    by_id = {a["id"]: a for a in atoms}
    inc_orders: dict[int, list[float]] = {a["id"]: [] for a in atoms}
    pi_nodes: set[int] = set()
    pi_edges: list[tuple[int, int]] = []

    for bond in bonds:
        a_id = int(bond["a"])
        b_id = int(bond["b"])
        order = float(bond["order"])
        inc_orders[a_id].append(order)
        inc_orders[b_id].append(order)
        if order == 1.5 and by_id[a_id]["element"] == "C" and by_id[b_id]["element"] == "C":
            pi_edges.append((a_id, b_id))

    has_higher_order = {
        aid for aid, ords in inc_orders.items() if any(o >= 2.0 for o in ords)
    }
    pi_adj: dict[int, list[int]] = {}
    for a_id, b_id in pi_edges:
        if a_id in has_higher_order or b_id in has_higher_order:
            continue
        pi_nodes.add(a_id)
        pi_nodes.add(b_id)
        pi_adj.setdefault(a_id, []).append(b_id)
        pi_adj.setdefault(b_id, []).append(a_id)

    for uid in pi_adj:
        pi_adj[uid].sort()

    ordered_nodes = sorted(pi_nodes)
    color: dict[int, int] = {}
    is_bipartite = True
    for start in ordered_nodes:
        if start in color:
            continue
        color[start] = 0
        q: deque[int] = deque([start])
        while q and is_bipartite:
            u = q.popleft()
            for v in pi_adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    q.append(v)
                elif color[v] == color[u]:
                    is_bipartite = False
                    break

    if not ordered_nodes:
        mate: dict[int, int] = {}
    elif is_bipartite:
        left_nodes = [u for u in ordered_nodes if color[u] == 0]
        mate = _hopcroft_karp_matching_dict(left_nodes, pi_adj)
    else:
        mate = _blossom_max_matching(ordered_nodes, pi_adj)

    unpaired: set[int] = {u for u in pi_nodes if u not in mate}
    for atom in atoms:
        aid = atom["id"]
        if (
            atom["element"] == "C"
            and aid not in pi_nodes
            and int(atom.get("formal_charge", 0)) == 0
        ):
            ords = inc_orders[aid]
            if len(ords) == 3 and abs(sum(ords) - 3.0) < 1e-6:
                unpaired.add(aid)

    return mate, unpaired


def write_sdf(model: dict, path) -> None:
    """Write chemically valid Kekule Molfile V3000 (.sdf) with formal charges and RAD=2."""
    report = assert_export_ready(model)
    if any(any(b.get("image", [0, 0, 0])) for b in model["bonds"]):
        raise ValueError("SDF cannot preserve periodic image bonds; use JSON")
    atoms = model["atoms"]
    bonds = model["bonds"]
    metadata = model.get("metadata", {})
    summary = electronic_summary(model)

    raw_name = "" if model.get("name") is None else str(model.get("name"))
    raw_geom = str(metadata.get("geometry_status", "unknown"))
    if "\n" in raw_name or "\r" in raw_name or "\n" in raw_geom or "\r" in raw_geom:
        raise ValueError("SDF single-line header fields cannot contain newlines")
    name = raw_name.strip() or "GO"
    geometry_status = raw_geom
    formal_charge = report["formal_charge"]
    electron_count = summary.get("electron_count")
    id_to_seq = {atom["id"]: idx for idx, atom in enumerate(atoms, 1)}

    mate, unpaired_pi = _kekule_matching_for_sdf(atoms, bonds)
    bond_order_map = {1.0: 1, 2.0: 2, 3.0: 3}

    lines = [
        f"{name}\n",
        "  GOBUILD 00000000003D\n",
        f"geometry_status={geometry_status}\n",
        "  0  0  0     0  0            999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        f"M  V30 COUNTS {len(atoms)} {len(bonds)} 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
    ]
    for idx, atom in enumerate(atoms, 1):
        x, y, z = atom["position"]
        q = int(atom.get("formal_charge", 0))
        chg_part = f" CHG={q}" if q != 0 else ""
        rad_part = " RAD=2" if atom["id"] in unpaired_pi else ""
        lines.append(
            f"M  V30 {idx} {atom['element']} {x:.6f} {y:.6f} {z:.6f} 0{chg_part}{rad_part}\n"
        )
    lines.append("M  V30 END ATOM\n")
    lines.append("M  V30 BEGIN BOND\n")
    for bidx, bond in enumerate(bonds, 1):
        a_id = int(bond["a"])
        b_id = int(bond["b"])
        order = float(bond["order"])
        if order == 1.5:
            order_int = 2 if mate.get(a_id) == b_id else 1
        else:
            order_int = bond_order_map[order]
        lines.append(f"M  V30 {bidx} {order_int} {id_to_seq[a_id]} {id_to_seq[b_id]}\n")
    lines.extend(
        [
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            ">  <formal_charge>\n",
            f"{formal_charge if formal_charge is not None else 'unknown'}\n",
            "\n",
            ">  <electron_count>\n",
            f"{electron_count if electron_count is not None else 'unknown'}\n",
            "\n",
            ">  <geometry_status>\n",
            f"{geometry_status}\n",
            "\n",
            ">  <unpaired_pi_sites>\n",
            f"{len(unpaired_pi)}\n",
            "\n",
            "$$$$\n",
        ]
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)


def _orthogonal_complement_basis(periodic_vecs: list[np.ndarray]) -> list[np.ndarray]:
    """Return orthonormal unit vectors spanning the orthogonal complement of periodic_vecs."""
    basis: list[np.ndarray] = []
    for vec in periodic_vecs:
        v = np.array(vec, dtype=float)
        for b in basis:
            v = v - float(np.dot(v, b)) * b
        nrm = float(np.linalg.norm(v))
        if nrm > 1e-10:
            basis.append(v / nrm)

    complement: list[np.ndarray] = []
    for trial in (
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    ):
        v = trial.copy()
        for b in basis + complement:
            v = v - float(np.dot(v, b)) * b
        nrm = float(np.linalg.norm(v))
        if nrm > 1e-8:
            complement.append(v / nrm)
            if len(basis) + len(complement) == 3:
                break
    return complement


def write_cif(model: dict, path) -> None:
    """Write standard P1 CIF for periodic or finite models with >= 20 A non-periodic clearance."""
    assert_export_ready(model)
    atoms = model["atoms"]
    pbc = [bool(x) for x in model["pbc"]]
    coords = np.asarray([a["position"] for a in atoms], dtype=float)
    orig_cell = (
        np.asarray(model["cell"], dtype=float) if model.get("cell") is not None else None
    )

    cell = np.zeros((3, 3), dtype=float)
    periodic_axes = [k for k in range(3) if pbc[k]]
    non_periodic_axes = [k for k in range(3) if not pbc[k]]

    if orig_cell is not None:
        for k in periodic_axes:
            cell[k] = orig_cell[k]

    if non_periodic_axes:
        if not periodic_axes:
            ortho_dirs = [
                np.array([1.0, 0.0, 0.0]),
                np.array([0.0, 1.0, 0.0]),
                np.array([0.0, 0.0, 1.0]),
            ]
        else:
            comp_vecs = _orthogonal_complement_basis([cell[k] for k in periodic_axes])
            # Match complement directions to non-periodic axes maximizing alignment
            if len(non_periodic_axes) == 1:
                ortho_dirs = [comp_vecs[0]]
            else:
                e0 = np.eye(3)[non_periodic_axes[0]]
                if abs(float(np.dot(comp_vecs[0], e0))) >= abs(float(np.dot(comp_vecs[1], e0))):
                    ortho_dirs = [comp_vecs[0], comp_vecs[1]]
                else:
                    ortho_dirs = [comp_vecs[1], comp_vecs[0]]

        for ax, u_dir in zip(non_periodic_axes, ortho_dirs, strict=True):
            # Orient positively along original cell vector or canonical standard axis if available
            ref_vec = orig_cell[ax] if orig_cell is not None else np.eye(3)[ax]
            if float(np.dot(u_dir, ref_vec)) < 0.0:
                u_dir = -u_dir
            proj = coords @ u_dir
            span = float(np.max(proj) - np.min(proj))
            length = max(span + 20.0, 20.0)
            cell[ax] = length * u_dir

        if float(np.linalg.det(cell)) < 0.0:
            cell[non_periodic_axes[-1]] = -cell[non_periodic_axes[-1]]

    va, vb, vc = cell[0], cell[1], cell[2]
    a_len = float(np.linalg.norm(va))
    b_len = float(np.linalg.norm(vb))
    c_len = float(np.linalg.norm(vc))
    alpha = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(vb, vc) / (b_len * c_len))))))
    beta = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(va, vc) / (a_len * c_len))))))
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(va, vb) / (a_len * b_len))))))

    inv_cell = np.linalg.inv(cell)
    frac = coords @ inv_cell
    for k in range(3):
        if pbc[k]:
            frac[:, k] = np.mod(frac[:, k], 1.0)
        else:
            fmin = float(np.min(frac[:, k]))
            fmax = float(np.max(frac[:, k]))
            frac[:, k] += 0.5 - 0.5 * (fmin + fmax)

    # Compute Cartesian coordinates in the standard crystallographic frame defined by
    # (_cell_length_*, _cell_angle_*) so _atom_site_Cartn_* and _atom_site_fract_* always agree
    # even for rigidly rotated ribbons or arbitrary cell orientations.
    al_rad, be_rad, ga_rad = math.radians(alpha), math.radians(beta), math.radians(gamma)
    cos_al, cos_be, cos_ga = math.cos(al_rad), math.cos(be_rad), math.cos(ga_rad)
    sin_ga = math.sin(ga_rad)
    va_std = np.array([a_len, 0.0, 0.0], dtype=float)
    vb_std = np.array([b_len * cos_ga, b_len * sin_ga, 0.0], dtype=float)
    cx_std = c_len * cos_be
    cy_std = c_len * (cos_al - cos_be * cos_ga) / sin_ga
    cz_std = math.sqrt(max(0.0, c_len * c_len - cx_std * cx_std - cy_std * cy_std))
    vc_std = np.array([cx_std, cy_std, cz_std], dtype=float)
    std_cell = np.vstack([va_std, vb_std, vc_std])
    cart = frac @ std_cell

    safe_name = re.sub(r"[^A-Za-z0-9_]", "_", str(model.get("name", "GO"))).strip("_") or "GO"
    lines = [
        f"data_{safe_name}\n",
        f"_cell_length_a    {a_len:.8f}\n",
        f"_cell_length_b    {b_len:.8f}\n",
        f"_cell_length_c    {c_len:.8f}\n",
        f"_cell_angle_alpha {alpha:.6f}\n",
        f"_cell_angle_beta  {beta:.6f}\n",
        f"_cell_angle_gamma {gamma:.6f}\n",
        "_symmetry_space_group_name_H-M 'P 1'\n",
        "_symmetry_Int_Tables_number 1\n",
        "\n",
        "loop_\n",
        "_atom_site_label\n",
        "_atom_site_type_symbol\n",
        "_atom_site_fract_x\n",
        "_atom_site_fract_y\n",
        "_atom_site_fract_z\n",
        "_atom_site_Cartn_x\n",
        "_atom_site_Cartn_y\n",
        "_atom_site_Cartn_z\n",
    ]
    for idx, atom in enumerate(atoms):
        label = f"{atom['element']}{atom['id']}"
        fx, fy, fz = frac[idx]
        cx, cy, cz = cart[idx]
        lines.append(
            f"{label:<8s} {atom['element']:<3s} {fx:14.9f} {fy:14.9f} {fz:14.9f} "
            f"{cx:14.8f} {cy:14.8f} {cz:14.8f}\n"
        )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)


def _tripos_atom_type(atom: dict, neighbors: dict, by_id: dict) -> str:
    element = atom["element"]
    ns = neighbors[atom["id"]]
    degree = len(ns)
    orders = [n[2] for n in ns]
    role = atom.get("metadata", {}).get("role")
    hybrid = atom.get("metadata", {}).get("hybridization")
    if element == "H":
        return "H"
    if element == "C":
        if hybrid == "sp3" or (degree == 4 and all(o == 1.0 for o in orders)):
            return "C.3"
        if any(o == 1.5 for o in orders):
            return "C.ar"
        if any(o == 2.0 for o in orders):
            return "C.2"
        if any(o == 3.0 for o in orders):
            return "C.1"
        if role == "sheet_carbon":
            return "C.ar"
        return "C.3" if degree == 4 else "C.2"
    if element == "O":
        if degree == 1:
            parent_id, _, order = ns[0]
            parent = by_id[parent_id]
            if parent["element"] == "C":
                parent_oxygens = [
                    by_id[nbr_id]
                    for nbr_id, _, _ in neighbors[parent_id]
                    if by_id[nbr_id]["element"] == "O" and len(neighbors[nbr_id]) == 1
                ]
                if len(parent_oxygens) >= 2 and any(
                    o.get("formal_charge", 0) == -1
                    or o.get("metadata", {}).get("group") == "COO"
                    or o.get("resname") == "CBD"
                    for o in parent_oxygens
                ):
                    return "O.co2"
            if order == 2.0:
                return "O.2"
        return "O.3"
    if element == "N":
        if atom.get("formal_charge", 0) == 1 or degree == 4:
            return "N.4"
        if any(o == 1.5 for o in orders):
            return "N.ar"
        if any(o == 2.0 for o in orders):
            return "N.2"
        return "N.3"
    if element == "S":
        return "S.2" if any(o == 2.0 for o in orders) else "S.3"
    return element


def write_mol2(model: dict, path) -> None:
    """Write Tripos MOL2 with contiguous residues and explicit bond orders; reject periodic image bonds."""
    assert_export_ready(model)
    if any(any(b.get("image", [0, 0, 0])) for b in model["bonds"]):
        raise ValueError("MOL2 cannot preserve periodic image bonds; use JSON")

    raw_atoms = model["atoms"]
    bonds = model["bonds"]
    neighbors = adjacency(model)
    by_id = {atom["id"]: atom for atom in raw_atoms}

    residue_order: list[tuple[int, str]] = []
    residue_groups: dict[tuple[int, str], list[dict]] = {}
    for atom in raw_atoms:
        key = (int(atom.get("resid", 1)), str(atom.get("resname", "GO")))
        if key not in residue_groups:
            residue_order.append(key)
            residue_groups[key] = [atom]
        else:
            residue_groups[key].append(atom)
    atoms = [atom for key in residue_order for atom in residue_groups[key]]
    id_to_seq = {atom["id"]: idx for idx, atom in enumerate(atoms, 1)}

    has_all_charges = all(atom.get("charge") is not None for atom in atoms)
    charge_type = "USER_CHARGES" if has_all_charges else "NO_CHARGES"
    mol_comment = (
        ""
        if has_all_charges
        else "Partial charges not assigned; missing charges written as 0.0000."
    )

    substructures: list[tuple[int, str, int]] = []
    substruct_map: dict[tuple[int, str], int] = {}
    for seq_idx, atom in enumerate(atoms, 1):
        key = (int(atom.get("resid", 1)), str(atom.get("resname", "GO")))
        if key not in substruct_map:
            substruct_map[key] = len(substructures) + 1
            substructures.append((key[0], key[1], seq_idx))

    bond_order_map = {1.0: "1", 1.5: "ar", 2.0: "2", 3.0: "3"}
    lines = [
        "@<TRIPOS>MOLECULE\n",
        f"{model.get('name', 'GO')}\n",
        f"{len(atoms)} {len(bonds)} {len(substructures)} 0 0\n",
        "SMALL\n",
        f"{charge_type}\n",
        "\n",
        f"{mol_comment}\n",
        "\n",
        "@<TRIPOS>ATOM\n",
    ]
    for seq_idx, atom in enumerate(atoms, 1):
        x, y, z = atom["position"]
        atom_name = str(atom.get("name") or f"{atom['element']}{atom['id']}")
        resname = str(atom.get("resname", "GO"))
        resid = int(atom.get("resid", 1))
        subst_id = substruct_map[(resid, resname)]
        tripos_type = _tripos_atom_type(atom, neighbors, by_id)
        charge = float(atom["charge"]) if atom.get("charge") is not None else 0.0
        lines.append(
            f"{seq_idx:7d} {atom_name:<8.8s} {x:10.4f} {y:10.4f} {z:10.4f} "
            f"{tripos_type:<6s} {subst_id:5d} {resname:<8.8s} {charge:9.4f}\n"
        )
    lines.append("@<TRIPOS>BOND\n")
    for bond_idx, bond in enumerate(bonds, 1):
        order_str = bond_order_map[float(bond["order"])]
        lines.append(
            f"{bond_idx:6d} {id_to_seq[bond['a']]:6d} {id_to_seq[bond['b']]:6d} {order_str}\n"
        )
    lines.append("@<TRIPOS>SUBSTRUCTURE\n")
    for subst_idx, (_, resname, root_atom) in enumerate(substructures, 1):
        lines.append(f"{subst_idx:6d} {resname:<8.8s} {root_atom:6d} RESIDUE 4 A {resname:<8.8s} 0 ROOT\n")

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)
