"""Direct periodic Amber/GROMACS/OpenMM topology writers for GO Workbench (Block 5).

Builds covalently periodic (and finite) ``gaff2_v3`` ParmEd / AmberParm / GROMACS /
OpenMM topologies directly from the authoritative JSON molecular graph and
provenance-verified ``GO_gaff2_v3.frcmod`` parameters without invoking ``tleap``.

Every periodic bond, angle, proper 1-4 dihedral, improper torsion, and nonbonded
exclusion is validated against Minimum Image Convention (MIC) unambiguity prior
to writing output files.
"""

from __future__ import annotations

import contextlib
import hashlib
import math
import os
import re
import tempfile
from collections import defaultdict, deque
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

from .ff_v3 import (
    assign_gaff2_v3,
    verify_v3_bundle,
)

ATOMIC_NUMBERS = {"C": 6, "H": 1, "O": 8}
ATOMIC_MASSES = {"C": 12.01, "H": 1.008, "O": 16.00}

_BOND_RE = re.compile(r"^([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s+(.+)$")
_ANGLE_RE = re.compile(
    r"^([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s+(.+)$"
)
_DIHE_RE = re.compile(
    r"^([A-Za-z0-9*X ]{1,2})-([A-Za-z0-9*X ]{1,2})-([A-Za-z0-9*X ]{1,2})-([A-Za-z0-9*X ]{1,2})\s+(.+)$"
)


def _cell_lengths_angles(cell: np.ndarray) -> tuple[list[float], list[float], list[float]]:
    """Return box lengths (A), angles (deg), and perpendicular face heights (A)."""
    a_vec, b_vec, c_vec = cell[0], cell[1], cell[2]
    la = float(np.linalg.norm(a_vec))
    lb = float(np.linalg.norm(b_vec))
    lc = float(np.linalg.norm(c_vec))

    def _ang(u: np.ndarray, lu: float, v: np.ndarray, lv: float) -> float:
        cos_val = float(np.clip(np.dot(u, v) / (lu * lv), -1.0, 1.0))
        return math.degrees(math.acos(cos_val))

    alpha = _ang(b_vec, lb, c_vec, lc)
    beta = _ang(a_vec, la, c_vec, lc)
    gamma = _ang(a_vec, la, b_vec, lb)

    vol = abs(float(np.linalg.det(cell)))
    heights = [
        vol / float(np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])))
        for i in range(3)
    ]
    return [la, lb, lc], [alpha, beta, gamma], heights


def _mic_vector(delta: np.ndarray, pbc: list[bool], cell: np.ndarray, inv_cell: np.ndarray) -> np.ndarray:
    """Apply Minimum Image Convention along periodic lattice axes."""
    if not any(pbc):
        return np.array(delta, dtype=float)
    frac = delta @ inv_cell
    for ax in range(3):
        if pbc[ax]:
            frac[..., ax] -= np.rint(frac[..., ax])
    return frac @ cell


def _parse_int_image(raw_image: Any, pbc: list[bool], context: str) -> tuple[int, int, int]:
    """Validate and normalize a 3-component lattice image shift."""
    if raw_image is None:
        raw_image = [0, 0, 0]
    if not isinstance(raw_image, (list, tuple)) or len(raw_image) != 3:
        raise ValueError(f"Invalid bond image {raw_image!r} on {context}; expected 3 integers")
    vals: list[int] = []
    for ax, v in enumerate(raw_image):
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)):
            raise ValueError(f"Bond image shift must contain integers on {context}, got {raw_image!r}")
        iv = int(v)
        if not pbc[ax] and iv != 0:
            raise ValueError(
                f"Nonzero bond image {raw_image!r} along non-periodic axis {ax} on {context}"
            )
        vals.append(iv)
    return (vals[0], vals[1], vals[2])


def load_v3_frcmod_parameters(
    reference_dir: str | Path | None = None,
    verified_bundle: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Load and verify ``GO_gaff2_v3.frcmod`` tables using bundle SHA-256 provenance."""
    bundle = (
        verified_bundle if verified_bundle is not None else verify_v3_bundle(reference_dir)
    )
    text = bundle["frcmod_bytes"].decode("utf-8")

    masses: dict[str, tuple[float, float]] = {}
    bonds: dict[tuple[str, str], tuple[float, float]] = {}
    angles: dict[tuple[str, str, str], tuple[float, float]] = {}
    dihe_specific: dict[tuple[str, str, str, str], list[tuple[float, float, float, int]]] = defaultdict(list)
    dihe_wildcard: dict[tuple[str, str], list[tuple[float, float, float, int]]] = defaultdict(list)
    impropers: list[tuple[str, str, str, str, float, float, int]] = []
    nonbon: dict[str, tuple[float, float]] = {}

    section: str | None = None
    for idx, raw_line in enumerate(text.splitlines(), 1):
        stripped = raw_line.strip()
        if idx == 1 and stripped not in {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}:
            continue
        if not stripped:
            continue
        if stripped in {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}:
            section = stripped
            continue
        if section is None or stripped.startswith((";", "#", "!")):
            continue

        if section == "MASS":
            parts = stripped.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed MASS entry in GO_gaff2_v3.frcmod: {raw_line!r}")
            masses[parts[0]] = (float(parts[1]), float(parts[2]) if len(parts) > 2 else 0.0)
        elif section == "BOND":
            m = _BOND_RE.match(stripped)
            if not m:
                raise ValueError(f"Malformed BOND entry in GO_gaff2_v3.frcmod: {raw_line!r}")
            t1, t2 = sorted((m.group(1).strip(), m.group(2).strip()))
            vals = m.group(3).split()
            bonds[(t1, t2)] = (float(vals[0]), float(vals[1]))
        elif section == "ANGLE":
            m = _ANGLE_RE.match(stripped)
            if not m:
                raise ValueError(f"Malformed ANGLE entry in GO_gaff2_v3.frcmod: {raw_line!r}")
            a, b, c = m.group(1).strip(), m.group(2).strip(), m.group(3).strip()
            vals = m.group(4).split()
            angles[(min(a, c), b, max(a, c))] = (float(vals[0]), float(vals[1]))
        elif section == "DIHE":
            m = _DIHE_RE.match(raw_line.rstrip())
            if not m:
                raise ValueError(f"Malformed DIHE entry in GO_gaff2_v3.frcmod: {raw_line!r}")
            t1, t2, t3, t4 = (m.group(k).strip() for k in range(1, 5))
            vals = m.group(5).split()
            idivf, pk, phase, pn_raw = float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])
            per = abs(int(round(pn_raw)))
            if per <= 0:
                raise ValueError(f"Unsupported dihedral periodicity {pn_raw!r} in GO_gaff2_v3.frcmod")
            term = (idivf, pk, phase, per)
            if t1 == "X" and t4 == "X":
                dihe_wildcard[(min(t2, t3), max(t2, t3))].append(term)
            else:
                fwd = (t1, t2, t3, t4)
                rev = (t4, t3, t2, t1)
                dihe_specific[min(fwd, rev)].append(term)
        elif section == "IMPROPER":
            m = _DIHE_RE.match(raw_line.rstrip())
            if not m:
                raise ValueError(f"Malformed IMPROPER entry in GO_gaff2_v3.frcmod: {raw_line!r}")
            t1, t2, t3, t4 = (m.group(k).strip() for k in range(1, 5))
            vals = m.group(5).split()
            pk, phase, pn_raw = float(vals[0]), float(vals[1]), float(vals[2])
            impropers.append((t1, t2, t3, t4, pk, phase, abs(int(round(pn_raw)))))
        elif section == "NONBON":
            if stripped in {"END", "RE"}:
                continue
            parts = stripped.split()
            if len(parts) < 3:
                raise ValueError(f"Malformed NONBON entry in GO_gaff2_v3.frcmod: {raw_line!r}")
            nonbon[parts[0]] = (float(parts[1]), float(parts[2]))

    return {
        "bundle": bundle,
        "masses": masses,
        "bonds": bonds,
        "angles": angles,
        "dihe_specific": dict(dihe_specific),
        "dihe_wildcard": dict(dihe_wildcard),
        "impropers": impropers,
        "nonbon": nonbon,
    }


def validate_periodic_model(model: dict[str, Any], *, require_periodic: bool = True) -> dict[str, Any]:
    """Validate periodic lattice cell, covalent graph, and 1-2/1-3/1-4 MIC unambiguity."""
    from .topology import _validate_graph

    _validate_graph(model)
    pbc_raw = model.get("pbc")
    if not isinstance(pbc_raw, (list, tuple)) or len(pbc_raw) != 3 or not all(
        isinstance(x, bool) for x in pbc_raw
    ):
        raise ValueError(f"Model 'pbc' must be a list of 3 booleans, got {pbc_raw!r}")
    pbc = [bool(x) for x in pbc_raw]
    is_periodic = any(pbc)
    if require_periodic and not is_periodic:
        raise ValueError(
            "Model is non-periodic (pbc=[False, False, False]); use standard finite export or supply a periodic model"
        )

    raw_cell = model.get("cell")
    if is_periodic and raw_cell is None:
        raise ValueError("Periodic model requires a valid 3x3 'cell' matrix")

    cell: np.ndarray | None = None
    inv_cell: np.ndarray | None = None
    lengths: list[float] | None = None
    angles: list[float] | None = None
    heights: list[float] | None = None
    min_periodic_height = math.inf

    if raw_cell is not None:
        cell = np.asarray(raw_cell, dtype=float)
        if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
            raise ValueError("Model 'cell' must be a finite 3x3 matrix")
        if float(np.linalg.det(cell)) <= 1e-6:
            raise ValueError("Model 'cell' must be right-handed with positive volume")
        if (
            cell[0, 0] <= 0.0
            or cell[1, 1] <= 0.0
            or cell[2, 2] <= 0.0
            or abs(cell[0, 1]) > 1e-8
            or abs(cell[0, 2]) > 1e-8
            or abs(cell[1, 2]) > 1e-8
        ):
            raise ValueError(
                "Periodic 'cell' must be in standard lower-triangular right-handed Cartesian form "
                "(a along +x, b in xy plane with by>0, cz>0)"
            )
        inv_cell = np.linalg.inv(cell)
        lengths, angles, heights = _cell_lengths_angles(cell)
        if is_periodic:
            min_periodic_height = min(heights[ax] for ax in range(3) if pbc[ax])

    atoms = model["atoms"]
    n_atoms = len(atoms)
    id_to_idx = {int(a["id"]): i for i, a in enumerate(atoms)}
    coords = np.asarray([a["position"] for a in atoms], dtype=float)

    adj: list[list[int]] = [[] for _ in range(n_atoms)]
    bond_image: dict[tuple[int, int], np.ndarray] = {}
    max_bond_mic_err = 0.0
    n_periodic_bonds = 0

    for b in model["bonds"]:
        i, j = id_to_idx[int(b["a"])], id_to_idx[int(b["b"])]
        img = _parse_int_image(b.get("image"), pbc, f"bond ({b['a']}, {b['b']})")
        if any(img):
            n_periodic_bonds += 1
        img_arr = np.asarray(img, dtype=int)
        if (i, j) in bond_image:
            raise ValueError(f"Duplicate bond entry between atoms {b['a']} and {b['b']}")
        bond_image[(i, j)] = img_arr
        bond_image[(j, i)] = -img_arr
        adj[i].append(j)
        adj[j].append(i)

        if cell is not None and inv_cell is not None:
            dr_exp = coords[j] + img_arr.astype(float) @ cell - coords[i]
            d_exp = float(np.linalg.norm(dr_exp))
            dr_mic = _mic_vector(dr_exp, pbc, cell, inv_cell)
            err = float(np.linalg.norm(dr_exp - dr_mic))
            max_bond_mic_err = max(max_bond_mic_err, err)
            if is_periodic and d_exp >= 0.5 * min_periodic_height:
                raise ValueError(
                    f"Periodic cell is too small for unambiguous bond ({b['a']}, {b['b']}): "
                    f"bond length {d_exp:.3f} A >= 0.5 * min box height ({0.5 * min_periodic_height:.3f} A)"
                )
            if d_exp <= 0.5 or err > 1e-5:
                raise ValueError(
                    f"Bond ({b['a']}, {b['b']}) vector with image {img} disagrees with MIC geometry "
                    f"(len={d_exp:.4f} A, mic_err={err:.2e} A)"
                )

    if is_periodic and cell is not None and inv_cell is not None:
        # Shortest bond-graph distances <= 3 across periodic lattice
        graph_dist: dict[tuple[int, int], int] = {}
        for start in range(n_atoms):
            dist_map = {start: 0}
            dq = deque([start])
            while dq:
                u = dq.popleft()
                du = dist_map[u]
                if du == 3:
                    continue
                for v in adj[u]:
                    if v not in dist_map:
                        dist_map[v] = du + 1
                        dq.append(v)
                        if start < v:
                            graph_dist[(start, v)] = dist_map[v]

        # Enforce unique lattice shift and strict MIC consistency across all 1-2, 1-3, 1-4 paths
        canonical_shift: dict[tuple[int, int], tuple[int, int, int]] = {}
        for (i, j), s in bond_image.items():
            if i < j:
                canonical_shift[(i, j)] = (int(s[0]), int(s[1]), int(s[2]))

        def _check_shift(u: int, v: int, s_uv: np.ndarray, kind: str) -> None:
            key = (u, v) if u < v else (v, u)
            s_canon = s_uv if u < v else -s_uv
            s_tup = (int(s_canon[0]), int(s_canon[1]), int(s_canon[2]))
            prev = canonical_shift.get(key)
            if prev is not None and prev != s_tup:
                raise ValueError(
                    f"Ambiguous periodic lattice winding across {kind} paths between atoms "
                    f"{atoms[key[0]]['id']} and {atoms[key[1]]['id']}: shifts {prev} vs {s_tup}"
                )
            canonical_shift.setdefault(key, s_tup)
            dr_exp = coords[key[1]] + s_canon.astype(float) @ cell - coords[key[0]]
            dr_mic = _mic_vector(coords[key[1]] - coords[key[0]], pbc, cell, inv_cell)
            err = float(np.linalg.norm(dr_exp - dr_mic))
            if err > 1e-4 or float(np.linalg.norm(dr_exp)) >= 0.5 * min_periodic_height:
                raise ValueError(
                    f"Periodic cell violates half-box MIC unambiguity on {kind} pair "
                    f"({atoms[key[0]]['id']}, {atoms[key[1]]['id']}): "
                    f"d={np.linalg.norm(dr_exp):.3f} A, half_height={0.5 * min_periodic_height:.3f} A"
                )

        # 1-3 angles (i - j - k)
        for j in range(n_atoms):
            nbs = adj[j]
            for a_i in range(len(nbs)):
                for a_k in range(a_i + 1, len(nbs)):
                    i, k = nbs[a_i], nbs[a_k]
                    if i == k:
                        raise ValueError("Degenerate periodic angle with identical endpoints")
                    s_ik = bond_image[(i, j)] + bond_image[(j, k)]
                    _check_shift(i, k, s_ik, "1-3 angle")

        # 1-4 dihedrals (i - j - k - m_idx, 4 distinct atoms)
        for (j, k), s_jk in bond_image.items():
            if j >= k:
                continue
            for i in adj[j]:
                if i == k:
                    continue
                s_ij = bond_image[(i, j)]
                for m_idx in adj[k]:
                    if m_idx in (i, j):
                        continue
                    pair = (min(i, m_idx), max(i, m_idx))
                    if graph_dist.get(pair, 99) <= 3:
                        s_im = s_ij + s_jk + bond_image[(k, m_idx)]
                        _check_shift(i, m_idx, s_im, "1-4 dihedral")

    return {
        "is_periodic": is_periodic,
        "pbc": pbc,
        "n_periodic_bonds": n_periodic_bonds,
        "box_lengths_A": lengths,
        "box_angles_deg": angles,
        "box_heights_nm": [0.1 * h for h in heights] if heights is not None else None,
        "min_box_height_nm": 0.1 * min(heights) if heights is not None else None,
        "min_periodic_height_nm": (0.1 * min_periodic_height) if math.isfinite(min_periodic_height) else None,
        "max_bond_mic_err_A": max_bond_mic_err,
    }


def _match_improper(
    c_idx: int,
    neighbors: list[int],
    atom_types: list[str],
    rules: list[tuple[str, str, str, str, float, float, int]],
) -> tuple[int, int, int, int, float, float, int] | None:
    """Match a degree-3 central atom deterministically to its single Amber improper rule."""
    if len(neighbors) != 3:
        return None
    tc = atom_types[c_idx]
    sorted_nbs = sorted(neighbors, key=lambda idx: (atom_types[idx], idx))
    for t1, t2, t3, t4, pk, phase, per in rules:
        if t3 != tc:
            continue
        for pos_d, nd in enumerate(sorted_nbs):
            if t4 not in ("X", atom_types[nd]):
                continue
            rem = [sorted_nbs[k] for k in range(3) if k != pos_d]
            na, nb = min(rem), max(rem)
            matches_ab = (t1 in ("X", atom_types[na]) and t2 in ("X", atom_types[nb]))
            matches_ba = (t1 in ("X", atom_types[nb]) and t2 in ("X", atom_types[na]))
            if matches_ab:
                return (na, nb, c_idx, nd, pk, phase, per)
            if matches_ba:
                return (nb, na, c_idx, nd, pk, phase, per)
    return None


def build_v3_parmed_structure(
    model: dict[str, Any],
    *,
    charge_model: str = "auto",
    custom_charges: str | Path | dict[str, Any] | list[float] | np.ndarray | None = None,
    reference_dir: str | Path | None = None,
    current_reference_dir: str | Path | None = None,
    allow_finite: bool = False,
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    """Build a parameterized ParmEd AmberParm directly from a GO graph and ``gaff2_v3`` frcmod."""
    import parmed as pmd

    cell_meta = validate_periodic_model(model, require_periodic=not allow_finite)
    frc = load_v3_frcmod_parameters(reference_dir=reference_dir)
    bundle = frc["bundle"]
    assigned = assign_gaff2_v3(
        model,
        reference_dir=reference_dir,
        verified_bundle=bundle,
        charge_model=charge_model,
        custom_charges=custom_charges,
        current_reference_dir=current_reference_dir,
    )

    from .topology import _residue_map

    residues = _residue_map(assigned)
    assigned["atoms"].sort(key=lambda atom: residues[(atom["resid"], atom["resname"])])

    atoms = assigned["atoms"]
    n_atoms = len(atoms)
    id_to_idx = {int(a["id"]): i for i, a in enumerate(atoms)}
    atom_types_list = [str(a["atom_type"]) for a in atoms]

    structure = pmd.Structure()
    pmd_atom_types: dict[str, pmd.AtomType] = {}
    pmd_atoms: list[pmd.Atom] = []
    for type_idx, atype in enumerate(sorted(set(atom_types_list)), 1):
        mass_val = frc["masses"].get(atype, (ATOMIC_MASSES["C"], 0.0))[0]
        el_guess = next(a["element"] for a in atoms if a["atom_type"] == atype)
        z_num = ATOMIC_NUMBERS[el_guess]
        if atype not in frc["nonbon"]:
            raise ValueError(f"Missing NONBON LJ parameters for v3 atom type {atype!r}")
        rmin_half, eps = frc["nonbon"][atype]
        at_obj = pmd.AtomType(atype, type_idx, mass_val, z_num)
        at_obj.set_lj_params(eps, rmin_half)
        pmd_atom_types[atype] = at_obj

    for atom in atoms:
        atype = str(atom["atom_type"])
        z_num = ATOMIC_NUMBERS[atom["element"]]
        mass_val = frc["masses"].get(atype, (ATOMIC_MASSES[atom["element"]], 0.0))[0]
        p_atom = pmd.Atom(
            name=str(atom["name"]),
            type=atype,
            atomic_number=z_num,
            mass=mass_val,
            charge=float(atom["charge"]),
        )
        p_atom.atom_type = pmd_atom_types[atype]
        structure.add_atom(p_atom, str(atom["resname"]), int(atom["resid"]))
        pmd_atoms.append(p_atom)

    # Bonds + adjacency
    adj: list[list[int]] = [[] for _ in range(n_atoms)]
    bond_type_cache: dict[tuple[str, str], pmd.BondType] = {}
    for b in assigned["bonds"]:
        i, j = id_to_idx[int(b["a"])], id_to_idx[int(b["b"])]
        adj[i].append(j)
        adj[j].append(i)
        bk = (min(atom_types_list[i], atom_types_list[j]), max(atom_types_list[i], atom_types_list[j]))
        if bk not in bond_type_cache:
            if bk not in frc["bonds"]:
                raise ValueError(f"Missing BOND parameters in GO_gaff2_v3.frcmod for {bk[0]}-{bk[1]}")
            kb, req = frc["bonds"][bk]
            bt = pmd.BondType(kb, req, list=structure.bond_types)
            structure.bond_types.append(bt)
            bond_type_cache[bk] = bt
        structure.bonds.append(pmd.Bond(pmd_atoms[i], pmd_atoms[j], type=bond_type_cache[bk]))

    # Angles (i < k around center j)
    angle_type_cache: dict[tuple[str, str, str], pmd.AngleType] = {}
    for j in range(n_atoms):
        nbs = sorted(adj[j])
        tj = atom_types_list[j]
        for a_i in range(len(nbs)):
            for a_k in range(a_i + 1, len(nbs)):
                i, k = nbs[a_i], nbs[a_k]
                ti, tk = atom_types_list[i], atom_types_list[k]
                ak = (min(ti, tk), tj, max(ti, tk))
                if ak not in angle_type_cache:
                    if ak not in frc["angles"]:
                        raise ValueError(
                            f"Missing ANGLE parameters in GO_gaff2_v3.frcmod for {ak[0]}-{ak[1]}-{ak[2]}"
                        )
                    ka, theteq = frc["angles"][ak]
                    at = pmd.AngleType(ka, theteq, list=structure.angle_types)
                    structure.angle_types.append(at)
                    angle_type_cache[ak] = at
                structure.angles.append(
                    pmd.Angle(pmd_atoms[i], pmd_atoms[j], pmd_atoms[k], type=angle_type_cache[ak])
                )

    # Shortest graph distance up to 3 (true 1-4 detection)
    dist3_pairs: set[tuple[int, int]] = set()
    for start in range(n_atoms):
        gdist = {start: 0}
        dq = deque([start])
        while dq:
            u = dq.popleft()
            du = gdist[u]
            if du == 3:
                continue
            for v in adj[u]:
                if v not in gdist:
                    nd = du + 1
                    gdist[v] = nd
                    dq.append(v)
                    if nd == 3 and start < v:
                        dist3_pairs.add((start, v))

    # Proper dihedrals (canonical simple 3-bond paths i-j-k-m_idx, 4 distinct atoms)
    seen_dihe_paths: set[tuple[int, int, int, int]] = set()
    active_14_assigned: set[tuple[int, int]] = set()
    dihe_type_cache: dict[tuple[float, int, float, float, float], pmd.DihedralType] = {}

    def _get_dihe_type(phi_k: float, per: int, phase: float, scee: float, scnb: float) -> pmd.DihedralType:
        dt_key = (phi_k, per, phase, scee, scnb)
        if dt_key not in dihe_type_cache:
            dt = pmd.DihedralType(phi_k, per, phase, scee=scee, scnb=scnb, list=structure.dihedral_types)
            structure.dihedral_types.append(dt)
            dihe_type_cache[dt_key] = dt
        return dihe_type_cache[dt_key]

    proper_candidates: list[tuple[int, int, int, int]] = []
    for j in range(n_atoms):
        for k in adj[j]:
            if j == k:
                continue
            for i in adj[j]:
                if i == k:
                    continue
                for m_idx in adj[k]:
                    if m_idx in (i, j):
                        continue
                    fwd = (i, j, k, m_idx)
                    rev = (m_idx, k, j, i)
                    canon = min(fwd, rev)
                    if canon not in seen_dihe_paths:
                        seen_dihe_paths.add(canon)
                        proper_candidates.append(canon)

    for i, j, k, m_idx in sorted(proper_candidates):
        t4_fwd = (atom_types_list[i], atom_types_list[j], atom_types_list[k], atom_types_list[m_idx])
        t4_rev = tuple(reversed(t4_fwd))
        terms = frc["dihe_specific"].get(min(t4_fwd, t4_rev))
        if terms is None:
            wk = (min(atom_types_list[j], atom_types_list[k]), max(atom_types_list[j], atom_types_list[k]))
            terms = frc["dihe_wildcard"].get(wk)
        if not terms:
            raise ValueError(
                f"Missing DIHE parameters in GO_gaff2_v3.frcmod for {'-'.join(t4_fwd)}"
            )
        end_pair = (min(i, m_idx), max(i, m_idx))
        needs_active_14 = end_pair in dist3_pairs and end_pair not in active_14_assigned
        if needs_active_14:
            active_14_assigned.add(end_pair)
        for term_idx, (idivf, pk, phase, per) in enumerate(terms):
            is_last_term = term_idx == len(terms) - 1
            ignore_end = not (needs_active_14 and is_last_term)
            dt = _get_dihe_type(pk / idivf, per, phase, 1.2, 2.0)
            structure.dihedrals.append(
                pmd.Dihedral(
                    pmd_atoms[i],
                    pmd_atoms[j],
                    pmd_atoms[k],
                    pmd_atoms[m_idx],
                    improper=False,
                    ignore_end=ignore_end,
                    type=dt,
                )
            )

    if active_14_assigned != dist3_pairs:
        raise ValueError(
            f"1-4 dihedral assignment incomplete: assigned {len(active_14_assigned)} of {len(dist3_pairs)} dist-3 pairs"
        )

    # Improper dihedrals (at most one per trigonal sp2 center atom3 == center)
    n_impropers = 0
    for c_idx in range(n_atoms):
        matched = _match_improper(c_idx, adj[c_idx], atom_types_list, frc["impropers"])
        if matched is None:
            continue
        na, nb, nc, nd, pk, phase, per = matched
        dt = _get_dihe_type(pk, per, phase, 0.0, 0.0)
        structure.dihedrals.append(
            pmd.Dihedral(
                pmd_atoms[na],
                pmd_atoms[nb],
                pmd_atoms[nc],
                pmd_atoms[nd],
                improper=True,
                ignore_end=True,
                type=dt,
            )
        )
        n_impropers += 1

    structure.coordinates = np.asarray([a["position"] for a in atoms], dtype=float)
    if cell_meta["box_lengths_A"] is not None and cell_meta["box_angles_deg"] is not None:
        structure.box = [*cell_meta["box_lengths_A"], *cell_meta["box_angles_deg"]]

    amber_parm = pmd.amber.AmberParm.from_structure(structure)
    amber_parm.rediscover_molecules(fix_broken=False)

    meta = {
        **cell_meta,
        "frcmod_sha256": bundle["frcmod_sha256"],
        "bundle_provenance_sha256": bundle["provenance_sha256"],
        "atom_count": len(amber_parm.atoms),
        "bond_count": len(amber_parm.bonds),
        "angle_count": len(amber_parm.angles),
        "proper_dihedral_count": sum(1 for d in amber_parm.dihedrals if not d.improper),
        "improper_dihedral_count": n_impropers,
        "dist3_14_pair_count": len(dist3_pairs),
        "topology_index_to_model_id": [int(a["id"]) for a in atoms],
        "model_id_to_topology_index": id_to_idx,
    }
    return amber_parm, assigned, meta


def _resolve_periodic_cutoff_nm(cell_meta: dict[str, Any], cutoff_nm: float | None) -> float:
    min_h = cell_meta.get("min_box_height_nm")
    if min_h is None:
        raise ValueError("Periodic box height unavailable for cutoff validation")
    max_allowed = 0.5 * float(min_h)
    if cutoff_nm is None:
        candidate = min(1.0, math.floor(0.45 * float(min_h) * 1000.0) / 1000.0)
        if candidate <= 0.0:
            raise ValueError(f"Periodic cell height ({min_h:.4f} nm) is too small")
        return float(candidate)
    val = float(cutoff_nm)
    if not math.isfinite(val) or val <= 0.0:
        raise ValueError("cutoff_nm must be positive and finite")
    if val >= max_allowed:
        raise ValueError(
            f"PME cutoff_nm={val} must be strictly less than half the periodic box height ({max_allowed:.6f} nm)"
        )
    return val


def _check_exported_mic_bonds(model: dict[str, Any], coords_ang: np.ndarray, id_to_idx: dict[int, int]) -> float:
    cell = np.asarray(model["cell"], dtype=float)
    inv_cell = np.linalg.inv(cell)
    pbc = [bool(x) for x in model["pbc"]]
    orig_xyz = np.asarray([a["position"] for a in model["atoms"]], dtype=float)
    max_err = 0.0
    for b in model["bonds"]:
        i, j = id_to_idx[int(b["a"])], id_to_idx[int(b["b"])]
        d_gro = float(np.linalg.norm(_mic_vector(coords_ang[j] - coords_ang[i], pbc, cell, inv_cell)))
        img = np.asarray(b.get("image", [0, 0, 0]), dtype=float)
        d_orig = float(np.linalg.norm(orig_xyz[j] + img @ cell - orig_xyz[i]))
        err = abs(d_gro - d_orig)
        max_err = max(max_err, err)
        if err > 0.02 or d_gro >= 2.2:
            raise ValueError(
                f"Exported GRO MIC bond check failed on ({b['a']}, {b['b']}): "
                f"gro_mic={d_gro:.4f} A, model={d_orig:.4f} A"
            )
    return max_err


def export_periodic_gromacs(
    model: dict[str, Any],
    output_prefix: str | Path,
    *,
    charge_model: str = "auto",
    charges_json: str | Path | None = None,
    reference_dir: str | Path | None = None,
    current_reference_dir: str | Path | None = None,
    cutoff_nm: float | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Export a covalently periodic GO model directly to GROMACS with ``periodic-molecules = yes``."""
    import parmed as pmd

    from .topology import (
        ENGINE_SOURCES,
        _json,
        _upgrade_gromacs_precision,
        _verify_parameter_roundtrip,
        _verify_supported_amber_terms,
    )

    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh prefix")
    if charges_json is not None and charge_model != "custom":
        raise ValueError(f"charges_json requires charge_model='custom' (got {charge_model!r})")

    prefix = Path(output_prefix).resolve()
    paths = {
        "top": prefix.with_suffix(".top"),
        "itp": prefix.with_suffix(".itp"),
        "gro": prefix.with_suffix(".gro"),
        "mdp": prefix.parent / f"{prefix.name}_periodic.mdp",
        "manifest": prefix.with_suffix(".manifest.json"),
        "provenance": prefix.with_suffix(".provenance.json"),
    }
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Periodic GROMACS export files already exist; choose a fresh prefix")

    amber_parm, assigned, meta = build_v3_parmed_structure(
        model,
        charge_model=charge_model,
        custom_charges=charges_json,
        reference_dir=reference_dir,
        current_reference_dir=current_reference_dir,
        allow_finite=False,
    )
    resolved_cutoff = _resolve_periodic_cutoff_nm(meta, cutoff_nm)
    _verify_supported_amber_terms(amber_parm, target="GROMACS", allow_explicit_exclusions=False)
    converted = pmd.gromacs.GromacsTopologyFile.from_structure(amber_parm)

    created_dirs: list[Path] = []
    cur_dir = prefix.parent
    while not cur_dir.exists():
        created_dirs.append(cur_dir)
        if cur_dir.parent == cur_dir:
            break
        cur_dir = cur_dir.parent
    prefix.parent.mkdir(parents=True, exist_ok=True)
    published: list[Path] = []
    try:
        with tempfile.TemporaryDirectory(dir=prefix.parent, prefix=f".{prefix.name}.tmp.") as stage_dir_str:
            stage_dir = Path(stage_dir_str)
            staged = {key: stage_dir / fp.name for key, fp in paths.items()}
            converted.write(str(staged["top"]), molfile=str(staged["itp"]), combine="all")
            top_text = staged["top"].read_text(encoding="utf-8").replace(
                f'"{staged["itp"]}"', f'"{paths["itp"].name}"'
            )
            header_note = "; GO Workbench direct periodic GROMACS topology (requires periodic-molecules = yes)\n"
            staged["top"].write_text(header_note + top_text, encoding="utf-8")

            bonds_upgraded = _upgrade_gromacs_precision(
                staged["top"], staged["itp"], amber_parm, converted
            )
            amber_parm.save(str(staged["gro"]), overwrite=False)
            roundtrip = pmd.load_file(str(staged["top"]), xyz=str(staged["gro"]))
            verification = _verify_parameter_roundtrip(
                amber_parm, roundtrip, top_path=staged["top"]
            )
            max_mic_bond_err = _check_exported_mic_bonds(
                assigned,
                np.asarray(roundtrip.coordinates, dtype=float),
                meta["model_id_to_topology_index"],
            )
            verification["max_gro_mic_bond_err_A"] = max_mic_bond_err
            verification["periodic_molecules"] = True

            mdp_text = "\n".join(
                [
                    "; Ready-to-use GROMACS energy-minimization MDP for covalently periodic GO",
                    "integrator               = steep",
                    "emtol                    = 100.0",
                    "emstep                   = 0.005",
                    "nsteps                   = 5000",
                    "pbc                      = xyz",
                    "periodic-molecules       = yes",
                    "cutoff-scheme            = Verlet",
                    "coulombtype              = PME",
                    f"rcoulomb                 = {resolved_cutoff:.4f}",
                    "vdwtype                  = Cut-off",
                    "vdw-modifier             = None",
                    f"rvdw                     = {resolved_cutoff:.4f}",
                    "DispCorr                 = EnerPres",
                    "",
                ]
            )
            staged["mdp"].write_text(mdp_text, encoding="utf-8")

            ff_meta = assigned.get("metadata", {}).get("force_field", {})
            prov = {
                "schema_version": 1,
                "format": "gromacs_periodic",
                "preset": ff_meta.get("preset", "gaff2_v3"),
                "charge_model": ff_meta.get("charge_model"),
                "frcmod_sha256": meta["frcmod_sha256"],
                "bundle_provenance_sha256": meta["bundle_provenance_sha256"],
                "itp_bond_precision_upgraded": bonds_upgraded,
                "cutoff_nm": resolved_cutoff,
                "periodic_metadata": {
                    k: meta[k]
                    for k in (
                        "pbc",
                        "n_periodic_bonds",
                        "box_lengths_A",
                        "box_angles_deg",
                        "box_heights_nm",
                        "max_bond_mic_err_A",
                    )
                },
                "outputs": {k: str(v) for k, v in paths.items() if k != "provenance"},
            }
            _json(staged["provenance"], prov)
            result = {
                **{k: str(v) for k, v in paths.items()},
                "status": "exported",
                "parameterized": True,
                "direct_periodic": True,
                "periodic_molecules": "yes",
                "cutoff_nm": resolved_cutoff,
                "roundtrip_verified": verification,
                "itp_bond_precision_upgraded": bonds_upgraded,
                "preset": ff_meta.get("preset", "gaff2_v3"),
                "charge_model": ff_meta.get("charge_model"),
                "atom_count": meta["atom_count"],
                "bond_count": meta["bond_count"],
                "angle_count": meta["angle_count"],
                "proper_dihedral_count": meta["proper_dihedral_count"],
                "improper_dihedral_count": meta["improper_dihedral_count"],
                "net_charge": float(sum(a["charge"] for a in assigned["atoms"])),
                "engine_sources": ENGINE_SOURCES,
            }
            _json(staged["manifest"], result)
            for key in ("top", "itp", "gro", "mdp", "provenance", "manifest"):
                os.link(staged[key], paths[key])
                published.append(paths[key])
    except Exception:
        for fp in published:
            with contextlib.suppress(FileNotFoundError):
                fp.unlink()
        for d in created_dirs:
            with contextlib.suppress(OSError):
                d.rmdir()
        raise
    return result


def export_periodic_openmm(
    model: dict[str, Any],
    output_xml: str | Path,
    *,
    nonbonded_method: str = "PME",
    cutoff_nm: float | None = None,
    charge_model: str = "auto",
    charges_json: str | Path | None = None,
    reference_dir: str | Path | None = None,
    current_reference_dir: str | Path | None = None,
    overwrite: bool = False,
    _force_pure_python: bool = False,
) -> dict[str, Any]:
    """Export a covalently periodic GO model directly to an OpenMM ``<System>`` XML."""
    from .topology import (
        ENGINE_SOURCES,
        _amber_parm_to_openmm_xml,
        _json,
    )

    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh output path")
    if nonbonded_method != "PME":
        raise ValueError(
            f"Periodic covalent OpenMM export requires nonbonded_method='PME' (got {nonbonded_method!r})"
        )
    if charges_json is not None and charge_model != "custom":
        raise ValueError(f"charges_json requires charge_model='custom' (got {charge_model!r})")

    xml_path = Path(output_xml).resolve()
    paths = {
        "xml": xml_path,
        "pdb": xml_path.with_suffix(".pdb"),
        "inpcrd": xml_path.with_suffix(".inpcrd"),
        "manifest": xml_path.with_suffix(".manifest.json"),
        "provenance": xml_path.with_suffix(".provenance.json"),
    }
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Periodic OpenMM export files already exist; choose a fresh output path")

    amber_parm, assigned, meta = build_v3_parmed_structure(
        model,
        charge_model=charge_model,
        custom_charges=charges_json,
        reference_dir=reference_dir,
        current_reference_dir=current_reference_dir,
        allow_finite=False,
    )
    resolved_cutoff = _resolve_periodic_cutoff_nm(meta, cutoff_nm)
    xml_text, counts = _amber_parm_to_openmm_xml(
        amber_parm,
        nonbonded_method="PME",
        cutoff_nm=resolved_cutoff,
        periodic_covalent=True,
    )

    serializer = "pure_python_periodic_amber_xml"
    if not _force_pure_python:
        try:
            import openmm

            system = openmm.XmlSerializer.deserialize(xml_text)
            if system.getNumParticles() != len(amber_parm.atoms):
                raise RuntimeError("Periodic OpenMM XML round-trip particle count mismatch")
            serializer = f"pure_python_periodic_amber_xml+openmm_{openmm.__version__}_verified"
        except ImportError:
            pass

    created_dirs: list[Path] = []
    cur_dir = xml_path.parent
    while not cur_dir.exists():
        created_dirs.append(cur_dir)
        if cur_dir.parent == cur_dir:
            break
        cur_dir = cur_dir.parent
    xml_path.parent.mkdir(parents=True, exist_ok=True)
    published: list[Path] = []
    try:
        with tempfile.TemporaryDirectory(dir=xml_path.parent, prefix=f".{xml_path.name}.tmp.") as stage_dir_str:
            stage_dir = Path(stage_dir_str)
            staged = {key: stage_dir / fp.name for key, fp in paths.items()}
            xml_bytes = xml_text.encode("utf-8")
            staged["xml"].write_bytes(xml_bytes)
            amber_parm.save(str(staged["pdb"]), overwrite=False)
            amber_parm.save(str(staged["inpcrd"]), overwrite=False)

            ff_meta = assigned.get("metadata", {}).get("force_field", {})
            prov = {
                "schema_version": 1,
                "format": "openmm_periodic_system_xml",
                "preset": ff_meta.get("preset", "gaff2_v3"),
                "charge_model": ff_meta.get("charge_model"),
                "frcmod_sha256": meta["frcmod_sha256"],
                "bundle_provenance_sha256": meta["bundle_provenance_sha256"],
                "nonbonded_method": "PME",
                "cutoff_nm": resolved_cutoff,
                "uses_periodic_bonded": True,
                "exceptions_use_periodic": True,
                "serializer": serializer,
                "counts": counts,
                "periodic_metadata": {
                    k: meta[k]
                    for k in ("pbc", "n_periodic_bonds", "box_lengths_A", "box_angles_deg", "box_heights_nm")
                },
                "outputs": {k: str(v) for k, v in paths.items() if k != "provenance"},
            }
            _json(staged["provenance"], prov)
            result = {
                **{k: str(v) for k, v in paths.items()},
                "status": "exported",
                "parameterized": True,
                "direct_periodic": True,
                "output": str(xml_path),
                "serializer": serializer,
                "nonbonded_method": "PME",
                "cutoff_nm": resolved_cutoff,
                "uses_periodic_bonded": True,
                "exceptions_use_periodic": True,
                "xml_sha256": hashlib.sha256(xml_bytes).hexdigest(),
                "preset": ff_meta.get("preset", "gaff2_v3"),
                "charge_model": ff_meta.get("charge_model"),
                "atom_count": counts["atoms"],
                "bond_count": counts["bonds"],
                "angle_count": counts["angles"],
                "torsion_count": counts["torsions"],
                "exception_count": counts["exceptions"],
                "net_charge": float(sum(a["charge"] for a in assigned["atoms"])),
                "engine_sources": ENGINE_SOURCES,
            }
            _json(staged["manifest"], result)
            for key in ("xml", "pdb", "inpcrd", "provenance", "manifest"):
                os.link(staged[key], paths[key])
                published.append(paths[key])
    except Exception:
        for fp in published:
            with contextlib.suppress(FileNotFoundError):
                fp.unlink()
        for d in created_dirs:
            with contextlib.suppress(OSError):
                d.rmdir()
        raise
    return result


def export_periodic_amber(
    model: dict[str, Any],
    output_dir: str | Path,
    *,
    charge_model: str = "auto",
    charges_json: str | Path | None = None,
    reference_dir: str | Path | None = None,
    current_reference_dir: str | Path | None = None,
    radii: str = "mbondi2",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Export a covalently periodic GO model directly to Amber ``prmtop`` and ``inpcrd``."""
    import parmed as pmd

    from .topology import (
        ENGINE_SOURCES,
        _json,
        _sha,
        _write_amber_provenance,
    )

    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh output directory")
    if radii not in ("bondi", "mbondi", "mbondi2", "mbondi3"):
        raise ValueError("Unsupported GB radii set")
    if charges_json is not None and charge_model != "custom":
        raise ValueError(f"charges_json requires charge_model='custom' (got {charge_model!r})")

    out_dir = Path(output_dir).resolve()
    if out_dir.exists():
        raise FileExistsError(f"Choose a fresh Amber output directory: {out_dir}")

    amber_parm, assigned, meta = build_v3_parmed_structure(
        model,
        charge_model=charge_model,
        custom_charges=charges_json,
        reference_dir=reference_dir,
        current_reference_dir=current_reference_dir,
        allow_finite=False,
    )
    pmd.tools.changeRadii(amber_parm, radii).execute()

    created_dirs: list[Path] = []
    cur_dir = out_dir
    while not cur_dir.exists():
        created_dirs.append(cur_dir)
        if cur_dir.parent == cur_dir:
            break
        cur_dir = cur_dir.parent
    out_dir.mkdir(parents=True)
    try:
        top_path = out_dir / "go.prmtop"
        crd_path = out_dir / "go.inpcrd"
        amber_parm.save(str(top_path), overwrite=False)
        amber_parm.save(str(crd_path), overwrite=False)

        reloaded = pmd.load_file(str(top_path), xyz=str(crd_path))
        if int(reloaded.ptr("IFBOX")) <= 0:
            raise ValueError("Periodic Amber prmtop missing box flag (IFBOX <= 0)")
        mols = list(reloaded.parm_data.get("ATOMS_PER_MOLECULE") or [])
        if mols != [len(reloaded.atoms)]:
            raise ValueError(
                f"Periodic sheet must remain a single macromolecule in prmtop, got ATOMS_PER_MOLECULE={mols}"
            )
        expected_q = np.asarray([float(a["charge"]) for a in assigned["atoms"]], dtype=float)
        reloaded_q = np.asarray([float(a.charge) for a in reloaded.atoms], dtype=float)
        if not np.allclose(expected_q, reloaded_q, atol=1e-6, rtol=0):
            raise ValueError("Periodic Amber prmtop charge round-trip drift beyond 1e-6 e")

        assigned_copy = deepcopy(assigned)
        _json(out_dir / "assigned_model.json", assigned_copy)
        top_sha = _sha(top_path)
        crd_sha = _sha(crd_path)
        prov_path = _write_amber_provenance(
            out_dir,
            preset=assigned.get("metadata", {}).get("force_field", {}).get("preset", "gaff2_v3"),
            frcmod_name="GO_gaff2_v3.frcmod",
            frcmod_sha256=meta["frcmod_sha256"],
            assigned=assigned,
            radii=radii,
            topology_sha256=top_sha,
            coordinates_sha256=crd_sha,
        )
        ff_meta = assigned.get("metadata", {}).get("force_field", {})
        result = {
            "status": "parameterized",
            "direct_periodic": True,
            "preset": ff_meta.get("preset", "gaff2_v3"),
            "charge_model": ff_meta.get("charge_model"),
            "parameterized": True,
            "topology": str(top_path),
            "coordinates": str(crd_path),
            "assigned_model": str(out_dir / "assigned_model.json"),
            "provenance": str(prov_path),
            "radii_set": radii,
            "ifbox": int(reloaded.ptr("IFBOX")),
            "atoms_per_molecule": mols,
            "atom_mapping": {str(i): aid for i, aid in enumerate(meta["topology_index_to_model_id"], 1)},
            "atom_count": meta["atom_count"],
            "bond_count": meta["bond_count"],
            "angle_count": meta["angle_count"],
            "proper_dihedral_count": meta["proper_dihedral_count"],
            "improper_dihedral_count": meta["improper_dihedral_count"],
            "net_charge": float(expected_q.sum()),
            "topology_sha256": top_sha,
            "coordinates_sha256": crd_sha,
            "periodic_note": (
                "Periodic covalent MD in Amber sander/pmemd must set iwrap=0 to avoid wrapping "
                "atoms within a single periodic sheet molecule; prefer GROMACS (periodic-molecules=yes) "
                "or OpenMM periodic export for periodic bonded sheets."
            ),
            "engine_sources": ENGINE_SOURCES,
        }
        _json(out_dir / "manifest.json", result)
        return result
    except Exception:
        for fp in out_dir.glob("*"):
            with contextlib.suppress(OSError):
                fp.unlink()
        for d in created_dirs:
            with contextlib.suppress(OSError):
                d.rmdir()
        raise
