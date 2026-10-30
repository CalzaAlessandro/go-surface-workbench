"""L0 transferable GO blocks: system RESP primary, explicit BCC projection control."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

from . import resp as atomic

TARGETS = {
    "GRT": 0.0,
    "OH": 0.0,
    "EPX": 0.0,
    "EPH": 0.0,
    "EPT": 0.0,
    "CAB": 0.0,
    "COOH": 0.0,
    "COO": -1.0,
}
STATUS = "PILOT_LIMITED_CONTEXT"
BLOCK_REGISTRY = {
    "GGG": {"block_type": "scaffold", "current_role_charges_e": {"scaffold_C": 0.0}},
    "GRT": {
        "block_type": "GRT",
        "current_role_charges_e": {"edge_C": -0.15, "edge_H": 0.15},
    },
    "HYD": {
        "block_type": "OH",
        "current_role_charges_e": {
            "hydroxyl_C": 0.25,
            "hydroxyl_O": -0.64,
            "hydroxyl_H": 0.39,
        },
    },
    "EPX": {
        "block_type": "EPX",
        "current_role_charges_e": {"epoxide_C": 0.185, "epoxide_O": -0.37},
    },
    "EPH": {
        "block_type": "EPH",
        "current_role_charges_e": {
            "EPH_C": 0.053,
            "EPH_CH": 0.137,
            "EPH_H": 0.12,
            "EPH_O": -0.31,
        },
    },
    "EPT": {
        "block_type": "EPT",
        "current_role_charges_e": {"EPT_C": 0.02, "EPT_H": 0.12, "EPT_O": -0.28},
    },
    "CAB": {
        "block_type": "CAB",
        "current_role_charges_e": {"carbonyl_C": 0.581, "carbonyl_O": -0.581},
    },
    "CBX": {
        "block_type": "COOH",
        "current_role_charges_e": {
            "COOH_anchor": -0.115,
            "COOH_C": 0.735,
            "COOH_O": -0.475,
            "COOH_OH": -0.605,
            "COOH_H": 0.46,
        },
    },
    "CBD": {
        "block_type": "COO",
        "current_role_charges_e": {"COO_anchor": -0.18, "COO_C": 0.76, "COO_O": -0.79},
    },
}


def _read(value):
    return (
        json.loads(Path(value).read_text())
        if isinstance(value, (str, Path))
        else copy.deepcopy(value)
    )


def _dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def partition_model(value):
    """Partition exact primary roles and intended bonds; never absorb FG into background."""
    model = _read(value)
    atomic._validate_model(model)
    atoms = {a["id"]: a for a in model["atoms"]}
    ns = {i: {} for i in atoms}
    for b in model["bonds"]:
        ns[b["a"]][b["b"]] = b["order"]
        ns[b["b"]][b["a"]] = b["order"]
    groups = []
    used = set()

    def neighbors(i, role):
        return sorted(j for j in ns[i] if atoms[j]["role"] == role)

    def one(i, role):
        matches = neighbors(i, role)
        if len(matches) != 1:
            raise ValueError(f"Incomplete/ambiguous {role} near atom {i}")
        return matches[0]

    def add(kind, ids, expected):
        if (
            Counter(atoms[i]["role"] for i in ids) != Counter(expected)
            or len(ids) != len(set(ids))
            or used.intersection(ids)
        ):
            raise ValueError(f"Overlapping or malformed {kind} block")
        group_formal = sum(atoms[i]["formal_charge"] for i in ids)
        if group_formal != TARGETS[kind]:
            raise ValueError(f"{kind} block formal charge mismatch")
        used.update(ids)
        groups.append(
            {
                "block_type": kind,
                "atom_ids": ids,
                "roles": [atoms[i]["role"] for i in ids],
                "target_charge_e": TARGETS[kind],
            }
        )

    for a in model["atoms"]:
        i = a["id"]
        role = a["role"]
        if role == "hydroxyl_O":
            c, h = one(i, "hydroxyl_C"), one(i, "hydroxyl_H")
            if (
                len(ns[i]) != 2
                or len(ns[h]) != 1
                or ns[i][h] != 1
                or ns[i][c] != 1
                or len(ns[c]) != 4
            ):
                raise ValueError("OH graph is not basal C-O-H")
            add("OH", [c, i, h], ["hydroxyl_C", "hydroxyl_O", "hydroxyl_H"])
        elif role == "epoxide_O":
            cs = neighbors(i, "epoxide_C")
            if (
                len(cs) != 2
                or len(ns[i]) != 2
                or cs[1] not in ns[cs[0]]
                or any(ns[i][c] != 1 or len(ns[c]) != 4 for c in cs)
            ):
                raise ValueError("EPX graph is not a C-C-O triangle")
            add("EPX", [*cs, i], ["epoxide_C", "epoxide_C", "epoxide_O"])
        elif role in ("EPH_O", "EPT_O"):
            kind = role[:3]
            cs = (
                neighbors(i, "EPH_C") + neighbors(i, "EPH_CH")
                if kind == "EPH"
                else neighbors(i, "EPT_C")
            )
            if (
                len(cs) != 2
                or len(ns[i]) != 2
                or cs[1] not in ns[cs[0]]
                or any(ns[i][c] != 1 for c in cs)
            ):
                raise ValueError("Terminal epoxide is not a C-C-O triangle")
            hs = []
            for c in cs:
                local = neighbors(c, kind + "_H")
                expected = 1 if atoms[c]["role"] in ("EPH_CH", "EPT_C") else 0
                if len(local) != expected or any(len(ns[h]) != 1 for h in local):
                    raise ValueError("Terminal epoxide hydrogen roles disagree")
                hs += local
            expected_roles = (
                ["EPH_C", "EPH_CH", "EPH_H", "EPH_O"]
                if kind == "EPH"
                else ["EPT_C", "EPT_C", "EPT_H", "EPT_H", "EPT_O"]
            )
            add(kind, [*cs, *hs, i], expected_roles)
        elif role == "carbonyl_O":
            c = one(i, "carbonyl_C")
            if len(ns[i]) != 1 or ns[i][c] != 2 or len(ns[c]) != 3:
                raise ValueError("CAB graph is not C=O")
            add("CAB", [c, i], ["carbonyl_C", "carbonyl_O"])
        elif role == "COOH_C":
            anchor, o, oh = one(i, "COOH_anchor"), one(i, "COOH_O"), one(i, "COOH_OH")
            h = one(oh, "COOH_H")
            if (
                len(ns[i]) != 3
                or ns[i][o] != 2
                or ns[i][oh] != 1
                or ns[i][anchor] != 1
                or len(ns[o]) != 1
                or len(ns[oh]) != 2
                or len(ns[h]) != 1
            ):
                raise ValueError("Malformed COOH graph")
            add(
                "COOH",
                [anchor, i, o, oh, h],
                ["COOH_anchor", "COOH_C", "COOH_O", "COOH_OH", "COOH_H"],
            )
        elif role == "COO_C":
            anchor = one(i, "COO_anchor")
            os = neighbors(i, "COO_O")
            if (
                len(os) != 2
                or len(ns[i]) != 3
                or sorted(ns[i][o] for o in os) != [1, 2]
                or any(len(ns[o]) != 1 for o in os)
            ):
                raise ValueError("Malformed COO graph/resonance representation")
            add("COO", [anchor, i, *os], ["COO_anchor", "COO_C", "COO_O", "COO_O"])
        elif role == "edge_H":
            c = one(i, "edge_C")
            if len(ns[i]) != 1 or len(ns[c]) != 3 or len(neighbors(c, "edge_H")) != 1:
                raise ValueError("Malformed GRT C-H cap")
            add("GRT", [c, i], ["edge_C", "edge_H"])
    scaffold = [a["id"] for a in model["atoms"] if a["role"] == "scaffold_C"]
    if any(
        atoms[i]["element"] != "C"
        or len(ns[i]) != 3
        or any(atoms[j]["element"] != "C" for j in ns[i])
        or atoms[i]["formal_charge"] != 0
        for i in scaffold
    ):
        raise ValueError("Background scaffold_C is not neutral three-carbon scaffold")
    if used.intersection(scaffold) or used | set(scaffold) != set(atoms):
        raise ValueError(
            f"Uncovered/overlapping chemistry: {sorted(set(atoms) - used - set(scaffold))}"
        )
    if sum(g["target_charge_e"] for g in groups) != model["formal_charge"]:
        raise ValueError("Block formal charge sums disagree with model")
    return {
        "model_id": model["model_id"],
        "groups": groups,
        "scaffold_atom_ids": scaffold,
        "formal_charge": model["formal_charge"],
        "full_coverage": True,
    }


def build_problem(values):
    models = [_read(v) for v in values]
    if not models or len({m["model_id"] for m in models}) != len(models):
        raise ValueError("Unique nonempty models required")
    if any(m["model_id"] == "pyrene_epoxide" for m in models):
        raise ValueError("pyrene_epoxide excluded: failed upstream QM bond gate")
    partitions = {m["model_id"]: partition_model(m) for m in models}
    if any(atomic.stage2_groups(m) for m in models):
        raise ValueError("L0 campaign stage2 supports no methyl/methylene groups")
    roles = sorted(
        {a["role"] for m in models for a in m["atoms"] if a["role"] != "scaffold_C"}
    )
    elements = {}
    for m in models:
        for a in m["atoms"]:
            if a["role"] in elements and elements[a["role"]] != a["element"]:
                raise ValueError("Shared role mixes elements")
            elements[a["role"]] = a["element"]
    constraints = []
    seen = {}
    for m in models:
        for g in partitions[m["model_id"]]["groups"]:
            coefficients = [g["roles"].count(role) for role in roles]
            if g["block_type"] in seen:
                if seen[g["block_type"]] != coefficients:
                    raise ValueError(
                        "Block occurrences have inconsistent role stoichiometry"
                    )
                continue
            seen[g["block_type"]] = coefficients
            constraints.append(
                {
                    "block_type": g["block_type"],
                    "coefficients": coefficients,
                    "target_charge_e": g["target_charge_e"],
                    "representative_model_id": m["model_id"],
                    "representative_atom_ids": g["atom_ids"],
                }
            )
    A = np.array([c["coefficients"] for c in constraints], float)
    b = np.array([c["target_charge_e"] for c in constraints])
    rank = int(np.linalg.matrix_rank(A))
    aug = int(np.linalg.matrix_rank(np.column_stack((A, b))))
    if rank != len(constraints) or aug != rank:
        raise ValueError("Dependent or infeasible shared block constraints")
    return {
        "status": STATUS,
        "model_order": [m["model_id"] for m in models],
        "models": models,
        "partitions": partitions,
        "block_registry": BLOCK_REGISTRY,
        "unfitted_terminal_blocks": [
            kind for kind in ("EPH", "EPT") if kind not in seen
        ],
        "unobserved_library_units": [
            unit
            for unit, entry in BLOCK_REGISTRY.items()
            if unit != "GGG" and entry["block_type"] not in seen
        ],
        "roles": roles,
        "role_elements": elements,
        "constraints": constraints,
        "constraint_rank": rank,
        "free_dimensions_after_constraints": len(roles) - rank,
        "fixed_roles": {"scaffold_C": 0.0},
        "scope": "L0 shared group roles; no distance/environment classes; EPX one mixed context only",
    }


def make_resp_input(problem, records, qwt=0.0005):
    counts = Counter(r["model_id"] for r in records)
    if set(counts) != set(problem["model_order"]) or set(counts.values()) != {3}:
        raise ValueError("Exactly three orientations per model required")
    if any(r["grid_points"] < 1 for r in records):
        raise ValueError("Positive ESP point counts required")
    lookup = {m["model_id"]: m for m in problem["models"]}
    nref = float(np.mean([r["grid_points"] for r in records]))
    molecules = []
    equivalences = []
    first_global = {}
    restraint_counts = Counter()
    for imol, r in enumerate(records, 1):
        local = {}
        rows = []
        for a in lookup[r["model_id"]]["atoms"]:
            role = a["role"]
            idx = a["id"]
            if role == "scaffold_C":
                vary = -1
            else:
                vary = local.get(role, 0)
                if role not in local:
                    local[role] = idx
                    if role in first_global:
                        equivalences.append([first_global[role], [imol, idx]])
                    else:
                        first_global[role] = [imol, idx]
                if a["element"] != "H":
                    restraint_counts[role] += 1
            rows.append(
                {
                    "atom_id": idx,
                    "role": role,
                    "atomic_number": atomic.ATOMIC_NUMBERS[a["element"]],
                    "ivary": vary,
                    "initial_charge_e": 0.0,
                }
            )
        molecules.append(
            {
                "model_id": r["model_id"],
                "header_charge": -99,
                "atoms": rows,
                "weight": float(f"{math.sqrt(nref / r['grid_points']):.5f}"),
                "grid_points": r["grid_points"],
            }
        )
    global_constraints = []
    for c in problem["constraints"]:
        imol = next(
            i
            for i, m in enumerate(molecules, 1)
            if m["model_id"] == c["representative_model_id"]
        )
        global_constraints.append(
            {
                "block_type": c["block_type"],
                "target_charge_e": c["target_charge_e"],
                "atom_pairs": [[imol, i] for i in c["representative_atom_ids"]],
            }
        )
    plan = {
        "molecules": molecules,
        "global_constraints": global_constraints,
        "intermolecular_equivalences": equivalences,
        "constraint_rank": problem["constraint_rank"],
        "mean_grid_points": nref,
        "qwt": qwt,
        "stage2": {"status": "noop_no_methyl_methylene"},
        "heavy_atom_copy_restraint_multiplicity": dict(restraint_counts),
        "weighting": "wtmol squared in RESP; equal3 orientations/model; objective divided by3 gives equal model orientation-mean MSE and physical atom-occurrence restraints",
    }
    lines = [
        "L0 MODULAR RESP PILOT LIMITED CONTEXT",
        " &cntrl",
        f" nmol = {len(molecules)},",
        " iqopt = 2,",
        " ihfree = 1,",
        " irstrnt = 1,",
        " ioutopt = 1,",
        f" qwt = {qwt},",
        " &end",
    ]
    for m in molecules:
        lines += [
            f"{m['weight']:10.5f}",
            m["model_id"],
            f"{-99:5d}{len(m['atoms']):5d}",
        ]
        lines += [f"{a['atomic_number']:5d}{a['ivary']:5d}" for a in m["atoms"]]
        lines.append("")
    for c in global_constraints:
        lines.append(f"{len(c['atom_pairs']):5d}{c['target_charge_e']:10.5f}")
        lines.append("".join(f"{imol:5d}{iatom:5d}" for imol, iatom in c["atom_pairs"]))
    lines.append("")
    for pair in equivalences:
        lines += [f"{2:5d}", "".join(f"{imol:5d}{iatom:5d}" for imol, iatom in pair)]
    lines.append("")
    return "\n".join(lines) + "\n", plan


def evaluate_library(problem, charges):
    charges = {k: float(v) for k, v in charges.items()}
    if (
        charges.get("scaffold_C") != 0
        or any(r not in charges for r in problem["roles"])
        or not np.isfinite(list(charges.values())).all()
    ):
        raise ValueError(
            "Library must cover finite shared roles with scaffold_C exactly zero"
        )
    groups = []
    models = []
    for m in problem["models"]:
        for g in problem["partitions"][m["model_id"]]["groups"]:
            total = sum(charges[r] for r in g["roles"])
            groups.append(
                {
                    "model_id": m["model_id"],
                    **g,
                    "charge_sum_e": total,
                    "charge_error_e": total - g["target_charge_e"],
                }
            )
        total = sum(charges[a["role"]] for a in m["atoms"])
        models.append(
            {
                "model_id": m["model_id"],
                "charge_sum_e": total,
                "formal_charge_e": m["formal_charge"],
                "charge_error_e": total - m["formal_charge"],
            }
        )
    return {"groups": groups, "models": models}


def current_library(problem):
    charges = {
        role: q
        for entry in BLOCK_REGISTRY.values()
        for role, q in entry["current_role_charges_e"].items()
    }
    charges = {
        role: q
        for role, q in charges.items()
        if role in problem["roles"] or role == "scaffold_C"
    }
    return {
        "method": "user_current_GAFF2_shared_block_control_same_zero_background",
        "shared_role_charges_e": charges,
        "evaluation": evaluate_library(problem, charges),
    }


def validate_fitted_charges(problem, plan, values):
    q = np.asarray(values, float)
    if (
        q.shape != (sum(len(m["atoms"]) for m in plan["molecules"]),)
        or not np.isfinite(q).all()
    ):
        raise ValueError("Invalid RESP charge count/finite values")
    shared = {"scaffold_C": 0.0}
    offset = 0
    for m in plan["molecules"]:
        for a, value in zip(m["atoms"], q[offset : offset + len(m["atoms"])]):
            role = a["role"]
            if role == "scaffold_C" and value != 0:
                raise ValueError("Frozen zero scaffold constraint failed")
            if role in shared and abs(shared[role] - value) > 2e-6:
                raise ValueError("Shared role equivalence failed")
            shared.setdefault(role, float(value))
        offset += len(m["atoms"])
    evaluation = evaluate_library(problem, shared)
    if any(
        abs(g["charge_error_e"]) > len(g["atom_ids"]) * 5.1e-7 + 2e-6
        for g in evaluation["groups"]
    ):
        raise ValueError("Block neutrality constraint failed")
    if any(
        abs(m["charge_error_e"])
        > len(
            next(v for v in problem["models"] if v["model_id"] == m["model_id"])[
                "atoms"
            ]
        )
        * 5.1e-7
        + 2e-6
        for m in evaluation["models"]
    ):
        raise ValueError("Model formal sum constraint failed")
    return {"shared_role_charges_e": shared, "evaluation": evaluation}


def project_am1_bcc(problem, raw):
    """Constrained L2 projection only; never substitutes for hyperbolic RESP."""
    raw = copy.deepcopy(raw)
    roles = problem["roles"]
    H = np.zeros(len(roles))
    z = H.copy()
    for m in problem["models"]:
        record = raw[m["model_id"]]
        atoms = record["atoms"]
        ids = [a["atom_id"] for a in atoms]
        if ids != [a["id"] for a in m["atoms"]]:
            raise ValueError("Raw AM1-BCC atom identity/order mismatch")
        q = np.array([a["charge_e"] for a in atoms], float)
        if not np.isfinite(q).all():
            raise ValueError("Nonfinite raw AM1-BCC charges")
        for atom, value in zip(m["atoms"], q):
            if atom["role"] == "scaffold_C":
                continue
            j = roles.index(atom["role"])
            H[j] += 1 / len(q)
            z[j] += value / len(q)
    if np.any(H <= 0):
        raise ValueError("Unobserved role in BCC projection")
    mu = z / H
    A = np.array([c["coefficients"] for c in problem["constraints"]], float)
    b = np.array([c["target_charge_e"] for c in problem["constraints"]])
    q = mu - (A.T @ np.linalg.solve((A / H) @ A.T, A @ mu - b)) / H
    shared = {"scaffold_C": 0.0, **dict(zip(roles, q.tolist()))}
    return {
        "status": STATUS,
        "method": "AM1_BCC_constrained_L2_block_projection_NOT_native_AM1_BCC",
        "objective": "equal-model mean squared atomic charge deviations; background fixed zero contributes constant loss",
        "raw_evidence": raw,
        "shared_role_charges_e": shared,
        "evaluation": evaluate_library(problem, shared),
    }


def _export(folder, problem, result):
    _dump(folder / "shared_library.json", result)
    with (folder / "shared_library.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["role", "charge_e"])
        writer.writerows(sorted(result["shared_role_charges_e"].items()))
    _dump(folder / "group_and_model_checks.json", result["evaluation"])


def _esp_geometry(model, coords):
    reference = np.asarray([a["position"] for a in model["atoms"]], float)
    if coords.shape != reference.shape:
        raise ValueError("ESP atom count differs from stable model")
    atomic._check_geometry(model, coords)
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    expected = np.linalg.norm(reference[:, None, :] - reference[None, :, :], axis=2)
    if not np.allclose(distances, expected, atol=5e-5, rtol=0):
        raise ValueError("ESP/model pair geometry mismatch")
    a, b = reference - reference.mean(axis=0), coords - coords.mean(axis=0)
    u, _, vh = np.linalg.svd(a.T @ b)
    rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vh)]) @ vh
    if not np.allclose(a @ rotation, b, atol=5e-5, rtol=0):
        raise ValueError("ESP geometry is reflection or changed atom identity")


def prepare_modular(spec, output_dir):
    entries = _read(spec)
    if isinstance(entries, dict):
        entries = entries["datasets"]
    problem = build_problem([e["model_path"] for e in entries])
    records = []
    texts = []
    evidence = []
    for entry, m in zip(entries, problem["models"]):
        paths = entry["esp_paths"]
        if len(paths) != 3 or len({str(Path(p).resolve()) for p in paths}) != 3:
            raise ValueError("Three distinct orientation ESP paths required")
        for path in paths:
            path = Path(path)
            data = atomic.read_esp(path)
            _esp_geometry(m, data[0] * atomic.BOHR_ANGSTROM)
            records.append({"model_id": m["model_id"], "grid_points": len(data[1])})
            texts.append(path.read_text().rstrip() + "\n")
            evidence.append(
                {
                    "model_id": m["model_id"],
                    "source_path": str(path.resolve()),
                    "sha256": atomic._digest(path),
                    "scope": entry.get("scope", "upstream_scope_not_declared"),
                }
            )
    cards, plan = make_resp_input(problem, records)
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=False)
    (folder / "resp.in").write_text(cards)
    (folder / "all.esp").write_text("".join(texts))
    nat = sum(len(m["atoms"]) for m in plan["molecules"])
    (folder / "qin.zero").write_text(
        "\n".join(
            "".join(f"{0.0:10.6f}" for _ in range(min(8, nat - i)))
            for i in range(0, nat, 8)
        )
        + "\n"
    )
    for i, (entry, text) in enumerate(zip(evidence, texts), 1):
        snapshot = folder / f"orientation_{i:03d}.esp"
        snapshot.write_text(text)
        entry["snapshot"] = snapshot.name
        entry["snapshot_sha256"] = atomic._digest(snapshot)
    for name, value in [
        ("problem.json", problem),
        ("plan.json", plan),
        ("inputs.json", evidence),
    ]:
        _dump(folder / name, value)
    return {
        "status": "prepared_not_fitted",
        "scope": STATUS,
        "output_dir": str(folder.resolve()),
        "constraint_rank": problem["constraint_rank"],
        "roles": problem["roles"],
        "atom_copies": nat,
        "stage2": plan["stage2"],
    }


def run_modular(output_dir, *, resp="resp"):
    folder = Path(output_dir).resolve()
    problem = _read(folder / "problem.json")
    plan = _read(folder / "plan.json")
    if (folder / "shared_library.json").exists():
        raise FileExistsError("Existing shared library result")
    program = atomic._executable(resp)
    atomic._command(
        [
            program["path"],
            "-i",
            "resp.in",
            "-o",
            "resp.out",
            "-p",
            "resp.punch",
            "-q",
            "qin.zero",
            "-e",
            "all.esp",
            "-t",
            "resp.charges",
            "-s",
            "resp.esout",
        ],
        folder,
        "system_resp",
    )
    for name in ("resp.out", "resp.punch", "resp.charges"):
        if not (folder / name).is_file() or not (folder / name).stat().st_size:
            raise ValueError("Missing system RESP output: " + name)
    text = (folder / "resp.out").read_text() + (folder / "system_resp.log").read_text()
    if re.search(
        r"\bNaN\b|\bInfinity\b|singular|not converg|no convergence|error|illegal|premature",
        text,
        re.IGNORECASE,
    ) or not re.search(r"Convergence in\s+\d+\s+iterations", text, re.IGNORECASE):
        raise ValueError("System RESP failure or missing convergence evidence")
    result = validate_fitted_charges(
        problem,
        plan,
        atomic._read_charges(
            folder / "resp.charges", sum(len(m["atoms"]) for m in plan["molecules"])
        ),
    )
    control = current_library(problem)
    _dump(folder / "current_library.json", control)
    lookup = {m["model_id"]: m for m in problem["models"]}
    metrics = {"modular": [], "current": []}
    for entry in _read(folder / "inputs.json"):
        data = atomic.read_esp(folder / entry["snapshot"])
        model = lookup[entry["model_id"]]
        for name, charges in (
            ("modular", result["shared_role_charges_e"]),
            ("current", control["shared_role_charges_e"]),
        ):
            q = np.array([charges[a["role"]] for a in model["atoms"]])
            metrics[name].append(
                {
                    "model_id": model["model_id"],
                    "snapshot": entry["snapshot"],
                    **atomic._metrics(*data, q),
                }
            )
    result.update(
        status=STATUS,
        method="system_Amber_RESP_L0_modular_hyperbolic",
        stage2=plan["stage2"],
        metrics=metrics,
        current_control=control,
        unfitted_terminal_blocks=problem["unfitted_terminal_blocks"],
        provenance={
            "resp": program,
            "code_sha256": atomic._digest(__file__),
            "inputs": _read(folder / "inputs.json"),
        },
        final_parameter_validation=False,
    )
    _export(folder, problem, result)
    return result


def fit_modular(spec, output_dir, *, resp="resp"):
    prepare_modular(spec, output_dir)
    return run_modular(output_dir, resp=resp)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "fit", "run", "project-bcc"])
    parser.add_argument("--dataset", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--resp", default="resp")
    parser.add_argument("--raw-bcc", type=Path)
    args = parser.parse_args(argv)
    if args.command == "run":
        result = run_modular(args.output_dir, resp=args.resp)
    elif args.command == "project-bcc":
        entries = _read(args.dataset)
        entries = entries["datasets"] if isinstance(entries, dict) else entries
        problem = build_problem([e["model_path"] for e in entries])
        result = project_am1_bcc(problem, _read(args.raw_bcc))
        args.output_dir.mkdir(parents=True, exist_ok=False)
        _export(args.output_dir, problem, result)
    else:
        result = (fit_modular if args.command == "fit" else prepare_modular)(
            args.dataset,
            args.output_dir,
            **({"resp": args.resp} if args.command == "fit" else {}),
        )
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
