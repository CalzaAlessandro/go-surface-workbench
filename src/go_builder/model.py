"""Explicit molecular graph validation and lossless JSON exchange.

Coordinates are Å, atom IDs are stable integers. PDB is a finite visualization
export only; JSON preserves image bonds and is the authoritative graph.
"""

from __future__ import annotations

import json
import math
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


def write_pdb(model: dict, path) -> None:
    """Write explicit finite CONECT records; never flatten periodic seam bonds."""
    validate_model(model)
    if any(any(b.get("image", [0, 0, 0])) for b in model["bonds"]):
        raise ValueError("PDB cannot preserve periodic image bonds; use JSON")
    if any(a["id"] > 99999 or a.get("resid", 1) > 9999 for a in model["atoms"]):
        raise ValueError("PDB identifier range exceeded; use JSON")
    lines = ["REMARK Unrelaxed constructed graph; authoritative connectivity is in JSON.\n"]
    for atom in model["atoms"]:
        x, y, z = atom["position"]
        if any(abs(v) >= 9999 for v in (x, y, z)):
            raise ValueError("PDB coordinate range exceeded")
        lines.append(
            f"HETATM{atom['id']:5d} {atom.get('name', atom['element']):<4.4s} "
            f"{atom.get('resname', 'GO'):>3.3s} A{atom.get('resid', 1):4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom['element']:>2s}\n"
        )
    for bond in model["bonds"]:
        lines.append(f"CONECT{bond['a']:5d}{bond['b']:5d}\n")
    lines.append("END\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("x") as handle:
        handle.writelines(lines)
