"""Finite-fragment Gaussian preparation and matched-geometry ESP scoring.

No quantum, RESP, Amber, or scheduler program is launched by this module.
Coordinates are Angstrom unless a serialized reference explicitly says bohr;
ESP is electrostatic potential in Hartree per elementary charge, not energy.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

BOHR_ANGSTROM = 0.529177210903
E_ANGSTROM_DEBYE = 4.80320471257
HARTREE_KCAL_MOL = 627.5094740631
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16, "Cl": 17, "Br": 35, "I": 53}


def _directory(path):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    return path


def _vector(value, label):
    try:
        a = np.asarray(value, dtype=float)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid {label}") from exc
    if a.shape != (3,) or not np.isfinite(a).all():
        raise ValueError(f"{label} must be a finite three-vector")
    return a


def _atoms(model, unit="angstrom"):
    if unit not in ("angstrom", "bohr"):
        raise ValueError("position_unit must be angstrom or bohr")
    result = {}
    for atom in model.get("atoms", []):
        atom_id = atom.get("id")
        if not isinstance(atom_id, int) or atom_id in result:
            raise ValueError("Atom IDs must be unique integers")
        element = atom.get("element")
        if element not in Z:
            raise ValueError(f"Unsupported QM element {element}")
        position = _vector(atom.get("position"), "atom position")
        if unit == "bohr":
            position = position * BOHR_ANGSTROM
        result[atom_id] = dict(atom, position=position, element=element)
    if not result:
        raise ValueError("Empty molecule")
    return result


def _capped_graph(model, atoms):
    neighbors = {i: [] for i in atoms}
    seen = set()
    for bond in model.get("bonds", []):
        if any(bond.get("image", [0, 0, 0])):
            raise ValueError("QM requires finite, unwrapped, H-capped geometry; image bonds remain")
        a, b = bond["a"], bond["b"]
        if a not in atoms or b not in atoms or a == b:
            raise ValueError("Invalid QM bond endpoint")
        key = tuple(sorted((a, b)))
        if key in seen:
            raise ValueError("Duplicate QM bond")
        seen.add(key)
        order = float(bond["order"])
        if order not in (1, 1.5, 2, 3):
            raise ValueError("Invalid QM bond order")
        neighbors[a].append(order)
        neighbors[b].append(order)
    for i, atom in atoms.items():
        orders = neighbors[i]
        valence = sum(orders)
        element = atom["element"]
        formal = atom.get("formal_charge", 0)
        if element == "C":
            aromatic = len(orders) == 3 and all(o in (1, 1.5) for o in orders) and 1.5 in orders
            valid = math.isclose(valence, 4) or aromatic
        elif element == "H":
            valid = math.isclose(valence, 1) or (len(atoms) == 1 and not orders)
        elif element == "O":
            valid = math.isclose(valence, 2 + formal) or (
                formal == -1 and len(orders) == 1 and orders[0] == 1.5
            )
        elif element == "N":
            valid = math.isclose(valence, 3 + formal)
        elif element in ("F", "Cl", "Br", "I"):
            valid = math.isclose(valence, 1)
        else:
            raise ValueError(
                f"Automatic capping validation unavailable for {element}; supply a supported finite fragment"
            )
        if not valid:
            raise ValueError(
                f"Atom {i} ({element}) is open/unsupported valence {valence}; explicitly cap the fragment"
            )


def prepare_qm(
    model: dict,
    output_dir,
    *,
    total_charge=None,
    multiplicity=1,
    method="HF",
    basis="6-31G*",
    nproc=4,
    memory="4GB",
) -> dict:
    """Prepare explicit-geometry Gaussian MK single-point and exploratory RESP.

    Accepts an already finite, chemically capped model. Does not cut a periodic
    sheet or guess broken-bond hydrogens. The default reference is a declared
    comparison protocol, not a reconstruction of CURRENT's unknown QM method.
    """
    out = _directory(output_dir)
    if any(model.get("pbc", [False, False, False])):
        raise ValueError(
            "Periodic sheets cannot be Gaussian finite references; build an H-capped fragment"
        )
    atoms = _atoms(model)
    _capped_graph(model, atoms)
    if total_charge is None:
        if any("formal_charge" not in a for a in atoms.values()):
            raise ValueError("Set integer total_charge explicitly or supply all formal charges")
        total_charge = sum(a["formal_charge"] for a in atoms.values())
    if (
        isinstance(total_charge, bool)
        or not isinstance(total_charge, (int, float))
        or not math.isfinite(total_charge)
        or int(total_charge) != total_charge
    ):
        raise ValueError("total_charge must be a finite integer")
    total_charge = int(total_charge)
    electrons = sum(Z[a["element"]] for a in atoms.values()) - total_charge
    if (
        not isinstance(multiplicity, int)
        or multiplicity < 1
        or electrons < 0
        or electrons < multiplicity - 1
        or (electrons - multiplicity + 1) % 2
    ):
        raise ValueError("Charge, electron count and spin multiplicity are inconsistent")
    if not isinstance(nproc, int) or nproc < 1:
        raise ValueError("nproc must be positive")
    if not re.fullmatch(r"[1-9][0-9]*(?:MB|GB)", memory, re.IGNORECASE):
        raise ValueError("memory must have a Gaussian unit, e.g. 4GB")
    if not re.fullmatch(r"[A-Za-z0-9()+*,.\-]+", method) or not re.fullmatch(
        r"[A-Za-z0-9()+*,.\-]+", basis
    ):
        raise ValueError("method/basis must be Gaussian tokens, without route directives")
    positions = [a["position"] for a in atoms.values()]
    for i, pos in enumerate(positions):
        if any(np.linalg.norm(pos - other) < 0.2 for other in positions[:i]):
            raise ValueError("Overlapping atom coordinates")
    lines = [
        "%chk=reference.chk",
        f"%mem={memory}",
        f"%nprocshared={nproc}",
        f"#P {method}/{basis} SCF=Tight NoSymm Pop=MK IOp(6/33=2,6/41=10,6/42=17)",
        "",
        "GO common-geometry reference; no optimization",
        "",
        f"{total_charge} {multiplicity}",
    ]
    lines += [
        f"{a['element']:2s} " + " ".join(f"{v:.10f}" for v in a["position"]) for a in atoms.values()
    ]
    resp = [
        "Exploratory single-stage RESP; equivalence constraints unset",
        " &cntrl",
        " nmol=1, ihfree=1, ioutopt=1, iqopt=1, qwt=0.0005,",
        " &end",
        " 1.0",
        "GO finite reference",
        f"{total_charge:5d}{len(atoms):5d}",
    ]
    resp += [f"{Z[a['element']]:5d}{0:5d}" for a in atoms.values()]
    manifest = {
        "schema_version": 1,
        "executed": False,
        "status": "inputs_prepared_only",
        "total_charge": total_charge,
        "multiplicity": multiplicity,
        "method": method,
        "basis": basis,
        "geometry": "single_point_fixed",
        "esp_scheme": "Gaussian MK",
        "nproc": nproc,
        "memory": memory,
        "atom_order_ids": list(atoms),
        "units": {"position": "angstrom", "charge": "e"},
        "limitations": [
            "Not the recovered CURRENT protocol.",
            "No QM or RESP result exists yet.",
            "Single-stage RESP is exploratory; define equivalences, conformations and two-stage constraints before production fitting.",
            "Geometry must be chemically reviewed; valence checks do not establish a stable electronic state.",
        ],
    }
    out.mkdir(parents=True)
    (out / "reference.gjf").write_text("\n".join(lines) + "\n\n")
    (out / "resp-stage1.in").write_text("\n".join(resp) + "\n\n")
    script = out / "prepare_resp.sh"
    script.write_text(
        '#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "$0")"\n'
        "for target in reference.esp resp-stage1.out resp-stage1.qout resp-stage1.punch resp-stage1.esout; do\n"
        '  if [[ -e "$target" ]]; then echo "Refusing overwrite: $target" >&2; exit 1; fi\ndone\n'
        "command -v espgen >/dev/null\ncommand -v resp >/dev/null\n"
        "test -s reference.log\n"
        "espgen -i reference.log -o reference.esp\n"
        "resp -i resp-stage1.in -o resp-stage1.out -e reference.esp -t resp-stage1.qout "
        "-p resp-stage1.punch -s resp-stage1.esout\n"
    )
    script.chmod(0o755)
    (out / "README.txt").write_text(
        "Prepared only. Review caps, protonation, integer charge and spin first.\n"
        "Run Gaussian externally, from this directory: g16 < reference.gjf > reference.log\n"
        "Inspect Gaussian normal termination, SCF convergence and ESP generation before fitting.\n"
        "Then ./prepare_resp.sh runs external AmberTools espgen/resp; it refuses output overwrites.\n"
        "Single-stage RESP has no equivalence constraints and is exploratory, not a validated final fit.\n"
        "For production RESP, define chemical equivalence and the stage-2 fixed/free atom mapping,\n"
        "consistent with the selected Amber protocol. Keep fitted hydrogen constraints documented.\n"
        "Do not optimize separately per charge candidate: use one common geometry for each comparison.\n"
        "Repeat on several reviewed conformations; do not average scores across different chemistry.\n"
        "Import espgen output with go_builder.qm.espgen_to_reference; potential units are Hartree/e.\n"
        "Use compare_esp with candidate atoms containing identical IDs/elements/Angstrom positions\n"
        "and source charges whose sum equals the reference integer total charge. No charge repair occurs.\n"
    )
    (out / "qm_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def compare_esp(
    reference_path, candidates: dict, output_dir, *, geometry_tolerance=1e-6, charge_tolerance=1e-6
) -> dict:
    """Score point-charge ESP on one common grid, same geometry and net charge.

    Reference schema: atoms[{id,element,position}], total_charge, position_unit,
    potential_unit='hartree/e', grid[{position,potential}], optional dipole,
    dipole_unit ('debye' or 'e*bohr' or 'e*angstrom'), dipole_origin (same
    position_unit). Candidate values are model dictionaries with atom charges.
    RMSE is unweighted; the reference grid therefore defines sampling/weighting.
    """
    out = _directory(output_dir)
    if not candidates:
        raise ValueError("Provide at least one matched charge candidate")
    if (
        not math.isfinite(geometry_tolerance)
        or geometry_tolerance < 0
        or not math.isfinite(charge_tolerance)
        or charge_tolerance < 0
    ):
        raise ValueError("Tolerances must be finite and nonnegative")
    ref_path = Path(reference_path)
    ref = json.loads(ref_path.read_text())
    if any(ref.get("pbc", [False, False, False])) or any(
        any(b.get("image", [0, 0, 0])) for b in ref.get("bonds", [])
    ):
        raise ValueError("Reference must describe a finite nonperiodic geometry")
    if ref.get("potential_unit") != "hartree/e":
        raise ValueError("Reference potential_unit must explicitly be hartree/e")
    unit = ref.get("position_unit")
    atoms = _atoms(ref, unit)
    ids = list(atoms)
    positions = np.array([atoms[i]["position"] for i in ids])
    total_charge = float(ref["total_charge"])
    if not math.isfinite(total_charge) or not math.isclose(
        total_charge, round(total_charge), abs_tol=1e-8
    ):
        raise ValueError("Reference total_charge must be a finite integer")
    origin = _vector(ref.get("dipole_origin", [0, 0, 0]), "dipole origin")
    if unit == "bohr":
        origin *= BOHR_ANGSTROM
    points, potentials = [], []
    for point in ref.get("grid", []):
        pos = _vector(point["position"], "grid position")
        points.append(pos * BOHR_ANGSTROM if unit == "bohr" else pos)
        potential = float(point["potential"])
        if not math.isfinite(potential):
            raise ValueError("Nonfinite reference ESP")
        potentials.append(potential)
    if not points:
        raise ValueError("Reference ESP grid is empty")
    potentials = np.asarray(potentials)
    ref_dipole = None
    if "dipole" in ref:
        factors = {
            "debye": 1.0,
            "e*angstrom": E_ANGSTROM_DEBYE,
            "e*bohr": E_ANGSTROM_DEBYE * BOHR_ANGSTROM,
        }
        if ref.get("dipole_unit") not in factors:
            raise ValueError("Explicit supported dipole_unit required")
        ref_dipole = _vector(ref["dipole"], "reference dipole") * factors[ref["dipole_unit"]]
    scores, rows = {}, []
    for name, candidate in candidates.items():
        if any(candidate.get("pbc", [False, False, False])) or any(
            any(b.get("image", [0, 0, 0])) for b in candidate.get("bonds", [])
        ):
            raise ValueError("ESP comparison requires finite nonperiodic candidates")
        cand = _atoms(candidate, candidate.get("position_unit", "angstrom"))
        if set(cand) != set(ids):
            raise ValueError(f"{name}: atom IDs differ from reference")
        charges = []
        for i in ids:
            if (
                cand[i]["element"] != atoms[i]["element"]
                or np.linalg.norm(cand[i]["position"] - atoms[i]["position"]) > geometry_tolerance
            ):
                raise ValueError(f"{name}: geometry/element mismatch at atom {i}")
            charge = float(cand[i]["charge"])
            if not math.isfinite(charge):
                raise ValueError(f"{name}: nonfinite charge")
            charges.append(charge)
        charges = np.asarray(charges)
        if abs(float(charges.sum()) - total_charge) > charge_tolerance:
            raise ValueError(
                f"{name}: total charge does not match reference; no automatic renormalization"
            )
        predicted = []
        # Chunk the grid to avoid constructing a large full distance tensor.
        for start in range(0, len(points), 1024):
            distances = np.linalg.norm(
                np.asarray(points[start : start + 1024])[:, None, :] - positions[None, :, :], axis=2
            )
            if np.any(distances < 1e-10):
                raise ValueError("ESP grid intersects an atom center")
            predicted.extend((BOHR_ANGSTROM / distances) @ charges)
        predicted = np.asarray(predicted)
        residuals = predicted - potentials
        rmse = float(np.sqrt(np.mean(residuals**2)))
        dipole = np.sum(charges[:, None] * (positions - origin), axis=0) * E_ANGSTROM_DEBYE
        denominator = float(np.mean(potentials**2))
        scores[str(name)] = {
            "charge_geometry_sha256": hashlib.sha256(
                json.dumps(
                    [
                        {
                            "id": i,
                            "element": cand[i]["element"],
                            "position_angstrom": cand[i]["position"].tolist(),
                            "charge_e": float(q),
                        }
                        for i, q in zip(ids, charges, strict=True)
                    ],
                    sort_keys=True,
                    allow_nan=False,
                ).encode()
            ).hexdigest(),
            "total_charge_e": float(charges.sum()),
            "rmse_hartree_per_e": rmse,
            "mae_hartree_per_e": float(np.mean(abs(residuals))),
            "max_abs_error_hartree_per_e": float(np.max(abs(residuals))),
            "relative_rmse": rmse / math.sqrt(denominator) if denominator else None,
            "rmse_kcal_mol_per_e": rmse * HARTREE_KCAL_MOL,
            "dipole_debye": dipole.tolist(),
            "dipole_error_debye": float(np.linalg.norm(dipole - ref_dipole))
            if ref_dipole is not None
            else None,
        }
        rows.extend(
            {
                "candidate": str(name),
                "point": i + 1,
                "reference_hartree_per_e": float(v),
                "predicted_hartree_per_e": float(p),
                "residual_hartree_per_e": float(e),
            }
            for i, (v, p, e) in enumerate(zip(potentials, predicted, residuals, strict=True))
        )
    result = {
        "schema_version": 1,
        "candidates": scores,
        "n_grid_points": len(points),
        "reference_sha256": hashlib.sha256(ref_path.read_bytes()).hexdigest(),
        "dipole_origin_angstrom": origin.tolist(),
        "geometry_tolerance_angstrom": geometry_tolerance,
        "charge_tolerance_e": charge_tolerance,
        "limits": [
            "Same-geometry finite-molecule point-charge ESP only; no polarization or periodic electrostatics.",
            "Charged-system dipoles depend on the declared origin.",
            "Grid/conformation selection controls the score; one benchmark cannot validate all GO chemistry.",
            "kcal/mol/e is potential converted to the energy of a unit positive test charge.",
        ],
    }
    out.mkdir(parents=True)
    (out / "esp_comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    import csv

    with (out / "esp_residuals.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return result


def espgen_to_reference(
    esp_path,
    model: dict,
    output_path,
    *,
    total_charge,
    geometry_tolerance=1e-5,
    dipole=None,
    dipole_unit="debye",
    dipole_origin=None,
) -> dict:
    """Import one-molecule Amber espgen/RESP ESP data, in bohr and atomic units.

    The model supplies ordered atom identities/elements; the ESP file supplies
    QM coordinates and potential grid. A coordinate mismatch rejects rotated,
    optimized, reordered or otherwise unmatched input rather than guessing a map.
    """
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(output_path)
    atoms = _atoms(model)
    if (
        not isinstance(total_charge, (int, float))
        or not math.isfinite(total_charge)
        or total_charge != int(total_charge)
    ):
        raise ValueError("Explicit integer reference total_charge required")
    raw_lines = [
        line
        for line in Path(esp_path).read_text().replace("D", "E").replace("d", "e").splitlines()
        if line.strip()
    ]
    lines = [line.split() for line in raw_lines]
    counts = []
    if raw_lines:
        # espgen uses adjacent fixed-width fields: five digits can touch.
        for width in (5, 7):
            try:
                counts.append((int(raw_lines[0][:width]), int(raw_lines[0][width : 2 * width])))
            except ValueError:
                pass
        try:
            counts.append(tuple(map(int, lines[0][:2])))
        except ValueError:
            pass
    matched = [
        pair
        for pair in counts
        if len(pair) == 2
        and pair[0] == len(atoms)
        and pair[1] > 0
        and len(lines) == 1 + pair[0] + pair[1]
    ]
    if not matched:
        raise ValueError("Amber ESP counts or model atom count mismatch")
    natoms, npoints = matched[0]
    ref_atoms = []
    for atom, fields in zip(atoms.values(), lines[1 : 1 + natoms], strict=True):
        if len(fields) != 3:
            raise ValueError(
                "Amber espgen atomic coordinate row must contain three bohr coordinates"
            )
        position = _vector(fields, "ESP atom coordinate")
        if np.linalg.norm(position * BOHR_ANGSTROM - atom["position"]) > geometry_tolerance:
            raise ValueError("ESP geometry/order differs from supplied model")
        ref_atoms.append(
            {"id": atom["id"], "element": atom["element"], "position": position.tolist()}
        )
    grid = []
    for fields in lines[1 + natoms :]:
        if len(fields) != 4:
            raise ValueError("Amber ESP grid row must contain potential followed by x,y,z")
        potential = float(fields[0])
        if not math.isfinite(potential):
            raise ValueError("Nonfinite ESP potential")
        grid.append(
            {
                "position": _vector(fields[1:], "ESP grid coordinate").tolist(),
                "potential": potential,
            }
        )
    result = {
        "schema_version": 1,
        "atoms": ref_atoms,
        "grid": grid,
        "total_charge": int(total_charge),
        "position_unit": "bohr",
        "potential_unit": "hartree/e",
        "source_esp_sha256": hashlib.sha256(Path(esp_path).read_bytes()).hexdigest(),
        "source_esp_path": str(Path(esp_path).resolve()),
    }
    if dipole is not None:
        if dipole_unit not in ("debye", "e*bohr", "e*angstrom"):
            raise ValueError("Unsupported reference dipole unit")
        result.update(dipole=_vector(dipole, "dipole").tolist(), dipole_unit=dipole_unit)
    if dipole_origin is not None:
        result["dipole_origin"] = (
            _vector(dipole_origin, "dipole origin Angstrom") / BOHR_ANGSTROM
        ).tolist()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
