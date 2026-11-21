"""Amber-first topology preparation and parameter-preserving ParmEd exchange.

Finite covalent graphs may use a periodic solvent box. Bonds carrying nonzero
image shifts are rejected: this exporter does not provide a validated periodic
bonded-force implementation for a covalently winding sheet.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

import numpy as np

from .ff import DEFAULT_REFERENCE, assign_current_fragments, normalize_libraries
from .ff_v3 import (
    RES_COLORS_EPX,
    assign_gaff2_v3,
    compute_l0_v2_charges,
    load_custom_charges,
    resolve_force_field_preset,
    verify_v3_bundle,
)

ENGINE_SOURCES = [
    "https://ambermd.org/FileFormats.php",
    "https://raw.githubusercontent.com/Amber-MD/AmberClassic/main/src/msander/ene.F90",
    "https://manual.gromacs.org/current/user-guide/mdp-options.html#periodic-molecules",
    "https://manual.gromacs.org/current/reference-manual/topologies/molecule-definition.html",
]


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def _periodic_edges(model):
    return any(any(b.get("image", [0, 0, 0])) for b in model["bonds"])


def _validate_graph(model):
    atoms = model.get("atoms", [])
    if not atoms:
        raise ValueError("Topology graph contains no atoms")
    ids = {a["id"] for a in atoms}
    if len(ids) != len(atoms):
        raise ValueError("Duplicate atom IDs")
    seen = set()
    for b in model["bonds"]:
        if b["a"] not in ids or b["b"] not in ids or b["a"] == b["b"]:
            raise ValueError("Invalid explicit bond endpoint")
        key = (min(b["a"], b["b"]), max(b["a"], b["b"]))
        if key in seen:
            raise ValueError("Repeated atom pair is unsupported for Amber topology")
        seen.add(key)
    for atom in atoms:
        if len(atom.get("position", [])) != 3 or not all(
            math.isfinite(float(x)) for x in atom["position"]
        ):
            raise ValueError("Finite explicit coordinates required")
        for key in ("name", "resname", "element"):
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", str(atom.get(key, ""))):
                raise ValueError(f"Unsafe or empty {key} in atom {atom['id']}")


def _residue_map(model):
    result = {}
    for atom in model["atoms"]:
        result.setdefault((atom["resid"], atom["resname"]), len(result) + 1)
    return result


def _write_amber_provenance(
    output_dir: Path,
    *,
    preset: str,
    frcmod_name: str,
    frcmod_sha256: str,
    assigned: dict,
    radii: str,
    topology_path: Path | None = None,
    coordinates_path: Path | None = None,
    topology_sha256: str | None = None,
    coordinates_sha256: str | None = None,
    junction_frcmod_name: str | None = None,
    junction_frcmod_sha256: str | None = None,
) -> Path:
    ff_meta = assigned.get("metadata", {}).get("force_field", {})
    files_prov = {
        "frcmod": frcmod_name,
        "frcmod_sha256": frcmod_sha256,
        "bundle_provenance_sha256": ff_meta.get("provenance_sha256"),
        "prmtop_sha256": (
            topology_sha256
            if topology_sha256 is not None
            else (_sha(topology_path) if topology_path else None)
        ),
        "inpcrd_sha256": (
            coordinates_sha256
            if coordinates_sha256 is not None
            else (_sha(coordinates_path) if coordinates_path else None)
        ),
    }
    if junction_frcmod_name is not None:
        files_prov["junction_frcmod"] = junction_frcmod_name
        files_prov["junction_frcmod_sha256"] = junction_frcmod_sha256
    prov = {
        "schema_version": 1,
        "preset": preset,
        "version": ff_meta.get("version", "1.0.0"),
        "tier": ff_meta.get("tier", "legacy"),
        "custom_combination": bool(ff_meta.get("custom_combination", False)),
        "bonded_model": ff_meta.get("bonded_model", ff_meta.get("name")),
        "charge_model": ff_meta.get("charge_model", "current_lib"),
        "charge_regime": ff_meta.get("charge_regime", "fixed_fragment_libraries"),
        "charge_info": ff_meta.get("charge_info"),
        "conventions": {
            "vdw_combining_rule": "Lorentz-Berthelot",
            "nonbonded_exclusions": "1-2 and 1-3 excluded",
            "scee_14": 1.2,
            "scnb_14": 2.0,
            "gb_radii_set": radii,
        },
        "files": files_prov,
        "atom_count": len(assigned["atoms"]),
        "bond_count": len(assigned["bonds"]),
        "net_charge_e": float(sum(a["charge"] for a in assigned["atoms"])),
    }
    prov_path = output_dir / "PROVENANCE.json"
    _json(prov_path, prov)
    return prov_path


def prepare_amber(
    model,
    output_dir,
    *,
    run=False,
    tleap="tleap",
    preset=None,
    charge_model="auto",
    charges_json=None,
    reference_dir=None,
    parameter_templates=None,
    radii="mbondi2",
    overwrite=False,
):
    """Prepare exact atom creation and explicit LEaP bonds, optionally execute.

    No bonds are reconstructed from coordinates. Existing output directories
    are protected even if overwrite=True (choose a fresh directory). Missing
    terms are fatal; inputs_prepared never means a usable parameterized system.
    """
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Choose a fresh Amber output directory: {output_dir}")
    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh output directory")
    _validate_graph(model)
    if _periodic_edges(model) or model.get("metadata", {}).get("covalently_periodic"):
        raise ValueError(
            "Amber export rejects periodic bond image/winding graphs; use a finite H-capped sheet. "
            "Periodic solvent around a finite sheet is a separate supported topology representation."
        )
    if radii not in ("bondi", "mbondi", "mbondi2", "mbondi3"):
        raise ValueError("Unsupported GB radii set")

    if preset is None:
        resolved_preset = "gaff2_current" if reference_dir is not None else "gaff2_v3"
    else:
        resolved_preset = resolve_force_field_preset(preset)
    if charges_json is not None and charge_model != "custom":
        raise ValueError(
            f"charges_json requires charge_model='custom' (got {charge_model!r})"
        )
    if parameter_templates is not None and (charge_model != "auto" or charges_json is not None):
        raise ValueError(
            "parameter_templates supplies dedicated charges and cannot be combined with explicit charge_model/charges_json overrides"
        )

    executable = shutil.which(str(tleap)) if run else None
    if run and executable is None:
        raise ValueError(f"tleap executable unavailable: {tleap}")

    if resolved_preset == "gaff2_v3":
        bundle = verify_v3_bundle(reference_dir)
        frcmod_name = "GO_gaff2_v3.frcmod"
        frcmod_sha = bundle["frcmod_sha256"]
        assigned = assign_gaff2_v3(
            model,
            reference_dir=reference_dir,
            verified_bundle=bundle,
            charge_model=charge_model,
            custom_charges=charges_json,
            parameter_templates=parameter_templates,
        )
        output_dir.mkdir(parents=True)
        v3_lib_out = output_dir / "v3_libraries"
        v3_lib_out.mkdir(parents=True)
        audit = {}
        ref_dir = Path(bundle["reference_dir"])
        manifest_files = bundle["provenance"]["files"]
        for lib_name, lib_raw in sorted(bundle["lib_bytes"].items()):
            dest = v3_lib_out / lib_name
            dest.write_bytes(lib_raw)
            res_key = lib_name.removesuffix("_gaff2_v3.lib")
            audit[res_key] = {
                "source": str(ref_dir / lib_name),
                "sha256": manifest_files[lib_name]["sha256"],
            }
        (output_dir / frcmod_name).write_bytes(bundle["frcmod_bytes"])
    else:
        cur_dir = Path(reference_dir or DEFAULT_REFERENCE)
        frcmod_name = "GO_gaff2.frcmod"
        source_frcmod = cur_dir / frcmod_name
        if not source_frcmod.is_file():
            raise ValueError("Current GO_gaff2.frcmod is missing")
        source_frcmod_bytes = source_frcmod.read_bytes()
        frcmod_sha = hashlib.sha256(source_frcmod_bytes).hexdigest()
        assigned = assign_current_fragments(
            model, cur_dir, parameter_templates=parameter_templates
        )
        if parameter_templates is not None:
            for atom in assigned["atoms"]:
                ameta = atom.setdefault("metadata", {})
                ameta["charge_source"] = "dedicated_parameter_templates"
                ameta.pop("v3_role", None)
                ameta.pop("epoxide_color", None)
            assigned.setdefault("metadata", {})["force_field"] = {
                "preset": "gaff2_current",
                "version": "1.0.0",
                "tier": "custom_combination",
                "custom_combination": True,
                "name": "GAFF2 + dedicated parameter_templates",
                "bonded_model": "dedicated_parameter_templates",
                "charge_model": "dedicated_parameter_templates",
                "charge_regime": "dedicated_parameter_templates",
                "charge_info": {
                    "charge_model": "dedicated_parameter_templates",
                    "regime": "dedicated_parameter_templates",
                    "provenance": parameter_templates.get("provenance"),
                },
            }
        else:
            ff_meta = assigned.setdefault("metadata", {}).setdefault("force_field", {})
            if charge_model in ("auto", "current", "current_lib"):
                q_info = {"charge_model": "current_lib", "regime": "fixed_fragment_libraries"}
                custom_combo = False
                q_model_name = "current_lib"
                if charge_model == "current_lib":
                    q_info["source_dir"] = str(cur_dir.resolve())
            elif charge_model in ("l0_v2", "L0_v2"):
                q_vec, q_info = compute_l0_v2_charges(assigned)
                for idx, atom in enumerate(assigned["atoms"]):
                    atom["charge"] = float(q_vec[idx])
                    atom.setdefault("metadata", {})["charge_source"] = "custom_combo:gaff2_current+L0_v2"
                custom_combo = True
                q_model_name = "L0_v2c_a030"
            elif charge_model == "custom":
                if charges_json is None:
                    raise ValueError("custom charges_json is required when charge_model='custom'")
                q_vec, q_info = load_custom_charges(charges_json, assigned)
                for idx, atom in enumerate(assigned["atoms"]):
                    atom["charge"] = float(q_vec[idx])
                    atom.setdefault("metadata", {})["charge_source"] = f"custom:{q_info['source']}"
                custom_combo = True
                q_model_name = "custom"
            else:
                raise ValueError(
                    f"Unsupported charge_model {charge_model!r} for gaff2_current"
                )
            ff_meta.update(
                {
                    "preset": "gaff2_current",
                    "version": "1.0.0",
                    "tier": "custom_combination" if custom_combo else "legacy",
                    "custom_combination": custom_combo,
                    "charge_model": q_model_name,
                    "charge_regime": q_info["regime"],
                    "charge_info": q_info,
                }
            )
        output_dir.mkdir(parents=True)
        audit = normalize_libraries(cur_dir, output_dir / "normalized_libraries")
        (output_dir / frcmod_name).write_bytes(source_frcmod_bytes)

    extra = []
    junction_frcmod_sha = None
    if parameter_templates:
        junction_bytes = Path(parameter_templates["frcmod"]).read_bytes()
        junction_frcmod_sha = hashlib.sha256(junction_bytes).hexdigest()
        (output_dir / "junction.frcmod").write_bytes(junction_bytes)
        extra = ["loadAmberParams junction.frcmod"]
    residues = _residue_map(assigned)
    assigned["atoms"].sort(key=lambda atom: residues[(atom["resid"], atom["resname"])])
    lines = [
        "source leaprc.gaff2",
        f"loadAmberParams {frcmod_name}",
        *extra,
        f"set default PBRadii {radii}",
        "mol = createUnit GO",
    ]
    for (_resid, name), index in residues.items():
        lines.extend([f"r{index} = createResidue {name}", f"add mol r{index}"])
    atom_map = {}
    used_names = {}
    for serial, atom in enumerate(assigned["atoms"], 1):
        index = residues[(atom["resid"], atom["resname"])]
        key = (index, atom["name"])
        if key in used_names:
            raise ValueError(f"Duplicate name in residue {key}")
        used_names[key] = True
        pos = " ".join(format(float(x), ".12g") for x in atom["position"])
        atom_map[atom["id"]] = f"mol.{index}.{atom['name']}"
        lines.extend(
            [
                f"a{serial} = createAtom {atom['name']} {atom['atom_type']} {atom['charge']:.12g}",
                f"set a{serial} element {atom['element']}",
                f"set a{serial} position {{ {pos} }}",
                f"add r{index} a{serial}",
            ]
        )
    for bond in assigned["bonds"]:
        order = {1: "S", 1.5: "A", 2: "D", 3: "T"}.get(bond.get("order", 1))
        if order is None:
            raise ValueError("Unsupported explicit bond order")
        lines.append(f'bond {atom_map[bond["a"]]} {atom_map[bond["b"]]} "{order}"')
    # Creating the finite solute does not automatically create a solvated system.
    # Retain requested box in the graph; users must choose solvation separately.
    lines.extend(
        ["check mol", "saveAmberParm mol go.prmtop go.inpcrd", "savePdb mol go.pdb", "quit"]
    )
    leap_input = output_dir / "tleap.in"
    leap_input.write_text("\n".join(lines) + "\n")
    _json(output_dir / "assigned_model.json", assigned)
    _json(output_dir / "library_audit.json", audit)
    prov_path = _write_amber_provenance(
        output_dir,
        preset=resolved_preset,
        frcmod_name=frcmod_name,
        frcmod_sha256=frcmod_sha,
        assigned=assigned,
        radii=radii,
        junction_frcmod_name="junction.frcmod" if parameter_templates else None,
        junction_frcmod_sha256=junction_frcmod_sha,
    )
    ff_meta = assigned.get("metadata", {}).get("force_field", {})
    result = {
        "status": "inputs_prepared",
        "preset": resolved_preset,
        "charge_model": ff_meta.get("charge_model"),
        "custom_combination": bool(ff_meta.get("custom_combination", False)),
        "parameterized": False,
        "leap_input": str(leap_input),
        "assigned_model": str(output_dir / "assigned_model.json"),
        "library_audit": str(output_dir / "library_audit.json"),
        "provenance": str(prov_path),
        "topology": None,
        "coordinates": None,
        "radii_set": radii,
        "atom_mapping": {str(i): a["id"] for i, a in enumerate(assigned["atoms"], 1)},
        "atom_count": len(assigned["atoms"]),
        "bond_count": len(assigned["bonds"]),
        "net_charge": sum(a["charge"] for a in assigned["atoms"]),
        "go_frcmod_sha256": frcmod_sha,
        "gbsa": gbsa_readiness(assigned),
        "periodic_policy": "finite solute; requested box is not solvation; nonzero bond images rejected",
        "engine_sources": ENGINE_SOURCES,
    }
    if junction_frcmod_sha is not None:
        result["junction_frcmod_sha256"] = junction_frcmod_sha
    if run:
        completed = subprocess.run(
            [executable, "-f", "tleap.in"],
            cwd=output_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        log = completed.stdout + "\n" + completed.stderr
        (output_dir / "tleap.log").write_text(log)
        result["returncode"] = completed.returncode
        result["status"] = "failed"
        _json(output_dir / "manifest.json", result)
        top = output_dir / "go.prmtop"
        crd = output_dir / "go.inpcrd"
        errors = re.search(
            r"FATAL|Could not find|does not have a type|Errors\s*=\s*[1-9]|Error!", log, re.I
        )
        if completed.returncode or errors or not top.is_file() or not crd.is_file():
            raise ValueError(
                f"LEaP parameterization failed (missing terms/junctions remain fatal); inspect {output_dir / 'tleap.log'}"
            )
        imported = import_amber(top, crd)
        actual = {frozenset((b["a"], b["b"])) for b in imported["bonds"]}
        serials = {a["id"]: i for i, a in enumerate(assigned["atoms"], 1)}
        expected = {frozenset((serials[b["a"]], serials[b["b"]])) for b in assigned["bonds"]}
        if len(imported["atoms"]) != len(assigned["atoms"]) or actual != expected:
            raise ValueError("LEaP output differs from authoritative atom/bond graph")
        if not np.allclose(
            [a["charge"] for a in imported["atoms"]],
            [a["charge"] for a in assigned["atoms"]],
            atol=1e-6,
        ):
            raise ValueError("LEaP altered assigned fragment charges")
        top_sha = imported["metadata"]["source_topology_sha256"]
        crd_sha = imported["metadata"]["source_coordinates_sha256"]
        _write_amber_provenance(
            output_dir,
            preset=resolved_preset,
            frcmod_name=frcmod_name,
            frcmod_sha256=frcmod_sha,
            assigned=assigned,
            radii=radii,
            topology_path=top,
            coordinates_path=crd,
            topology_sha256=top_sha,
            coordinates_sha256=crd_sha,
            junction_frcmod_name="junction.frcmod" if parameter_templates else None,
            junction_frcmod_sha256=junction_frcmod_sha,
        )
        result.update(
            status="parameterized",
            parameterized=True,
            topology=str(top),
            coordinates=str(crd),
            topology_sha256=top_sha,
            gbsa=gbsa_readiness(imported),
        )
    _json(output_dir / "manifest.json", result)
    return result


def _load_amber_with_digests(topology, coordinates=None):
    import parmed as pmd

    top_path = Path(topology).resolve()
    top_bytes = top_path.read_bytes()
    top_sha = hashlib.sha256(top_bytes).hexdigest()
    crd_bytes = None
    crd_sha = None
    crd_name = None
    if coordinates is not None:
        crd_path = Path(coordinates).resolve()
        crd_bytes = crd_path.read_bytes()
        crd_sha = hashlib.sha256(crd_bytes).hexdigest()
        crd_name = crd_path.name
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_root = Path(tmp_dir)
        tmp_top = tmp_root / (top_path.name or "input.prmtop")
        tmp_top.write_bytes(top_bytes)
        load_kwargs = {}
        if crd_bytes is not None:
            safe_crd_name = (
                crd_name
                if crd_name and crd_name != tmp_top.name
                else f"coords_{crd_name or 'input.inpcrd'}"
            )
            tmp_crd = tmp_root / safe_crd_name
            tmp_crd.write_bytes(crd_bytes)
            load_kwargs["xyz"] = str(tmp_crd)
        structure = pmd.load_file(str(tmp_top), **load_kwargs)
    if not isinstance(structure, pmd.amber.AmberParm):
        raise ValueError("An existing fully parameterized Amber prmtop is required")
    if structure.coordinates is None:
        raise ValueError("Amber coordinates are required; coordinates are never fabricated")
    _check_parameterized(structure)
    return structure, top_sha, crd_sha


def _load_amber(topology, coordinates=None):
    return _load_amber_with_digests(topology, coordinates)[0]


def _check_parameterized(structure):
    for bond in structure.bonds:
        if (
            bond.type is None
            or not math.isfinite(float(bond.type.k))
            or not math.isfinite(float(bond.type.req))
        ):
            raise ValueError("Unparameterized or nonfinite bonds in topology")
    for angle in structure.angles:
        if (
            angle.type is None
            or not math.isfinite(float(angle.type.k))
            or not math.isfinite(float(angle.type.theteq))
        ):
            raise ValueError("Unparameterized or nonfinite angles in topology")
    for dihedral in structure.dihedrals:
        types = dihedral.type if isinstance(dihedral.type, list) else [dihedral.type]
        if not types:
            raise ValueError("Unparameterized dihedrals in topology")
        for dt in types:
            if dt is None:
                raise ValueError("Unparameterized dihedrals in topology")
            phi_k = float(dt.phi_k)
            per = float(dt.per)
            phase = float(dt.phase)
            scee = float(getattr(dt, "scee", 0.0))
            scnb = float(getattr(dt, "scnb", 0.0))
            if not (
                math.isfinite(phi_k)
                and math.isfinite(per)
                and math.isfinite(phase)
                and math.isfinite(scee)
                and math.isfinite(scnb)
            ):
                raise ValueError("Nonfinite dihedral parameter in topology")
            if per < 0.0 or not per.is_integer():
                raise ValueError(
                    f"Negative or non-integer dihedral periodicity {dt.per!r} in topology"
                )
            if per == 0.0 and abs(phi_k) > 0.0:
                raise ValueError(
                    f"Zero dihedral periodicity with nonzero force constant ({phi_k}) in topology"
                )
            if scee < 0.0 or scnb < 0.0:
                raise ValueError(
                    f"Negative 1-4 scaling factor in topology (scee={scee}, scnb={scnb})"
                )
    for improper in structure.impropers:
        if (
            improper.type is None
            or not math.isfinite(float(improper.type.psi_k))
            or not math.isfinite(float(improper.type.psi_eq))
        ):
            raise ValueError("Unparameterized or nonfinite impropers in topology")
    parm_data = getattr(structure, "parm_data", None) or {}
    for flag in ("DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PHASE"):
        for raw_val in parm_data.get(flag) or []:
            if not math.isfinite(float(raw_val)):
                raise ValueError(f"Nonfinite {flag} parameter in topology")
    raw_pers = parm_data.get("DIHEDRAL_PERIODICITY") or []
    raw_fks = parm_data.get("DIHEDRAL_FORCE_CONSTANT") or []
    for idx, raw_per in enumerate(raw_pers):
        fper = float(raw_per)
        if not math.isfinite(fper) or fper < 0.0 or not fper.is_integer():
            raise ValueError(
                f"Negative, non-integer, or nonfinite dihedral periodicity {raw_per!r} in topology"
            )
        if fper == 0.0 and idx < len(raw_fks) and abs(float(raw_fks[idx])) > 0.0:
            raise ValueError(
                f"Zero dihedral periodicity with nonzero force constant ({raw_fks[idx]}) in topology"
            )
    for flag in ("SCEE_SCALE_FACTOR", "SCNB_SCALE_FACTOR"):
        for raw_sc in parm_data.get(flag) or []:
            fsc = float(raw_sc)
            if not math.isfinite(fsc) or fsc < 0.0:
                raise ValueError(
                    f"Negative or nonfinite 1-4 scaling factor {raw_sc!r} in {flag}"
                )
    if any(
        not math.isfinite(a.charge)
        or not math.isfinite(a.mass)
        or a.mass <= 0
        or not math.isfinite(a.epsilon)
        or not math.isfinite(a.rmin)
        for a in structure.atoms
    ):
        raise ValueError("Missing or invalid atom charge/mass/Lennard-Jones parameters")
    if not np.isfinite(structure.coordinates).all():
        raise ValueError("Nonfinite imported coordinates")
    if structure.box is not None:
        box_arr = np.asarray(structure.box, dtype=float)
        if (
            box_arr.shape[0] < 6
            or not np.isfinite(box_arr[:6]).all()
            or np.any(box_arr[:3] <= 0.0)
            or np.any((box_arr[3:6] <= 0.0) | (box_arr[3:6] >= 180.0))
        ):
            raise ValueError(
                "Invalid nonfinite or nonpositive periodic box dimensions in topology"
            )
        if getattr(structure, "box_vectors", None) is not None:
            import parmed as pmd

            vecs = np.asarray(
                structure.box_vectors.value_in_unit(pmd.unit.nanometers), dtype=float
            )
            if (
                vecs.shape != (3, 3)
                or not np.isfinite(vecs).all()
                or abs(float(np.linalg.det(vecs))) <= 1e-12
            ):
                raise ValueError("Invalid or degenerate periodic box vectors in topology")


def annotate_imported_epoxides(model):
    """Annotate intact basal epoxide motifs (EPX/EP2/EP3/EP4) from explicit bonds.

    This recognizes known named basal epoxide templates and their complete local
    graph. It never assigns bonds from distances, guesses arbitrary oxygen
    chemistry, or treats edge epoxides as supported basal graft sites.
    """
    model = deepcopy(model)
    by_id = {a["id"]: a for a in model["atoms"]}
    neighbors = {i: set() for i in by_id}
    for bond in model["bonds"]:
        neighbors[bond["a"]].add(bond["b"])
        neighbors[bond["b"]].add(bond["a"])
    grouped = {}
    for atom in model["atoms"]:
        if atom["resname"] in RES_COLORS_EPX:
            grouped.setdefault((atom["resid"], atom["resname"]), []).append(atom)
    metadata = model.setdefault("metadata", {})
    existing = list(metadata.get("groups", []))
    represented = {frozenset(g.get("site_ids", [])) for g in existing}
    next_group = max((g["id"] for g in existing), default=0) + 1
    for _, atoms in sorted(grouped.items()):
        carbons = sorted(a["id"] for a in atoms if a["element"] == "C")
        oxygens = [a["id"] for a in atoms if a["element"] == "O"]
        if len(carbons) != 2 or len(oxygens) != 1:
            continue
        a, b = carbons
        oxygen = oxygens[0]
        if frozenset(carbons) in represented:
            continue
        if neighbors[oxygen] != set(carbons) or b not in neighbors[a]:
            continue
        if any(
            len(neighbors[c]) != 4 or sum(by_id[n]["element"] == "C" for n in neighbors[c]) != 3
            for c in carbons
        ):
            continue
        existing.append(
            {
                "id": next_group,
                "type": "epoxide",
                "site_ids": carbons,
                "atom_ids": [*carbons, oxygen],
                "oxygen_id": oxygen,
                "face": 1,
                "annotation_source": "named_EPX_and_explicit_C_C_O_triangle",
            }
        )
        next_group += 1
    metadata["groups"] = existing
    metadata.setdefault("events", [])
    return model


def import_amber(topology, coordinates=None, **options):
    """Read explicit Amber bonds and parameters via ParmEd; never guess bonds.

    Amber prmtop stores no bond-image shifts: import preserves its finite,
    unwrapped covalent interpretation. Supplied wrapped/broken coordinates need
    separate trajectory processing, never graph reconstruction here.
    """
    if options:
        raise ValueError(f"Unknown import options: {sorted(options)}")
    import parmed as pmd

    structure, top_sha, crd_sha = _load_amber_with_digests(topology, coordinates)
    atoms = []
    xyz = structure.coordinates
    for i, atom in enumerate(structure.atoms):
        element = pmd.periodic_table.Element[atom.atomic_number] if atom.atomic_number > 0 else "EP"
        atoms.append(
            {
                "id": i + 1,
                "name": atom.name,
                "element": element,
                "resname": atom.residue.name,
                "resid": atom.residue.idx + 1,
                "position": [float(v) for v in xyz[i]],
                "atom_type": atom.type,
                "charge": float(atom.charge),
                "formal_charge": int(atom.formal_charge or 0),
                "gb_radius": float(atom.solvent_radius),
                "gb_screen": float(atom.screen),
                "metadata": {
                    "amber_atom_index": i + 1,
                    "mass": float(atom.mass),
                    "epsilon": float(atom.epsilon),
                    "rmin": float(atom.rmin),
                },
            }
        )
    bonds = [
        {
            "a": min(b.atom1.idx, b.atom2.idx) + 1,
            "b": max(b.atom1.idx, b.atom2.idx) + 1,
            "order": float(b.order or 1),
            "image": [0, 0, 0],
            "metadata": {"k": float(b.type.k), "req": float(b.type.req)},
        }
        for b in structure.bonds
    ]
    box = structure.box
    cell = (
        None
        if box is None
        else np.asarray(structure.box_vectors.value_in_unit(pmd.unit.angstroms)).tolist()
    )
    stretched = any(
        np.linalg.norm(xyz[b.atom1.idx] - xyz[b.atom2.idx]) > max(3.0, 3 * b.type.req)
        for b in structure.bonds
    )
    result = {
        "schema_version": 1,
        "name": Path(topology).stem,
        "atoms": atoms,
        "bonds": bonds,
        "cell": cell,
        "pbc": [box is not None] * 3,
        "metadata": {
            "representation": "amber_import",
            "bond_orders_unknown": True,
            "formal_charges_unknown": True,
            "parameterized": True,
            "source_topology": str(Path(topology).resolve()),
            "source_topology_sha256": top_sha,
            "source_coordinates": str(Path(coordinates).resolve()) if coordinates else None,
            "source_coordinates_sha256": crd_sha,
            "radii_set": " ".join(structure.parm_data.get("RADIUS_SET", [])),
            "bond_image_policy": "Amber endpoints imported verbatim; no bond images in prmtop; no distance inference",
            "wrapped_or_stretched_bonds": stretched,
        },
    }
    return annotate_imported_epoxides(result)


def _verify_parameter_roundtrip(original, converted, *, top_path=None):
    """Verify numerical atoms, valence terms and 1-4 scaling after conversion."""
    import parmed as pmd

    is_gmx = isinstance(converted, pmd.gromacs.GromacsTopologyFile)

    def atom_values(structure):
        return [(a.charge, a.mass, a.epsilon, a.rmin) for a in structure.atoms]

    def canonical(indices):
        return min(tuple(indices), tuple(reversed(indices)))

    def valence(structure, field):
        result = []
        for term in getattr(structure, field):
            if field == "bonds":
                ids = canonical([term.atom1.idx, term.atom2.idx])
                values = (term.type.k, term.type.req)
            elif field == "angles":
                ids = canonical([term.atom1.idx, term.atom2.idx, term.atom3.idx])
                values = (term.type.k, term.type.theteq)
            elif field == "dihedrals":
                ids = canonical([term.atom1.idx, term.atom2.idx, term.atom3.idx, term.atom4.idx])
                types = term.type if isinstance(term.type, list) else [term.type]
                result.extend(
                    (
                        (*ids, int(term.improper)),
                        (t.phi_k, t.per, t.phase)
                        if is_gmx
                        else (
                            t.phi_k,
                            t.per,
                            t.phase,
                            0.0 if (term.improper or term.ignore_end) else t.scee,
                            0.0 if (term.improper or term.ignore_end) else t.scnb,
                        ),
                    )
                    for t in types
                )
                continue
            elif field == "impropers":
                ids = canonical([term.atom1.idx, term.atom2.idx, term.atom3.idx, term.atom4.idx])
                values = (term.type.psi_k, term.type.psi_eq)
            elif field == "adjusts":
                ids = canonical([term.atom1.idx, term.atom2.idx])
                values = (term.type.rmin, term.type.epsilon, term.type.chgscale)
            result.append((ids, values))
        return sorted(result)

    if len(original.atoms) != len(converted.atoms) or not np.allclose(
        atom_values(original), atom_values(converted), atol=2e-6, rtol=2e-6
    ):
        raise ValueError("GROMACS round-trip atom charge/mass/Lennard-Jones mismatch")
    if bool(getattr(converted, "unknown_functional", False)):
        raise ValueError("GROMACS round-trip encountered unknown_functional in topology")
    for extra_field in (
        "rb_torsions",
        "urey_bradleys",
        "trigonal_angles",
        "out_of_plane_bends",
        "pi_torsions",
        "stretch_bends",
        "torsion_torsions",
        "chiral_frames",
        "multipole_frames",
    ):
        if len(getattr(original, extra_field, ())) != len(getattr(converted, extra_field, ())):
            raise ValueError(f"GROMACS round-trip {extra_field} count mismatch")
    if len(getattr(original, "cmaps", ())) != len(getattr(converted, "cmaps", ())):
        raise ValueError("GROMACS round-trip CMAP torsion count mismatch")
    if (hasattr(original, "has_NBFIX") and original.has_NBFIX()) != (
        hasattr(converted, "has_NBFIX") and converted.has_NBFIX()
    ):
        raise ValueError("GROMACS round-trip NBFIX pair override mismatch")
    summary = {"atoms": len(original.atoms)}
    for field in ("bonds", "angles", "dihedrals", "impropers"):
        before, after = valence(original, field), valence(converted, field)
        if [x[0] for x in before] != [x[0] for x in after] or (
            before
            and not np.allclose([x[1] for x in before], [x[1] for x in after], atol=2e-6, rtol=2e-6)
        ):
            raise ValueError(f"GROMACS round-trip {field} parameter/scaling mismatch")
        summary["torsion_terms" if field == "dihedrals" else field] = len(before)

    orig_eps = sum(1 for a in original.atoms if isinstance(a, pmd.ExtraPoint))
    conv_eps = sum(1 for a in converted.atoms if isinstance(a, pmd.ExtraPoint))
    if orig_eps != conv_eps:
        raise ValueError("GROMACS round-trip extra-point count mismatch")

    if is_gmx:
        orig_12_13, orig_pairs_14, orig_extra_excl, orig_unexcl_14, orig_zero_excl, _orig_g = (
            _extract_amber_exclusions(original)
        )
        if orig_unexcl_14:
            raise ValueError(
                f"GROMACS round-trip failed: source Amber topology leaves {len(orig_unexcl_14)} 1-4 pair(s) unexcluded"
            )
        if top_path is not None:
            _verify_gromacs_topology_sections(Path(top_path), converted=converted)
        else:
            for itp in getattr(converted, "itps", ()) or ():
                itp_p = Path(itp)
                if itp_p.is_absolute() and itp_p.is_file():
                    _verify_gromacs_topology_sections(itp_p, converted=converted)
        orig_comb = str(getattr(original, "combining_rule", "lorentz") or "lorentz").lower()
        conv_comb = str(getattr(converted, "combining_rule", "lorentz") or "lorentz").lower()
        if orig_comb != conv_comb:
            raise ValueError(
                f"GROMACS round-trip combining_rule mismatch ({orig_comb!r} vs {conv_comb!r})"
            )
        expected_comb_rule = (
            2 if orig_comb == "lorentz" else 3 if orig_comb == "geometric" else None
        )
        if expected_comb_rule is None:
            raise ValueError(f"GROMACS round-trip unsupported combining_rule={orig_comb!r}")
        defaults = getattr(converted, "defaults", None)
        if defaults is not None:
            nbfunc = getattr(defaults, "nbfunc", None)
            if nbfunc is not None and int(nbfunc) != 1:
                raise ValueError(
                    f"GROMACS round-trip [ defaults ] nbfunc={nbfunc} mismatch (expected 1)"
                )
            comb_rule = getattr(defaults, "comb_rule", None)
            if comb_rule is not None and int(comb_rule) != expected_comb_rule:
                raise ValueError(
                    f"GROMACS round-trip [ defaults ] comb_rule={comb_rule} mismatch "
                    f"(expected {expected_comb_rule})"
                )
            gen_pairs = str(getattr(defaults, "gen_pairs", "yes")).strip().lower()
            if gen_pairs not in {"yes", "no"}:
                raise ValueError(
                    f"GROMACS round-trip invalid [ defaults ] gen_pairs={getattr(defaults, 'gen_pairs', None)!r}"
                )
            for fudge_attr in ("fudgeLJ", "fudgeQQ"):
                fudge_val = getattr(defaults, fudge_attr, None)
                if (
                    fudge_val is None
                    or not math.isfinite(float(fudge_val))
                    or float(fudge_val) < 0.0
                ):
                    raise ValueError(
                        f"GROMACS round-trip invalid [ defaults ] {fudge_attr}={fudge_val!r}"
                    )
        expected_adjusts = {
            pair: (
                float(original.atoms[pair[0]].rmin) + float(original.atoms[pair[1]].rmin),
                math.sqrt(
                    float(original.atoms[pair[0]].epsilon)
                    * float(original.atoms[pair[1]].epsilon)
                )
                / float(scnb),
                1.0 / float(scee),
            )
            for pair, (scee, scnb) in orig_pairs_14.items()
        }
        converted_adjusts: dict[tuple[int, int], tuple[float, float, float]] = {}
        for adj in converted.adjusts:
            pair = tuple(sorted((int(adj.atom1.idx), int(adj.atom2.idx))))
            if pair in converted_adjusts:
                raise ValueError(f"GROMACS round-trip duplicate 1-4 adjust pair {pair}")
            converted_adjusts[pair] = (
                float(adj.type.rmin),
                float(adj.type.epsilon),
                float(adj.type.chgscale),
            )
        if sorted(converted_adjusts) != sorted(expected_adjusts):
            raise ValueError(
                f"GROMACS round-trip 1-4 pair identity mismatch "
                f"(expected {len(expected_adjusts)} pairs, got {len(converted_adjusts)})"
            )
        if expected_adjusts:
            exp_arr = [expected_adjusts[k] for k in sorted(expected_adjusts)]
            got_arr = [converted_adjusts[k] for k in sorted(converted_adjusts)]
            if not np.allclose(exp_arr, got_arr, atol=2e-6, rtol=2e-6):
                raise ValueError("GROMACS round-trip 1-4 pair LJ/electrostatic parameter mismatch")
        conv_explicit_excl = _molecule_instance_explicit_exclusions(converted)
        if conv_explicit_excl != orig_extra_excl:
            raise ValueError("GROMACS round-trip explicit nonbonded exclusion mismatch")
        mol_nrexcls: set[int] = set()
        for mol_entry in (getattr(converted, "molecules", None) or {}).values():
            if isinstance(mol_entry, tuple) and len(mol_entry) >= 2:
                mol_obj = mol_entry[0]
                raw_nrexcl = getattr(mol_obj, "nrexcl", None)
                if len(mol_entry) == 2:
                    tuple_nrexcl = int(mol_entry[1])
                    if raw_nrexcl is not None and int(raw_nrexcl) != tuple_nrexcl:
                        raise ValueError(
                            f"GROMACS round-trip inconsistent moleculetype nrexcl ({raw_nrexcl} vs {tuple_nrexcl})"
                        )
                    mol_nrexcls.add(tuple_nrexcl)
                else:
                    mol_nrexcls.add(int(raw_nrexcl) if raw_nrexcl is not None else int(mol_entry[2]))
        if not mol_nrexcls:
            mol_nrexcls = {int(getattr(converted, "nrexcl", 3))}
        if len(mol_nrexcls) != 1:
            raise ValueError(
                f"GROMACS round-trip mixed nrexcl across moleculetypes: {sorted(mol_nrexcls)}"
            )
        nrexcl = next(iter(mol_nrexcls))
        if nrexcl < 0:
            raise ValueError(f"GROMACS round-trip invalid nrexcl={nrexcl}")
        conv_g = _bond_graph_pairs_upto(converted, nrexcl)
        conv_graph_excl = set().union(*conv_g.values()) if conv_g else set()
        unexcluded_adjusts = set(converted_adjusts) - (conv_graph_excl | conv_explicit_excl)
        if unexcluded_adjusts:
            raise ValueError(
                f"GROMACS round-trip nrexcl={nrexcl} leaves {len(unexcluded_adjusts)} 1-4 adjust pair(s) unexcluded"
            )
        conv_zero_excl = (conv_graph_excl | conv_explicit_excl) - set(converted_adjusts)
        if conv_zero_excl != orig_zero_excl:
            raise ValueError(
                f"GROMACS round-trip effective nonbonded exclusion mismatch "
                f"(expected {len(orig_zero_excl)} zero-interaction pairs, got {len(conv_zero_excl)})"
            )
        summary["pairs_14"] = len(expected_adjusts)
    orig_box_nm = _structure_box_vectors_nm(original)
    if orig_box_nm is not None:
        conv_box_nm = _structure_box_vectors_nm(converted)
        if conv_box_nm is None:
            raise ValueError(
                "GROMACS round-trip periodic box mismatch: missing periodic box in exported GRO"
            )
        orig_red = np.asarray(_reduce_openmm_box_vectors(orig_box_nm), dtype=float)
        conv_red = np.asarray(_reduce_openmm_box_vectors(conv_box_nm), dtype=float)
        if abs(float(np.linalg.det(orig_red))) < 1e-9 or abs(float(np.linalg.det(conv_red))) < 1e-9:
            raise ValueError("GROMACS round-trip periodic box mismatch: degenerate cell vectors")
        direct_ok = np.allclose(orig_red, conv_red, atol=1e-3, rtol=1e-5)
        if not direct_ok:
            m_mat = conv_red @ np.linalg.inv(orig_red)
            m_int = np.rint(m_mat)
            orig_metric = orig_red @ orig_red.T
            conv_metric = conv_red @ conv_red.T
            unimodular_ok = (
                np.allclose(m_mat, m_int, atol=1e-3, rtol=1e-4)
                and abs(abs(float(np.linalg.det(m_int))) - 1.0) < 1e-3
                and np.allclose(m_int @ orig_metric @ m_int.T, conv_metric, atol=2e-3, rtol=1e-4)
            )
            if not unimodular_ok:
                raise ValueError(
                    "GROMACS round-trip periodic box mismatch between source Amber and exported GRO"
                )
        summary["periodic_box_verified"] = True
    else:
        summary["periodic_box_verified"] = False
    summary["coordinate_tolerance_angstrom"] = 0.006
    if not np.allclose(original.coordinates, converted.coordinates, atol=0.006, rtol=0):
        raise ValueError("GROMACS round-trip coordinate mismatch beyond GRO precision")
    return summary


def _atom_signature(atom) -> tuple[str, str, str]:
    res = getattr(atom, "residue", None)
    return (
        str(getattr(atom, "name", "")),
        str(getattr(res, "name", "")) if res is not None else "",
        str(getattr(atom, "type", "")),
    )


def _structure_box_vectors_nm(structure) -> np.ndarray | None:
    box = getattr(structure, "box", None)
    box_vectors = getattr(structure, "box_vectors", None)
    if box is None or box_vectors is None:
        return None
    import parmed as pmd

    vecs_nm = np.asarray(box_vectors.value_in_unit(pmd.unit.nanometer), dtype=float)
    if vecs_nm.shape != (3, 3) or not np.all(np.isfinite(vecs_nm)):
        raise ValueError("GROMACS round-trip periodic box mismatch: non-finite box vectors")
    return vecs_nm


def _molecule_instance_explicit_exclusions(converted) -> set[tuple[int, int]]:
    explicit_excl = {
        tuple(sorted((int(a.idx), int(p.idx))))
        for a in converted.atoms
        for p in (getattr(a, "exclusion_partners", ()) or ())
        if int(p.idx) != int(a.idx)
    }
    molecules = getattr(converted, "molecules", None) or {}
    if not molecules:
        return explicit_excl

    mol_specs: list[tuple[int, list[tuple[str, str, str]], set[tuple[int, int]]]] = []
    for mol_entry in molecules.values():
        mol_obj = mol_entry[0] if isinstance(mol_entry, tuple) and mol_entry else mol_entry
        mol_atoms = list(getattr(mol_obj, "atoms", ()) or ())
        if not mol_atoms:
            continue
        idx_by_id = {id(a): idx for idx, a in enumerate(mol_atoms)}
        mol_pairs: set[tuple[int, int]] = set()
        for idx, atom in enumerate(mol_atoms):
            for partner in (getattr(atom, "exclusion_partners", ()) or ()):
                p_idx = idx_by_id.get(id(partner), int(partner.idx))
                if p_idx != idx:
                    mol_pairs.add((min(idx, p_idx), max(idx, p_idx)))
        mol_specs.append((len(mol_atoms), [_atom_signature(a) for a in mol_atoms], mol_pairs))

    if not any(pairs for _, _, pairs in mol_specs):
        return explicit_excl

    conv_sigs = [_atom_signature(a) for a in converted.atoms]
    total_atoms = len(conv_sigs)
    if len(mol_specs) == 1:
        mol_len, mol_sig, mol_pairs = mol_specs[0]
        if mol_len <= 0 or total_atoms % mol_len != 0:
            raise ValueError(
                "GROMACS round-trip unmappable moleculetype explicit exclusions"
            )
        for offset in range(0, total_atoms, mol_len):
            if conv_sigs[offset : offset + mol_len] != mol_sig:
                raise ValueError(
                    "GROMACS round-trip unmappable moleculetype explicit exclusions"
                )
            for u, v in mol_pairs:
                explicit_excl.add((offset + u, offset + v))
        return explicit_excl

    cursor = 0
    while cursor < total_atoms:
        matches = [
            (mol_len, mol_pairs)
            for mol_len, mol_sig, mol_pairs in mol_specs
            if cursor + mol_len <= total_atoms and conv_sigs[cursor : cursor + mol_len] == mol_sig
        ]
        if len(matches) != 1:
            raise ValueError(
                "GROMACS round-trip unmappable multi-moleculetype explicit exclusions"
            )
        mol_len, mol_pairs = matches[0]
        for u, v in mol_pairs:
            explicit_excl.add((cursor + u, cursor + v))
        cursor += mol_len
    return explicit_excl


_ALLOWED_GROMACS_SECTIONS = frozenset(
    {
        "defaults",
        "atomtypes",
        "bondtypes",
        "pairtypes",
        "angletypes",
        "dihedraltypes",
        "nonbond_params",
        "moleculetype",
        "atoms",
        "bonds",
        "pairs",
        "angles",
        "dihedrals",
        "exclusions",
        "system",
        "molecules",
        "settles",
    }
)


def _verify_gromacs_topology_sections(top_path: Path, *, converted) -> None:
    import parmed as pmd

    molecules = getattr(converted, "molecules", None) or {}
    mol_by_name = {
        str(mol_name): mol_entry[0]
        for mol_name, mol_entry in molecules.items()
        if isinstance(mol_entry, tuple) and mol_entry
    }
    _moltype_reset_sections = frozenset(
        {
            "defaults",
            "atomtypes",
            "bondtypes",
            "pairtypes",
            "angletypes",
            "dihedraltypes",
            "nonbond_params",
            "moleculetype",
            "system",
            "molecules",
        }
    )
    _max_section_tokens = {
        "atoms": 8,
        "bonds": 5,
        "angles": 6,
        "dihedrals": 8,
        "pairs": 5,
    }

    def _is_eligible_settle_water(mol_name: str | None) -> bool:
        if mol_name is None or mol_name not in mol_by_name:
            return False
        mol_obj = mol_by_name[mol_name]
        real_atoms = [a for a in mol_obj.atoms if not isinstance(a, pmd.ExtraPoint)]
        return (
            len(real_atoms) == 3
            and len(mol_obj.atoms) >= 3
            and [int(a.atomic_number) for a in mol_obj.atoms[:3]] == [8, 1, 1]
        )

    def _validate_settle_row(mol_name: str | None, toks: list[str], filename: str) -> None:
        if not _is_eligible_settle_water(mol_name):
            raise ValueError(
                f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
            )
        assert mol_name is not None
        mol_obj = mol_by_name[mol_name]
        try:
            if len(toks) != 4:
                raise ValueError
            row_atom = int(toks[0])
            funct = int(toks[1])
            doh_nm = float(toks[2])
            dhh_nm = float(toks[3])
        except ValueError as exc:
            raise ValueError(
                f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
            ) from exc
        if (
            row_atom != 1
            or funct != 1
            or not (math.isfinite(doh_nm) and doh_nm > 0.0)
            or not (math.isfinite(dhh_nm) and dhh_nm > 0.0)
        ):
            raise ValueError(
                f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
            )
        idx_by_id = {id(a): idx for idx, a in enumerate(mol_obj.atoms)}
        bond_reqs: dict[tuple[int, int], float] = {}
        for bond in getattr(mol_obj, "bonds", ()) or ():
            if bond.type is None:
                continue
            i1 = idx_by_id.get(id(bond.atom1), int(bond.atom1.idx))
            i2 = idx_by_id.get(id(bond.atom2), int(bond.atom2.idx))
            bond_reqs[(min(i1, i2), max(i1, i2))] = float(bond.type.req)
        r01 = bond_reqs.get((0, 1))
        r02 = bond_reqs.get((0, 2))
        if (
            r01 is None
            or r02 is None
            or not (math.isfinite(r01) and r01 > 0.0 and math.isfinite(r02) and r02 > 0.0)
            or abs(r01 - r02) > 1e-4
        ):
            raise ValueError(
                f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
            )
        doh_a = doh_nm * 10.0
        if abs(doh_a - r01) > 2e-3 or abs(doh_a - r02) > 2e-3:
            raise ValueError(
                f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
            )
        dhh_a = dhh_nm * 10.0
        r12 = bond_reqs.get((1, 2))
        if r12 is not None and math.isfinite(r12) and r12 > 0.0:
            expected_hh = r12
        else:
            theta_deg = None
            for ang in getattr(mol_obj, "angles", ()) or ():
                if ang.type is None:
                    continue
                i1 = idx_by_id.get(id(ang.atom1), int(ang.atom1.idx))
                i2 = idx_by_id.get(id(ang.atom2), int(ang.atom2.idx))
                i3 = idx_by_id.get(id(ang.atom3), int(ang.atom3.idx))
                if i2 == 0 and {i1, i3} == {1, 2}:
                    theta_deg = float(ang.type.theteq)
                    break
            if theta_deg is None or not (math.isfinite(theta_deg) and 0.0 < theta_deg < 180.0):
                raise ValueError(
                    f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
                )
            expected_hh = 2.0 * r01 * math.sin(math.radians(theta_deg) / 2.0)
        if abs(dhh_a - expected_hh) > 2e-3:
            raise ValueError(
                f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {filename}"
            )

    def _scan_with_defines(defines: frozenset[str]) -> tuple[set[str], set[str]]:
        active_stack: set[Path] = set()
        include_count = [0]
        global_cond_depth = [0]
        seen_moltypes: set[str] = set()
        flex_cond_moltypes: set[str] = set()
        rigid_settle_moltypes: set[str] = set()
        state: dict[str, object] = {
            "section": None,
            "moltype": None,
            "moltype_has_settles": False,
            "in_settles": False,
            "settles_rows": 0,
            "settles_file": None,
        }

        def _finalize_settles() -> None:
            if state["in_settles"] and int(state["settles_rows"]) != 1:  # type: ignore[arg-type]
                raise ValueError(
                    f"GROMACS round-trip invalid [settles] row count in {state['settles_file']}"
                )
            state["in_settles"] = False
            state["settles_rows"] = 0
            state["settles_file"] = None

        def _check_conditional_scope(sec_name: str | None, filename: str) -> None:
            if global_cond_depth[0] <= 0:
                return
            if "FLEXIBLE" not in defines:
                if sec_name != "settles":
                    raise ValueError(
                        f"GROMACS round-trip unverified non-SETTLE directive/data inside rigid FLEXIBLE conditional in {filename}"
                    )
                return
            mol_name = state["moltype"] if isinstance(state["moltype"], str) else None
            if sec_name not in {"bonds", "angles"} or not _is_eligible_settle_water(mol_name):
                raise ValueError(
                    f"GROMACS round-trip unverified non-water/non-bonded directive/data inside FLEXIBLE conditional in {filename}"
                )
            assert mol_name is not None
            flex_cond_moltypes.add(mol_name)

        def _scan(cur: Path) -> None:
            resolved = cur.resolve()
            if resolved in active_stack:
                raise ValueError(f"GROMACS round-trip cyclic #include detected: {resolved.name}")
            if not resolved.is_file():
                raise ValueError(f"GROMACS round-trip missing topology file {resolved}")
            if include_count[0] >= 1024:
                raise ValueError("GROMACS round-trip exceeded maximum #include traversal depth")
            include_count[0] += 1
            active_stack.add(resolved)
            cond_stack: list[dict[str, bool]] = []
            try:
                for raw_line in resolved.read_text(encoding="utf-8").splitlines():
                    stripped_raw = raw_line.rstrip()
                    code = raw_line.split(";", 1)[0].strip()
                    if stripped_raw.endswith("\\") or code.endswith("\\"):
                        raise ValueError(
                            f"GROMACS round-trip unsupported line continuation (\\) in {resolved.name}"
                        )
                    if not code:
                        continue
                    is_active = all(frame["active"] for frame in cond_stack)
                    if code.startswith("#"):
                        if re.match(r"^#\s*(?:ifdef|ifndef)(?=$|[^A-Za-z0-9_])", code):
                            m_if = re.match(r"^#\s*(ifdef|ifndef)\s+FLEXIBLE\s*$", code)
                            if not m_if:
                                raise ValueError(
                                    f"GROMACS round-trip unsupported preprocessor directive in {resolved.name}: {code!r}"
                                )
                            if global_cond_depth[0] >= 1:
                                raise ValueError(
                                    f"GROMACS round-trip nested preprocessor conditional in {resolved.name}: {code!r}"
                                )
                            cond_true = (
                                ("FLEXIBLE" in defines)
                                if m_if.group(1) == "ifdef"
                                else ("FLEXIBLE" not in defines)
                            )
                            cond_stack.append({"active": cond_true, "saw_else": False})
                            global_cond_depth[0] += 1
                            continue
                        if re.match(r"^#\s*else(?=$|[^A-Za-z0-9_])", code):
                            if not re.match(r"^#\s*else\s*$", code):
                                raise ValueError(
                                    f"GROMACS round-trip unsupported preprocessor directive in {resolved.name}: {code!r}"
                                )
                            if not cond_stack or cond_stack[-1]["saw_else"]:
                                raise ValueError(
                                    f"GROMACS round-trip unexpected #else directive in {resolved.name}"
                                )
                            cond_stack[-1]["saw_else"] = True
                            cond_stack[-1]["active"] = not cond_stack[-1]["active"]
                            continue
                        if re.match(r"^#\s*endif(?=$|[^A-Za-z0-9_])", code):
                            if not re.match(r"^#\s*endif\s*$", code):
                                raise ValueError(
                                    f"GROMACS round-trip unsupported preprocessor directive in {resolved.name}: {code!r}"
                                )
                            if not cond_stack:
                                raise ValueError(
                                    f"GROMACS round-trip unexpected #endif directive in {resolved.name}"
                                )
                            cond_stack.pop()
                            global_cond_depth[0] -= 1
                            continue
                        if re.match(r"^#\s*include(?=$|[^A-Za-z0-9_])", code):
                            inc_match = re.match(
                                r'^#\s*include\s*(?:"([^"]+)"|<([^>]+)>)\s*$', code
                            )
                            if not inc_match:
                                raise ValueError(
                                    f"GROMACS round-trip malformed #include directive in {resolved.name}: {code!r}"
                                )
                            if is_active:
                                inc_target = inc_match.group(1) or inc_match.group(2)
                                inc_path = (resolved.parent / inc_target).resolve()
                                if not inc_path.is_file():
                                    raise ValueError(
                                        f"GROMACS round-trip missing #include target {inc_target!r} in {resolved.name}"
                                    )
                                _scan(inc_path)
                            continue
                        raise ValueError(
                            f"GROMACS round-trip unsupported preprocessor directive in {resolved.name}: {code!r}"
                        )
                    if not is_active:
                        continue
                    sec_match = re.match(r"^\[\s*([^\]]+?)\s*\]$", code)
                    if sec_match:
                        _finalize_settles()
                        sec = sec_match.group(1).strip().lower()
                        state["section"] = sec
                        if sec in _moltype_reset_sections:
                            state["moltype"] = None
                            state["moltype_has_settles"] = False
                        if sec not in _ALLOWED_GROMACS_SECTIONS:
                            raise ValueError(
                                f"GROMACS round-trip unsupported/unverified directive [{sec}] in {resolved.name}"
                            )
                        _check_conditional_scope(sec, resolved.name)
                        if sec == "settles":
                            mol_name = state["moltype"] if isinstance(state["moltype"], str) else None
                            if mol_name is None or mol_name not in mol_by_name:
                                raise ValueError(
                                    f"GROMACS round-trip unsupported/unverified directive [settles] in {resolved.name}"
                                )
                            mol_obj = mol_by_name[mol_name]
                            real_atoms = [a for a in mol_obj.atoms if not isinstance(a, pmd.ExtraPoint)]
                            if (
                                len(real_atoms) != 3
                                or sorted(int(a.atomic_number) for a in real_atoms) != [1, 1, 8]
                            ):
                                raise ValueError(
                                    f"GROMACS round-trip unsupported/unverified directive [settles] in {resolved.name}"
                                )
                            if not _is_eligible_settle_water(mol_name):
                                raise ValueError(
                                    f"GROMACS round-trip SETTLE water geometry/ordering mismatch in {resolved.name}"
                                )
                            if state["moltype_has_settles"]:
                                raise ValueError(
                                    f"GROMACS round-trip duplicate [settles] directive in moleculetype {mol_name!r} of {resolved.name}"
                                )
                            state["moltype_has_settles"] = True
                            rigid_settle_moltypes.add(mol_name)
                            state["in_settles"] = True
                            state["settles_rows"] = 0
                            state["settles_file"] = resolved.name
                        continue
                    if code.startswith("["):
                        raise ValueError(
                            f"GROMACS round-trip malformed section header in {resolved.name}: {code!r}"
                        )
                    toks = code.split()
                    sec = state["section"] if isinstance(state["section"], str) else None
                    _check_conditional_scope(sec, resolved.name)
                    if sec == "moleculetype":
                        if state["moltype"] is not None:
                            raise ValueError(
                                f"GROMACS round-trip unexpected multiple entries in [moleculetype] of {resolved.name}"
                            )
                        mol_name = toks[0]
                        if mol_name in seen_moltypes:
                            raise ValueError(
                                f"GROMACS round-trip duplicate [moleculetype] {mol_name!r} in {resolved.name}"
                            )
                        seen_moltypes.add(mol_name)
                        state["moltype"] = mol_name
                    elif sec in _max_section_tokens and len(toks) > _max_section_tokens[sec]:
                        raise ValueError(
                            f"GROMACS round-trip unsupported B-state columns in [{sec}] of {resolved.name}"
                        )
                    elif sec == "settles":
                        rows = int(state["settles_rows"]) + 1  # type: ignore[arg-type]
                        state["settles_rows"] = rows
                        if rows > 1:
                            raise ValueError(
                                f"GROMACS round-trip invalid [settles] row count in {resolved.name}"
                            )
                        mol_name = state["moltype"] if isinstance(state["moltype"], str) else None
                        _validate_settle_row(mol_name, toks, resolved.name)
                if cond_stack:
                    raise ValueError(
                        f"GROMACS round-trip unterminated preprocessor conditional in {resolved.name}"
                    )
            finally:
                active_stack.discard(resolved)

        _scan(Path(top_path))
        _finalize_settles()
        return flex_cond_moltypes, rigid_settle_moltypes

    _, rigid_settle_moltypes = _scan_with_defines(frozenset())
    flex_cond_moltypes, _ = _scan_with_defines(frozenset({"FLEXIBLE"}))
    missing_rigid = sorted(flex_cond_moltypes - rigid_settle_moltypes)
    if missing_rigid:
        raise ValueError(
            f"GROMACS round-trip FLEXIBLE-conditional water missing rigid [settles] for {missing_rigid}"
        )


def _link_upstream_provenance(topology, coordinates, *, top_sha: str | None = None, crd_sha: str | None = None):
    top_path = Path(topology).resolve()
    crd_path = Path(coordinates).resolve()
    top_sha = top_sha if top_sha is not None else _sha(top_path)
    crd_sha = crd_sha if crd_sha is not None else _sha(crd_path)
    upstream_path = top_path.parent / "PROVENANCE.json"
    upstream_data = None
    upstream_sha = None
    verified = False
    if upstream_path.is_file():
        try:
            raw_bytes = upstream_path.read_bytes()
            candidate = json.loads(raw_bytes.decode("utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            raw_bytes = None
            candidate = None
        files_meta = candidate.get("files") if isinstance(candidate, dict) else None
        if (
            raw_bytes is not None
            and isinstance(files_meta, dict)
            and files_meta.get("prmtop_sha256") == top_sha
            and files_meta.get("inpcrd_sha256") == crd_sha
        ):
            upstream_data = candidate
            upstream_sha = hashlib.sha256(raw_bytes).hexdigest()
            verified = True
    return {
        "source_topology": str(top_path),
        "source_topology_sha256": top_sha,
        "source_coordinates": str(crd_path),
        "source_coordinates_sha256": crd_sha,
        "amber_provenance": str(upstream_path) if verified else None,
        "amber_provenance_sha256": upstream_sha if verified else None,
        "amber_provenance_verified": verified,
        "preset": upstream_data.get("preset") if verified else None,
        "charge_model": upstream_data.get("charge_model") if verified else None,
        "charge_regime": upstream_data.get("charge_regime") if verified else None,
        "custom_combination": upstream_data.get("custom_combination") if verified else None,
    }


def _upgrade_gromacs_precision(top_path: Path, itp_path: Path, structure, converted) -> int:
    """Rewrite ParmEd GROMACS force-field numbers at full float64 precision.

    ParmEd formats .top/.itp values with 4-8 digits (e.g. fudgeQQ=0.83333333,
    6-8 decimal LJ sigma/eps and partial charges, and 4-5 digit bond force
    constants), introducing up to ~0.08 kJ/mol energy drift and >1e-3 kJ/(mol*nm)
    force drift on unrelaxed contacts relative to Amber/OpenMM.
    """
    scale_sig = (2.0 ** (5.0 / 6.0)) * 0.1
    type_lj: dict[str, tuple[float, float]] = {}
    for s_atom, c_atom in zip(structure.atoms, converted.atoms, strict=True):
        rmin = float(s_atom.rmin)
        sig = 0.0 if rmin == 0.0 else rmin * scale_sig
        eps = float(s_atom.epsilon) * 4.184
        type_lj.setdefault(str(c_atom.type), (sig, eps))

    scee_vals = [
        float(dt.scee)
        for d in structure.dihedrals
        if not getattr(d, "improper", False) and not getattr(d, "ignore_end", False)
        for dt in (d.type if isinstance(d.type, list) else [d.type])
        if dt is not None and float(getattr(dt, "scee", 0.0)) > 0.0
    ]
    scnb_vals = [
        float(dt.scnb)
        for d in structure.dihedrals
        if not getattr(d, "improper", False) and not getattr(d, "ignore_end", False)
        for dt in (d.type if isinstance(d.type, list) else [d.type])
        if dt is not None and float(getattr(dt, "scnb", 0.0)) > 0.0
    ]
    fudge_qq = 1.0 / scee_vals[0] if scee_vals and all(abs(x - scee_vals[0]) < 1e-9 for x in scee_vals) else None
    fudge_lj = 1.0 / scnb_vals[0] if scnb_vals and all(abs(x - scnb_vals[0]) < 1e-9 for x in scnb_vals) else None

    top_section = None
    top_out: list[str] = []
    for raw_line in top_path.read_text(encoding="utf-8").splitlines(keepends=True):
        stripped = raw_line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            top_section = stripped[1:-1].strip().lower()
            top_out.append(raw_line)
            continue
        if stripped and not stripped.startswith((";", "#", "*")):
            body, sep, comment = raw_line.rstrip("\r\n").partition(";")
            parts = body.split()
            suffix = f" ; {comment.strip()}" if sep else ""
            if top_section == "defaults" and len(parts) >= 5:
                flj = fudge_lj if fudge_lj is not None else float(parts[3])
                fqq = fudge_qq if fudge_qq is not None else float(parts[4])
                top_out.append(f"{parts[0]} {parts[1]} {parts[2]} {flj:.15g} {fqq:.15g}{suffix}\n")
                continue
            if top_section == "atomtypes" and len(parts) >= 7 and parts[0] in type_lj:
                sig, eps = type_lj[parts[0]]
                top_out.append(
                    f"{parts[0]:<10s} {parts[1]:>5s} {parts[2]:>11s} {parts[3]:>11s} {parts[4]:>2s} "
                    f"{sig:21.15e} {eps:21.15e}{suffix}\n"
                )
                continue
        top_out.append(raw_line)
    top_path.write_text("".join(top_out), encoding="utf-8")

    bond_lookup: dict[tuple[int, int], tuple[float, float]] = {}
    for bond in structure.bonds:
        if bond.type is None:
            continue
        i1, i2 = bond.atom1.idx + 1, bond.atom2.idx + 1
        key = (min(i1, i2), max(i1, i2))
        bond_lookup[key] = (float(bond.type.req) * 0.1, float(bond.type.k) * 836.8)

    updated_bonds = 0
    section = None
    out_lines: list[str] = []
    for raw_line in itp_path.read_text(encoding="utf-8").splitlines(keepends=True):
        stripped = raw_line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            out_lines.append(raw_line)
            continue
        if section == "atoms" and stripped and not stripped.startswith((";", "#")):
            body, sep, comment = raw_line.rstrip("\r\n").partition(";")
            parts = body.split()
            if len(parts) >= 8:
                try:
                    nr = int(parts[0])
                except ValueError:
                    out_lines.append(raw_line)
                    continue
                if 1 <= nr <= len(structure.atoms):
                    s_atom = structure.atoms[nr - 1]
                    suffix = f" ; {comment.strip()}" if sep else ""
                    extra = (" " + " ".join(parts[8:])) if len(parts) > 8 else ""
                    out_lines.append(
                        f"{nr:6d} {parts[1]:>10s} {parts[2]:>6s} {parts[3]:>6s} "
                        f"{parts[4]:>6s} {parts[5]:>6s} {float(s_atom.charge):22.15e} "
                        f"{float(s_atom.mass):15.10f}{extra}{suffix}\n"
                    )
                    continue
        if section == "bonds" and stripped and not stripped.startswith((";", "#")):
            body, sep, comment = raw_line.rstrip("\r\n").partition(";")
            parts = body.split()
            if len(parts) >= 5 and parts[2] == "1":
                try:
                    ai, aj = int(parts[0]), int(parts[1])
                except ValueError:
                    out_lines.append(raw_line)
                    continue
                key = (min(ai, aj), max(ai, aj))
                if key in bond_lookup:
                    c0, c1 = bond_lookup[key]
                    suffix = f" ; {comment.strip()}" if sep else ""
                    out_lines.append(f"{ai:6d} {aj:6d} {1:5d} {c0:21.15e} {c1:21.15e}{suffix}\n")
                    updated_bonds += 1
                    continue
        out_lines.append(raw_line)
    itp_path.write_text("".join(out_lines), encoding="utf-8")
    return updated_bonds


def export_gromacs(topology, coordinates, output_prefix, *, overwrite=False, model=None):
    """Convert the same complete Amber parameterization with ParmEd.

    Winding graph export is deliberately gated even though GROMACS documents
    periodic-molecules=yes: preserving arbitrary image/angle/dihedral semantics
    from this graph through an Amber intermediate is not validated here.
    """
    if model and (_periodic_edges(model) or model.get("metadata", {}).get("covalently_periodic")):
        raise ValueError(
            "Periodic winding GROMACS graph export not implemented/validated; periodic-molecules=yes alone is insufficient"
        )
    prefix = Path(output_prefix).resolve()
    paths = {
        key: prefix.with_suffix(suffix)
        for key, suffix in [
            ("top", ".top"),
            ("itp", ".itp"),
            ("gro", ".gro"),
            ("manifest", ".manifest.json"),
            ("provenance", ".provenance.json"),
        ]
    }
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("GROMACS export files already exist; choose a fresh prefix")
    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh prefix")
    structure, top_sha, crd_sha = _load_amber_with_digests(topology, coordinates)
    upstream = _link_upstream_provenance(
        topology, coordinates, top_sha=top_sha, crd_sha=crd_sha
    )
    _verify_supported_amber_terms(
        structure, target="GROMACS", allow_explicit_exclusions=False
    )
    import parmed as pmd

    converted = pmd.gromacs.GromacsTopologyFile.from_structure(structure)
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
        with tempfile.TemporaryDirectory(
            dir=prefix.parent, prefix=f".{prefix.name}.tmp."
        ) as stage_dir_str:
            stage_dir = Path(stage_dir_str)
            staged = {key: stage_dir / fp.name for key, fp in paths.items()}
            converted.write(str(staged["top"]), molfile=str(staged["itp"]), combine="all")
            text = staged["top"].read_text().replace(
                f'"{staged["itp"]}"', f'"{paths["itp"].name}"'
            )
            staged["top"].write_text(text)
            bonds_upgraded = _upgrade_gromacs_precision(
                staged["top"], staged["itp"], structure, converted
            )
            structure.save(str(staged["gro"]), overwrite=False)
            roundtrip = pmd.load_file(str(staged["top"]), xyz=str(staged["gro"]))
            verification = _verify_parameter_roundtrip(
                structure, roundtrip, top_path=staged["top"]
            )
            prov = {
                "schema_version": 1,
                "format": "gromacs",
                **upstream,
                "itp_bond_precision_upgraded": bonds_upgraded,
                "outputs": {k: str(v) for k, v in paths.items() if k != "provenance"},
            }
            _json(staged["provenance"], prov)
            result = {key: str(value) for key, value in paths.items()}
            result.update(
                status="exported",
                parameterized=True,
                roundtrip_verified=verification,
                source_topology_sha256=upstream["source_topology_sha256"],
                itp_bond_precision_upgraded=bonds_upgraded,
                preset=upstream["preset"],
                charge_model=upstream["charge_model"],
                atom_count=len(structure.atoms),
                bond_count=len(structure.bonds),
                periodic_molecules="no (finite covalent topology only)",
                gbsa_note="GB radii/screens remain in source Amber topology; standard GROMACS export does not transfer these fields",
                engine_sources=ENGINE_SOURCES,
            )
            _json(staged["manifest"], result)
            for key in ("top", "itp", "gro", "provenance", "manifest"):
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


def _fmt_xml_float(value: float) -> str:
    val = float(value)
    if val == 0.0:
        return "0"
    return format(val, ".16g")


def _validate_14_scaling_factor(scee: float, scnb: float, *, context: str) -> tuple[float, float]:
    """Validate finite, positive, standard Amber 1-4 scaling factors (scee=1.2, scnb=2.0)."""
    f_ee = float(scee)
    f_nb = float(scnb)
    if not (math.isfinite(f_ee) and math.isfinite(f_nb)) or f_ee <= 0.0 or f_nb <= 0.0:
        raise ValueError(
            f"1-4 scaling factors on {context} must be positive and finite "
            f"(got scee={scee!r}, scnb={scnb!r})"
        )
    if not (
        math.isclose(f_ee, 1.2, rel_tol=1e-6, abs_tol=1e-9)
        and math.isclose(f_nb, 2.0, rel_tol=1e-6, abs_tol=1e-9)
    ):
        raise ValueError(
            f"Non-standard Amber 1-4 scaling factors on {context} "
            f"(scee={f_ee!r}, scnb={f_nb!r}; exporter requires scee=1.2, scnb=2.0)"
        )
    return f_ee, f_nb


def _bond_graph_pairs_upto(structure, max_dist: int = 3) -> dict[int, set[tuple[int, int]]]:
    """Return canonical (i < j) atom pairs keyed by shortest bond-graph distance 1..max_dist."""
    if max_dist < 1:
        return {}
    from collections import deque

    n_atoms = len(structure.atoms)
    adj: list[list[int]] = [[] for _ in range(n_atoms)]
    for bond in structure.bonds:
        i, j = int(bond.atom1.idx), int(bond.atom2.idx)
        if i != j:
            adj[i].append(j)
            adj[j].append(i)
    by_dist: dict[int, set[tuple[int, int]]] = {d: set() for d in range(1, max_dist + 1)}
    for start in range(n_atoms):
        gdist = {start: 0}
        dq = deque([start])
        while dq:
            u = dq.popleft()
            du = gdist[u]
            if du == max_dist:
                continue
            for v in adj[u]:
                if v not in gdist:
                    nd = du + 1
                    gdist[v] = nd
                    dq.append(v)
                    if start < v:
                        by_dist[nd].add((start, v))
    return by_dist


def _extract_amber_exclusions(
    structure,
) -> tuple[
    set[tuple[int, int]],
    dict[tuple[int, int], tuple[float, float]],
    set[tuple[int, int]],
    set[tuple[int, int]],
    set[tuple[int, int]],
    dict[int, set[tuple[int, int]]],
]:
    """Return 1-2/1-3 pairs, scaled 1-4 pairs, extra/unexcluded/zero-excl sets, and graph distances."""
    n_atoms = len(structure.atoms)
    excluded_12_13: set[tuple[int, int]] = set()
    for bond in structure.bonds:
        i, j = bond.atom1.idx, bond.atom2.idx
        if i != j:
            excluded_12_13.add((min(i, j), max(i, j)))
    for angle in structure.angles:
        i, j = angle.atom1.idx, angle.atom3.idx
        if i != j:
            excluded_12_13.add((min(i, j), max(i, j)))

    graph_d = _bond_graph_pairs_upto(structure, 3)
    graph_12 = graph_d[1]
    graph_13 = graph_d[2]
    graph_14 = graph_d[3]

    explicit_partner_pairs: set[tuple[int, int]] = set()
    amber_excluded: set[tuple[int, int]] = set()
    for atom in structure.atoms:
        i = atom.idx
        for partner in getattr(atom, "exclusion_partners", ()) or ():
            j = partner.idx
            if j != i:
                pair = (min(i, j), max(i, j))
                explicit_partner_pairs.add(pair)
                amber_excluded.add(pair)

    parm_data = getattr(structure, "parm_data", None) or {}
    num_excl = parm_data.get("NUMBER_EXCLUDED_ATOMS")
    excl_list = parm_data.get("EXCLUDED_ATOMS_LIST")
    has_amber_excl_list = num_excl is not None and excl_list is not None
    if has_amber_excl_list:
        pos = 0
        for i, raw_cnt in enumerate(num_excl):
            if i >= n_atoms:
                break
            cnt = int(raw_cnt)
            if cnt < 0 or pos + cnt > len(excl_list):
                raise ValueError(
                    "Malformed NUMBER_EXCLUDED_ATOMS / EXCLUDED_ATOMS_LIST in Amber topology"
                )
            for k in range(cnt):
                j1 = int(excl_list[pos + k])
                if j1 <= 0:
                    continue
                if j1 > n_atoms:
                    raise ValueError(
                        f"Out-of-range excluded atom index {j1} in EXCLUDED_ATOMS_LIST"
                    )
                j = j1 - 1
                if j != i:
                    amber_excluded.add((min(i, j), max(i, j)))
            pos += cnt
        missing_bonded = excluded_12_13 - amber_excluded
        if missing_bonded:
            raise ValueError(
                f"Amber EXCLUDED_ATOMS_LIST omits {len(missing_bonded)} 1-2/1-3 bonded pair(s); "
                f"non-excluded bonded pairs are unsupported"
            )

    for flag in ("DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"):
        raw_dih = parm_data.get(flag)
        if raw_dih is None:
            continue
        n_raw = len(raw_dih)
        if n_raw % 5 != 0:
            raise ValueError(f"Malformed {flag} section length ({n_raw}) in Amber topology")
        for k in range(0, n_raw, 5):
            if int(raw_dih[k + 2]) == 0 or int(raw_dih[k + 3]) == 0:
                raise ValueError(
                    f"Amber {flag} entry {k // 5} places atom index 0 in sign-flagged slot 3/4 "
                    f"({[int(x) for x in raw_dih[k : k + 5]]}), which has ambiguous 1-4/improper semantics"
                )

    pairs_14: dict[tuple[int, int], tuple[float, float]] = {}
    dihedral_14_ends: set[tuple[int, int]] = set()
    for dihedral in structure.dihedrals:
        if getattr(dihedral, "improper", False):
            continue
        i, j = dihedral.atom1.idx, dihedral.atom4.idx
        pair = (min(i, j), max(i, j))
        is_active = not getattr(dihedral, "ignore_end", False)
        if is_active:
            if (
                i == j
                or pair in excluded_12_13
                or pair in graph_12
                or pair in graph_13
                or pair in explicit_partner_pairs
            ):
                raise ValueError(
                    f"Amber topology contains active 1-4 dihedral on bonded/excluded pair {(i, j)}, unsupported by export"
                )
        elif (
            i != j
            and pair not in excluded_12_13
            and pair not in graph_12
            and pair not in graph_13
        ):
            dihedral_14_ends.add(pair)
        if not is_active:
            continue
        dt = dihedral.type[0] if isinstance(dihedral.type, list) else dihedral.type
        if dt is None:
            raise ValueError(f"Untyped dihedral on 1-4 pair {pair}")
        scee, scnb = _validate_14_scaling_factor(
            getattr(dt, "scee", 0.0),
            getattr(dt, "scnb", 0.0),
            context=f"1-4 dihedral pair {pair}",
        )
        if pair in pairs_14:
            raise ValueError(
                f"Amber topology contains duplicate active 1-4 dihedral interaction on atom pair {pair}, unsupported by export"
            )
        pairs_14[pair] = (scee, scnb)

    if has_amber_excl_list:
        missing_14_excl = set(pairs_14) - amber_excluded
        if missing_14_excl:
            raise ValueError(
                f"Amber EXCLUDED_ATOMS_LIST omits {len(missing_14_excl)} active 1-4 pair(s)"
            )
        unexcluded_14 = (graph_14 | dihedral_14_ends) - amber_excluded - set(pairs_14)
        orig_zero_excl = (amber_excluded | explicit_partner_pairs) - set(pairs_14)
    else:
        unexcluded_14 = set()
        orig_zero_excl = (excluded_12_13 | explicit_partner_pairs) - set(pairs_14)

    extra_exclusions = (amber_excluded | explicit_partner_pairs) - excluded_12_13 - set(pairs_14)
    return excluded_12_13, pairs_14, extra_exclusions, unexcluded_14, orig_zero_excl, graph_d


def _verify_standard_lorentz_berthelot_lj(structure, *, target: str = "OpenMM XML") -> None:
    """Reject NBFIX, non-Lorentz mixing, 12-6-4 C4 terms, and 10-12 H-bond terms."""
    comb = str(getattr(structure, "combining_rule", "lorentz") or "lorentz").lower()
    if comb != "lorentz":
        raise ValueError(
            f"Unsupported Amber LJ combining_rule {comb!r}; {target} export requires Lorentz-Berthelot"
        )
    if hasattr(structure, "has_NBFIX") and structure.has_NBFIX():
        raise ValueError(
            f"Amber topology contains pair-specific Lennard-Jones overrides (NBFIX); "
            f"{target} export cannot represent NBFIX"
        )
    if hasattr(structure, "has_1012") and structure.has_1012():
        raise ValueError(
            f"Amber topology contains 10-12 HBOND terms unsupported by {target} export"
        )
    parm_data = getattr(structure, "parm_data", None) or {}
    for flag in (
        "LENNARD_JONES_ACOEF",
        "LENNARD_JONES_BCOEF",
        "LENNARD_JONES_14_ACOEF",
        "LENNARD_JONES_14_BCOEF",
        "HBOND_ACOEF",
        "HBOND_BCOEF",
    ):
        vals = parm_data.get(flag)
        if vals and any(not math.isfinite(float(x)) for x in vals):
            raise ValueError(f"Nonfinite entry in Amber {flag} table unsupported by {target} export")
    ccoef = parm_data.get("LENNARD_JONES_CCOEF")
    if ccoef and any((not math.isfinite(float(x))) or abs(float(x)) > 0.0 for x in ccoef):
        raise ValueError(
            f"Amber topology contains nonzero or nonfinite 12-6-4 LENNARD_JONES_CCOEF terms unsupported by {target} export"
        )
    pointers = parm_data.get("POINTERS")
    nb_index = parm_data.get("NONBONDED_PARM_INDEX")
    acoef = parm_data.get("LENNARD_JONES_ACOEF")
    bcoef = parm_data.get("LENNARD_JONES_BCOEF")
    c14 = parm_data.get("LENNARD_JONES_14_CCOEF")
    if c14 and any((not math.isfinite(float(x))) or abs(float(x)) > 0.0 for x in c14):
        raise ValueError(
            f"Amber/Chamber topology contains nonzero or nonfinite 12-6-4 LENNARD_JONES_14_CCOEF terms unsupported by {target} export"
        )
    a14 = parm_data.get("LENNARD_JONES_14_ACOEF")
    b14 = parm_data.get("LENNARD_JONES_14_BCOEF")
    if a14 is not None or b14 is not None:
        if (
            acoef is None
            or bcoef is None
            or a14 is None
            or b14 is None
            or len(a14) != len(acoef)
            or len(b14) != len(bcoef)
            or not all(
                math.isclose(float(x14), float(x), rel_tol=1e-6, abs_tol=1e-9)
                for x14, x in zip(a14, acoef, strict=True)
            )
            or not all(
                math.isclose(float(x14), float(x), rel_tol=1e-6, abs_tol=1e-9)
                for x14, x in zip(b14, bcoef, strict=True)
            )
        ):
            raise ValueError(
                f"Amber/Chamber topology contains distinct 1-4 Lennard-Jones tables "
                f"(LENNARD_JONES_14_ACOEF/BCOEF) unsupported by {target} export"
            )
    for atom in structure.atoms:
        raw_r14 = getattr(atom, "rmin_14", None)
        raw_e14 = getattr(atom, "epsilon_14", None)
        if raw_r14 is not None or raw_e14 is not None:
            r14 = float(atom.rmin if raw_r14 is None else raw_r14)
            e14 = float(atom.epsilon if raw_e14 is None else raw_e14)
            if not (math.isfinite(r14) and math.isfinite(e14)) or not (
                math.isclose(r14, float(atom.rmin), rel_tol=1e-6, abs_tol=1e-9)
                and math.isclose(e14, float(atom.epsilon), rel_tol=1e-6, abs_tol=1e-9)
            ):
                raise ValueError(
                    f"Amber/Chamber topology contains distinct per-atom 1-4 Lennard-Jones parameters "
                    f"(rmin_14/epsilon_14) unsupported by {target} export"
                )
    if pointers is not None and nb_index is not None and acoef is not None and bcoef is not None and len(pointers) > 1:
        ntypes = int(pointers[1])
        hba = parm_data.get("HBOND_ACOEF") or []
        hbb = parm_data.get("HBOND_BCOEF") or []
        nb_map: dict[int, tuple[float, float]] = {}
        for atom in structure.atoms:
            nb_i = int(getattr(atom, "nb_idx", 0) or 0)
            if nb_i <= 0:
                continue
            r_eps = (float(atom.rmin), float(atom.epsilon))
            if nb_i in nb_map:
                prev = nb_map[nb_i]
                if not (math.isclose(prev[0], r_eps[0], rel_tol=1e-6, abs_tol=1e-9) and math.isclose(prev[1], r_eps[1], rel_tol=1e-6, abs_tol=1e-9)):
                    raise ValueError(f"Inconsistent LJ parameters within Amber nb_idx {nb_i}")
            else:
                nb_map[nb_i] = r_eps
        active = sorted(nb_map)
        for idx_pos, nb_i in enumerate(active):
            r_i, e_i = nb_map[nb_i]
            for nb_j in active[idx_pos:]:
                r_j, e_j = nb_map[nb_j]
                eps_lb = math.sqrt(e_i * e_j)
                ptr = int(nb_index[ntypes * (nb_i - 1) + (nb_j - 1)])
                ptr_sym = int(nb_index[ntypes * (nb_j - 1) + (nb_i - 1)])
                if ptr != ptr_sym:
                    raise ValueError(
                        f"Asymmetric Amber NONBONDED_PARM_INDEX entry for types ({nb_i}, {nb_j}) unsupported by {target} export"
                    )
                if ptr < 0:
                    hb_i = -ptr - 1
                    ha = float(hba[hb_i]) if hb_i < len(hba) else 0.0
                    hb = float(hbb[hb_i]) if hb_i < len(hbb) else 0.0
                    if abs(ha) > 0.0 or abs(hb) > 0.0:
                        raise ValueError(
                            f"Amber topology contains nonzero 10-12 HBOND terms unsupported by {target} export"
                        )
                    if eps_lb > 0.0:
                        raise ValueError(
                            f"Amber topology disables 12-6 LJ via 10-12 HBOND index on an active LJ pair; "
                            f"{target} export would silently introduce a nonzero LJ interaction"
                        )
                    continue
                if ptr == 0:
                    if eps_lb > 0.0:
                        raise ValueError(
                            f"Amber topology omits 12-6 LJ pair entry on an active LJ pair; "
                            f"{target} export would silently introduce a nonzero LJ interaction"
                        )
                    continue
                r_lb = r_i + r_j
                a_lb = eps_lb * (r_lb**12)
                b_lb = 2.0 * eps_lb * (r_lb**6)
                a_p = float(acoef[ptr - 1])
                b_p = float(bcoef[ptr - 1])
                if not (
                    math.isclose(a_p, a_lb, rel_tol=1e-4, abs_tol=1e-6)
                    and math.isclose(b_p, b_lb, rel_tol=1e-4, abs_tol=1e-6)
                ):
                    raise ValueError(
                        f"Amber topology nonbonded pair coefficients deviate from Lorentz-Berthelot "
                        f"combining rules (geometric mixing or pair-specific LJ overrides unsupported by {target} export)"
                    )


def _verify_supported_amber_terms(
    structure, *, target: str, allow_explicit_exclusions: bool
) -> None:
    """Reject Amber topology terms that cannot be represented faithfully by the exporter."""
    import parmed as pmd

    _verify_standard_lorentz_berthelot_lj(structure, target=target)
    parm_data = getattr(structure, "parm_data", None) or {}
    if any(int(x) != 0 for x in (parm_data.get("IPOL") or [])):
        raise ValueError(
            f"Amber topology enables polarization (IPOL != 0), unsupported by {target} export"
        )
    pointers = parm_data.get("POINTERS") or []
    ptr_dict = getattr(structure, "pointers", None) or {}
    if (
        int(ptr_dict.get("IFPERT", 0) or 0) != 0
        or any(len(pointers) > idx and int(pointers[idx]) != 0 for idx in range(20, 27))
        or any(
            str(k).startswith("PERT_") and bool(np.asarray(v).size > 0)
            for k, v in parm_data.items()
        )
    ):
        raise ValueError(
            f"Amber topology enables perturbation (IFPERT != 0), unsupported by {target} export"
        )
    if (
        (len(pointers) > 29 and int(pointers[29]) != 0)
        or bool(parm_data.get("CAP_INFO"))
        or bool(parm_data.get("CAP_INFO2"))
    ):
        raise ValueError(
            f"Amber topology contains solvent cap (CAP_INFO) unsupported by {target} export"
        )
    if any(isinstance(a, pmd.ExtraPoint) for a in structure.atoms):
        raise ValueError(
            f"Amber topology contains virtual sites / ExtraPoints unsupported by {target} export"
        )
    cmap_count = [int(x) for x in (parm_data.get("CMAP_COUNT") or [])]
    charmm_cmap_count = [int(x) for x in (parm_data.get("CHARMM_CMAP_COUNT") or [])]
    if (
        getattr(structure, "cmaps", None)
        or any(x > 0 for x in cmap_count)
        or any(x > 0 for x in charmm_cmap_count)
        or bool(parm_data.get("CMAP_INDEX"))
        or bool(parm_data.get("CHARMM_CMAP_INDEX"))
    ):
        raise ValueError(
            f"Amber topology contains CMAP torsion terms unsupported by {target} export"
        )
    unsupported_attrs = (
        "urey_bradleys",
        "impropers",
        "rb_torsions",
        "trigonal_angles",
        "out_of_plane_bends",
        "pi_torsions",
        "stretch_bends",
        "torsion_torsions",
        "multipole_frames",
        "chiral_frames",
        "drude_pairs",
        "adjusts",
    )
    for attr in unsupported_attrs:
        if getattr(structure, attr, None):
            raise ValueError(
                f"Amber topology contains {attr} terms unsupported by {target} export"
            )
    excluded_12_13, pairs_14, extra_exclusions, unexcluded_14, orig_zero_excl, graph_d = (
        _extract_amber_exclusions(structure)
    )
    if not allow_explicit_exclusions and extra_exclusions:
        raise ValueError(
            f"Amber topology contains {len(extra_exclusions)} explicit nonbonded exclusion(s) "
            f"unsupported by {target} export"
        )
    if not allow_explicit_exclusions and unexcluded_14:
        raise ValueError(
            f"Amber topology leaves {len(unexcluded_14)} 1-4 pair(s) unexcluded; "
            f"{target} nrexcl=3 export cannot disable 1-4 scaling on unexcluded pairs"
        )
    if not allow_explicit_exclusions:
        g123 = graph_d[1] | graph_d[2] | graph_d[3]
        missing_graph_excl = (g123 - set(pairs_14)) - orig_zero_excl
        if missing_graph_excl:
            raise ValueError(
                f"Amber topology leaves {len(missing_graph_excl)} bond-graph <=3 pair(s) unexcluded; "
                f"{target} nrexcl=3 export automatically excludes all 1-2/1-3/1-4 graph pairs"
            )
        extra_beyond_3 = orig_zero_excl - g123
        if extra_beyond_3:
            raise ValueError(
                f"Amber topology contains {len(extra_beyond_3)} exclusion(s) beyond bond-graph distance 3, "
                f"unsupported by {target} nrexcl=3 export"
            )
        non_top_14 = set(pairs_14) - graph_d[3]
        if non_top_14:
            raise ValueError(
                f"Amber topology contains {len(non_top_14)} 1-4 pair(s) not at bond-graph distance 3, "
                f"unsupported by {target} export"
            )


def _reduce_openmm_box_vectors(vecs_nm: np.ndarray) -> list[list[float]]:
    """Convert periodic box vectors to OpenMM reduced lower-triangular right-handed form."""
    arr = np.asarray(vecs_nm, dtype=float)
    if arr.shape != (3, 3) or not np.isfinite(arr).all():
        raise ValueError("Invalid periodic box vectors in topology")
    if abs(float(np.linalg.det(arr))) <= 1e-12:
        raise ValueError("Degenerate zero-volume periodic box vectors in topology")

    # Preserve metric tensor (lengths and angles) in OpenMM lower-triangular right-handed orientation.
    if (
        abs(arr[0, 1]) > 1e-10
        or abs(arr[0, 2]) > 1e-10
        or abs(arr[1, 2]) > 1e-10
        or arr[0, 0] <= 0.0
        or arr[1, 1] <= 0.0
        or arr[2, 2] <= 0.0
    ):
        gram = arr @ arr.T
        low = np.linalg.cholesky(gram)
    else:
        low = arr.copy()
        low[0, 1] = 0.0
        low[0, 2] = 0.0
        low[1, 2] = 0.0

    a, b, c = low[0].copy(), low[1].copy(), low[2].copy()
    for _ in range(4):
        c = c - b * round(float(c[1] / b[1]))
        c = c - a * round(float(c[0] / a[0]))
        b = b - a * round(float(b[0] / a[0]))
    a[1] = a[2] = b[2] = 0.0
    if not (
        a[0] > 0.0
        and b[1] > 0.0
        and c[2] > 0.0
        and a[0] >= 2.0 * abs(b[0]) - 1e-9
        and a[0] >= 2.0 * abs(c[0]) - 1e-9
        and b[1] >= 2.0 * abs(c[1]) - 1e-9
    ):
        raise ValueError("Periodic box vectors cannot be brought into OpenMM reduced form")
    return [[float(x) for x in a], [float(x) for x in b], [float(x) for x in c]]


def _amber_parm_to_openmm_xml(
    structure,
    *,
    nonbonded_method: str = "NoCutoff",
    cutoff_nm: float = 1.0,
) -> tuple[str, dict[str, int]]:
    """Serialize an AmberParm directly into an OpenMM <System> XML string."""
    method_codes = {"NoCutoff": 0, "CutoffNonPeriodic": 1, "PME": 4}
    if nonbonded_method not in method_codes:
        raise ValueError(f"Unsupported nonbonded_method {nonbonded_method!r}")
    if not math.isfinite(float(cutoff_nm)) or float(cutoff_nm) <= 0.0:
        raise ValueError("cutoff_nm must be positive and finite")
    _verify_supported_amber_terms(
        structure, target="OpenMM XML", allow_explicit_exclusions=True
    )
    if nonbonded_method == "PME" and structure.box is None:
        raise ValueError("PME requires periodic box vectors in the source Amber topology/coordinates")

    # Match OpenMM AmberPrmtopFile.createSystem: bonded forces and 1-4 exceptions use direct distances ("0")
    use_pbc = "0"
    if structure.box is not None and getattr(structure, "box_vectors", None) is not None:
        import parmed as pmd

        vecs = structure.box_vectors.value_in_unit(pmd.unit.nanometers)
        box_rows = _reduce_openmm_box_vectors(np.asarray(vecs, dtype=float))
        if nonbonded_method == "PME":
            mat = np.asarray(box_rows, dtype=float)
            vol = abs(float(np.linalg.det(mat)))
            heights = [
                vol / float(np.linalg.norm(np.cross(mat[(i + 1) % 3], mat[(i + 2) % 3])))
                for i in range(3)
            ]
            max_allowed = 0.5 * min(box_rows[0][0], box_rows[1][1], box_rows[2][2], *heights)
            if float(cutoff_nm) >= max_allowed:
                raise ValueError(
                    f"PME cutoff_nm={cutoff_nm} must be strictly less than half the periodic box "
                    f"dimension/height ({max_allowed:.6f} nm)"
                )
    else:
        box_rows = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]

    root = ET.Element("System", {"openmmVersion": "8.0", "type": "System", "version": "1"})
    box_el = ET.SubElement(root, "PeriodicBoxVectors")
    for label, row in zip(("A", "B", "C"), box_rows, strict=True):
        ET.SubElement(
            box_el,
            label,
            {"x": _fmt_xml_float(row[0]), "y": _fmt_xml_float(row[1]), "z": _fmt_xml_float(row[2])},
        )

    particles_el = ET.SubElement(root, "Particles")
    for atom in structure.atoms:
        ET.SubElement(particles_el, "Particle", {"mass": _fmt_xml_float(atom.mass)})
    ET.SubElement(root, "Constraints")

    forces_el = ET.SubElement(root, "Forces")

    bond_force = ET.SubElement(
        forces_el,
        "Force",
        {
            "forceGroup": "0",
            "name": "HarmonicBondForce",
            "type": "HarmonicBondForce",
            "usesPeriodic": use_pbc,
            "version": "2",
        },
    )
    bonds_el = ET.SubElement(bond_force, "Bonds")
    for bond in structure.bonds:
        if bond.type is None:
            raise ValueError("Untyped bond encountered during OpenMM XML export")
        ET.SubElement(
            bonds_el,
            "Bond",
            {
                "d": _fmt_xml_float(float(bond.type.req) * 0.1),
                "k": _fmt_xml_float(float(bond.type.k) * 836.8),
                "p1": str(bond.atom1.idx),
                "p2": str(bond.atom2.idx),
            },
        )

    angle_force = ET.SubElement(
        forces_el,
        "Force",
        {
            "forceGroup": "0",
            "name": "HarmonicAngleForce",
            "type": "HarmonicAngleForce",
            "usesPeriodic": use_pbc,
            "version": "2",
        },
    )
    angles_el = ET.SubElement(angle_force, "Angles")
    for angle in structure.angles:
        if angle.type is None:
            raise ValueError("Untyped angle encountered during OpenMM XML export")
        ET.SubElement(
            angles_el,
            "Angle",
            {
                "a": _fmt_xml_float(math.radians(float(angle.type.theteq))),
                "k": _fmt_xml_float(float(angle.type.k) * 8.368),
                "p1": str(angle.atom1.idx),
                "p2": str(angle.atom2.idx),
                "p3": str(angle.atom3.idx),
            },
        )

    torsion_force = ET.SubElement(
        forces_el,
        "Force",
        {
            "forceGroup": "0",
            "name": "PeriodicTorsionForce",
            "type": "PeriodicTorsionForce",
            "usesPeriodic": use_pbc,
            "version": "2",
        },
    )
    torsions_el = ET.SubElement(torsion_force, "Torsions")
    torsion_count = 0
    for dihedral in structure.dihedrals:
        types = dihedral.type if isinstance(dihedral.type, list) else [dihedral.type]
        for dt in types:
            if dt is None:
                raise ValueError("Untyped dihedral encountered during OpenMM XML export")
            per_f = float(dt.per)
            if not math.isfinite(per_f) or per_f < 0.0 or not per_f.is_integer():
                raise ValueError(f"Invalid torsion periodicity {dt.per!r} during OpenMM XML export")
            per = int(per_f)
            phi_k = float(dt.phi_k)
            phase_f = float(dt.phase)
            if not math.isfinite(phi_k) or not math.isfinite(phase_f):
                raise ValueError("Nonfinite torsion parameters on dihedral during OpenMM XML export")
            for scale_name in ("scee", "scnb"):
                raw_scale = getattr(dt, scale_name, 0.0)
                if raw_scale is not None and not math.isfinite(float(raw_scale)):
                    raise ValueError(
                        f"Nonfinite torsion {scale_name}={raw_scale!r} during OpenMM XML export"
                    )
            if per == 0 and phi_k == 0.0:
                continue
            if per == 0:
                raise ValueError(
                    f"Invalid torsion periodicity 0 with nonzero force constant {phi_k!r} during OpenMM XML export"
                )
            ET.SubElement(
                torsions_el,
                "Torsion",
                {
                    "k": _fmt_xml_float(phi_k * 4.184),
                    "p1": str(dihedral.atom1.idx),
                    "p2": str(dihedral.atom2.idx),
                    "p3": str(dihedral.atom3.idx),
                    "p4": str(dihedral.atom4.idx),
                    "periodicity": str(per),
                    "phase": _fmt_xml_float(math.radians(float(dt.phase))),
                },
            )
            torsion_count += 1

    nb_force = ET.SubElement(
        forces_el,
        "Force",
        {
            "alpha": "0",
            "cutoff": _fmt_xml_float(cutoff_nm),
            "dispersionCorrection": "1",
            "ewaldTolerance": "0.0005",
            "exceptionsUsePeriodic": use_pbc,
            "forceGroup": "0",
            "includeDirectSpace": "1",
            "ljAlpha": "0",
            "ljnx": "0",
            "ljny": "0",
            "ljnz": "0",
            "method": str(method_codes[nonbonded_method]),
            "name": "NonbondedForce",
            "nx": "0",
            "ny": "0",
            "nz": "0",
            "recipForceGroup": "-1",
            "rfDielectric": "78.3",
            "switchingDistance": "-1",
            "type": "NonbondedForce",
            "useSwitchingFunction": "0",
            "version": "4",
        },
    )
    ET.SubElement(nb_force, "GlobalParameters")
    ET.SubElement(nb_force, "ParticleOffsets")
    ET.SubElement(nb_force, "ExceptionOffsets")

    nb_particles = ET.SubElement(nb_force, "Particles")
    q_list: list[float] = []
    sig_list: list[float] = []
    eps_list: list[float] = []
    scale_sig = (2.0 ** (5.0 / 6.0)) * 0.1
    for atom in structure.atoms:
        q = float(atom.charge)
        eps = float(atom.epsilon) * 4.184
        rmin = float(atom.rmin)
        sig = 0.0 if rmin == 0.0 else rmin * scale_sig
        q_list.append(q)
        sig_list.append(sig)
        eps_list.append(eps)
        ET.SubElement(
            nb_particles,
            "Particle",
            {
                "eps": _fmt_xml_float(eps),
                "q": _fmt_xml_float(q),
                "sig": _fmt_xml_float(sig),
            },
        )

    excluded_12_13, pairs_14, extra_exclusions, _, _, _ = _extract_amber_exclusions(structure)

    exceptions_el = ET.SubElement(nb_force, "Exceptions")
    all_exceptions: dict[tuple[int, int], tuple[float, float, float]] = {}
    for pair in excluded_12_13 | extra_exclusions:
        all_exceptions[pair] = (0.0, 0.0, 0.5)
    for pair, (scee, scnb) in pairs_14.items():
        i, j = pair
        q12 = (q_list[i] * q_list[j]) / scee
        sig12 = 0.5 * (sig_list[i] + sig_list[j])
        eps12 = math.sqrt(eps_list[i] * eps_list[j]) / scnb
        if eps12 == 0.0 and sig12 == 0.0:
            sig12 = 0.5
        all_exceptions[pair] = (eps12, q12, sig12)

    for (i, j), (eps12, q12, sig12) in sorted(all_exceptions.items()):
        ET.SubElement(
            exceptions_el,
            "Exception",
            {
                "eps": _fmt_xml_float(eps12),
                "p1": str(i),
                "p2": str(j),
                "q": _fmt_xml_float(q12),
                "sig": _fmt_xml_float(sig12),
            },
        )

    ET.indent(root, space="  ")
    xml_text = '<?xml version="1.0" ?>\n' + ET.tostring(root, encoding="unicode") + "\n"
    counts = {
        "atoms": len(structure.atoms),
        "bonds": len(structure.bonds),
        "angles": len(structure.angles),
        "torsions": torsion_count,
        "exceptions": len(all_exceptions),
        "excluded_12_13": len(excluded_12_13),
        "explicit_exclusions": len(extra_exclusions),
        "scaled_14": len(pairs_14),
    }
    return xml_text, counts


def export_openmm(
    topology,
    coordinates,
    output_xml,
    *,
    nonbonded_method="NoCutoff",
    cutoff_nm=1.0,
    overwrite=False,
    model=None,
    _force_pure_python=False,
):
    """Export an Amber prmtop/inpcrd topology into a standalone OpenMM System XML."""
    if model and (_periodic_edges(model) or model.get("metadata", {}).get("covalently_periodic")):
        raise ValueError(
            "Periodic winding OpenMM covalent export not implemented/validated"
        )
    xml_path = Path(output_xml).resolve()
    paths = {
        "xml": xml_path,
        "manifest": xml_path.with_suffix(".manifest.json"),
        "provenance": xml_path.with_suffix(".provenance.json"),
    }
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("OpenMM export files already exist; choose a fresh output path")
    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh output path")
    structure, top_sha, crd_sha = _load_amber_with_digests(topology, coordinates)
    if not structure.atoms:
        raise ValueError("Amber topology contains no atoms")
    upstream = _link_upstream_provenance(
        topology, coordinates, top_sha=top_sha, crd_sha=crd_sha
    )

    xml_text, counts = _amber_parm_to_openmm_xml(
        structure,
        nonbonded_method=nonbonded_method,
        cutoff_nm=cutoff_nm,
    )
    serializer = "pure_python_amber_xml"
    if not _force_pure_python:
        try:
            import openmm

            system = openmm.XmlSerializer.deserialize(xml_text)
            if system.getNumParticles() != len(structure.atoms):
                raise RuntimeError("OpenMM XML round-trip particle count mismatch")
            serializer = f"pure_python_amber_xml+openmm_{openmm.__version__}_verified"
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
        with tempfile.TemporaryDirectory(
            dir=xml_path.parent, prefix=f".{xml_path.name}.tmp."
        ) as stage_dir_str:
            stage_dir = Path(stage_dir_str)
            staged = {key: stage_dir / fp.name for key, fp in paths.items()}
            xml_bytes = xml_text.encode("utf-8")
            staged["xml"].write_bytes(xml_bytes)
            prov = {
                "schema_version": 1,
                "format": "openmm_system_xml",
                **upstream,
                "nonbonded_method": nonbonded_method,
                "cutoff_nm": float(cutoff_nm),
                "serializer": serializer,
                "counts": counts,
                "outputs": {k: str(v) for k, v in paths.items() if k != "provenance"},
            }
            _json(staged["provenance"], prov)
            result = {
                "status": "exported",
                "parameterized": True,
                "output": str(xml_path),
                "xml": str(xml_path),
                "manifest": str(paths["manifest"]),
                "provenance": str(paths["provenance"]),
                "serializer": serializer,
                "nonbonded_method": nonbonded_method,
                "cutoff_nm": float(cutoff_nm),
                "source_topology_sha256": upstream["source_topology_sha256"],
                "xml_sha256": hashlib.sha256(xml_bytes).hexdigest(),
                "preset": upstream["preset"],
                "charge_model": upstream["charge_model"],
                "atom_count": counts["atoms"],
                "bond_count": counts["bonds"],
                "angle_count": counts["angles"],
                "torsion_count": counts["torsions"],
                "exception_count": counts["exceptions"],
                "periodic_policy": "finite solute unless solvated/boxed PME prmtop supplied",
                "gbsa_note": "Standard explicit NonbondedForce XML does not include implicit GBSA force; source prmtop holds GB radii",
                "engine_sources": ENGINE_SOURCES,
            }
            _json(staged["manifest"], result)
            for key in ("xml", "provenance", "manifest"):
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


def gbsa_readiness(model):
    """Diagnose numerical GB prerequisites and explicit noncovalent partitions.

    ready means topology/radius prerequisites, never scientific validation of a
    particular GO charge model. Binding readiness additionally requires ligand
    selection and no bond across the receptor/ligand partition.
    """
    blockers = []
    warnings = []
    atoms = model.get("atoms", [])
    metadata = model.get("metadata", {})
    if not atoms:
        blockers.append("No atoms")
    if _periodic_edges(model) or metadata.get("covalently_periodic"):
        blockers.append("Covalently periodic sheet is not a finite GBSA solute")
    if (
        metadata.get("requires_parameterization")
        or any(a.get("metadata", {}).get("junction") for a in atoms)
        and not metadata.get("junction_template")
    ):
        blockers.append("New covalent junction lacks dedicated parameterization")
    if not metadata.get("parameterized"):
        blockers.append("No verified fully parameterized topology")
    if any(a.get("charge") is None or not math.isfinite(float(a["charge"])) for a in atoms):
        blockers.append("Missing finite partial charge(s)")
    if any(not a.get("atom_type") for a in atoms):
        blockers.append("Missing force-field atom type(s)")
    if any(
        not isinstance(a.get("gb_radius"), (float, int))
        or a["gb_radius"] <= 0
        or not math.isfinite(a["gb_radius"])
        or not isinstance(a.get("gb_screen"), (float, int))
        or a["gb_screen"] <= 0
        or not math.isfinite(a["gb_screen"])
        for a in atoms
    ):
        blockers.append("Missing positive GB radii/screens; set compatible radii in Amber")
    if metadata.get("wrapped_or_stretched_bonds"):
        blockers.append("Wrapped or stretched bonded coordinates must be made whole before GBSA")
    mapping = _residue_map(model)
    solvent = {"WAT", "HOH", "SOL", "Na+", "Cl-", "NA", "CL", "K+"}
    ligand_resids = set(metadata.get("ligand_resids", []))
    receptor_resids = set(metadata.get("receptor_resids", []))
    receptor = []
    ligand = []
    solvent_indices = []
    for (resid, name), index in mapping.items():
        components = {
            a.get("metadata", {}).get("component")
            for a in atoms
            if a["resid"] == resid and a["resname"] == name
        }
        if name in solvent or "solvent" in components or "ion" in components:
            solvent_indices.append(index)
        elif resid in ligand_resids or "ligand" in components:
            ligand.append(index)
        elif not receptor_resids or resid in receptor_resids or "receptor" in components:
            receptor.append(index)

    def mask(values):
        if not values:
            return None
        values = sorted(set(values))
        ranges = []
        start = end = values[0]
        for value in values[1:]:
            if value == end + 1:
                end = value
            else:
                ranges.append(str(start) if start == end else f"{start}-{end}")
                start = end = value
        ranges.append(str(start) if start == end else f"{start}-{end}")
        return ":" + ",".join(ranges)

    if not receptor:
        blockers.append("No receptor atoms selected")
    atom_component = {a["id"]: mapping[(a["resid"], a["resname"])] for a in atoms}
    cross_bonds = any(
        (atom_component[b["a"]] in receptor and atom_component[b["b"]] in ligand)
        or (atom_component[b["b"]] in receptor and atom_component[b["a"]] in ligand)
        for b in model.get("bonds", [])
    )
    if cross_bonds:
        warnings.append(
            "Receptor and ligand are covalently linked: ordinary binding MM/GBSA partition is invalid"
        )
    if not ligand:
        warnings.append("No ligand partition supplied; receptor preparation only")
    warnings.append("Transferred GO fragment charges require whole-sheet electrostatic validation")
    return {
        "ready": not blockers,
        "binding_ready": not blockers and bool(ligand) and not cross_bonds,
        "blockers": blockers,
        "warnings": warnings,
        "radii_set": metadata.get("radii_set"),
        "masks": {
            "receptor": mask(receptor),
            "ligand": mask(ligand),
            "complex": mask(sorted(receptor + ligand)),
            "strip": mask(solvent_indices),
        },
        "net_charge": sum(a["charge"] for a in atoms)
        if atoms and all(a.get("charge") is not None for a in atoms)
        else None,
    }
