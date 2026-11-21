"""Amber-first topology preparation and parameter-preserving ParmEd exchange.

Finite covalent graphs may use a periodic solvent box. Bonds carrying nonzero
image shifts are rejected: this exporter does not provide a validated periodic
bonded-force implementation for a covalently winding sheet.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

import numpy as np

from .ff import DEFAULT_REFERENCE, assign_current_fragments, normalize_libraries
from .ff_v3 import (
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
) -> Path:
    ff_meta = assigned.get("metadata", {}).get("force_field", {})
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
        "files": {
            "frcmod": frcmod_name,
            "frcmod_sha256": frcmod_sha256,
            "bundle_provenance_sha256": ff_meta.get("provenance_sha256"),
            "prmtop_sha256": _sha(topology_path) if topology_path else None,
            "inpcrd_sha256": _sha(coordinates_path) if coordinates_path else None,
        },
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
        source_frcmod = Path(bundle["frcmod_path"])
        assigned = assign_gaff2_v3(
            model,
            reference_dir=reference_dir,
            charge_model=charge_model,
            custom_charges=charges_json,
            parameter_templates=parameter_templates,
        )
        output_dir.mkdir(parents=True)
        v3_lib_out = output_dir / "v3_libraries"
        v3_lib_out.mkdir(parents=True)
        audit = {}
        ref_dir = Path(bundle["reference_dir"])
        for lib_file in sorted(ref_dir.glob("*_gaff2_v3.lib")):
            dest = v3_lib_out / lib_file.name
            shutil.copyfile(lib_file, dest)
            audit[lib_file.stem.split("_")[0]] = {
                "source": str(lib_file),
                "sha256": _sha(dest),
            }
        shutil.copyfile(source_frcmod, output_dir / frcmod_name)
    else:
        cur_dir = Path(reference_dir or DEFAULT_REFERENCE)
        frcmod_name = "GO_gaff2.frcmod"
        source_frcmod = cur_dir / frcmod_name
        if not source_frcmod.is_file():
            raise ValueError("Current GO_gaff2.frcmod is missing")
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
        shutil.copyfile(source_frcmod, output_dir / frcmod_name)

    extra = []
    if parameter_templates:
        source = Path(parameter_templates["frcmod"])
        shutil.copyfile(source, output_dir / "junction.frcmod")
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
    frcmod_sha = _sha(source_frcmod)
    prov_path = _write_amber_provenance(
        output_dir,
        preset=resolved_preset,
        frcmod_name=frcmod_name,
        frcmod_sha256=frcmod_sha,
        assigned=assigned,
        radii=radii,
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
        _write_amber_provenance(
            output_dir,
            preset=resolved_preset,
            frcmod_name=frcmod_name,
            frcmod_sha256=frcmod_sha,
            assigned=assigned,
            radii=radii,
            topology_path=top,
            coordinates_path=crd,
        )
        result.update(
            status="parameterized",
            parameterized=True,
            topology=str(top),
            coordinates=str(crd),
            topology_sha256=_sha(top),
            gbsa=gbsa_readiness(imported),
        )
    _json(output_dir / "manifest.json", result)
    return result


def _load_amber(topology, coordinates=None):
    import parmed as pmd

    structure = pmd.load_file(
        str(topology), **({"xyz": str(coordinates)} if coordinates is not None else {})
    )
    if not isinstance(structure, pmd.amber.AmberParm):
        raise ValueError("An existing fully parameterized Amber prmtop is required")
    if structure.coordinates is None:
        raise ValueError("Amber coordinates are required; coordinates are never fabricated")
    _check_parameterized(structure)
    return structure


def _check_parameterized(structure):
    for field in ("bonds", "angles", "dihedrals", "impropers"):
        if any(term.type is None for term in getattr(structure, field)):
            raise ValueError(f"Unparameterized {field} in topology")
    if any(
        not math.isfinite(a.charge)
        or a.mass <= 0
        or not math.isfinite(a.epsilon)
        or not math.isfinite(a.rmin)
        for a in structure.atoms
    ):
        raise ValueError("Missing or invalid atom charge/mass/Lennard-Jones parameters")
    if not np.isfinite(structure.coordinates).all():
        raise ValueError("Nonfinite imported coordinates")


def annotate_imported_epoxides(model):
    """Annotate intact current EPX basal motifs from existing explicit bonds.

    This recognizes a known named template and its complete local graph. It
    never assigns bonds from distances, guesses arbitrary oxygen chemistry,
    or treats edge epoxides as supported basal graft sites.
    """
    model = deepcopy(model)
    by_id = {a["id"]: a for a in model["atoms"]}
    neighbors = {i: set() for i in by_id}
    for bond in model["bonds"]:
        neighbors[bond["a"]].add(bond["b"])
        neighbors[bond["b"]].add(bond["a"])
    grouped = {}
    for atom in model["atoms"]:
        if atom["resname"] == "EPX":
            grouped.setdefault(atom["resid"], []).append(atom)
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

    structure = _load_amber(topology, coordinates)
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
            "source_topology_sha256": _sha(topology),
            "source_coordinates": str(Path(coordinates).resolve()) if coordinates else None,
            "radii_set": " ".join(structure.parm_data.get("RADIUS_SET", [])),
            "bond_image_policy": "Amber endpoints imported verbatim; no bond images in prmtop; no distance inference",
            "wrapped_or_stretched_bonds": stretched,
        },
    }
    return annotate_imported_epoxides(result)


def _verify_parameter_roundtrip(original, converted):
    """Verify numerical atoms, valence terms and 1-4 scaling after conversion."""

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
                        (
                            t.phi_k,
                            t.per,
                            t.phase,
                            0.0 if term.improper else t.scee,
                            0.0 if term.improper else t.scnb,
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
    summary = {"atoms": len(original.atoms)}
    for field in ("bonds", "angles", "dihedrals", "impropers"):
        before, after = valence(original, field), valence(converted, field)
        if [x[0] for x in before] != [x[0] for x in after] or (
            before
            and not np.allclose([x[1] for x in before], [x[1] for x in after], atol=2e-6, rtol=2e-6)
        ):
            raise ValueError(f"GROMACS round-trip {field} parameter/scaling mismatch")
        summary["torsion_terms" if field == "dihedrals" else field] = len(before)
    # ParmEd expands standard Amber 1-4 interactions into explicit GROMACS pairs;
    # proper torsion scee/scnb above verify their retained scaling factors.
    summary["coordinate_tolerance_angstrom"] = 0.006
    if not np.allclose(original.coordinates, converted.coordinates, atol=0.006, rtol=0):
        raise ValueError("GROMACS round-trip coordinate mismatch beyond GRO precision")
    return summary


def _link_upstream_provenance(topology, coordinates):
    top_path = Path(topology).resolve()
    crd_path = Path(coordinates).resolve()
    top_sha = _sha(top_path)
    crd_sha = _sha(crd_path)
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
    structure = _load_amber(topology, coordinates)
    upstream = _link_upstream_provenance(topology, coordinates)
    import parmed as pmd

    converted = pmd.gromacs.GromacsTopologyFile.from_structure(structure)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    converted.write(str(paths["top"]), molfile=str(paths["itp"]), combine="all")
    text = paths["top"].read_text().replace(f'"{paths["itp"]}"', f'"{paths["itp"].name}"')
    paths["top"].write_text(text)
    bonds_upgraded = _upgrade_gromacs_precision(paths["top"], paths["itp"], structure, converted)
    structure.save(str(paths["gro"]), overwrite=False)
    roundtrip = pmd.load_file(str(paths["top"]), xyz=str(paths["gro"]))
    verification = _verify_parameter_roundtrip(structure, roundtrip)
    prov = {
        "schema_version": 1,
        "format": "gromacs",
        **upstream,
        "itp_bond_precision_upgraded": bonds_upgraded,
        "outputs": {k: str(v) for k, v in paths.items() if k != "provenance"},
    }
    _json(paths["provenance"], prov)
    result = {key: str(value) for key, value in paths.items()}
    result.update(
        status="exported",
        parameterized=True,
        roundtrip_verified=verification,
        source_topology_sha256=_sha(topology),
        itp_bond_precision_upgraded=bonds_upgraded,
        preset=upstream["preset"],
        charge_model=upstream["charge_model"],
        atom_count=len(structure.atoms),
        bond_count=len(structure.bonds),
        periodic_molecules="no (finite covalent topology only)",
        gbsa_note="GB radii/screens remain in source Amber topology; standard GROMACS export does not transfer these fields",
        engine_sources=ENGINE_SOURCES,
    )
    _json(paths["manifest"], result)
    return result


def _fmt_xml_float(value: float) -> str:
    val = float(value)
    if val == 0.0:
        return "0"
    return format(val, ".16g")


def _verify_standard_lorentz_berthelot_lj(structure) -> None:
    """Reject NBFIX, non-Lorentz mixing, 12-6-4 C4 terms, and 10-12 H-bond terms."""
    comb = str(getattr(structure, "combining_rule", "lorentz") or "lorentz").lower()
    if comb != "lorentz":
        raise ValueError(
            f"Unsupported Amber LJ combining_rule {comb!r}; standalone NonbondedForce XML requires Lorentz-Berthelot"
        )
    if hasattr(structure, "has_NBFIX") and structure.has_NBFIX():
        raise ValueError(
            "Amber topology contains pair-specific Lennard-Jones overrides (NBFIX); "
            "standalone NonbondedForce OpenMM XML export cannot represent NBFIX"
        )
    parm_data = getattr(structure, "parm_data", None) or {}
    ccoef = parm_data.get("LENNARD_JONES_CCOEF")
    if ccoef and any(abs(float(x)) > 0.0 for x in ccoef):
        raise ValueError(
            "Amber topology contains nonzero 12-6-4 LENNARD_JONES_CCOEF terms unsupported by standalone NonbondedForce XML"
        )
    pointers = parm_data.get("POINTERS")
    nb_index = parm_data.get("NONBONDED_PARM_INDEX")
    acoef = parm_data.get("LENNARD_JONES_ACOEF")
    bcoef = parm_data.get("LENNARD_JONES_BCOEF")
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
                if ptr < 0:
                    hb_i = -ptr - 1
                    ha = float(hba[hb_i]) if hb_i < len(hba) else 0.0
                    hb = float(hbb[hb_i]) if hb_i < len(hbb) else 0.0
                    if abs(ha) > 0.0 or abs(hb) > 0.0:
                        raise ValueError(
                            "Amber topology contains nonzero 10-12 HBOND terms unsupported by standalone NonbondedForce XML"
                        )
                    if eps_lb > 0.0:
                        raise ValueError(
                            "Amber topology disables 12-6 LJ via 10-12 HBOND index on an active LJ pair; "
                            "standalone NonbondedForce XML would silently introduce a nonzero LJ interaction"
                        )
                    continue
                if ptr == 0:
                    if eps_lb > 0.0:
                        raise ValueError(
                            "Amber topology omits 12-6 LJ pair entry on an active LJ pair; "
                            "standalone NonbondedForce XML would silently introduce a nonzero LJ interaction"
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
                        "Amber topology nonbonded pair coefficients deviate from Lorentz-Berthelot "
                        "combining rules (geometric mixing or pair-specific LJ overrides unsupported "
                        "in standalone NonbondedForce XML)"
                    )


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
    if getattr(structure, "impropers", None):
        raise ValueError("Unsupported CHARMM-style impropers in AmberParm")
    _verify_standard_lorentz_berthelot_lj(structure)
    if nonbonded_method == "PME" and structure.box is None:
        raise ValueError("PME requires periodic box vectors in the source Amber topology/coordinates")

    # Match OpenMM AmberPrmtopFile.createSystem: bonded forces and 1-4 exceptions use direct distances ("0")
    use_pbc = "0"
    if structure.box is not None and getattr(structure, "box_vectors", None) is not None:
        import parmed as pmd

        vecs = structure.box_vectors.value_in_unit(pmd.unit.nanometers)
        box_rows = [[float(c) for c in row] for row in vecs]
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
            per = int(round(abs(float(dt.per))))
            phi_k = float(dt.phi_k)
            if per == 0 and phi_k == 0.0:
                continue
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

    excluded_12_13: set[tuple[int, int]] = set()
    for bond in structure.bonds:
        i, j = bond.atom1.idx, bond.atom2.idx
        excluded_12_13.add((min(i, j), max(i, j)))
    for angle in structure.angles:
        i, j = angle.atom1.idx, angle.atom3.idx
        excluded_12_13.add((min(i, j), max(i, j)))

    pairs_14: dict[tuple[int, int], tuple[float, float]] = {}
    for dihedral in structure.dihedrals:
        if getattr(dihedral, "improper", False) or getattr(dihedral, "ignore_end", False):
            continue
        i, j = dihedral.atom1.idx, dihedral.atom4.idx
        if i == j:
            continue
        pair = (min(i, j), max(i, j))
        if pair in excluded_12_13:
            continue
        dt = dihedral.type[0] if isinstance(dihedral.type, list) else dihedral.type
        scee = float(dt.scee) if dt is not None and float(getattr(dt, "scee", 0.0)) > 0.0 else 1.2
        scnb = float(dt.scnb) if dt is not None and float(getattr(dt, "scnb", 0.0)) > 0.0 else 2.0
        pairs_14.setdefault(pair, (scee, scnb))

    exceptions_el = ET.SubElement(nb_force, "Exceptions")
    all_exceptions: dict[tuple[int, int], tuple[float, float, float]] = {}
    for pair in excluded_12_13:
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
    structure = _load_amber(topology, coordinates)
    if not structure.atoms:
        raise ValueError("Amber topology contains no atoms")
    upstream = _link_upstream_provenance(topology, coordinates)

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

    xml_path.parent.mkdir(parents=True, exist_ok=True)
    xml_path.write_text(xml_text, encoding="utf-8")
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
    _json(paths["provenance"], prov)
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
        "source_topology_sha256": _sha(topology),
        "xml_sha256": _sha(xml_path),
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
    _json(paths["manifest"], result)
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
