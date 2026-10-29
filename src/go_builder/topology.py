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
from copy import deepcopy
from pathlib import Path

import numpy as np

from .ff import DEFAULT_REFERENCE, assign_current_fragments, normalize_libraries

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


def prepare_amber(
    model,
    output_dir,
    *,
    run=False,
    tleap="tleap",
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
    reference_dir = Path(reference_dir or DEFAULT_REFERENCE)
    assigned = assign_current_fragments(
        model, reference_dir, parameter_templates=parameter_templates
    )
    executable = shutil.which(str(tleap)) if run else None
    if run and executable is None:
        raise ValueError(f"tleap executable unavailable: {tleap}")
    source_frcmod = reference_dir / "GO_gaff2.frcmod"
    if not source_frcmod.is_file():
        raise ValueError("Current GO_gaff2.frcmod is missing")
    output_dir.mkdir(parents=True)
    audit = normalize_libraries(reference_dir, output_dir / "normalized_libraries")
    shutil.copyfile(source_frcmod, output_dir / "GO_gaff2.frcmod")
    extra = []
    if parameter_templates:
        source = Path(parameter_templates["frcmod"])
        shutil.copyfile(source, output_dir / "junction.frcmod")
        extra = ["loadAmberParams junction.frcmod"]
    residues = _residue_map(assigned)
    assigned["atoms"].sort(key=lambda atom: residues[(atom["resid"], atom["resname"])])
    lines = [
        "source leaprc.gaff2",
        "loadAmberParams GO_gaff2.frcmod",
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
    result = {
        "status": "inputs_prepared",
        "parameterized": False,
        "leap_input": str(leap_input),
        "assigned_model": str(output_dir / "assigned_model.json"),
        "library_audit": str(output_dir / "library_audit.json"),
        "topology": None,
        "coordinates": None,
        "radii_set": radii,
        "atom_mapping": {str(i): a["id"] for i, a in enumerate(assigned["atoms"], 1)},
        "atom_count": len(assigned["atoms"]),
        "bond_count": len(assigned["bonds"]),
        "net_charge": sum(a["charge"] for a in assigned["atoms"]),
        "go_frcmod_sha256": _sha(source_frcmod),
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
            atol=1e-7,
        ):
            raise ValueError("LEaP altered assigned fragment charges")
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
        ]
    }
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("GROMACS export files already exist; choose a fresh prefix")
    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh prefix")
    structure = _load_amber(topology, coordinates)
    import parmed as pmd

    converted = pmd.gromacs.GromacsTopologyFile.from_structure(structure)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    converted.write(str(paths["top"]), molfile=str(paths["itp"]), combine="all")
    text = paths["top"].read_text().replace(f'"{paths["itp"]}"', f'"{paths["itp"].name}"')
    paths["top"].write_text(text)
    structure.save(str(paths["gro"]), overwrite=False)
    roundtrip = pmd.load_file(str(paths["top"]), xyz=str(paths["gro"]))
    verification = _verify_parameter_roundtrip(structure, roundtrip)
    result = {key: str(value) for key, value in paths.items()}
    result.update(
        status="exported",
        parameterized=True,
        roundtrip_verified=verification,
        source_topology_sha256=_sha(topology),
        atom_count=len(structure.atoms),
        bond_count=len(structure.bonds),
        periodic_molecules="no (finite covalent topology only)",
        gbsa_note="GB radii/screens remain in source Amber topology; standard GROMACS export does not transfer these fields",
        engine_sources=ENGINE_SOURCES,
    )
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
