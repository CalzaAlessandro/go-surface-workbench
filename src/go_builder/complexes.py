"""Assemble complete, independently parameterized docking components.

All poses share one reference frame; checked atom mappings recover Amber order.
No fitting, bond inference, solvation, charge assignment or MD is performed.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import parmed as pmd
from parmed.tools.changeradii import ChRad

from .model import write_pdb
from .topology import _check_parameterized, _verify_parameter_roundtrip


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _mask(first, last):
    return f":{first}" if first == last else f":{first}-{last}"


# Element identities from the public Vina AutoDock4 atom type definitions.
# Closure pseudo-atoms G*/CG* and hydrated dummy W are deliberately unsupported.
_AD_ELEMENTS = {
    "H": 1,
    "HD": 1,
    "C": 6,
    "A": 6,
    "N": 7,
    "NA": 7,
    "O": 8,
    "OA": 8,
    "P": 15,
    "S": 16,
    "SA": 16,
    "F": 9,
    "Cl": 17,
    "Br": 35,
    "I": 53,
    "Si": 14,
    "At": 85,
    "Mg": 12,
    "Mn": 25,
    "Zn": 30,
    "Ca": 20,
    "Fe": 26,
}


def _read_pdbqt(path, structure):
    if len(structure.residues) != 1 or len({a.name for a in structure.atoms}) != len(
        structure.atoms
    ):
        raise ValueError(
            "PDBQT mapping requires one residue with globally unique template atom names"
        )
    lines = path.read_text().splitlines()
    models = [i for i, line in enumerate(lines) if line.startswith("MODEL ")]
    if models:
        start = models[0] + 1
        stop = next(
            (i for i in range(start, len(lines)) if lines[i].startswith(("ENDMDL", "MODEL "))),
            len(lines),
        )
        lines = lines[start:stop]
    records = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if len(records) != len(structure.atoms):
        raise ValueError(
            "PDBQT first pose must contain every template atom, including all hydrogens"
        )
    by_name = {a.name: a for a in structure.atoms}
    coords = np.empty((len(structure.atoms), 3))
    seen = set()
    residues = set()
    for line in records:
        name = line[12:16].strip()
        if name not in by_name or name in seen:
            raise ValueError(
                "PDBQT atom names are unknown or ambiguous; explicit template identity required"
            )
        seen.add(name)
        residues.add((line[17:20], line[21:27]))
        if line[16:17] != " ":
            raise ValueError("PDBQT alternate atom locations are unsupported")
        ad_type = line.split()[-1] if line.split() else ""
        if _AD_ELEMENTS.get(ad_type) != by_name[name].atomic_number:
            raise ValueError(f"PDBQT AutoDock type/element identity mismatch for atom {name}")
        try:
            coords[by_name[name].idx] = [float(line[k : k + 8]) for k in (30, 38, 46)]
        except ValueError as exc:
            raise ValueError("Invalid PDBQT pose coordinates") from exc
    if len(residues) != 1 or not np.isfinite(coords).all():
        raise ValueError("PDBQT requires one residue and finite coordinates for every atom")
    return coords, "pdbqt_first_model_named"


def _plain_graph(numbers, edges, selected=None):
    from rdkit import Chem

    selected = list(range(len(numbers))) if selected is None else list(selected)
    local = {index: i for i, index in enumerate(selected)}
    rw = Chem.RWMol()
    for index in selected:
        atom = Chem.Atom(int(numbers[index]))
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    for a, b in edges:
        if a in local and b in local:
            rw.AddBond(local[a], local[b], Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol


def _ff_equivalent(first, second):
    if first.atomic_number != second.atomic_number or first.type != second.type:
        return False
    fields = ("charge", "mass", "epsilon", "rmin", "solvent_radius", "screen")
    return np.allclose(
        [getattr(first, k) for k in fields],
        [getattr(second, k) for k in fields],
        atol=1e-8,
        rtol=1e-8,
    )


def _verify_tetrahedral_stereo(structure, coords, graph):
    from rdkit import Chem

    centers = [
        item.centeredOn
        for item in Chem.FindPotentialStereo(graph)
        if str(item.type) == "Atom_Tetrahedral"
        and structure.atoms[item.centeredOn].atomic_number in (6, 14)
    ]
    if not centers:
        return
    reference = structure.coordinates
    if (
        reference is None
        or np.asarray(reference).shape != coords.shape
        or not np.isfinite(reference).all()
    ):
        raise ValueError(
            "SDF tetrahedral stereochemistry requires explicit ordered template reference coordinates"
        )
    for center in centers:
        neighbors = sorted(atom.GetIdx() for atom in graph.GetAtomWithIdx(center).GetNeighbors())
        if len(neighbors) != 4:
            raise ValueError("Unsupported incomplete tetrahedral center; all hydrogens required")

        def sign(xyz, indices=neighbors):
            points = xyz[indices]
            vectors = points[:3] - points[3]
            norms = np.linalg.norm(vectors, axis=1)
            if np.any(norms < 1e-6):
                raise ValueError("Degenerate stereocenter geometry cannot be checked")
            volume = float(np.linalg.det(vectors / norms[:, None]))
            if abs(volume) < 1e-3:
                raise ValueError("Nearly planar stereocenter geometry cannot be checked")
            return volume > 0

        if sign(np.asarray(reference)) != sign(coords):
            raise ValueError(
                f"SDF tetrahedral stereo inversion at template atom {structure.atoms[center].name}"
            )


def _read_sdf(path, structure):
    from rdkit import Chem

    if path.suffix.lower() == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False, strictParsing=True)
        mol = next(iter(supplier), None)
    else:
        mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False, strictParsing=True)
    if mol is None or mol.GetNumConformers() != 1:
        raise ValueError("Cannot read first complete SDF/MOL pose with explicit coordinates")
    if mol.GetNumAtoms() != len(structure.atoms):
        raise ValueError(
            "SDF pose atom count differs from template; all hydrogens/protonation must be explicit"
        )
    pose_xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if not np.isfinite(pose_xyz).all():
        raise ValueError("SDF pose contains nonfinite coordinates")
    numbers = [a.atomic_number for a in structure.atoms]
    pose_numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
    if sorted(numbers) != sorted(pose_numbers):
        raise ValueError("SDF element/hydrogen composition differs from template")
    edges = {tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in structure.bonds}
    pose_edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    if len(edges) != len(pose_edges):
        raise ValueError("SDF connectivity graph differs from template")
    graph = _plain_graph(numbers, edges)
    mappings = []  # Each mapping lists pose indices in template order.
    if mol.HasProp("GO_BUILDER_ATOM_NAMES"):
        try:
            names = json.loads(mol.GetProp("GO_BUILDER_ATOM_NAMES"))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "GO_BUILDER_ATOM_NAMES must be a JSON array in SDF atom order"
            ) from exc
        template_names = [a.name for a in structure.atoms]
        if (
            not isinstance(names, list)
            or len(names) != len(numbers)
            or len(set(template_names)) != len(numbers)
            or any(not isinstance(name, str) for name in names)
            or set(names) != set(template_names)
        ):
            raise ValueError(
                "Explicit SDF atom name map must cover unique complete template atom names"
            )
        mappings = [tuple(names.index(name) for name in template_names)]
    else:
        heavy = [i for i, z in enumerate(numbers) if z != 1]
        pose_heavy = [i for i, z in enumerate(pose_numbers) if z != 1]
        if not heavy:
            raise ValueError("Hydrogen-only SDF requires explicit atom name map")
        template_heavy = _plain_graph(numbers, edges, heavy)
        source_heavy = _plain_graph(pose_numbers, pose_edges, pose_heavy)
        candidates = source_heavy.GetSubstructMatches(
            template_heavy, uniquify=False, maxMatches=4097
        )
        if len(candidates) >= 4097:
            raise ValueError(
                "SDF heavy-atom symmetry search limit reached; provide explicit atom name map"
            )

        def hydrogen_groups(nums, links):
            neighbors = {i: [] for i in range(len(nums))}
            for a, b in links:
                neighbors[a].append(b)
                neighbors[b].append(a)
            groups = {i: [] for i, z in enumerate(nums) if z != 1}
            for i, z in enumerate(nums):
                if z != 1:
                    continue
                if len(neighbors[i]) != 1 or nums[neighbors[i][0]] == 1:
                    raise ValueError(
                        "SDF/template hydrogen connectivity is incomplete or unsupported"
                    )
                groups[neighbors[i][0]].append(i)
            return groups

        h_groups = hydrogen_groups(numbers, edges)
        pose_h_groups = hydrogen_groups(pose_numbers, pose_edges)
        for candidate in candidates:
            mapping = [-1] * len(numbers)
            valid = True
            for local, index in enumerate(heavy):
                pose_index = pose_heavy[candidate[local]]
                mapping[index] = pose_index
                target_h, source_h = sorted(h_groups[index]), sorted(pose_h_groups[pose_index])
                if len(target_h) != len(source_h):
                    valid = False
                    break
                if target_h and any(
                    not _ff_equivalent(structure.atoms[target_h[0]], structure.atoms[h])
                    for h in target_h[1:]
                ):
                    raise ValueError(
                        "SDF hydrogen symmetry has different RESP/type/mass; explicit atom name map required"
                    )
                for a, b in zip(target_h, source_h, strict=True):
                    mapping[a] = b
            if valid:
                mappings.append(tuple(mapping))
        if not mappings:
            raise ValueError("SDF explicit graph/hydrogen protonation does not match template")
    validated = []
    for mapping in mappings:
        if len(set(mapping)) != len(numbers) or any(
            numbers[i] != pose_numbers[j] for i, j in enumerate(mapping)
        ):
            continue
        inverse = {pose: i for i, pose in enumerate(mapping)}
        if {tuple(sorted((inverse[a], inverse[b]))) for a, b in pose_edges} != edges:
            continue
        validated.append(mapping)
    if not validated:
        raise ValueError("SDF named/automatic connectivity graph mapping differs from template")
    # Any graph symmetry that exchanges unequal Amber parameters is ambiguous.
    first_inverse = {pose: i for i, pose in enumerate(validated[0])}
    for mapping in validated[1:]:
        for i, pose in enumerate(mapping):
            if not _ff_equivalent(structure.atoms[i], structure.atoms[first_inverse[pose]]):
                raise ValueError(
                    "SDF heavy-atom symmetry exchanges unequal charge/type/mass; explicit atom name map required"
                )
    mapping = min(validated)
    coords = pose_xyz[list(mapping)]
    _verify_tetrahedral_stereo(structure, coords, graph)
    return coords, "sdf_first_pose_explicit_graph"


def _read_pose(path, structure):
    path = Path(path)
    if path.suffix.lower() == ".pdbqt":
        return _read_pdbqt(path, structure)
    if path.suffix.lower() in (".sdf", ".mol"):
        return _read_sdf(path, structure)
    if path.suffix.lower() in (".pdb", ".ent"):
        lines = path.read_text().splitlines()
        if sum(line.startswith("MODEL ") for line in lines) > 1:
            raise ValueError("PDB pose must contain one complete model")
        records = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
        if len(records) != len(structure.atoms):
            raise ValueError("PDB pose atom count differs from complete parameterized topology")
        coords = []
        previous_residue = None
        residue_index = -1
        for index, (line, atom) in enumerate(zip(records, structure.atoms, strict=True), 1):
            if len(line) < 78 or line[16:17] not in ("", " "):
                raise ValueError(
                    "PDB pose requires explicit element fields and no alternate locations"
                )
            element = pmd.periodic_table.Element[atom.atomic_number]
            if line[12:16].strip() != atom.name or line[76:78].strip().upper() != element.upper():
                raise ValueError(
                    f"PDB pose atom identity/name/element/order mismatch at atom {index}"
                )
            residue_key = (line[21:22], line[22:26], line[26:27])
            if residue_key != previous_residue:
                residue_index += 1
                previous_residue = residue_key
            if residue_index != atom.residue.idx or line[17:20].strip() != atom.residue.name[:3]:
                raise ValueError(f"PDB pose residue identity/order mismatch at atom {index}")
            try:
                coords.append([float(line[start : start + 8]) for start in (30, 38, 46)])
            except ValueError as exc:
                raise ValueError(f"Invalid PDB pose coordinates at atom {index}") from exc
        return np.asarray(coords, dtype=float), "pdb_checked"
    if path.suffix.lower() not in (".rst7", ".rst", ".restrt", ".inpcrd", ".ncrst", ".crd"):
        raise ValueError(
            "Unknown pose format; use a complete ordered Amber restart or strictly matched PDB"
        )
    try:
        restart = pmd.amber.Rst7.open(str(path))
    except (ValueError, TypeError, IndexError, OSError, pmd.exceptions.AmberError) as exc:
        raise ValueError(f"Cannot read complete Amber pose: {path}") from exc
    if restart.natom != len(structure.atoms):
        raise ValueError("Amber pose atom count differs from complete parameterized topology")
    coords = np.asarray(restart.coordinates, dtype=float).reshape((-1, 3))
    return coords, "amber_restart_topology_order"


def _load_component(topology, coordinates, reference_coordinates=None):
    topology, coordinates = Path(topology), Path(coordinates)
    try:
        structure = pmd.load_file(str(topology))
    except pmd.exceptions.AmberError as exc:
        raise ValueError(f"Invalid/corrupt Amber topology {topology}: {exc}") from exc
    if type(structure) is not pmd.amber.AmberParm:
        raise ValueError(
            "Complete standard Amber prmtop inputs are required (no unparameterized or CHAMBER input)"
        )
    if reference_coordinates is not None:
        if Path(reference_coordinates).suffix.lower() in (".sdf", ".mol", ".pdbqt"):
            raise ValueError(
                "Template stereo reference must be an ordered Amber restart or checked PDB"
            )
        reference_xyz, _ = _read_pose(reference_coordinates, structure)
        structure.coordinates = reference_xyz
    xyz, coordinate_format = _read_pose(coordinates, structure)
    if xyz.shape != (len(structure.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Pose must contain one finite coordinate triplet per topology atom")
    structure.coordinates = xyz
    _check_parameterized(structure)
    if any(
        res.name.upper() in ("WAT", "HOH", "SOL", "TIP3", "TIP4", "TIP5")
        for res in structure.residues
    ):
        raise ValueError(
            "Supply dry receptor and ligand inputs; explicit solvent is not silently stripped"
        )
    if any(atom.atomic_number <= 0 for atom in structure.atoms):
        raise ValueError("Extra-point/virtual-site complexes are not supported by this assembler")
    if (
        structure.combining_rule != "lorentz"
        or structure.has_NBFIX()
        or any(index < 0 for index in structure.parm_data["NONBONDED_PARM_INDEX"])
    ):
        raise ValueError(
            "Unsupported pair-specific Lennard-Jones or mixing protocol; unsafe merge refused"
        )
    if any(
        getattr(structure, name)
        for name in ("cmaps", "urey_bradleys", "impropers", "rb_torsions", "adjusts")
    ):
        raise ValueError(
            "Nonstandard Amber bonded/nonbonded terms require a separately validated merge"
        )
    radius_label = " ".join(structure.parm_data.get("RADIUS_SET", [])).strip()
    # The explanatory text in "modified Bondi radii (mbondi)" also contains
    # "Bondi". Read the explicit protocol tag, not the first matching word.
    match = re.search(r"\((mbondi3|mbondi2|mbondi|bondi)\)\s*$", radius_label.lower())
    if match is None:
        match = re.fullmatch(r"(mbondi3|mbondi2|mbondi|bondi)", radius_label.lower())
    if not match:
        raise ValueError(
            f"Unknown GB radii protocol {radius_label!r}; assign explicit compatible radii first"
        )
    radii = match.group(1)
    # Structure copies contain atoms and connectivity but no Amber flag table.
    # ChRad supports these directly; an AmberParm copy lacks RADIUS_SET and
    # incorrectly appears to have a complete flag table.
    check = structure.copy(pmd.Structure)
    ChRad(check, radii)
    original_values = np.asarray([(a.solvent_radius, a.screen) for a in structure.atoms])
    expected_values = np.asarray([(a.solvent_radius, a.screen) for a in check.atoms])
    if (
        not np.isfinite(original_values).all()
        or np.any(original_values <= 0)
        or not np.allclose(original_values, expected_values, atol=1e-7, rtol=0)
    ):
        raise ValueError("Atom GB radii/screens do not match declared radii protocol")
    active_scaling = set()
    for term in structure.dihedrals:
        if term.improper or term.ignore_end:
            continue
        values = term.type if isinstance(term.type, list) else [term.type]
        for value in values:
            if not all(math.isfinite(x) and x > 0 for x in (value.scee, value.scnb)):
                raise ValueError("Invalid proper-torsion SCEE/SCNB protocol")
            active_scaling.add((round(value.scee, 7), round(value.scnb, 7)))
    if len(active_scaling) > 1:
        raise ValueError("Mixed proper-torsion SCEE/SCNB protocols require explicit validation")
    if any(
        np.linalg.norm(xyz[b.atom1.idx] - xyz[b.atom2.idx]) > max(3.0, 3 * b.type.req)
        for b in structure.bonds
    ):
        raise ValueError(
            "Wrapped/stretched intracomponent bonds: provide a complete made-whole docking pose"
        )
    structure.box = None
    return structure, {
        "source_topology": str(topology.resolve()),
        "source_coordinates": str(coordinates.resolve()),
        "topology_sha256": _hash(topology),
        "coordinates_sha256": _hash(coordinates),
        "coordinate_format": coordinate_format,
        "reference_coordinates": str(Path(reference_coordinates).resolve())
        if reference_coordinates
        else None,
        "reference_coordinates_sha256": _hash(reference_coordinates)
        if reference_coordinates
        else None,
        "radii": radii,
        "radii_label": radius_label,
        "scee_scnb": list(next(iter(active_scaling))) if active_scaling else None,
        "atom_count": len(structure.atoms),
        "residue_count": len(structure.residues),
        "net_charge": float(sum(a.charge for a in structure.atoms)),
    }


def _compatible(parts, information):
    if len({entry["radii"] for entry in information}) != 1:
        raise ValueError("Receptor and ligand GB radii protocols differ")
    protocols = {tuple(entry["scee_scnb"]) for entry in information if entry["scee_scnb"]}
    if len(protocols) > 1:
        raise ValueError("Receptor and ligand SCEE/SCNB 1-4 protocols differ")
    types = {}
    for structure in parts:
        for atom in structure.atoms:
            values = (
                atom.atomic_number,
                atom.mass,
                atom.epsilon,
                atom.rmin,
                atom.epsilon_14,
                atom.rmin_14,
            )
            if atom.type in types and not np.allclose(
                types[atom.type], values, atol=1e-7, rtol=1e-7
            ):
                raise ValueError(
                    f"Incompatible same-name Lennard-Jones type {atom.type!r}; no automatic renaming/guessing"
                )
            types[atom.type] = values
    contacts = []
    for i, first in enumerate(parts):
        for j, second in enumerate(parts[i + 1 :], i + 1):
            # Bounded-memory blocks; distances diagnose severe overlap only.
            minimum = float("inf")
            closest = None
            for start in range(0, len(first.atoms), 256):
                delta = (
                    first.coordinates[start : start + 256, None, :] - second.coordinates[None, :, :]
                )
                squared = np.sum(delta * delta, axis=2)
                local = np.unravel_index(np.argmin(squared), squared.shape)
                value = float(np.sqrt(squared[local]))
                if value < minimum:
                    minimum = value
                    closest = [start + int(local[0]) + 1, int(local[1]) + 1]
            if minimum < 0.6:
                raise ValueError(
                    f"Severe cross-component overlap/contact ({minimum:.3f} Å): supply complete nonoverlapping simultaneous poses"
                )
            contacts.append(
                {
                    "components": [i, j],
                    "minimum_distance_angstrom": minimum,
                    "local_atom_indices": closest,
                }
            )
    return contacts


def _exclusions(structure):
    return {
        tuple(sorted((atom.idx, partner)))
        for atom in structure.atoms
        for partner in atom.nonbonded_exclusions(only_greater=True)
    }


def _active_pairs(structure):
    result = {}
    for term in structure.dihedrals:
        if term.improper or term.ignore_end:
            continue
        values = term.type if isinstance(term.type, list) else [term.type]
        key = tuple(sorted((term.atom1.idx, term.atom4.idx)))
        result[key] = tuple((round(t.scee, 7), round(t.scnb, 7)) for t in values)
    return result


def _combine(parts):
    merged = parts[0].copy(pmd.Structure)
    for part in parts[1:]:
        merged += part.copy(pmd.Structure)
    offset = 0
    for part in parts:
        # Structure concatenation does not promise preservation of unusual
        # explicit exclusions. Carry the source exclusion list explicitly.
        for a, b in _exclusions(part):
            merged.atoms[offset + a].exclude(merged.atoms[offset + b])
        offset += len(part.atoms)
    merged.box = None
    return merged


def _pdb_model(structure):
    return {
        "schema_version": 1,
        "name": "parameterized_docking_complex",
        "atoms": [
            {
                "id": atom.idx + 1,
                "name": atom.name,
                "element": pmd.periodic_table.Element[atom.atomic_number],
                "resid": atom.residue.idx + 1,
                "resname": atom.residue.name,
                "position": list(map(float, structure.coordinates[atom.idx])),
                "charge": float(atom.charge),
                "atom_type": atom.type,
                "formal_charge": 0,
                "metadata": {},
            }
            for atom in structure.atoms
        ],
        "bonds": [
            {
                "a": min(b.atom1.idx, b.atom2.idx) + 1,
                "b": max(b.atom1.idx, b.atom2.idx) + 1,
                "order": 1.0,
                "image": [0, 0, 0],
            }
            for b in structure.bonds
        ],
        "cell": None,
        "pbc": [False] * 3,
        "metadata": {
            "representation": "amber_import",
            "bond_orders_unknown": True,
            "formal_charges_unknown": True,
            "parameterized": True,
        },
    }


def _write_component(structure, prefix, radii_label):
    parm = pmd.amber.AmberParm.from_structure(structure.copy(pmd.Structure))
    parm.box = None
    parm.parm_data["RADIUS_SET"] = [radii_label]
    top, crd, pdb = (Path(str(prefix) + suffix) for suffix in (".parm7", ".rst7", ".pdb"))
    parm.save(str(top), overwrite=False)
    parm.save(str(crd), overwrite=False)
    reread = pmd.load_file(str(top), xyz=str(crd))
    verification = _verify_parameter_roundtrip(structure, reread)
    if _exclusions(structure) != _exclusions(reread) or _active_pairs(structure) != _active_pairs(
        reread
    ):
        raise ValueError("Amber assembly altered nonbonded exclusions or active 1-4 scaling")
    before = [(a.solvent_radius, a.screen) for a in structure.atoms]
    after = [(a.solvent_radius, a.screen) for a in reread.atoms]
    if not np.allclose(before, after, atol=1e-7, rtol=0):
        raise ValueError("Amber assembly altered GB radii/screens")
    write_pdb(_pdb_model(reread), pdb)
    # The shared writer includes a generic construction remark; replace it for
    # these already parameterized inputs without changing ATOM/CONECT records.
    text = pdb.read_text().splitlines(keepends=True)
    text[0] = "REMARK Parameterized docking components; topology is authoritative.\n"
    pdb.write_text("".join(text))
    return {"topology": str(top), "coordinates": str(crd), "pdb": str(pdb)}, verification


def assemble_complex(
    receptor_topology,
    receptor_coordinates,
    ligand_topologies,
    ligand_coordinates,
    output_dir,
    *,
    ligand_reference_coordinates=None,
):
    """Merge one complete receptor and one or more simultaneous ligand poses.

    A receptor's existing covalent grafts remain receptor atoms. Ligands are
    noncovalent independent components. This writes dry Amber inputs and masks;
    it does not establish the scientific validity of a charge model or pose.
    ligand_reference_coordinates optionally supplies one ordered restart/PDB
    per ligand for SDF tetrahedral stereo verification; None entries are allowed.
    """
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Choose a fresh complex output directory: {output_dir}")
    if not isinstance(ligand_topologies, (list, tuple)) or not isinstance(
        ligand_coordinates, (list, tuple)
    ):
        raise ValueError("Ligand topologies and coordinates must be explicit lists")
    if not ligand_topologies or len(ligand_topologies) != len(ligand_coordinates):
        raise ValueError(
            "Supply one coordinate file for each ligand; list lengths/counts must agree"
        )
    pairs = [
        (receptor_topology, receptor_coordinates),
        *zip(ligand_topologies, ligand_coordinates, strict=True),
    ]
    if ligand_reference_coordinates is None:
        references = [None] * len(ligand_topologies)
    elif not isinstance(ligand_reference_coordinates, (list, tuple)) or len(
        ligand_reference_coordinates
    ) != len(ligand_topologies):
        raise ValueError("Supply one optional template reference coordinate entry per ligand")
    else:
        references = list(ligand_reference_coordinates)
    loaded = [
        _load_component(top, xyz, reference)
        for (top, xyz), reference in zip(pairs, [None, *references], strict=True)
    ]
    parts = [item[0] for item in loaded]
    information = [item[1] for item in loaded]
    contacts = _compatible(parts, information)
    complex_structure = _combine(parts)
    ligands_structure = _combine(parts[1:])
    count = len(parts[0].atoms)
    expected = {tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in complex_structure.bonds}
    if any(a < count <= b for a, b in expected):
        raise ValueError("Unexpected receptor-ligand covalent bond; assembly aborted")
    total_residues = len(complex_structure.residues)
    receptor_residues = len(parts[0].residues)
    ligand_masks = []
    next_residue = receptor_residues + 1
    for part in parts[1:]:
        ligand_masks.append(_mask(next_residue, next_residue + len(part.residues) - 1))
        next_residue += len(part.residues)
    masks = {
        "receptor": _mask(1, receptor_residues),
        "ligand": _mask(receptor_residues + 1, total_residues),
        "complex": _mask(1, total_residues),
        "ligands": ligand_masks,
    }
    # Validate PDB representability before creating any output.
    if len(complex_structure.atoms) > 99999 or total_residues > 9999:
        raise ValueError("Explicit PDB identifier limits exceeded; use a smaller assembly")
    if any(np.max(np.abs(part.coordinates)) >= 9999 for part in parts):
        raise ValueError("Docking coordinates exceed PDB representation limits")
    output_dir.mkdir(parents=True)
    outputs = {}
    verification = {}
    for label, structure, basename in [
        ("complex", complex_structure, "compl"),
        ("receptor", parts[0], "rec"),
        ("ligand", ligands_structure, "lig"),
    ]:
        outputs[label], verification[label] = _write_component(
            structure, output_dir / basename, information[0]["radii_label"]
        )
    individual = []
    for i, part in enumerate(parts[1:], 1):
        files, checked = _write_component(
            part, output_dir / f"lig_{i:03d}", information[0]["radii_label"]
        )
        individual.append(
            {
                **information[i],
                **files,
                "complex_mask": ligand_masks[i - 1],
                "atom_range_complex": [
                    sum(len(p.atoms) for p in parts[:i]) + 1,
                    sum(len(p.atoms) for p in parts[: i + 1]),
                ],
                "roundtrip_verified": checked,
            }
        )
    net_charge = float(sum(a.charge for a in complex_structure.atoms))
    warnings = [
        "Poses and GO/graft charge models require scientific validation; assembly is technical preparation only.",
        "Multiple ligands represent simultaneous occupancy; combined MM/GBSA is not a separate affinity for each ligand.",
    ]
    if abs(net_charge - round(net_charge)) > 1e-4:
        warnings.append(
            "Input partial charges have nonintegral total; preserved without renormalization."
        )
    result = {
        "schema_version": 1,
        "status": "assembled",
        "parameterized": True,
        "complex_topology": outputs["complex"]["topology"],
        "complex_coordinates": outputs["complex"]["coordinates"],
        "complex_pdb": outputs["complex"]["pdb"],
        "receptor_topology": outputs["receptor"]["topology"],
        "receptor_coordinates": outputs["receptor"]["coordinates"],
        "receptor_pdb": outputs["receptor"]["pdb"],
        "ligand_topology": outputs["ligand"]["topology"],
        "ligand_coordinates": outputs["ligand"]["coordinates"],
        "ligand_pdb": outputs["ligand"]["pdb"],
        "receptor": information[0],
        "ligands": individual,
        "masks": masks,
        "atom_count": len(complex_structure.atoms),
        "bond_count": len(complex_structure.bonds),
        "net_charge": net_charge,
        "radii_set": information[0]["radii"],
        "contacts": contacts,
        "cross_component_bonds": 0,
        "roundtrip_verified": verification,
        "gbsa": {
            "ready": True,
            "binding_ready": True,
            "blockers": [],
            "warnings": warnings,
            "masks": masks,
            "radii_set": information[0]["radii"],
            "interpretation": "receptor versus joint ligand group"
            if len(individual) > 1
            else "receptor versus one ligand",
        },
        "policy": "complete topology atom order; same coordinate frame; no alignment, bond inference, solvation or MD",
    }
    (output_dir / "masks.json").write_text(json.dumps(masks, indent=2) + "\n")
    (output_dir / "manifest.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
