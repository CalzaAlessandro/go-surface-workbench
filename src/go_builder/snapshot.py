"""Prepare a complete docking receptor from an actual Amber conformation.

Periodic solvent wrapping is removed through existing covalent bonds. Atoms
are never projected onto a plane, rotated, retyped or connected by proximity.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import deque
from pathlib import Path

import numpy as np
import parmed as pmd

from .complexes import _read_pose, _write_component
from .model import save_model
from .topology import _check_parameterized, import_amber


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _cell(box):
    if box is None:
        return None
    box = np.asarray(box, dtype=float)
    if box.shape != (6,) or not np.isfinite(box).all() or np.any(box[:3] <= 0):
        raise ValueError(
            "Snapshot cell must contain three positive lengths and three finite angles"
        )
    if np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
        raise ValueError("Invalid snapshot cell angles")
    vectors = pmd.geometry.box_lengths_and_angles_to_vectors(*box)
    if hasattr(vectors, "value_in_unit"):
        vectors = vectors.value_in_unit(pmd.unit.angstroms)
    cell = np.asarray(vectors, dtype=float)
    if cell.shape != (3, 3) or abs(np.linalg.det(cell)) < 1e-8:
        raise ValueError("Singular snapshot cell")
    return cell


def make_whole(coordinates, bonds, cell):
    """Unwrap finite components by their explicit zero-based covalent edges.

    Conflicting cycle translations indicate a covalently winding graph or
    inconsistent input. Such a graph cannot be exported as a finite receptor.
    """
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("Finite N×3 coordinates required")
    if cell is None:
        return xyz.copy()
    cell = np.asarray(cell, dtype=float)
    if cell.shape != (3, 3) or not np.isfinite(cell).all() or abs(np.linalg.det(cell)) < 1e-8:
        raise ValueError("A finite nonsingular cell is required for periodic unwrapping")
    inverse = np.linalg.inv(cell)
    offsets = np.asarray(list(itertools.product((-1, 0, 1), repeat=3)))
    neighbors = [[] for _ in xyz]
    for a, b in bonds:
        if a == b or not 0 <= a < len(xyz) or not 0 <= b < len(xyz):
            raise ValueError("Invalid explicit bond endpoint")
        neighbors[a].append(b)
        neighbors[b].append(a)
    result = np.empty_like(xyz)
    visited = np.zeros(len(xyz), dtype=bool)
    for root in range(len(xyz)):
        if visited[root]:
            continue
        visited[root] = True
        result[root] = xyz[root]
        queue = deque([root])
        while queue:
            a = queue.popleft()
            for b in neighbors[a]:
                delta = xyz[b] - xyz[a]
                translations = np.rint(delta @ inverse) + offsets
                candidates = delta - translations @ cell
                vector = candidates[np.argmin(np.sum(candidates * candidates, axis=1))]
                target = result[a] + vector
                if visited[b]:
                    if np.linalg.norm(target - result[b]) > 1e-3:
                        raise ValueError(
                            "Covalently winding/infinite periodic graph or inconsistent cell: cannot make a finite docking receptor"
                        )
                else:
                    visited[b] = True
                    result[b] = target
                    queue.append(b)
    return result


def cover_boxes(coordinates, *, box_size=(20, 20, 20), overlap=5.0, padding=5.0, max_boxes=4096):
    """Axis-aligned overlapping 3D boxes cover the whole actual receptor extent.

    The grid is geometric coverage, not a solvent-accessibility calculation.
    Fold interiors and pockets are retained as possible search regions.
    """
    xyz = np.asarray(coordinates, dtype=float)
    size = np.asarray(box_size, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("Finite nonempty N×3 receptor coordinates required")
    if size.shape != (3,) or not np.isfinite(size).all() or np.any(size <= 0):
        raise ValueError("box_size requires three finite positive lengths")
    if not math.isfinite(overlap) or overlap < 0 or overlap >= min(size):
        raise ValueError("overlap must be nonnegative and smaller than every box dimension")
    if not math.isfinite(padding) or padding < 0:
        raise ValueError("padding must be finite and nonnegative")
    if type(max_boxes) is not int or max_boxes < 1:
        raise ValueError("max_boxes must be a positive integer")
    lower, upper = xyz.min(axis=0) - padding, xyz.max(axis=0) + padding
    counts = [
        max(1, math.ceil(max(0.0, float(width - length)) / (length - overlap)) + 1)
        for width, length in zip(upper - lower, size, strict=True)
    ]
    if math.prod(counts) > max_boxes:
        raise ValueError(
            f"Requested 3D coverage requires {math.prod(counts)} boxes; increase box size or max_boxes explicitly"
        )
    axes = [
        np.asarray([(lo + hi) / 2]) if n == 1 else np.linspace(lo + length / 2, hi - length / 2, n)
        for lo, hi, length, n in zip(lower, upper, size, counts, strict=True)
    ]
    return [
        {
            "id": f"box_{i:03d}",
            "center": list(map(float, center)),
            "size": size.tolist(),
            "config": f"boxes/box_{i:03d}.txt",
        }
        for i, center in enumerate(itertools.product(*axes), 1)
    ]


def _frame(topology, coordinate_path, frame):
    if type(frame) is not int or frame < 1:
        raise ValueError("frame is a one-based positive integer")
    path = Path(coordinate_path)
    if path.suffix.lower() in (".nc", ".netcdf"):
        from netCDF4 import Dataset

        with Dataset(path, "r") as nc:
            if "coordinates" not in nc.variables:
                raise ValueError("Amber NetCDF trajectory lacks coordinates")
            data = nc.variables["coordinates"]
            if len(data.shape) != 3 or tuple(data.shape[1:]) != (len(topology.atoms), 3):
                raise ValueError("NetCDF atom count differs from topology")
            if frame > data.shape[0]:
                raise ValueError("Requested frame exceeds trajectory length")
            units = str(getattr(data, "units", "angstrom")).lower()
            if units not in ("angstrom", "angstroms", "a"):
                raise ValueError(f"Unsupported NetCDF coordinate units: {units}")
            xyz = np.asarray(data[frame - 1], dtype=float)
            if ("cell_lengths" in nc.variables) != ("cell_angles" in nc.variables):
                raise ValueError("NetCDF cell requires lengths and angles together")
            box = None
            if "cell_lengths" in nc.variables:
                box = np.r_[
                    nc.variables["cell_lengths"][frame - 1], nc.variables["cell_angles"][frame - 1]
                ].astype(float)
            time = float(nc.variables["time"][frame - 1]) if "time" in nc.variables else None
            if time is not None and not math.isfinite(time):
                raise ValueError("Nonfinite trajectory frame time")
        return xyz, box, time, "amber_netcdf"
    if frame != 1:
        raise ValueError("A restart/PDB contains one frame; use frame=1")
    xyz, kind = _read_pose(path, topology)
    if kind == "pdb_checked":
        cells = [line for line in path.read_text().splitlines() if line.startswith("CRYST1")]
        box = (
            None
            if not cells
            else [
                float(cells[0][a:b])
                for a, b in [(6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54)]
            ]
        )
        return xyz, box, None, kind
    restart = pmd.amber.Rst7.open(str(path))
    time = float(restart.time) if restart.time is not None else None
    return xyz, restart.box, time, kind


def prepare_docking_receptor(
    topology,
    coordinates,
    output_dir,
    *,
    receptor_mask=":*",
    frame=1,
    box_size=(20, 20, 20),
    overlap=5.0,
    padding=5.0,
    max_boxes=4096,
):
    """Extract, make whole and export the chosen native or equilibrated frame."""
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Choose a fresh docking receptor directory: {output_dir}")
    source = pmd.load_file(str(topology))
    if type(source) is not pmd.amber.AmberParm:
        raise ValueError("A complete standard Amber topology is required")
    xyz, box, time, kind = _frame(source, coordinates, frame)
    if xyz.shape != (len(source.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Snapshot must contain finite coordinates for every topology atom")
    source.coordinates = xyz
    _check_parameterized(source)
    try:
        selection = np.asarray(pmd.amber.AmberMask(source, receptor_mask).Selection(), dtype=bool)
    except (ValueError, TypeError, pmd.exceptions.MaskError) as exc:
        raise ValueError(f"Invalid receptor mask: {receptor_mask}") from exc
    if not selection.any():
        raise ValueError("Receptor mask selects no atoms")
    if any(selection[b.atom1.idx] != selection[b.atom2.idx] for b in source.bonds):
        raise ValueError(
            "Receptor mask cuts a covalent bond; include the complete GO and its grafts"
        )
    selected = source[selection.tolist()]
    if any(
        res.name.upper() in ("WAT", "HOH", "SOL", "TIP3", "TIP4", "TIP5")
        for res in selected.residues
    ):
        raise ValueError("Select a dry receptor mask; solvent is not a docking surface")
    raw = np.asarray(selected.coordinates).copy()
    selected.coordinates = make_whole(
        raw, [(b.atom1.idx, b.atom2.idx) for b in selected.bonds], _cell(box)
    )
    for bond in selected.bonds:
        length = np.linalg.norm(
            selected.coordinates[bond.atom1.idx] - selected.coordinates[bond.atom2.idx]
        )
        if not 0.3 < length < max(3.0, 3 * bond.type.req):
            raise ValueError(
                "Invalid/stretched receptor bond after unwrapping; correct snapshot/cell before docking"
            )
    selected.box = None
    heavy = np.asarray([a.atomic_number > 1 for a in selected.atoms])
    if not heavy.any():
        raise ValueError("Receptor has no heavy atoms")
    boxes = cover_boxes(
        selected.coordinates[heavy],
        box_size=box_size,
        overlap=overlap,
        padding=padding,
        max_boxes=max_boxes,
    )
    radius_label = " ".join(source.parm_data.get("RADIUS_SET", []))
    output_dir.mkdir(parents=True)
    files, verification = _write_component(selected, output_dir / "receptor", radius_label)
    imported = import_amber(files["topology"], files["coordinates"])
    imported["metadata"]["snapshot"] = {
        "source_frame": frame,
        "source_time_ps": time,
        "flattened": False,
    }
    save_model(imported, output_dir / "receptor.json")
    (output_dir / "boxes").mkdir()
    for entry in boxes:
        values = [
            *(f"center_{axis} = {v:.10f}" for axis, v in zip("xyz", entry["center"], strict=True)),
            *(f"size_{axis} = {v:.10f}" for axis, v in zip("xyz", entry["size"], strict=True)),
        ]
        (output_dir / entry["config"]).write_text("\n".join(values) + "\n")
    result = {
        "schema_version": 1,
        "status": "receptor_prepared",
        "topology": files["topology"],
        "coordinates": files["coordinates"],
        "receptor_pdb": "receptor.pdb",
        "receptor_pdb_sha256": _hash(files["pdb"]),
        "model": str(output_dir / "receptor.json"),
        "boxes": boxes,
        "boxes_file": str(output_dir / "boxes.json"),
        "source_topology": str(Path(topology).resolve()),
        "source_coordinates": str(Path(coordinates).resolve()),
        "source_topology_sha256": _hash(topology),
        "source_coordinates_sha256": _hash(coordinates),
        "source_frame": frame,
        "source_time_ps": time,
        "source_coordinate_format": kind,
        "receptor_mask": receptor_mask,
        "source_atom_ids": (np.flatnonzero(selection) + 1).tolist(),
        "source_cell": None if box is None else np.asarray(box).tolist(),
        "atom_count": len(selected.atoms),
        "heavy_atom_count": int(heavy.sum()),
        "net_charge": float(sum(a.charge for a in selected.atoms)),
        "flattened": False,
        "coordinate_policy": "make whole through explicit bonds; retain root image; no alignment/projection",
        "max_unwrapping_displacement_angstrom": float(
            np.linalg.norm(selected.coordinates - raw, axis=1).max()
        ),
        "heavy_atom_extent_angstrom": np.ptp(selected.coordinates[heavy], axis=0).tolist(),
        "roundtrip_verified": verification,
        "coverage": "full 3D padded bounding volume; not a solvent-accessibility or reaction model",
        "docking_executed": False,
    }
    (output_dir / "boxes.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
