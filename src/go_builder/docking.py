"""Prepare and run AutoDock Vina boxes tied to a specific receptor snapshot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem

_CONFIG_KEYS = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
_BOX_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_VINA_ELEMENT = {
    "A": "C",
    "C": "C",
    "N": "N",
    "NA": "N",
    "OA": "O",
    "O": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "H": "H",
    "HD": "H",
    "HS": "H",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _read_pdb_atoms(path: Path, *, pdbqt: bool) -> list[tuple[str, tuple[float, float, float]]]:
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise ValueError(f"Malformed atom record in {path} at line {line_number}")
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError as exc:
            raise ValueError(f"Invalid coordinates in {path} at line {line_number}") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Non-finite coordinates in {path} at line {line_number}")
        if pdbqt:
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Missing AutoDock atom type in {path} at line {line_number}")
            atom_type = fields[-1].upper()
            element = _VINA_ELEMENT.get(atom_type)
            if element is None:
                raise ValueError(f"Unsupported AutoDock atom type {atom_type!r} in {path}")
        else:
            element = line[76:78].strip()
            if not element:
                raise ValueError(
                    f"Reference PDB requires explicit element columns at line {line_number}"
                )
            element = element[0].upper() + element[1:].lower()
        atoms.append((element, xyz))
    if not atoms:
        raise ValueError(f"No atom records found in {path}")
    return atoms


def _validate_receptor_identity(reference: Path, receptor: Path, *, pdbqt: bool = True) -> None:
    expected = [
        (element, xyz)
        for element, xyz in _read_pdb_atoms(reference, pdbqt=False)
        if element not in {"H", "D"}
    ]
    observed = [
        (element, xyz)
        for element, xyz in _read_pdb_atoms(receptor, pdbqt=pdbqt)
        if element not in {"H", "D"}
    ]
    if len(expected) != len(observed):
        raise ValueError("receptor PDBQT heavy-atom count does not match reference snapshot")
    # Vectorized per-reference distance checks avoid Python-level all-pairs work.
    observed_coordinates = np.asarray([xyz for _, xyz in observed], dtype=float)
    observed_elements = np.asarray([element for element, _ in observed], dtype=str)
    candidates = []
    tolerance_squared = 0.0200001**2
    for element, xyz in expected:
        same_element = np.flatnonzero(observed_elements == element)
        deltas = observed_coordinates[same_element] - np.asarray(xyz)
        squared_distances = np.einsum("ij,ij->i", deltas, deltas)
        candidates.append(same_element[squared_distances <= tolerance_squared].tolist())
    if any(not row for row in candidates):
        raise ValueError(
            "receptor PDBQT heavy atoms do not match reference snapshot coordinates/elements"
        )
    left_match = [-1] * len(expected)
    right_match = [-1] * len(observed)
    # Iterative augmenting paths preserve one-to-one matching without recursion.
    for root in range(len(expected)):
        if left_match[root] != -1:
            continue
        queue = [root]
        queue_index = 0
        visited_left = {root}
        parent_right: dict[int, int] = {}
        free_right = -1
        while queue_index < len(queue) and free_right == -1:
            left = queue[queue_index]
            queue_index += 1
            for right in candidates[left]:
                if right in parent_right:
                    continue
                parent_right[right] = left
                matched_left = right_match[right]
                if matched_left == -1:
                    free_right = right
                    break
                if matched_left not in visited_left:
                    visited_left.add(matched_left)
                    queue.append(matched_left)
        if free_right == -1:
            raise ValueError(
                "receptor PDBQT heavy-atom set is truncated or differs from reference snapshot"
            )
        right = free_right
        while right != -1:
            left = parent_right[right]
            previous_right = left_match[left]
            left_match[left] = right
            right_match[right] = left
            right = previous_right


def _config_values(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        match = re.fullmatch(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^\s#]+)\s*", line)
        if not match:
            raise ValueError(f"Invalid box config line {line_number} in {path}")
        key, raw = match.groups()
        if key not in _CONFIG_KEYS or key in values:
            raise ValueError(f"Box config has unexpected or duplicate key {key!r}")
        try:
            value = float(raw)
        except ValueError as exc:
            raise ValueError(f"Box config value for {key} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"Box config value for {key} must be finite")
        values[key] = value
    if set(values) != set(_CONFIG_KEYS):
        raise ValueError("Box config must contain exactly the six center/size values")
    return values


def _load_boxes(boxes_path: Path) -> tuple[dict, Path, list[dict]]:
    try:
        document = json.loads(boxes_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read boxes JSON: {exc}") from exc
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ValueError("Unsupported boxes schema; expected schema_version 1")
    receptor_rel = document.get("receptor_pdb")
    if not isinstance(receptor_rel, str) or not receptor_rel:
        raise ValueError("boxes JSON must name receptor_pdb relative to itself")
    if Path(receptor_rel).is_absolute():
        raise ValueError("receptor_pdb must be relative to boxes JSON")
    reference = (boxes_path.parent / receptor_rel).resolve(strict=True)
    try:
        reference.relative_to(boxes_path.parent.resolve())
    except ValueError as exc:
        raise ValueError("receptor_pdb must remain under boxes JSON directory") from exc
    expected_hash = document.get("receptor_pdb_sha256")
    if not isinstance(expected_hash, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_hash):
        raise ValueError("boxes JSON must contain a SHA-256 for receptor_pdb")
    if _sha256(reference).lower() != expected_hash.lower():
        raise ValueError("reference receptor PDB SHA-256 mismatch")
    boxes = document.get("boxes")
    if not isinstance(boxes, list) or not boxes:
        raise ValueError("boxes JSON must contain at least one box")
    ids: set[str] = set()
    parsed: list[dict] = []
    for box in boxes:
        if not isinstance(box, dict):
            raise ValueError("each box must be an object")
        box_id = box.get("id")
        if not isinstance(box_id, str) or not _BOX_ID.fullmatch(box_id) or box_id in ids:
            raise ValueError("box ids must be unique simple filename-safe identifiers")
        ids.add(box_id)
        center, size = box.get("center"), box.get("size")
        if (
            not isinstance(center, list)
            or len(center) != 3
            or not isinstance(size, list)
            or len(size) != 3
        ):
            raise ValueError(f"box {box_id} must have three center and size values")
        values = [*center, *size]
        if any(
            isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
            for v in values
        ):
            raise ValueError(f"box {box_id} center and size must be finite numbers")
        if any(v <= 0 for v in size):
            raise ValueError(f"box {box_id} size values must be positive")
        config_rel = box.get("config")
        if not isinstance(config_rel, str) or not config_rel:
            raise ValueError(f"box {box_id} must name a config file")
        if Path(config_rel) != Path("boxes") / f"{box_id}.txt":
            raise ValueError(f"box {box_id} config must be boxes/{box_id}.txt")
        config_path = (boxes_path.parent / config_rel).resolve(strict=True)
        try:
            config_path.relative_to(boxes_path.parent.resolve())
        except ValueError as exc:
            raise ValueError(f"box {box_id} config must remain under boxes JSON directory") from exc
        supplied = _config_values(config_path)
        target = dict(zip(_CONFIG_KEYS, values, strict=True))
        if any(
            not math.isclose(supplied[key], float(target[key]), rel_tol=0.0, abs_tol=1e-8)
            for key in _CONFIG_KEYS
        ):
            raise ValueError(f"box {box_id} config values do not match boxes JSON")
        parsed.append(
            {"id": box_id, "center": [float(v) for v in center], "size": [float(v) for v in size]}
        )
    return document, reference, parsed


def _config_text(box: dict) -> str:
    values = [*box["center"], *box["size"]]
    return "".join(
        f"{key} = {value:.8g}\n" for key, value in zip(_CONFIG_KEYS, values, strict=True)
    )


def _resolve_binary(executable: str, label: str) -> str:
    if not isinstance(executable, str) or not executable or "\0" in executable:
        raise ValueError(f"{label} must be a nonempty executable name or path")
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(f"{label} executable not found: {executable}")
    return resolved


def _read_sdf(path: Path, *, exactly_one: bool) -> list:
    try:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True, strictParsing=True)
        molecules = list(supplier)
    except Exception as exc:
        raise ValueError(f"Cannot parse ligand SDF {path}: {exc}") from exc
    if exactly_one and (len(molecules) != 1 or molecules[0] is None):
        raise ValueError("GNINA SDF ligand must contain exactly one valid molecule")
    if not molecules or molecules[0] is None:
        raise ValueError(f"SDF contains no valid molecule: {path}")
    molecule = molecules[0]
    if molecule.GetNumAtoms() == 0 or molecule.GetNumConformers() == 0:
        raise ValueError(f"SDF molecule must have atoms and coordinates: {path}")
    conformer = molecule.GetConformer()
    if not conformer.Is3D():
        raise ValueError(f"SDF molecule must have 3D coordinates: {path}")
    coordinates = np.asarray(
        [
            (position.x, position.y, position.z)
            for index in range(molecule.GetNumAtoms())
            for position in [conformer.GetAtomPosition(index)]
        ],
        dtype=float,
    )
    if coordinates.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(coordinates).all():
        raise ValueError(f"SDF coordinates must be finite 3D coordinates: {path}")
    return molecules


def _input_heavy_atom_count(path: Path, file_format: str) -> int:
    if file_format == "sdf":
        molecule = _read_sdf(path, exactly_one=True)[0]
        return sum(atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms())
    atoms = _read_pdb_atoms(path, pdbqt=True)
    return sum(element not in {"H", "D"} for element, _ in atoms)


def _gnina_score(path: Path) -> float:
    molecule = _read_sdf(path, exactly_one=False)[0]
    if not molecule.HasProp("minimizedAffinity"):
        raise ValueError(f"No minimizedAffinity property in first GNINA pose of {path}")
    try:
        score = float(molecule.GetProp("minimizedAffinity"))
    except ValueError as exc:
        raise ValueError(f"Invalid minimizedAffinity property in {path}") from exc
    if not math.isfinite(score):
        raise ValueError(f"Non-finite minimizedAffinity property in {path}")
    return score


def _validate_pose_heavy_atom_count(path: Path, engine: str, expected_count: int) -> None:
    if engine == "gnina":
        atoms = [
            atom
            for atom in _read_sdf(path, exactly_one=False)[0].GetAtoms()
            if atom.GetAtomicNum() > 1
        ]
    else:
        atoms = [
            element for element, _ in _read_pdb_atoms(path, pdbqt=True) if element not in {"H", "D"}
        ]
    if len(atoms) != expected_count:
        raise ValueError(f"{engine} output heavy-atom count does not match ligand input")


def _command(
    vina: str,
    receptor: Path,
    ligand: Path,
    config: Path,
    output: Path,
    cpus: int,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
) -> list[str]:
    return [
        vina,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--config",
        str(config),
        "--cpu",
        str(cpus),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--seed",
        str(seed),
        "--out",
        str(output),
    ]


def _gnina_command(
    gnina: str,
    receptor: Path,
    ligand: Path,
    config: Path,
    output: Path,
    cpus: int,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
) -> list[str]:
    return [
        gnina,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--config",
        str(config),
        "--cpu",
        str(cpus),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--seed",
        str(seed),
        "--scoring",
        "vina",
        "--cnn_scoring",
        "none",
        "--no_gpu",
        "--pose_sort_order",
        "Energy",
        "--out",
        str(output),
    ]


def _best_score(path: Path) -> float:
    for line in path.read_text().splitlines():
        match = re.match(r"\s*REMARK\s+VINA RESULT:\s*([-+0-9.eE]+)(?:\s|$)", line)
        if match:
            score = float(match.group(1))
            if math.isfinite(score):
                return score
    raise ValueError(f"No finite REMARK VINA RESULT score in {path}")


def dock_boxes(
    boxes_path: str | Path,
    receptor_pdbqt: str | Path,
    ligand_pdbqt: str | Path,
    output_dir: str | Path,
    *,
    vina: str = "vina",
    cpus: int = 1,
    exhaustiveness: int = 16,
    num_modes: int = 10,
    runs: int = 3,
    seed: int = 1,
    run: bool = False,
    engine: str = "vina",
    gnina: str = "gnina",
) -> dict:
    """Prepare reviewable Vina/GNINA jobs or run boxes and seeds sequentially.

    Receptor identity is checked against the hashed snapshot PDB by one-to-one
    element and coordinate matching. GNINA accepts a receptor PDB or PDBQT and
    a one-record graph-bearing ligand SDF or PDBQT. No atom typing or charges
    are converted by this function.
    """
    cpus = _positive_int(cpus, "cpus")
    exhaustiveness = _positive_int(exhaustiveness, "exhaustiveness")
    num_modes = _positive_int(num_modes, "num_modes")
    runs = _positive_int(runs, "runs")
    seed = _positive_int(seed, "seed")
    if not isinstance(run, bool):
        raise ValueError("run must be a boolean")
    if engine not in {"vina", "gnina"}:
        raise ValueError("engine must be 'vina' or 'gnina'")
    boxes_path = Path(boxes_path).resolve(strict=True)
    receptor = Path(receptor_pdbqt).resolve(strict=True)
    ligand = Path(ligand_pdbqt).resolve(strict=True)
    if not receptor.is_file() or not ligand.is_file():
        raise ValueError("receptor and ligand inputs must be files")
    document, reference, boxes = _load_boxes(boxes_path)
    receptor_format = receptor.suffix.lower().lstrip(".")
    ligand_format = ligand.suffix.lower().lstrip(".")
    if engine == "vina":
        if receptor_format != "pdbqt" or ligand_format != "pdbqt":
            raise ValueError("Vina requires receptor and ligand PDBQT inputs")
    else:
        if receptor_format not in {"pdb", "pdbqt"}:
            raise ValueError("GNINA receptor must be a PDB or PDBQT file")
        if ligand_format not in {"sdf", "pdbqt"}:
            raise ValueError("GNINA ligand must be an SDF or PDBQT file")
    _validate_receptor_identity(reference, receptor, pdbqt=receptor_format == "pdbqt")
    ligand_heavy_atoms = _input_heavy_atom_count(ligand, ligand_format)
    if ligand_heavy_atoms < 1:
        raise ValueError("ligand input must contain at least one heavy atom")
    if engine == "vina" and (not isinstance(vina, str) or not vina or "\0" in vina):
        raise ValueError("vina must be a nonempty executable name or path")
    if engine == "gnina" and (not isinstance(gnina, str) or not gnina or "\0" in gnina):
        raise ValueError("gnina must be a nonempty executable name or path")
    out_dir = Path(output_dir).resolve()
    if out_dir.exists():
        if not out_dir.is_dir() or any(out_dir.iterdir()):
            raise ValueError("output_dir must be new or an existing empty directory")
    else:
        out_dir.mkdir(parents=True)
    configs: dict[str, Path] = {}
    for box in boxes:
        config_path = out_dir / "boxes" / f"{box['id']}.txt"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        config_path.write_text(_config_text(box))
        configs[box["id"]] = config_path
    commands: list[dict] = []
    for box in boxes:
        for run_index in range(runs):
            run_seed = seed + len(commands)
            extension = "pdbqt" if engine == "vina" else "sdf"
            result_path = out_dir / f"{box['id']}_run_{run_index + 1:03d}.{extension}"
            command_builder = _command if engine == "vina" else _gnina_command
            executable = vina if engine == "vina" else gnina
            argv = command_builder(
                executable,
                receptor,
                ligand,
                configs[box["id"]],
                result_path,
                cpus,
                exhaustiveness,
                num_modes,
                run_seed,
            )
            commands.append(
                {
                    "box_id": box["id"],
                    "run": run_index + 1,
                    "seed": run_seed,
                    "argv": argv,
                    "output": str(result_path),
                }
            )
    (out_dir / "commands.json").write_text(json.dumps({"commands": commands}, indent=2) + "\n")
    script_lines = ["#!/usr/bin/env bash", "set -euo pipefail"]
    script_lines.extend(shlex.join(item["argv"]) for item in commands)
    script = out_dir / "dock_boxes.sh"
    script.write_text("\n".join(script_lines) + "\n")
    script.chmod(0o755)
    manifest = {
        "status": "inputs_prepared",
        "docking_completed": False,
        "engine": engine,
        "score_field": "vina_score" if engine == "vina" else "minimizedAffinity",
        "schema_version": document["schema_version"],
        "boxes_file": str(boxes_path),
        "reference_receptor_pdb": str(reference),
        "reference_receptor_pdb_sha256": _sha256(reference),
        "receptor_input": str(receptor),
        "receptor_input_format": receptor_format,
        "receptor_input_sha256": _sha256(receptor),
        "ligand_input": str(ligand),
        "ligand_input_format": ligand_format,
        "ligand_input_sha256": _sha256(ligand),
        "ligand_heavy_atoms": ligand_heavy_atoms,
        "source_frame": document.get("source_frame"),
        "box_ids": [box["id"] for box in boxes],
        "cpus": cpus,
        "exhaustiveness": exhaustiveness,
        "num_modes": num_modes,
        "runs_per_box": runs,
        "initial_seed": seed,
        "commands": len(commands),
    }
    manifest_path = out_dir / "manifest.json"
    if engine == "vina":
        manifest["receptor_pdbqt"] = str(receptor)
        manifest["ligand_pdbqt"] = str(ligand)
    else:
        manifest["receptor_pdbqt"] = None
        manifest["ligand_pdbqt"] = None
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    if not run:
        return manifest

    manifest["status"] = "docking_running"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    scored: list[dict] = []
    engine_label = "Vina" if engine == "vina" else "GNINA"
    try:
        executable = _resolve_binary(
            vina if engine == "vina" else gnina, "Vina" if engine == "vina" else "GNINA"
        )
        for item in commands:
            item["argv"][0] = executable
        (out_dir / "commands.json").write_text(json.dumps({"commands": commands}, indent=2) + "\n")
        for item in commands:
            stem = f"{item['box_id']}_run_{item['run']:03d}"
            result = subprocess.run(item["argv"], capture_output=True, text=True, check=False)
            (out_dir / f"{stem}.stdout.txt").write_text(result.stdout)
            (out_dir / f"{stem}.stderr.txt").write_text(result.stderr)
            output = Path(item["output"])
            if result.returncode != 0:
                raise RuntimeError(
                    f"{engine_label} failed for {stem} with exit code {result.returncode}"
                )
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(f"{engine_label} produced no nonempty output for {stem}")
            _validate_pose_heavy_atom_count(output, engine, ligand_heavy_atoms)
            score = _best_score(output) if engine == "vina" else _gnina_score(output)
            scored.append({**item, "score": score})
        best_per_box = []
        for box in boxes:
            best = min(
                (item for item in scored if item["box_id"] == box["id"]),
                key=lambda item: item["score"],
            )
            best_per_box.append(
                {
                    "box_id": box["id"],
                    "score": best["score"],
                    "best_pose": item_relpath(best["output"], out_dir),
                    "run": best["run"],
                    "seed": best["seed"],
                }
            )
        global_best = min(best_per_box, key=lambda item: item["score"])
        with (out_dir / "summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=("box_id", "score", "best_pose", "run", "seed")
            )
            writer.writeheader()
            writer.writerows(best_per_box)
            writer.writerow({**global_best, "box_id": "GLOBAL"})
        manifest.update(
            status="docking_completed",
            docking_completed=True,
            best_per_box=best_per_box,
            global_best=global_best,
        )
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        return manifest
    except Exception as exc:
        manifest.update(status="failed", docking_completed=False, error=str(exc))
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        raise


def item_relpath(path: str, root: Path) -> str:
    """Return a result path relative to the docking output directory."""
    return str(Path(path).relative_to(root))
