"""Ordered command-line entry points with one scientific worker for both backends."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

from .errors import PlacementError
from .execution import submit_slurm

ROOT = Path(__file__).resolve().parents[2]
STAGES = (
    "fragment",
    "build",
    "functionalize",
    "amber",
    "gromacs",
    "audit",
    "qm",
    "esp-import",
    "esp",
    "import",
    "check",
    "complex",
    "snapshot",
    "dock",
    "alignment",
)


def parser_for(stage: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"GO Builder: {stage}. Coordinates in Å; Amber/Gaussian remain external."
    )
    if stage == "fragment":
        source = parser.add_mutually_exclusive_group(required=True)
        source.add_argument("--input", type=Path)
        source.add_argument("--smiles")
        parser.add_argument(
            "--inspect",
            action="store_true",
            help="Show one-based atom map without generating geometry; use --json for full map",
        )
        parser.add_argument("--attachment-atom", type=int)
        parser.add_argument("--transfer-hydrogen", type=int)
        parser.add_argument("--count", type=int, default=1)
        parser.add_argument("--seed", type=int, default=0)
        parser.add_argument("--generate-3d", action="store_true")
        parser.add_argument("--allow-undefined-stereo", action="store_true")
        parser.add_argument("--no-optimize", action="store_true")
        parser.add_argument("--output", type=Path)
    elif stage == "build":
        parser.add_argument("--recipe", required=True, type=Path)
        parser.add_argument("--output", required=True, type=Path)
        parser.add_argument("--pdb", type=Path)
    elif stage == "functionalize":
        parser.add_argument("--model", required=True, type=Path)
        parser.add_argument("--recipe", required=True, type=Path)
        parser.add_argument("--output", required=True, type=Path)
        parser.add_argument("--pdb", type=Path)
        parser.add_argument(
            "--count",
            type=int,
            help="Number of grafts; random sites by default, fraction calculated",
        )
        parser.add_argument(
            "--id",
            nargs="+",
            type=int,
            help="Optional explicit epoxide group IDs; must match count if supplied",
        )
        parser.add_argument("--seed", type=int, help="Reproducible site/orientation selection seed")
    elif stage == "amber":
        parser.add_argument("--model", required=True, type=Path)
        parser.add_argument("--output-dir", required=True, type=Path)
        parser.add_argument(
            "--run", action="store_true", help="Execute external tleap after preparing its inputs"
        )
        parser.add_argument("--tleap", default="tleap")
        parser.add_argument("--reference-dir", type=Path)
        parser.add_argument("--parameter-templates", type=Path)
        parser.add_argument(
            "--radii", choices=("bondi", "mbondi", "mbondi2", "mbondi3"), default="mbondi2"
        )
    elif stage == "gromacs":
        parser.add_argument("--topology", required=True, type=Path)
        parser.add_argument("--coordinates", required=True, type=Path)
        parser.add_argument("--output-prefix", required=True, type=Path)
    elif stage == "audit":
        parser.add_argument("--output-dir", required=True, type=Path)
        parser.add_argument("--current-dir", type=Path)
        parser.add_argument("--old-dir", type=Path)
        parser.add_argument("--opls-dir", type=Path)
        parser.add_argument("--current-topology", type=Path)
        parser.add_argument("--current-go-residue-count", type=int, default=494)
    elif stage == "qm":
        parser.add_argument("--model", required=True, type=Path)
        parser.add_argument("--output-dir", required=True, type=Path)
        parser.add_argument("--total-charge", type=int)
        parser.add_argument("--multiplicity", type=int, default=1)
        parser.add_argument("--method", default="HF")
        parser.add_argument("--basis", default="6-31G*")
        parser.add_argument("--nproc", type=int, default=4)
        parser.add_argument("--qm-memory", default="4GB")
    elif stage == "esp-import":
        parser.add_argument("--esp-file", required=True, type=Path)
        parser.add_argument("--model", required=True, type=Path)
        parser.add_argument("--output", required=True, type=Path)
        parser.add_argument("--total-charge", type=int, required=True)
    elif stage == "esp":
        parser.add_argument("--reference", required=True, type=Path)
        parser.add_argument(
            "--candidates",
            required=True,
            type=Path,
            help="JSON mapping candidate names to atom/charge models",
        )
        parser.add_argument("--output-dir", required=True, type=Path)
    elif stage == "import":
        parser.add_argument("--topology", required=True, type=Path)
        parser.add_argument("--coordinates", type=Path)
        parser.add_argument("--output", required=True, type=Path)
        parser.add_argument("--pdb", type=Path)
    elif stage == "complex":
        parser.add_argument("--receptor-topology", required=True, type=Path)
        parser.add_argument("--receptor-coordinates", required=True, type=Path)
        parser.add_argument("--ligand-topology", required=True, type=Path, action="append")
        parser.add_argument("--ligand-coordinates", required=True, type=Path, action="append")
        parser.add_argument(
            "--ligand-reference-coordinates",
            type=Path,
            action="append",
            help="One ordered pre-docking reference per ligand, for SDF stereo validation",
        )
        parser.add_argument("--output-dir", required=True, type=Path)
    elif stage == "snapshot":
        parser.add_argument("--topology", required=True, type=Path)
        parser.add_argument("--coordinates", required=True, type=Path)
        parser.add_argument("--receptor-mask", default=":*")
        parser.add_argument("--frame", type=int, default=1)
        parser.add_argument("--box-size", nargs=3, type=float, default=[20.0, 20.0, 20.0])
        parser.add_argument("--overlap", type=float, default=5.0)
        parser.add_argument("--padding", type=float, default=5.0)
        parser.add_argument("--max-boxes", type=int, default=4096)
        parser.add_argument("--output-dir", required=True, type=Path)
    elif stage == "dock":
        parser.add_argument("--boxes", required=True, type=Path)
        parser.add_argument(
            "--receptor-pdbqt", "--receptor", dest="receptor_pdbqt", required=True, type=Path
        )
        parser.add_argument(
            "--ligand-pdbqt", "--ligand", dest="ligand_pdbqt", required=True, type=Path
        )
        parser.add_argument("--output-dir", required=True, type=Path)
        parser.add_argument("--engine", choices=("vina", "gnina"), default="vina")
        parser.add_argument("--vina", default="vina")
        parser.add_argument("--gnina", default="gnina")
        parser.add_argument("--exhaustiveness", type=int, default=16)
        parser.add_argument("--num-modes", type=int, default=10)
        parser.add_argument("--runs", type=int, default=3)
        parser.add_argument("--seed", type=int, default=1)
        parser.add_argument("--run", action="store_true")
    elif stage == "alignment":
        parser.add_argument("--paths", nargs="+", type=Path, required=True)
        parser.add_argument("--backup-dir", type=Path, required=True)
        parser.add_argument("--default-mask")
        parser.add_argument("--apply", action="store_true")
    else:
        parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--slurm", action="store_true", help="Submit the same command to Slurm")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show plan; with --slurm write worker without submission",
    )
    parser.add_argument("--submission-dir", type=Path, default=Path("go_jobs"))
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--memory", default="2G")
    parser.add_argument("--time", default="01:00:00")
    parser.add_argument("--partition")
    parser.add_argument("--account")
    parser.add_argument(
        "--setup", type=Path, help="Source a site software setup script before the worker"
    )
    parser.add_argument(
        "--python", dest="worker_python", help="Python interpreter available on the worker"
    )
    parser.add_argument("--json", action="store_true", help="Print machine-readable report")
    return parser


def _read(path: Path) -> dict:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _worker_argv(argv: list[str]) -> list[str]:
    result = []
    skip = False
    for token in argv:
        if skip:
            skip = False
        elif token in ("--slurm", "--dry-run"):
            continue
        elif token in ("--setup", "--python"):
            skip = True
        elif token.startswith(("--setup=", "--python=")):
            continue
        else:
            result.append(token)
    return result


def _report(stage: str, result: dict, machine: bool, error: bool = False) -> None:
    if machine:
        print(json.dumps(result, indent=2, default=str))
    else:
        stream = sys.stderr if error else sys.stdout
        timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S")
        level = "WARNING" if result.get("status") == "INCOMPLETE" else "ERROR" if error else "INFO"
        summary = {
            key: result[key]
            for key in (
                "status",
                "message",
                "output",
                "output_dir",
                "atoms",
                "atom_count",
                "bonds",
                "bond_count",
                "parameterized",
                "job_id",
                "worker",
                "receipt",
                "topology",
                "coordinates",
                "complex_topology",
                "n_grid_points",
                "composition",
                "group_counts",
                "formal_charge_e",
                "geometry_status",
                "requires_parameterization",
                "electron_parity",
                "functionalization",
                "requested",
                "achieved",
                "capacity",
                "attempts",
                "suggestions",
            )
            if result.get(key) is not None
        }
        summary.setdefault("status", "COMPLETED")
        if "datasets" in result:
            summary["datasets"] = len(result["datasets"])
        if "candidates" in result:
            summary["candidates"] = list(result["candidates"])
        if "gbsa" in result:
            summary["gbsa_ready"] = result["gbsa"]["ready"]
        print(
            f"[{timestamp}] [GO-{stage.upper()}] [{level}] "
            + json.dumps(summary, default=str, ensure_ascii=False),
            file=stream,
        )


def _dispatch(stage: str, args: argparse.Namespace) -> dict:
    if stage == "alignment":
        from .alignment import patch_inputs

        return patch_inputs(
            args.paths, args.backup_dir, apply=args.apply, default_mask=args.default_mask
        )
    if stage == "snapshot":
        from .snapshot import prepare_docking_receptor

        result = prepare_docking_receptor(
            args.topology,
            args.coordinates,
            args.output_dir,
            receptor_mask=args.receptor_mask,
            frame=args.frame,
            box_size=args.box_size,
            overlap=args.overlap,
            padding=args.padding,
            max_boxes=args.max_boxes,
        )
        result["output_dir"] = str(args.output_dir)
        return result
    if stage == "dock":
        from .docking import dock_boxes

        result = dock_boxes(
            args.boxes,
            args.receptor_pdbqt,
            args.ligand_pdbqt,
            args.output_dir,
            vina=args.vina,
            cpus=args.cpus,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            runs=args.runs,
            seed=args.seed,
            run=args.run,
            engine=args.engine,
            gnina=args.gnina,
        )
        result["output_dir"] = str(args.output_dir)
        return result
    if stage == "complex":
        from .complexes import assemble_complex

        result = assemble_complex(
            args.receptor_topology,
            args.receptor_coordinates,
            args.ligand_topology,
            args.ligand_coordinates,
            args.output_dir,
            ligand_reference_coordinates=args.ligand_reference_coordinates,
        )
        result["output_dir"] = str(args.output_dir)
        return result
    if stage == "fragment":
        from .fragments import inspect_fragment, prepare_fragment

        if args.output and args.output.exists():
            raise FileExistsError(args.output)
        if args.inspect:
            result = inspect_fragment(
                args.input, smiles=args.smiles, allow_undefined_stereo=args.allow_undefined_stereo
            )
        else:
            if args.attachment_atom is None or args.output is None:
                raise ValueError(
                    "Fragment preparation requires --attachment-atom and --output; use --inspect --json to choose the atom"
                )
            result = prepare_fragment(
                args.input,
                smiles=args.smiles,
                attachment_atom=args.attachment_atom,
                transfer_hydrogen=args.transfer_hydrogen,
                count=args.count,
                seed=args.seed,
                generate_3d=args.generate_3d,
                allow_undefined_stereo=args.allow_undefined_stereo,
                optimize=not args.no_optimize,
            )
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x") as handle:
                json.dump(result, handle, indent=2, allow_nan=False)
                handle.write("\n")
        if not args.inspect:
            return {
                "status": "PREPARED",
                "output": str(args.output),
                "atom_count": len(result["fragment"]["atoms"]),
                "attachment_atom": args.attachment_atom,
            }
        return result
    if stage in ("build", "functionalize", "import"):
        from .model import save_model, write_pdb

        # Check both destinations before writing either one.
        for path in (args.output, args.pdb):
            if path is not None and path.exists():
                raise FileExistsError(path)
        if stage == "build":
            from .build import build_sheet

            model = build_sheet(_read(args.recipe))
        elif stage == "functionalize":
            from .functionalize import open_epoxides

            recipe = _read(args.recipe)
            if args.id is not None:
                recipe["id"] = args.id
                # Explicit CLI sites define loading if --count is omitted,
                # instead of inheriting the fragment recipe's default count.
                if args.count is None:
                    recipe.pop("count", None)
            if args.count is not None:
                recipe["count"] = args.count
                if args.id is None:
                    recipe.pop("id", None)
            if args.seed is not None:
                recipe["seed"] = args.seed
            model = open_epoxides(_read(args.model), recipe)
        else:
            from .topology import import_amber

            model = import_amber(args.topology, args.coordinates)
        save_model(model, args.output)
        if args.pdb:
            write_pdb(model, args.pdb)
        metadata = model.get("metadata", {})
        return {
            "status": "BUILT" if stage == "build" else "PREPARED",
            "output": str(args.output),
            "atoms": len(model["atoms"]),
            "bonds": len(model["bonds"]),
            "composition": dict(sorted(Counter(a["element"] for a in model["atoms"]).items())),
            "group_counts": dict(
                sorted(Counter(g["type"] for g in metadata.get("groups", [])).items())
            ),
            "formal_charge_e": (
                None
                if metadata.get("formal_charges_unknown")
                else sum(a.get("formal_charge", 0) for a in model["atoms"])
            ),
            "geometry_status": metadata.get("geometry_status", "imported_coordinates"),
            "electron_parity": metadata.get("electron_parity", "unknown"),
            "requires_parameterization": metadata.get("requires_parameterization", False)
            or any(a.get("charge") is None or a.get("atom_type") is None for a in model["atoms"]),
            **(
                {"functionalization": metadata["functionalization_history"][-1]}
                if metadata.get("functionalization_history")
                else {}
            ),
        }
    if stage == "amber":
        from .topology import prepare_amber

        templates = _read(args.parameter_templates) if args.parameter_templates else None
        return prepare_amber(
            _read(args.model),
            args.output_dir,
            run=args.run,
            tleap=args.tleap,
            reference_dir=args.reference_dir,
            parameter_templates=templates,
            radii=args.radii,
        )
    if stage == "gromacs":
        from .topology import export_gromacs

        return export_gromacs(args.topology, args.coordinates, args.output_prefix)
    if stage == "audit":
        from .charges import audit_charge_sources

        result = audit_charge_sources(
            args.output_dir,
            current_dir=args.current_dir,
            old_dir=args.old_dir,
            opls_dir=args.opls_dir,
            current_topology=args.current_topology,
            current_go_residue_count=args.current_go_residue_count,
        )
        result["output_dir"] = str(args.output_dir)
        return result
    if stage == "qm":
        from .qm import prepare_qm

        result = prepare_qm(
            _read(args.model),
            args.output_dir,
            total_charge=args.total_charge,
            multiplicity=args.multiplicity,
            method=args.method,
            basis=args.basis,
            nproc=args.nproc,
            memory=args.qm_memory,
        )
        result["output_dir"] = str(args.output_dir)
        return result
    if stage == "esp-import":
        from .qm import espgen_to_reference

        result = espgen_to_reference(
            args.esp_file, _read(args.model), args.output, total_charge=args.total_charge
        )
        return {
            "status": "IMPORTED",
            "output": str(args.output),
            "atoms": len(result["atoms"]),
            "grid_points": len(result["grid"]),
        }
    if stage == "esp":
        from .qm import compare_esp

        candidates = _read(args.candidates)
        return compare_esp(args.reference, candidates, args.output_dir)
    from .model import validate_model
    from .topology import gbsa_readiness

    model = _read(args.model)
    return {"status": "CHECKED", "graph": validate_model(model), "gbsa": gbsa_readiness(model)}


def main(stage: str | None = None, argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if stage is None:
        top = argparse.ArgumentParser(
            description="GO Builder: explicit graph, anti opening, Amber and charge validation"
        )
        top.add_argument("stage", choices=STAGES)
        if not argv or argv[0] in ("-h", "--help"):
            top.print_help()
            return 0
        first = top.parse_args(argv[:1])
        stage = first.stage
        argv = argv[1:]
    args = parser_for(stage).parse_args(argv)
    command = [
        args.worker_python or sys.executable,
        str(ROOT / "go.py"),
        stage,
        *_worker_argv(argv),
    ]
    try:
        if args.slurm:
            result = submit_slurm(
                command,
                Path.cwd(),
                args.submission_dir,
                cpus=args.cpus,
                memory=args.memory,
                time=args.time,
                partition=args.partition,
                account=args.account,
                setup=args.setup,
                dry_run=args.dry_run,
            )
        elif args.dry_run:
            result = {"status": "PLANNED", "command": command, "backend": "local"}
        elif args.setup or args.worker_python:
            if args.setup and not args.setup.is_file():
                raise FileNotFoundError(args.setup)
            shell = (
                ("source " + shlex.quote(str(args.setup.resolve())) + "\n") if args.setup else ""
            )
            completed = subprocess.run(
                ["bash", "-c", "set -euo pipefail\n" + shell + "exec " + shlex.join(command)]
            )
            return completed.returncode
        else:
            result = _dispatch(stage, args)
        _report(stage, result, args.json)
        return 0
    except PlacementError as error:
        _report(stage, error.as_report(), args.json, error=True)
        return 1
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError, KeyError) as error:
        _report(stage, {"status": "ERROR", "message": str(error)}, args.json, error=True)
        return 1
