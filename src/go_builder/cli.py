"""Ordered command-line entry points with one scientific worker for both backends."""

from __future__ import annotations

import argparse
import json
import math
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
    "openmm",
    "lammps",
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
        from .presets import PRESET_NAMES

        src_group = parser.add_argument_group("Source / Config")
        source = src_group.add_mutually_exclusive_group(required=False)
        source.add_argument(
            "--config",
            "--recipe",
            dest="recipe",
            type=Path,
            metavar="PATH",
            help="Path to .toml or .json config/recipe file",
        )
        source.add_argument(
            "--preset",
            choices=PRESET_NAMES,
            help="Built-in preset or illustrative template name",
        )
        source.add_argument(
            "--list-presets",
            action="store_true",
            help="List available presets and templates and exit",
        )
        source.add_argument(
            "--init-config",
            type=Path,
            metavar="PATH",
            help="Write a complete annotated config template (.toml or .json) and exit",
        )
        src_group.add_argument(
            "--write-config",
            type=Path,
            metavar="PATH",
            help="Save the resolved recipe as a .toml or .json config file",
        )

        geom_group = parser.add_argument_group("Sheet Geometry")
        geom_group.add_argument(
            "--size",
            nargs=2,
            type=float,
            metavar=("LX", "LY"),
            help="Target physical sheet dimensions in nm (sets lx_nm, ly_nm)",
        )
        geom_group.add_argument(
            "--diam",
            type=float,
            metavar="D",
            help="Nanoflake diameter in nm (for hexagonal or circular shape)",
        )
        geom_group.add_argument(
            "--nx",
            type=int,
            help="Honeycomb cell count along x (>= 2; default: 6 if no size/diam)",
        )
        geom_group.add_argument(
            "--ny",
            type=int,
            help="Honeycomb cell count along y (>= 2; default: 6 if no size/diam)",
        )
        geom_group.add_argument(
            "--shape",
            choices=("rectangular", "hexagonal", "circular"),
            help="Sheet shape (default: rectangular)",
        )
        geom_group.add_argument(
            "--edge",
            dest="edge_type",
            choices=("zigzag", "armchair"),
            help="Hexagonal nanoflake edge termination (default: zigzag)",
        )
        geom_group.add_argument(
            "--pbc",
            action="store_true",
            default=None,
            help="Enable 2D periodic boundary conditions [True, True, False] (default: False)",
        )
        geom_group.add_argument(
            "--name",
            type=str,
            help="Structure label in metadata and exports (default: GO)",
        )
        geom_group.add_argument(
            "--seed",
            type=int,
            help="Random seed >= 0 (default: 0, or overrides config/preset seed)",
        )

        holes_group = parser.add_argument_group("Holes / Pores")
        holes_group.add_argument(
            "--holes",
            type=str,
            metavar="SPEC",
            help="Vacancy hole profile: DIST:PCT (e.g. mixed:8, small:5) or DIST:Nx (e.g. medium:3x)",
        )

        comp_group = parser.add_argument_group("Composition & Groups")
        comp_group.add_argument(
            "--oc",
            dest="o_over_c",
            type=float,
            help="Target atomic O/C ratio (mutually exclusive with --co)",
        )
        comp_group.add_argument(
            "--o-over-c",
            dest="o_over_c",
            type=float,
            help=argparse.SUPPRESS,
        )
        comp_group.add_argument(
            "--co",
            dest="c_over_o",
            type=float,
            help="Target atomic C/O ratio (mutually exclusive with --oc)",
        )
        comp_group.add_argument(
            "--c-over-o",
            dest="c_over_o",
            type=float,
            help=argparse.SUPPRESS,
        )
        comp_group.add_argument(
            "--ph",
            type=float,
            help="Target pH for carboxyl protonation via Henderson-Hasselbalch (pKa=4.5)",
        )
        comp_group.add_argument(
            "--priority",
            dest="group_priority",
            nargs="+",
            metavar="G",
            help="Ordered group priority list (e.g. OH epoxide COOH COO)",
        )
        comp_group.add_argument(
            "--group-priority",
            dest="group_priority",
            nargs="+",
            help=argparse.SUPPRESS,
        )
        comp_group.add_argument(
            "--comp-mode",
            dest="composition_mode",
            choices=("balanced", "random"),
            help="Composition resolution mode for O/C or C/O targets (default: balanced)",
        )
        comp_group.add_argument(
            "--composition-mode",
            dest="composition_mode",
            choices=("balanced", "random"),
            help=argparse.SUPPRESS,
        )
        comp_group.add_argument(
            "--groups",
            type=str,
            metavar="SPEC",
            help="Explicit group counts or percentages (e.g. OH=10,epoxide=8,COO=2 or OH=15%%,epoxide=10%%)",
        )

        place_group = parser.add_argument_group("Placement & Electronic Filter")
        place_group.add_argument(
            "--place",
            dest="placement",
            choices=("random", "island", "irregular_island", "dispersed"),
            help="Functional group spatial placement mode (default: random)",
        )
        place_group.add_argument(
            "--placement",
            dest="placement",
            choices=("random", "island", "irregular_island", "dispersed"),
            help=argparse.SUPPRESS,
        )
        place_group.add_argument(
            "--islands",
            type=int,
            metavar="N",
            help="Number of island nuclei when using --place island or irregular_island (default: 1)",
        )
        place_group.add_argument(
            "--island-style",
            dest="island_style",
            choices=("compact", "irregular"),
            help="Island morphology when using --place island (default: compact; forced irregular for irregular_island)",
        )
        place_group.add_argument(
            "--branching",
            type=float,
            metavar="F",
            help="Branching intensity in [0.0, 1.0] for irregular islands (default: 0.65)",
        )
        place_group.add_argument(
            "--clar",
            action=argparse.BooleanOptionalAction,
            default=None,
            help=(
                "Enforce Clar closed-shell Kekule pi-system balance (even bipartite sublattice parity, "
                "perfect Kekule matching, nonzero Huckel gap) across all placement modes (default: ON; disable with --no-clar)"
            ),
        )
        place_group.add_argument(
            "--min-huckel-gap",
            dest="min_huckel_gap",
            type=float,
            metavar="BETA",
            help="Minimum Huckel HOMO-LUMO pi gap (|beta| units) when --clar is active (default: 0.0)",
        )
        place_group.add_argument(
            "--faces",
            choices=("both", "up", "down"),
            help="Basal functionalization faces (default: both)",
        )

        out_group = parser.add_argument_group("Output & Log")
        out_group.add_argument(
            "--out",
            dest="output",
            type=Path,
            metavar="PATH",
            help="Output GO model JSON path",
        )
        out_group.add_argument(
            "--output",
            dest="output",
            type=Path,
            help=argparse.SUPPRESS,
        )
        out_group.add_argument("--pdb", type=Path, metavar="PATH", help="Export finite PDB (.pdb)")
        out_group.add_argument("--xyz", type=Path, metavar="PATH", help="Export non-periodic XYZ (.xyz)")
        out_group.add_argument(
            "--extxyz",
            type=Path,
            metavar="PATH",
            help="Export Extended XYZ (.extxyz, periodic or finite)",
        )
        out_group.add_argument("--mol2", type=Path, metavar="PATH", help="Export finite Tripos MOL2 (.mol2)")
        out_group.add_argument("--sdf", type=Path, metavar="PATH", help="Export finite Molfile V3000 (.sdf)")
        out_group.add_argument("--cif", type=Path, metavar="PATH", help="Export P1 CIF (.cif)")
        out_group.add_argument(
            "--log",
            type=Path,
            metavar="PATH",
            help="Write full build report, ISO timestamp, and effective config to file",
        )
        out_group.add_argument(
            "--dry-run",
            action="store_true",
            help="Show a local creation plan without building",
        )
        out_group.add_argument(
            "--json",
            action="store_true",
            help="Print machine-readable JSON report",
        )
        return parser
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
        parser.add_argument("--model", type=Path)
        parser.add_argument("--output-dir", type=Path)
        parser.add_argument(
            "--preset",
            choices=("gaff2_v3", "gaff2_current"),
            default=None,
            help="AMBER force-field preset (default: gaff2_v3, or gaff2_current if --reference-dir is set)",
        )
        parser.add_argument(
            "--list-presets",
            action="store_true",
            help="List packaged AMBER force-field presets and exit",
        )
        parser.add_argument(
            "--charge-model",
            choices=("auto", "l0_v2", "current_lib", "custom"),
            default="auto",
            help="Charge model override (auto selects preset default)",
        )
        parser.add_argument(
            "--charges-json",
            type=Path,
            help="Explicit per-atom charges JSON for --charge-model custom",
        )
        parser.add_argument(
            "--run", action="store_true", help="Execute external tleap after preparing its inputs"
        )
        parser.add_argument("--tleap", default="tleap")
        parser.add_argument("--reference-dir", type=Path)
        parser.add_argument("--parameter-templates", type=Path)
        parser.add_argument(
            "--radii", choices=("bondi", "mbondi", "mbondi2", "mbondi3"), default="mbondi2"
        )
        parser.add_argument(
            "--direct-periodic",
            action="store_true",
            help="Directly build periodic Amber prmtop/inpcrd from gaff2_v3 frcmod without tleap",
        )
    elif stage == "gromacs":
        parser.add_argument("--topology", type=Path)
        parser.add_argument("--coordinates", type=Path)
        parser.add_argument("--model", type=Path)
        parser.add_argument("--output-prefix", required=True, type=Path)
        parser.add_argument(
            "--charge-model",
            choices=("auto", "l0_v2", "current_lib", "custom"),
            default="auto",
        )
        parser.add_argument("--charges-json", type=Path)
        parser.add_argument("--reference-dir", type=Path)
        parser.add_argument("--cutoff-nm", type=float, default=None)
    elif stage == "openmm":
        parser.add_argument("--topology", type=Path)
        parser.add_argument("--coordinates", type=Path)
        parser.add_argument("--model", type=Path)
        parser.add_argument("--output", required=True, type=Path)
        parser.add_argument(
            "--nonbonded-method",
            choices=("NoCutoff", "CutoffNonPeriodic", "PME"),
            default=None,
        )
        parser.add_argument("--cutoff-nm", type=float, default=None)
        parser.add_argument(
            "--charge-model",
            choices=("auto", "l0_v2", "current_lib", "custom"),
            default="auto",
        )
        parser.add_argument("--charges-json", type=Path)
        parser.add_argument("--reference-dir", type=Path)
    elif stage == "lammps":
        parser.add_argument("--topology", type=Path)
        parser.add_argument("--coordinates", type=Path)
        parser.add_argument("--model", type=Path)
        parser.add_argument("--output-prefix", required=True, type=Path)
        parser.add_argument("--cutoff-ang", type=float, default=None)
        parser.add_argument(
            "--charge-model",
            choices=("auto", "l0_v2", "current_lib", "custom"),
            default="auto",
        )
        parser.add_argument("--charges-json", type=Path)
        parser.add_argument("--reference-dir", type=Path)
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
    elif stage == "check":
        parser.add_argument("--model", required=True, type=Path)
        parser.add_argument(
            "--strict",
            action="store_true",
            help="Return exit code 1 if geometry integrity checks fail",
        )
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


def _validate_write_config_path(path: Path | None) -> None:
    if path is None:
        return
    suffix = path.suffix.lower()
    if suffix not in (".toml", ".json"):
        raise ValueError(
            f"Unsupported config extension {suffix!r}; use .toml or .json"
        )


def _validate_output_paths(paths: list[Path | None]) -> None:
    active = [p for p in paths if p is not None]
    resolved = [p.resolve() for p in active]
    if len(set(resolved)) != len(resolved):
        raise ValueError("All output paths must be distinct")
    for path in active:
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)


def _parse_holes_spec(spec: str) -> dict:
    parts = [p.strip() for p in spec.split(":")]
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(
            f"Invalid --holes spec {spec!r}; expected DIST:PCT (e.g. mixed:8) or DIST:Nx (e.g. medium:3x)"
        )
    dist, rhs = parts
    allowed_dists = ("small", "medium", "large", "mixed", "random")
    if dist not in allowed_dists:
        raise ValueError(
            f"Invalid --holes distribution {dist!r}; choose from {', '.join(allowed_dists)}"
        )
    if rhs.lower().endswith("x"):
        count_str = rhs[:-1].strip()
        try:
            count = int(count_str)
        except ValueError as exc:
            raise ValueError(f"Invalid --holes count in {spec!r}") from exc
        if count < 1:
            raise ValueError("--holes count must be >= 1")
        return {"distribution": dist, "count": count}
    pct_str = rhs[:-1].strip() if rhs.endswith("%") else rhs
    try:
        pct = float(pct_str)
    except ValueError as exc:
        raise ValueError(f"Invalid --holes percentage in {spec!r}") from exc
    if not (0.0 < pct <= 50.0):
        raise ValueError("--holes percent_vacancies must be in (0, 50]")
    return {"distribution": dist, "percent_vacancies": pct}


def _parse_groups_spec(spec: str) -> dict:
    allowed_groups = {
        "OH",
        "epoxide",
        "carbonyl",
        "COOH",
        "COO",
        "epoxide_eph",
        "epoxide_ept",
    }
    items = [item.strip() for item in spec.split(",") if item.strip()]
    if not items:
        raise ValueError("--groups requires at least one GROUP=VALUE entry")
    parsed: dict = {}
    for item in items:
        if "=" not in item:
            raise ValueError(
                f"Invalid --groups entry {item!r}; expected GROUP=COUNT or GROUP=PCT%"
            )
        group, raw_val = (p.strip() for p in item.split("=", 1))
        if group not in allowed_groups:
            raise ValueError(
                f"Unknown group {group!r} in --groups; choose from {', '.join(sorted(allowed_groups))}"
            )
        if group in parsed:
            raise ValueError(f"Duplicate group {group!r} in --groups")
        if raw_val.endswith("%"):
            try:
                pct = float(raw_val[:-1].strip())
            except ValueError as exc:
                raise ValueError(f"Invalid percentage for {group!r} in --groups") from exc
            if not (0.0 <= pct <= 100.0):
                raise ValueError(f"Percentage for {group!r} must be between 0 and 100")
            parsed[group] = {"percent": pct}
        else:
            try:
                cnt = int(raw_val)
            except ValueError as exc:
                raise ValueError(f"Invalid integer count for {group!r} in --groups") from exc
            if cnt < 0:
                raise ValueError(f"Count for {group!r} must be >= 0")
            parsed[group] = cnt
    return parsed


def _check_build_cli_conflicts(args: argparse.Namespace) -> None:
    if args.c_over_o is not None and args.o_over_c is not None:
        raise ValueError("--co/--c-over-o and --oc/--o-over-c are mutually exclusive")
    if args.size is not None and args.diam is not None:
        raise ValueError("--size and --diam are mutually exclusive")
    if args.size is not None and (args.nx is not None or args.ny is not None):
        raise ValueError("--size and --nx/--ny are mutually exclusive")
    if args.diam is not None and (args.nx is not None or args.ny is not None):
        raise ValueError("--diam and --nx/--ny are mutually exclusive")
    if args.islands is not None:
        if args.islands < 1:
            raise ValueError("--islands must be >= 1")
        if args.placement is not None and args.placement not in ("island", "irregular_island"):
            raise ValueError("--islands requires --place island or irregular_island")
    if getattr(args, "island_style", None) is not None:
        if args.placement is not None and args.placement not in ("island", "irregular_island"):
            raise ValueError("--island-style requires --place island or irregular_island")
        if args.placement == "irregular_island" and args.island_style != "irregular":
            raise ValueError("--island-style compact conflicts with --place irregular_island")
    if getattr(args, "branching", None) is not None:
        if not math.isfinite(args.branching) or not 0.0 <= args.branching <= 1.0:
            raise ValueError("--branching must be between 0 and 1")
        if args.placement is not None and args.placement not in ("island", "irregular_island"):
            raise ValueError("--branching requires --place island or irregular_island")
        if getattr(args, "island_style", None) == "compact":
            raise ValueError("--branching requires --place irregular_island or --island-style irregular")
    if getattr(args, "min_huckel_gap", None) is not None:
        if getattr(args, "clar", None) is False:
            raise ValueError("--min-huckel-gap cannot be combined with --no-clar")
        if not math.isfinite(args.min_huckel_gap) or args.min_huckel_gap < 0.0:
            raise ValueError("--min-huckel-gap must be >= 0")


def _resolve_build_recipe(args: argparse.Namespace) -> tuple[dict, dict | None]:
    _check_build_cli_conflicts(args)
    provenance = None
    cli_overrides: dict = {}
    if args.preset:
        from .presets import resolve_preset

        recipe, provenance = resolve_preset(args.preset, seed=args.seed)
        if args.seed is not None:
            cli_overrides["seed"] = args.seed
    elif args.recipe:
        from .config import load_config

        recipe = load_config(args.recipe)
        if args.seed is not None:
            recipe["seed"] = args.seed
            cli_overrides["seed"] = args.seed
    else:
        recipe = {}
        if args.seed is not None:
            recipe["seed"] = args.seed
            cli_overrides["seed"] = args.seed

    if args.name is not None:
        recipe["name"] = args.name
        cli_overrides["name"] = args.name
    if args.shape is not None:
        recipe["shape"] = args.shape
        cli_overrides["shape"] = args.shape
        if args.shape != "hexagonal" and args.edge_type is None:
            recipe.pop("edge_type", None)
    if args.edge_type is not None:
        recipe["edge_type"] = args.edge_type
        cli_overrides["edge_type"] = args.edge_type
    if args.pbc:
        recipe["pbc"] = [True, True, False]
        cli_overrides["pbc"] = [True, True, False]
    if args.faces is not None:
        recipe["faces"] = args.faces
        cli_overrides["faces"] = args.faces

    if args.size is not None:
        lx, ly = float(args.size[0]), float(args.size[1])
        if lx <= 0.0 or ly <= 0.0:
            raise ValueError("--size dimensions must be positive")
        for key in ("nx", "ny", "diameter_nm", "radius_nm"):
            recipe.pop(key, None)
        recipe["lx_nm"] = lx
        recipe["ly_nm"] = ly
        cli_overrides["lx_nm"] = lx
        cli_overrides["ly_nm"] = ly
    elif args.diam is not None:
        d = float(args.diam)
        if d <= 0.0:
            raise ValueError("--diam must be positive")
        for key in ("nx", "ny", "lx_nm", "ly_nm", "radius_nm"):
            recipe.pop(key, None)
        recipe["diameter_nm"] = d
        cli_overrides["diameter_nm"] = d
    elif args.nx is not None or args.ny is not None:
        for key in ("lx_nm", "ly_nm", "diameter_nm", "radius_nm"):
            recipe.pop(key, None)
        if args.nx is not None:
            recipe["nx"] = args.nx
            cli_overrides["nx"] = args.nx
        if args.ny is not None:
            recipe["ny"] = args.ny
            cli_overrides["ny"] = args.ny

    if args.holes is not None:
        parsed_holes = _parse_holes_spec(args.holes)
        recipe["holes"] = parsed_holes
        cli_overrides["holes"] = parsed_holes

    island_style = getattr(args, "island_style", None)
    branching = getattr(args, "branching", None)
    if (
        args.placement is not None
        or args.islands is not None
        or island_style is not None
        or branching is not None
    ):
        existing_placement = recipe.get("placement")
        has_island_opts = (
            args.islands is not None or island_style is not None or branching is not None
        )
        if args.placement in ("random", "dispersed"):
            recipe["placement"] = args.placement
        elif has_island_opts or args.placement in ("island", "irregular_island"):
            if isinstance(existing_placement, dict) and existing_placement.get("mode") in (
                "island",
                "irregular_island",
            ):
                updated = dict(existing_placement)
            elif existing_placement in ("island", "irregular_island"):
                updated = {"mode": existing_placement}
            else:
                updated = {}
            target_mode = (
                args.placement
                if args.placement in ("island", "irregular_island")
                else updated.get("mode", "island")
            )
            if args.placement is None and existing_placement in ("random", "dispersed"):
                raise ValueError(
                    f"Cannot combine island flags with placement mode {existing_placement!r}"
                )
            updated["mode"] = target_mode
            if args.islands is not None:
                updated["n_islands"] = args.islands
            if island_style is not None:
                updated["style"] = island_style
            if branching is not None:
                eff_style = updated.get(
                    "style",
                    "irregular" if target_mode == "irregular_island" else "compact",
                )
                if eff_style != "irregular":
                    raise ValueError(
                        "--branching requires --place irregular_island or --island-style irregular"
                    )
                updated["branching"] = branching
            if (
                len(updated) == 1
                and target_mode in ("island", "irregular_island")
                and args.placement is not None
            ):
                recipe["placement"] = target_mode
            else:
                recipe["placement"] = updated
        cli_overrides["placement"] = recipe["placement"]

    use_clar = getattr(args, "clar", None)
    min_huckel_gap = getattr(args, "min_huckel_gap", None)
    if use_clar is not None or min_huckel_gap is not None:
        existing_clar = recipe.get("clar")
        clar_obj = dict(existing_clar) if isinstance(existing_clar, dict) else {}
        if use_clar is False:
            clar_obj["enabled"] = False
        elif use_clar is True or min_huckel_gap is not None:
            clar_obj["enabled"] = True
        if min_huckel_gap is not None:
            clar_obj["min_huckel_gap"] = float(min_huckel_gap)
        recipe["clar"] = clar_obj
        cli_overrides["clar"] = clar_obj
    elif "clar" not in recipe or recipe["clar"] is None:
        from .presets import LAYERONE_PRESET_NAMES

        if args.preset not in LAYERONE_PRESET_NAMES:
            recipe["clar"] = {"enabled": True}

    has_comp_cli = any(
        val is not None
        for val in (
            args.c_over_o,
            args.o_over_c,
            args.ph,
            args.composition_mode,
            args.group_priority,
        )
    )
    if args.groups is not None:
        if has_comp_cli:
            raise ValueError(
                "Cannot combine --groups with --oc/--co/--ph/--priority/--comp-mode"
            )
        if (
            recipe.get("composition") is not None
            or recipe.get("c_over_o") is not None
            or recipe.get("o_over_c") is not None
        ):
            raise ValueError(
                "Cannot combine --groups with a config/recipe that specifies composition"
            )
        parsed_groups = _parse_groups_spec(args.groups)
        recipe["groups"] = parsed_groups
        cli_overrides["groups"] = parsed_groups

    if has_comp_cli:
        if recipe.get("groups"):
            raise ValueError(
                "Cannot combine --c-over-o/--o-over-c composition flags with a recipe that already specifies groups"
            )
        if "composition" in recipe and recipe["composition"] is not None:
            if not isinstance(recipe["composition"], dict):
                raise ValueError("Recipe composition must be a JSON/TOML object")
            comp_cfg = dict(recipe["composition"])
        else:
            comp_cfg = {}
        if "c_over_o" in recipe:
            comp_cfg.setdefault("c_over_o", recipe.pop("c_over_o"))
        if "o_over_c" in recipe:
            comp_cfg.setdefault("o_over_c", recipe.pop("o_over_c"))
        comp_ovr: dict = {}
        if args.c_over_o is not None:
            comp_cfg.pop("O_over_C", None)
            comp_cfg.pop("o_over_c", None)
            comp_cfg.pop("C_over_O", None)
            comp_cfg["c_over_o"] = args.c_over_o
            comp_ovr["c_over_o"] = args.c_over_o
        if args.o_over_c is not None:
            comp_cfg.pop("C_over_O", None)
            comp_cfg.pop("c_over_o", None)
            comp_cfg.pop("O_over_C", None)
            comp_cfg["o_over_c"] = args.o_over_c
            comp_ovr["o_over_c"] = args.o_over_c
        if args.ph is not None:
            if (
                "deprotonation_fraction" in comp_cfg
                and comp_cfg["deprotonation_fraction"] is not None
            ):
                raise ValueError(
                    "Cannot specify --ph when recipe composition already sets deprotonation_fraction"
                )
            comp_cfg["ph"] = args.ph
            comp_ovr["ph"] = args.ph
        if args.composition_mode is not None:
            comp_cfg["mode"] = args.composition_mode
            comp_ovr["mode"] = args.composition_mode
        if args.group_priority is not None:
            if (
                comp_cfg.get("preferences") is not None
                or comp_cfg.get("group_preferences") is not None
            ):
                raise ValueError(
                    "Cannot specify --group-priority when recipe composition already sets preferences"
                )
            comp_cfg["priority"] = list(args.group_priority)
            comp_ovr["priority"] = list(args.group_priority)
        if not any(
            comp_cfg.get(k) is not None
            for k in ("C_over_O", "c_over_o", "O_over_C", "o_over_c")
        ):
            raise ValueError(
                "--ph, --comp-mode, and --priority require --co, --oc, or a config with composition"
            )
        recipe["composition"] = comp_cfg
        if comp_ovr:
            cli_overrides["composition"] = comp_ovr

    if provenance is not None and cli_overrides:
        provenance["cli_overrides"] = cli_overrides

    return recipe, provenance


def _edge_motifs_summary(model: dict) -> dict:
    morphology = model.get("metadata", {}).get("morphology")
    if isinstance(morphology, dict) and "armchair_edge_pairs" in morphology:
        return {
            "perimeter_carbons": morphology.get("perimeter_carbons", 0),
            "armchair_edge_pairs": morphology.get("armchair_edge_pairs", 0),
            "zigzag_edge_carbons": morphology.get("zigzag_edge_carbons", 0),
        }
    atoms = model["atoms"]
    scaffold_ids = {
        a["id"]
        for a in atoms
        if a["element"] == "C" and a.get("metadata", {}).get("role") == "sheet_carbon"
    }
    if not scaffold_ids:
        scaffold_ids = {a["id"] for a in atoms if a["element"] == "C"}
    scaffold_adj: dict[int, list[int]] = {cid: [] for cid in scaffold_ids}
    for bond in model["bonds"]:
        a_id, b_id = bond["a"], bond["b"]
        if a_id in scaffold_adj and b_id in scaffold_adj:
            scaffold_adj[a_id].append(b_id)
            scaffold_adj[b_id].append(a_id)
    edge_set = {cid for cid, nbrs in scaffold_adj.items() if len(nbrs) == 2}
    armchair_pairs = 0
    zigzag_carbons = 0
    seen: set[int] = set()
    for cid in sorted(edge_set):
        if cid in seen:
            continue
        comp = []
        stack = [cid]
        seen.add(cid)
        while stack:
            curr = stack.pop()
            comp.append(curr)
            for nbr in scaffold_adj[curr]:
                if nbr in edge_set and nbr not in seen:
                    seen.add(nbr)
                    stack.append(nbr)
        armchair_pairs += len(comp) // 2
        zigzag_carbons += len(comp) % 2
    return {
        "perimeter_carbons": len(edge_set),
        "armchair_edge_pairs": armchair_pairs,
        "zigzag_edge_carbons": zigzag_carbons,
    }


def _build_display_context(model: dict, recipe: dict, written_files: dict[str, str]) -> dict:
    from .model import check_geometry_integrity, compute_chemistry_summary

    metadata = model.get("metadata", {})
    dims = metadata.get("dimensions") or {}
    morph = metadata.get("morphology") or {}
    nx = dims.get("nx", recipe.get("nx", 6))
    ny = dims.get("ny", recipe.get("ny", 6))
    residues = len({(int(a.get("resid", 1)), str(a.get("resname", "GO"))) for a in model["atoms"]})
    return {
        "name": model.get("name", recipe.get("name", "GO")),
        "shape": morph.get("shape", recipe.get("shape", "rectangular")),
        "edge_type": morph.get("edge_type", recipe.get("edge_type")),
        "seed": metadata.get("seed", recipe.get("seed", 0)),
        "geometry_status": metadata.get("geometry_status", "unknown"),
        "nx": nx,
        "ny": ny,
        "atoms": len(model["atoms"]),
        "bonds": len(model["bonds"]),
        "residues": residues,
        "edge_motifs": _edge_motifs_summary(model),
        "hole_summary": metadata.get("hole_summary"),
        "group_counts": dict(
            sorted(Counter(g["type"] for g in metadata.get("groups", [])).items())
        ),
        "placement": metadata.get("placement") or recipe.get("placement", "random"),
        "faces": recipe.get("faces", "both"),
        "chemistry": compute_chemistry_summary(model),
        "geometry": check_geometry_integrity(model),
        "written_files": dict(written_files),
        "effective_recipe": recipe,
    }


def _fmt_num(value: float | None, fmt: str = ".3f", fallback: str = "n/a") -> str:
    if value is None:
        return fallback
    return format(value, fmt)


def _format_chemistry_and_qm_lines(
    chem: dict,
    geom: dict,
    group_counts: dict | None = None,
) -> list[str]:
    lines: list[str] = []
    c_over_o = _fmt_num(chem.get("c_over_o"), ".3f", "inf (no O)")
    o_over_c = _fmt_num(chem.get("o_over_c"), ".3f", "0.000")
    q_val = chem.get("formal_charge")
    q_str = "unknown" if q_val is None else f"{q_val:+d} e" if q_val != 0 else "0 e"
    lines.append("Composition & Chemistry")
    lines.append(
        f"  Formula        : {chem.get('formula', 'n/a')}   |   C/O: {c_over_o}   |   "
        f"O/C: {o_over_c}   |   Charge: {q_str}"
    )
    if group_counts is not None:
        nonzero = {k: v for k, v in group_counts.items() if v}
        grp_str = (
            ", ".join(f"{k}={v}" for k, v in nonzero.items())
            if nonzero
            else "none (pristine scaffold)"
        )
        lines.append(f"  Groups         : {grp_str}")
    sp2_pct = 100.0 * float(chem.get("sp2_fraction", 0.0))
    sp3_pct = 100.0 * float(chem.get("sp3_fraction", 0.0))
    dom = chem.get("sp2_domains", {})
    lines.append(
        f"  Hybridization  : sp2 {sp2_pct:.1f}% ({chem.get('sp2_carbons', 0)} C) / "
        f"sp3 {sp3_pct:.1f}% ({chem.get('sp3_carbons', 0)} C)"
    )
    lines.append(
        f"  sp2 Domains    : count={dom.get('count', 0)}, "
        f"largest={dom.get('largest_domain_carbons', 0)} C, "
        f"mean={float(dom.get('mean_domain_carbons', 0.0)):.1f} C"
    )
    envs = chem.get("chemical_environments", {})
    ec = envs.get("counts", {})
    ef = envs.get("fractions", {})
    cbx = envs.get("carboxyl_breakdown", {})
    lines.append(
        "  C Environments : "
        f"C-C sp2={ec.get('C-C_sp2', 0)} ({100.0 * ef.get('C-C_sp2', 0.0):.1f}%), "
        f"C-O={ec.get('C-O_sp3', 0)} ({100.0 * ef.get('C-O_sp3', 0.0):.1f}%), "
        f"C=O={ec.get('C=O_carbonyl', 0)} ({100.0 * ef.get('C=O_carbonyl', 0.0):.1f}%), "
        f"O-C=O={ec.get('O-C=O_carboxyl', 0)} ({100.0 * ef.get('O-C=O_carboxyl', 0.0):.1f}%)"
    )
    anchor_cy = cbx.get("anchor_CY", 0)
    if anchor_cy:
        lines.append(
            f"                   (anchor CY={anchor_cy} included in C-C sp2; "
            "structural graph environments, not XPS intensities)"
        )
    else:
        lines.append("                   (structural graph environments, not XPS intensities)")
    lines.append("")

    qm = chem.get("qm_diagnostics", {})
    gate = "PASS" if geom.get("passed") else "FAIL"
    lines.append("Integrity Gate & QM Diagnostics")
    min_b = _fmt_num(geom.get("min_bond_length"), ".3f")
    max_b = _fmt_num(geom.get("max_bond_length"), ".3f")
    lines.append(
        f"  Integrity Gate : {gate}   |   Clashes: {geom.get('clash_count', 0)}   |   "
        f"Bonds: [{min_b}, {max_b}] A ({'valid' if geom.get('bond_lengths_valid') else 'INVALID'})"
    )
    min_h = _fmt_num(geom.get("min_nonbonded_h_distance"), ".3f")
    min_hv = _fmt_num(geom.get("min_nonbonded_heavy_distance"), ".3f")
    cutoff = geom.get("nonbonded_minima_cutoff")
    cutoff_suffix = (
        f" (within {float(cutoff):.1f} A cutoff)"
        if geom.get("nonbonded_minima_cutoff_limited") and cutoff is not None
        else ""
    )
    lines.append(f"  Min Nonbonded  : H={min_h} A   |   Heavy={min_hv} A{cutoff_suffix}")
    imb = qm.get("sp2_sublattice_imbalance")
    sc = qm.get("sp2_sublattice_counts")
    imb_str = (
        f"|N_A - N_B|={imb} (N_A={sc['A']}, N_B={sc['B']})"
        if imb is not None and sc is not None
        else "n/a"
    )
    lines.append(
        f"  Electrons      : {qm.get('electron_count', geom.get('electron_count', 'n/a'))} "
        f"({qm.get('electron_parity', geom.get('electron_parity', 'unknown'))})   |   "
        f"Sublattice: {imb_str}"
    )
    for warn in geom.get("warnings") or []:
        lines.append(f"  Warning        : {warn}")
    return lines


def format_build_report(ctx: dict) -> str:
    chem = ctx["chemistry"]
    geom = ctx["geometry"]
    dims = chem.get("dimensions_nm", {})
    ext = dims.get("extent_nm", [0.0, 0.0, 0.0])
    cell = dims.get("cell_nm")
    pbc = dims.get("periodic", [False, False, False])
    motifs = ctx.get("edge_motifs", {})

    shape_str = str(ctx.get("shape", "rectangular"))
    if shape_str == "hexagonal" and ctx.get("edge_type"):
        shape_str = f"{shape_str} ({ctx['edge_type']})"

    lines = [
        "=== GO Builder — Build Summary ===",
        f"Name: {ctx['name']}   |   Shape: {shape_str}   |   Seed: {ctx['seed']}   |   Status: {ctx['geometry_status']}",
        "",
        "Geometry & Dimensions",
        f"  Grid (nx x ny) : {ctx.get('nx', '?')} x {ctx.get('ny', '?')}   |   PBC: {pbc}",
        f"  Extent (nm)    : {ext[0]:.3f} x {ext[1]:.3f} x {ext[2]:.3f}",
        (
            f"  Cell (nm)      : {cell[0]:.3f} x {cell[1]:.3f} x {cell[2]:.3f}"
            if cell is not None
            else "  Cell (nm)      : none"
        ),
        f"  Topology       : {ctx['atoms']} atoms, {ctx['bonds']} bonds, {ctx['residues']} residues",
        (
            f"  Edge Motifs    : armchair_pairs={motifs.get('armchair_edge_pairs', 0)}, "
            f"zigzag_carbons={motifs.get('zigzag_edge_carbons', 0)}, "
            f"perimeter_C={motifs.get('perimeter_carbons', 0)}"
        ),
        "",
    ]

    holes = ctx.get("hole_summary")
    if holes:
        lines.extend(
            [
                "Pores / Holes",
                f"  Distribution   : {holes.get('distribution')}   |   "
                f"Excavation events: {holes.get('excavation_events', holes.get('holes_created', 0))}   |   "
                f"Pruned tips: {holes.get('pruned_tip_carbons', 0)} C",
                f"  Components     : {holes.get('final_pore_count', holes.get('final_pore_components', 0))} physical openings "
                f"({holes.get('internal_pores', holes.get('internal_pore_count', 0))} internal pores, "
                f"{holes.get('outer_edge_notches', holes.get('outer_edge_notch_count', 0))} outer-edge notches; "
                f"{holes.get('deleted_atom_clusters', holes.get('final_pore_components', 0))} deleted clusters)",
                f"  Vacancy        : {holes.get('removed_carbons', 0)}/{holes.get('initial_carbons', 0)} C "
                f"({float(holes.get('actual_vacancy_percent', 0.0)):.2f}%)",
                f"  Accessible Area (2D bare-C scaffold estimate): {float(holes.get('pore_accessible_area_angstrom2', 0.0)):.2f} A^2 "
                f"({float(holes.get('pore_accessible_area_nm2', 0.0)):.4f} nm^2, "
                f"probe r={float(holes.get('probe_radius_angstrom', 1.4)):.2f} A)",
                "",
            ]
        )

    lines.extend(_format_chemistry_and_qm_lines(chem, geom, ctx.get("group_counts", {})))
    lines.append("")

    placement = ctx.get("placement", "random")
    if isinstance(placement, dict):
        pmode = placement.get("mode", "random")
        extra = []
        if pmode in ("island", "irregular_island"):
            if placement.get("style"):
                extra.append(f"style={placement['style']}")
            if placement.get("branching") is not None:
                extra.append(f"branching={placement['branching']}")
            extra.append(f"islands={placement.get('n_islands', 1)}")
            if placement.get("compactness") is not None:
                extra.append(f"compactness={placement['compactness']}")
            if placement.get("edge_mode"):
                extra.append(f"edge_mode={placement['edge_mode']}")
        nn_d = placement.get("mean_basal_group_nn_distance")
        if nn_d is not None:
            extra.append(f"mean_basal_nn={nn_d:.3f} A")
        place_str = f"{pmode} ({', '.join(extra)})" if extra else pmode
    else:
        place_str = str(placement)
    lines.extend(
        [
            "Placement",
            f"  Mode           : {place_str}   |   Faces: {ctx.get('faces', 'both')}",
            "",
        ]
    )

    label_map = {
        "output": "JSON model",
        "pdb": "PDB",
        "xyz": "XYZ",
        "extxyz": "Extended XYZ",
        "mol2": "MOL2",
        "sdf": "SDF (V3000)",
        "cif": "CIF (P1)",
        "write_config": "Config",
        "log": "Build log",
    }
    lines.append("Written Files")
    for key, path in ctx.get("written_files", {}).items():
        lines.append(f"  {label_map.get(key, key):<14s} : {path}")
    return "\n".join(lines)


def format_check_report(result: dict) -> str:
    chem = result.get("chemistry", {})
    geom = result.get("geometry", {})
    graph = result.get("graph", {})
    gbsa = result.get("gbsa", {})
    dims = chem.get("dimensions_nm", {})
    ext = dims.get("extent_nm", [0.0, 0.0, 0.0])
    cell = dims.get("cell_nm")
    pbc = dims.get("periodic", [False, False, False])

    lines = [
        "=== GO Builder — Check Summary ===",
        f"Name: {result.get('model_name', 'GO')}   |   "
        f"Graph: {'VALID' if graph.get('valid') else 'INVALID'}   |   "
        f"Geometry: {'PASS' if geom.get('passed') else 'FAIL'}   |   "
        f"GBSA Ready: {gbsa.get('ready')}",
        "",
        "Geometry & Dimensions",
        f"  Extent (nm)    : {ext[0]:.3f} x {ext[1]:.3f} x {ext[2]:.3f}   |   PBC: {pbc}",
        (
            f"  Cell (nm)      : {cell[0]:.3f} x {cell[1]:.3f} x {cell[2]:.3f}"
            if cell is not None
            else "  Cell (nm)      : none"
        ),
        f"  Topology       : {geom.get('atom_count', graph.get('atom_count', 0))} atoms, "
        f"{geom.get('bond_count', graph.get('bond_count', 0))} bonds",
        "",
    ]
    lines.extend(_format_chemistry_and_qm_lines(chem, geom))
    return "\n".join(lines)


def _format_build_log(report_text: str, recipe: dict) -> str:
    iso_ts = datetime.now().astimezone().isoformat(timespec="seconds")
    return (
        f"{report_text}\n\n"
        f"Timestamp (ISO): {iso_ts}\n\n"
        "Effective Configuration (JSON):\n"
        f"{json.dumps(recipe, indent=2, ensure_ascii=False, sort_keys=False)}\n"
    )


def _report(stage: str, result: dict, machine: bool, error: bool = False) -> None:
    display_ctx = result.pop("_display", None)
    if machine:
        print(json.dumps(result, indent=2, default=str))
        return
    if not error and stage == "build" and display_ctx is not None:
        print(format_build_report(display_ctx))
        return
    if not error and stage == "check" and "chemistry" in result:
        print(format_check_report(result))
        return
    if not error and stage == "build" and result.get("status") == "AVAILABLE_PRESETS":
        print("Available presets and illustrative templates:")
        for p in result["presets"]:
            tag = p.get("family", "Preset")
            enabled = "buildable" if p.get("build_enabled") else "disabled"
            print(f"  - {p['name']:<20s} [{tag}, {enabled}, {p.get('readiness')}]")
        return
    if not error and stage == "amber" and result.get("status") == "AVAILABLE_FORCE_FIELD_PRESETS":
        print("Available AMBER force-field presets:")
        for p in result["presets"]:
            default_tag = " [DEFAULT]" if p.get("default") else ""
            print(f"  - {p['name']:<16s} ({p.get('tier', 'preset')}, v{p.get('version', '?')}){default_tag}")
            print(f"      Bonded : {p.get('bonded_model')}")
            print(f"      Charges: {p.get('charge_model')} | {p.get('status')}")
        return
    if not error and stage == "build" and result.get("status") in ("CONFIG_INITIALIZED", "CONFIG_WRITTEN"):
        print(f"[{result['status']}] Wrote configuration to {result['output']}")
        return

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
            "xml",
            "provenance",
            "preset",
            "charge_model",
            "custom_combination",
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
    if "presets" in result:
        summary["presets"] = result["presets"]
    if "gbsa" in result:
        summary["gbsa_ready"] = result["gbsa"]["ready"]
    if "geometry" in result:
        summary["geometry_passed"] = result["geometry"]["passed"]
        summary["clash_count"] = result["geometry"]["clash_count"]
        if "electron_count" in result["geometry"]:
            summary["electron_count"] = result["geometry"]["electron_count"]
        if "electron_parity" in result["geometry"]:
            summary["electron_parity"] = result["geometry"]["electron_parity"]
        if result["geometry"].get("warnings"):
            summary["warnings"] = result["geometry"]["warnings"]
    print(
        f"[{timestamp}] [GO-{stage.upper()}] [{level}] "
        + json.dumps(summary, default=str, ensure_ascii=False),
        file=stream,
    )


def _update_preset_provenance(provenance: dict, recipe: dict, model: dict) -> dict:
    from .model import compute_chemistry_summary

    meta = model.get("metadata", {})
    dims = meta.get("dimensions") or {}
    chem = compute_chemistry_summary(model)
    cc = dict(provenance.get("constructed_composition", {}))
    for k in ("nx", "ny", "lx_nm", "ly_nm", "diameter_nm", "radius_nm"):
        cc.pop(k, None)
    for k in ("nx", "ny", "lx_nm", "ly_nm", "diameter_nm", "radius_nm", "shape", "edge_type", "faces"):
        if k in recipe:
            cc[k] = recipe[k]
        elif k in ("nx", "ny") and k in dims:
            cc[k] = dims[k]
    if "holes" in recipe:
        cc["holes"] = recipe["holes"]
    if "placement" in recipe:
        cc["placement"] = recipe["placement"]
    if "composition" in recipe:
        cc["composition"] = recipe["composition"]
    if "groups" in recipe:
        cc["groups"] = recipe["groups"]
    if "resolved_group_counts" in meta:
        cc["resolved_group_counts"] = meta["resolved_group_counts"]
    elif isinstance(meta.get("composition"), dict) and "group_counts" in meta["composition"]:
        cc["resolved_group_counts"] = meta["composition"]["group_counts"]
    provenance["constructed_composition"] = cc

    realized_groups = dict(
        sorted(Counter(g["type"] for g in meta.get("groups", [])).items())
    )
    provenance["realized_composition"] = {
        "formula": chem.get("formula"),
        "atom_counts": chem.get("atom_counts"),
        "formal_charge": chem.get("formal_charge"),
        "c_over_o": chem.get("c_over_o"),
        "o_over_c": chem.get("o_over_c"),
        "group_counts": realized_groups,
        "extent_nm": chem.get("dimensions_nm", {}).get("extent_nm"),
        "scaffold_diameter_nm": chem.get("dimensions_nm", {}).get("scaffold_diameter_nm"),
        "seed": meta.get("seed"),
    }
    return provenance


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
        from .model import (
            assert_export_ready,
            save_model,
            write_cif,
            write_extxyz,
            write_mol2,
            write_pdb,
            write_sdf,
            write_xyz,
        )

        xyz_path = getattr(args, "xyz", None)
        extxyz_path = getattr(args, "extxyz", None)
        mol2_path = getattr(args, "mol2", None)
        sdf_path = getattr(args, "sdf", None)
        cif_path = getattr(args, "cif", None)
        log_path = getattr(args, "log", None)
        write_cfg_path = getattr(args, "write_config", None)

        # Check all destinations for distinctness, extension validity, and non-existence before writing any file.
        _validate_write_config_path(write_cfg_path)
        _validate_output_paths(
            [
                args.output,
                args.pdb,
                xyz_path,
                extxyz_path,
                mol2_path,
                sdf_path,
                cif_path,
                log_path,
                write_cfg_path,
            ]
        )

        recipe: dict = {}
        if stage == "build":
            from .build import build_sheet

            recipe, provenance = _resolve_build_recipe(args)
            if provenance is not None:
                readiness = provenance["source_info"]
                if not readiness.get("build_enabled", False):
                    raise ValueError(
                        f"Preset {args.preset} is disabled: {readiness.get('readiness', 'validation pending')}; "
                        "placement validation and force-field parameterization selection are pending. "
                        "No model was written."
                    )
            model = build_sheet(recipe)
            if provenance is not None:
                model.setdefault("metadata", {})["preset"] = _update_preset_provenance(
                    provenance, recipe, model
                )
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

        has_image_bonds = any(any(b.get("image", [0, 0, 0])) for b in model["bonds"])
        if any((args.pdb, xyz_path, extxyz_path, mol2_path, sdf_path, cif_path)):
            assert_export_ready(model, allow_undercoordinated=(stage == "import"))
        if args.pdb and has_image_bonds:
            raise ValueError("PDB cannot preserve periodic image bonds; use JSON")
        if xyz_path and has_image_bonds:
            raise ValueError("XYZ cannot preserve periodic image bonds; use write_extxyz or JSON")
        if mol2_path and has_image_bonds:
            raise ValueError("MOL2 cannot preserve periodic image bonds; use JSON")
        if sdf_path and has_image_bonds:
            raise ValueError("SDF cannot preserve periodic image bonds; use JSON")

        save_model(model, args.output)
        if args.pdb:
            write_pdb(model, args.pdb, allow_undercoordinated=(stage == "import"))
        if xyz_path:
            write_xyz(model, xyz_path)
        if extxyz_path:
            write_extxyz(model, extxyz_path)
        if mol2_path:
            write_mol2(model, mol2_path)
        if sdf_path:
            write_sdf(model, sdf_path)
        if cif_path:
            write_cif(model, cif_path)
        if write_cfg_path is not None:
            from .config import dump_config

            dump_config(recipe, write_cfg_path)

        written_files: dict[str, str] = {"output": str(args.output)}
        if args.pdb:
            written_files["pdb"] = str(args.pdb)
        if xyz_path:
            written_files["xyz"] = str(xyz_path)
        if extxyz_path:
            written_files["extxyz"] = str(extxyz_path)
        if mol2_path:
            written_files["mol2"] = str(mol2_path)
        if sdf_path:
            written_files["sdf"] = str(sdf_path)
        if cif_path:
            written_files["cif"] = str(cif_path)
        if write_cfg_path is not None:
            written_files["write_config"] = str(write_cfg_path)
        if log_path is not None:
            written_files["log"] = str(log_path)

        display_ctx = None
        if stage == "build":
            display_ctx = _build_display_context(model, recipe, written_files)
            if log_path is not None:
                log_path.parent.mkdir(parents=True, exist_ok=True)
                report_text = format_build_report(display_ctx)
                with log_path.open("x", encoding="utf-8") as handle:
                    handle.write(_format_build_log(report_text, recipe))

        metadata = model.get("metadata", {})
        result = {
            "status": "BUILT" if stage == "build" else "PREPARED",
            "output": str(args.output),
            **({"pdb": str(args.pdb)} if args.pdb else {}),
            **({"xyz": str(xyz_path)} if xyz_path else {}),
            **({"extxyz": str(extxyz_path)} if extxyz_path else {}),
            **({"mol2": str(mol2_path)} if mol2_path else {}),
            **({"sdf": str(sdf_path)} if sdf_path else {}),
            **({"cif": str(cif_path)} if cif_path else {}),
            **({"write_config": str(write_cfg_path)} if write_cfg_path else {}),
            **({"log": str(log_path)} if log_path else {}),
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
            **(
                {"hole_summary": metadata["hole_summary"]}
                if metadata.get("hole_summary") is not None
                else {}
            ),
            **(
                {"composition_report": metadata["composition"]}
                if metadata.get("composition") is not None
                else {}
            ),
        }
        if display_ctx is not None:
            result["_display"] = display_ctx
        return result
    if stage == "amber":
        from .ff_v3 import list_force_field_presets
        from .topology import prepare_amber

        if args.list_presets:
            return {
                "status": "AVAILABLE_FORCE_FIELD_PRESETS",
                "presets": list_force_field_presets(),
            }
        if args.direct_periodic:
            from .periodic_topology import export_periodic_amber

            return export_periodic_amber(
                _read(args.model),
                args.output_dir,
                charge_model=args.charge_model,
                charges_json=args.charges_json,
                reference_dir=args.reference_dir,
                radii=args.radii,
            )
        templates = _read(args.parameter_templates) if args.parameter_templates else None
        return prepare_amber(
            _read(args.model),
            args.output_dir,
            run=args.run,
            tleap=args.tleap,
            preset=args.preset,
            charge_model=args.charge_model,
            charges_json=args.charges_json,
            reference_dir=args.reference_dir,
            parameter_templates=templates,
            radii=args.radii,
        )
    if stage == "gromacs":
        if args.model is not None:
            from .periodic_topology import export_periodic_gromacs

            return export_periodic_gromacs(
                _read(args.model),
                args.output_prefix,
                charge_model=args.charge_model,
                charges_json=args.charges_json,
                reference_dir=args.reference_dir,
                cutoff_nm=args.cutoff_nm,
            )
        from .topology import export_gromacs

        return export_gromacs(args.topology, args.coordinates, args.output_prefix)
    if stage == "openmm":
        if args.model is not None:
            from .periodic_topology import export_periodic_openmm

            return export_periodic_openmm(
                _read(args.model),
                args.output,
                nonbonded_method=args.nonbonded_method or "PME",
                cutoff_nm=args.cutoff_nm,
                charge_model=args.charge_model,
                charges_json=args.charges_json,
                reference_dir=args.reference_dir,
            )
        from .topology import export_openmm

        return export_openmm(
            args.topology,
            args.coordinates,
            args.output,
            nonbonded_method=args.nonbonded_method or "NoCutoff",
            cutoff_nm=1.0 if args.cutoff_nm is None else args.cutoff_nm,
        )
    if stage == "lammps":
        if args.model is not None:
            from .lammps_export import export_periodic_lammps

            return export_periodic_lammps(
                _read(args.model),
                args.output_prefix,
                cutoff_ang=args.cutoff_ang,
                charge_model=args.charge_model,
                charges_json=args.charges_json,
                reference_dir=args.reference_dir,
            )
        from .lammps_export import export_lammps

        return export_lammps(
            args.topology,
            args.coordinates,
            args.output_prefix,
            cutoff_ang=120.0 if args.cutoff_ang is None else args.cutoff_ang,
        )
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
    from .model import check_geometry_integrity, compute_chemistry_summary, validate_model
    from .topology import gbsa_readiness

    model = _read(args.model)
    return {
        "status": "CHECKED",
        "model_name": model.get("name", args.model.stem),
        "graph": validate_model(model),
        "geometry": check_geometry_integrity(model),
        "chemistry": compute_chemistry_summary(model),
        "gbsa": gbsa_readiness(model),
    }


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
    parser = parser_for(stage)
    args = parser.parse_args(argv)

    write_config_only = False
    if stage == "build":
        has_export_outputs = any(
            getattr(args, attr, None) is not None
            for attr in ("output", "pdb", "xyz", "extxyz", "mol2", "sdf", "cif", "log")
        )
        if args.list_presets:
            if args.write_config is not None:
                parser.error("--list-presets cannot be combined with --write-config")
        elif args.init_config is not None:
            if args.write_config is not None or has_export_outputs:
                parser.error("--init-config cannot be combined with build output or --write-config")
        elif args.write_config is not None and not has_export_outputs:
            write_config_only = True
        elif args.output is None:
            parser.error(
                "--out/--output is required for building (or use --list-presets / --init-config / --write-config)"
            )
    elif stage == "amber":
        if not args.list_presets and (args.model is None or args.output_dir is None):
            parser.error("--model and --output-dir are required (or use --list-presets)")
        if args.direct_periodic and args.parameter_templates is not None:
            parser.error("--direct-periodic cannot be combined with --parameter-templates")
        if args.direct_periodic and args.run:
            parser.error("--direct-periodic does not use external tleap (--run)")
    elif stage in ("gromacs", "openmm", "lammps"):
        has_amber_pair = args.topology is not None or args.coordinates is not None
        if args.model is not None and has_amber_pair:
            parser.error("--model cannot be combined with --topology or --coordinates")
        if args.model is None and (args.topology is None or args.coordinates is None):
            parser.error("supply either --model or both --topology and --coordinates")
        if args.model is None and (
            args.charge_model != "auto"
            or args.charges_json is not None
            or args.reference_dir is not None
            or (stage == "gromacs" and args.cutoff_nm is not None)
        ):
            parser.error(
                "--charge-model, --charges-json, --reference-dir, and periodic options require --model"
            )

    worker_py = getattr(args, "worker_python", None) or sys.executable
    packaged_resources = Path(__file__).resolve().parent / "_resources"
    go_script = ROOT / "go.py"
    if not packaged_resources.exists() and go_script.is_file():
        command = [worker_py, str(go_script), stage, *_worker_argv(argv)]
    else:
        command = [worker_py, "-m", "go_builder", stage, *_worker_argv(argv)]
    try:
        if stage == "build":
            if args.list_presets:
                from .presets import list_presets

                result = {"status": "AVAILABLE_PRESETS", "presets": list_presets()}
            elif args.init_config is not None:
                from .config import (
                    render_annotated_json_template,
                    render_annotated_toml_template,
                )

                dest = args.init_config
                suffix = dest.suffix.lower()
                if suffix == ".toml":
                    content = render_annotated_toml_template()
                elif suffix == ".json":
                    content = render_annotated_json_template()
                else:
                    raise ValueError(
                        f"Unsupported config template extension {suffix!r}; use .toml or .json"
                    )
                _validate_output_paths([dest])
                dest.parent.mkdir(parents=True, exist_ok=True)
                with dest.open("x", encoding="utf-8") as handle:
                    handle.write(content)
                result = {
                    "status": "CONFIG_INITIALIZED",
                    "output": str(dest),
                    "format": suffix.lstrip("."),
                }
            elif args.dry_run:
                _check_build_cli_conflicts(args)
                _validate_write_config_path(getattr(args, "write_config", None))
                result = {"status": "PLANNED", "command": command, "backend": "local"}
            elif write_config_only:
                from .config import dump_config

                _validate_write_config_path(args.write_config)
                _validate_output_paths([args.write_config])
                recipe, _ = _resolve_build_recipe(args)
                dump_config(recipe, args.write_config)
                result = {
                    "status": "CONFIG_WRITTEN",
                    "output": str(args.write_config),
                    **({"preset": args.preset} if args.preset else {}),
                    "recipe": recipe,
                }
            else:
                result = _dispatch(stage, args)
        elif args.slurm:
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
        if stage == "check" and getattr(args, "strict", False):
            if not result.get("geometry", {}).get("passed", True):
                return 1
        return 0
    except PlacementError as error:
        _report(stage, error.as_report(), args.json, error=True)
        return 1
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError, KeyError) as error:
        _report(stage, {"status": "ERROR", "message": str(error)}, args.json, error=True)
        return 1
