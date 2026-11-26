"""Native LAMMPS (.data atom_style full + .in) exporter for GO Workbench (Block B4).

Exports both covalently periodic and finite GO models (directly from ``gaff2_v3``
or from existing AMBER ``prmtop``/``inpcrd`` files) into LAMMPS ``units real`` /
``atom_style full`` data files with ready-to-run input templates and SHA-256
provenance manifests.

Conventions matched against LAMMPS ``src/MOLECULE/improper_cvff.cpp``,
``dihedral_style fourier``, and ``special_bonds amber``:
- Harmonic bonds/angles use AMBER/real units directly (``K`` in kcal/(mol*A^2)
  and kcal/(mol*rad^2); ``r0`` in A; ``theta0`` in degrees).
- Proper dihedrals sharing an ordered quartet ``(i, j, k, l)`` are grouped into
  one ``dihedral_style fourier`` entry ``m K_1 n_1 d_1 ...``.
- Improper torsions keep native AMBER/ParmEd quartet order ``(i, j, k, l)``
  which matches ``improper_style cvff`` (``K d n`` with ``d=+1`` for 0 deg,
  ``d=-1`` for 180 deg) verbatim.
"""

from __future__ import annotations

import contextlib
import hashlib
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from .periodic_topology import _resolve_periodic_cutoff_nm, build_v3_parmed_structure
from .topology import (
    ENGINE_SOURCES,
    _json,
    _link_upstream_provenance,
    _load_amber_with_digests,
    _periodic_edges,
    _reduce_openmm_box_vectors,
    _verify_supported_amber_terms,
)

LAMMPS_SOURCES = [
    *ENGINE_SOURCES,
    "https://docs.lammps.org/atom_style.html",
    "https://docs.lammps.org/dihedral_fourier.html",
    "https://docs.lammps.org/improper_cvff.html",
    "https://docs.lammps.org/special_bonds.html",
]


def _fmt_f(val: float) -> str:
    v = float(val)
    if v == 0.0:
        v = 0.0
    return f"{v:.10g}"


def _rmin2_to_sigma_ang(rmin_half: float) -> float:
    val = float(rmin_half)
    if val == 0.0:
        return 0.0
    return float(2.0 * val * (2.0 ** (-1.0 / 6.0)))


def _align_box_and_coords(
    structure: Any, *, is_periodic: bool
) -> tuple[dict[str, float], np.ndarray, bool]:
    """Return LAMMPS box bounds/tilts (Angstrom), aligned coordinates, and triclinic flag."""
    coords = np.asarray(structure.coordinates, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or not np.isfinite(coords).all():
        raise ValueError("Invalid nonfinite atom coordinates in topology")

    if is_periodic:
        if getattr(structure, "box_vectors", None) is None:
            raise ValueError("Periodic LAMMPS export requires periodic box vectors")
        import parmed as pmd

        vecs_ang = np.asarray(
            structure.box_vectors.value_in_unit(pmd.unit.angstroms), dtype=float
        )
        # Convert reduced nm vectors cleanly back to A while preserving exact lower-triangular layout
        reduced_nm = np.asarray(_reduce_openmm_box_vectors(vecs_ang * 0.1), dtype=float)
        cell_lmp = reduced_nm * 10.0

        # Rigid orthogonal alignment if input box is not already aligned to lower-triangular frame
        gram = vecs_ang @ vecs_ang.T
        low_orig = np.linalg.cholesky(gram)
        # If vecs_ang is already diagonal/lower-triangular, preserve coordinates directly
        if (
            np.allclose(vecs_ang, np.tril(vecs_ang), atol=1e-8)
            and vecs_ang[0, 0] > 0.0
            and vecs_ang[1, 1] > 0.0
            and vecs_ang[2, 2] > 0.0
        ):
            coords_out = coords.copy()
        else:
            rot = np.linalg.solve(vecs_ang, low_orig)
            coords_out = coords @ rot

        lx = float(cell_lmp[0, 0])
        ly = float(cell_lmp[1, 1])
        lz = float(cell_lmp[2, 2])
        xy = float(cell_lmp[1, 0])
        xz = float(cell_lmp[2, 0])
        yz = float(cell_lmp[2, 1])
        triclinic = any(abs(v) > 1e-9 for v in (xy, xz, yz))
        box = {
            "xlo": 0.0,
            "xhi": lx,
            "ylo": 0.0,
            "yhi": ly,
            "zlo": 0.0,
            "zhi": lz,
            "xy": xy,
            "xz": xz,
            "yz": yz,
        }
        return box, coords_out, triclinic

    # Non-periodic finite system: use tight padded orthorhombic box around coordinates
    mins = coords.min(axis=0) - 20.0
    maxs = coords.max(axis=0) + 20.0
    box = {
        "xlo": float(mins[0]),
        "xhi": float(maxs[0]),
        "ylo": float(mins[1]),
        "yhi": float(maxs[1]),
        "zlo": float(mins[2]),
        "zhi": float(maxs[2]),
        "xy": 0.0,
        "xz": 0.0,
        "yz": 0.0,
    }
    return box, coords.copy(), False


def _infer_element_label(atom: Any) -> str:
    el = str(getattr(atom, "element_name", "") or "").strip()
    if el and el != "EP":
        return el
    mass = float(getattr(atom, "mass", 0.0) or 0.0)
    if mass < 3.0:
        return "H"
    if abs(mass - 12.01) < 1.5:
        return "C"
    if abs(mass - 16.00) < 1.5:
        return "O"
    return "?"


def _serialize_amber_parm_to_lammps(
    structure: Any,
    *,
    is_periodic: bool,
    cutoff_ang: float,
) -> tuple[str, str, dict[str, Any]]:
    """Serialize a ParmEd AmberParm into LAMMPS ``data`` and ``in`` file contents."""
    _verify_supported_amber_terms(
        structure, target="LAMMPS", allow_explicit_exclusions=False
    )
    if not structure.atoms:
        raise ValueError("Topology contains no atoms")
    if not math.isfinite(float(cutoff_ang)) or float(cutoff_ang) <= 0.0:
        raise ValueError("cutoff_ang must be positive and finite")

    box, coords, triclinic = _align_box_and_coords(structure, is_periodic=is_periodic)
    if is_periodic:
        max_h = 0.5 * min(
            box["xhi"] - box["xlo"],
            box["yhi"] - box["ylo"],
            box["zhi"] - box["zlo"],
        )
        if float(cutoff_ang) >= max_h:
            raise ValueError(
                f"Periodic LAMMPS cutoff_ang={cutoff_ang} must be strictly less than half the box length ({max_h:.4f} A)"
            )

    # Atom types deduplicated by (type_name, mass, eps, sigma)
    atom_type_id_map: dict[tuple[str, float, float, float], int] = {}
    atom_type_records: list[tuple[int, str, float, float, float, str]] = []
    atom_type_indices: list[int] = []
    for atom in structure.atoms:
        atype = str(getattr(atom, "type", "") or atom.name)
        elem_lbl = _infer_element_label(atom)
        mass = float(atom.mass)
        eps = float(atom.epsilon)
        sig = _rmin2_to_sigma_ang(float(atom.rmin))
        if eps == 0.0 and sig == 0.0:
            sig = 0.5
        sig_key = (atype, round(mass, 8), round(eps, 10), round(sig, 10))
        if sig_key not in atom_type_id_map:
            tid = len(atom_type_records) + 1
            atom_type_id_map[sig_key] = tid
            atom_type_records.append((tid, atype, mass, eps, sig, elem_lbl))
        atom_type_indices.append(atom_type_id_map[sig_key])

    # Bonds (harmonic: K r0)
    bond_type_map: dict[tuple[float, float], int] = {}
    bond_type_records: list[tuple[int, float, float, str]] = []
    bond_records: list[tuple[int, int, int, int]] = []
    for idx, b in enumerate(structure.bonds, 1):
        if b.type is None:
            raise ValueError("Untyped bond encountered during LAMMPS export")
        k = float(b.type.k)
        req = float(b.type.req)
        bkey = (round(k, 10), round(req, 10))
        if bkey not in bond_type_map:
            btid = len(bond_type_records) + 1
            bond_type_map[bkey] = btid
            lbl = f"{b.atom1.type}-{b.atom2.type}"
            bond_type_records.append((btid, k, req, lbl))
        bond_records.append((idx, bond_type_map[bkey], b.atom1.idx + 1, b.atom2.idx + 1))

    # Angles (harmonic: K theta0)
    angle_type_map: dict[tuple[float, float], int] = {}
    angle_type_records: list[tuple[int, float, float, str]] = []
    angle_records: list[tuple[int, int, int, int, int]] = []
    for idx, ang in enumerate(structure.angles, 1):
        if ang.type is None:
            raise ValueError("Untyped angle encountered during LAMMPS export")
        k = float(ang.type.k)
        theteq = float(ang.type.theteq)
        akey = (round(k, 10), round(theteq, 10))
        if akey not in angle_type_map:
            atid = len(angle_type_records) + 1
            angle_type_map[akey] = atid
            lbl = f"{ang.atom1.type}-{ang.atom2.type}-{ang.atom3.type}"
            angle_type_records.append((atid, k, theteq, lbl))
        angle_records.append(
            (idx, angle_type_map[akey], ang.atom1.idx + 1, ang.atom2.idx + 1, ang.atom3.idx + 1)
        )

    # Proper & improper dihedrals separated cleanly
    proper_groups: dict[tuple[int, int, int, int], list[tuple[float, int, float]]] = defaultdict(list)
    proper_labels: dict[tuple[int, int, int, int], str] = {}
    improper_raw: list[tuple[tuple[int, int, int, int], float, int, int, str]] = []

    for dih in structure.dihedrals:
        types = dih.type if isinstance(dih.type, list) else [dih.type]
        quad = (dih.atom1.idx + 1, dih.atom2.idx + 1, dih.atom3.idx + 1, dih.atom4.idx + 1)
        lbl = f"{dih.atom1.type}-{dih.atom2.type}-{dih.atom3.type}-{dih.atom4.type}"
        for dt in types:
            if dt is None:
                raise ValueError("Untyped dihedral encountered during LAMMPS export")
            per_f = float(dt.per)
            if not math.isfinite(per_f) or per_f < 0.0 or not per_f.is_integer():
                raise ValueError(f"Invalid torsion periodicity {dt.per!r} during LAMMPS export")
            per = int(per_f)
            phi_k = float(dt.phi_k)
            phase = float(dt.phase)
            if not math.isfinite(phi_k) or not math.isfinite(phase):
                raise ValueError("Nonfinite torsion parameters during LAMMPS export")
            if per == 0 and phi_k == 0.0:
                continue
            if per == 0:
                raise ValueError(
                    f"Zero dihedral periodicity with nonzero force constant ({phi_k}) unsupported"
                )
            if dih.improper:
                phase_mod = round(phase % 360.0, 4)
                if abs(phase_mod - 180.0) < 1e-3:
                    d_sign = -1
                elif phase_mod < 1e-3 or abs(phase_mod - 360.0) < 1e-3:
                    d_sign = 1
                else:
                    raise ValueError(
                        f"improper_style cvff supports phase 0 or 180 deg only (got {phase})"
                    )
                improper_raw.append((quad, phi_k, d_sign, per, lbl))
            else:
                proper_groups[quad].append((phi_k, per, phase))
                proper_labels.setdefault(quad, lbl)

    # Multi-term proper dihedrals grouped per quartet for LAMMPS dihedral_style fourier
    dihe_type_map: dict[tuple[tuple[float, int, float], ...], int] = {}
    dihe_type_records: list[tuple[int, tuple[tuple[float, int, float], ...], str]] = []
    dihe_records: list[tuple[int, int, int, int, int, int]] = []
    for idx, (quad, terms) in enumerate(proper_groups.items(), 1):
        sig = tuple((round(k, 10), per, round(ph, 8)) for k, per, ph in terms)
        if sig not in dihe_type_map:
            dtid = len(dihe_type_records) + 1
            dihe_type_map[sig] = dtid
            dihe_type_records.append((dtid, sig, proper_labels[quad]))
        dihe_records.append((idx, dihe_type_map[sig], *quad))

    # Impropers (cvff: K d n, exact Amber quartet order)
    imp_type_map: dict[tuple[float, int, int], int] = {}
    imp_type_records: list[tuple[int, float, int, int, str]] = []
    imp_records: list[tuple[int, int, int, int, int, int]] = []
    for idx, (quad, k, d_sign, per, lbl) in enumerate(improper_raw, 1):
        ikey = (round(k, 10), d_sign, per)
        if ikey not in imp_type_map:
            itid = len(imp_type_records) + 1
            imp_type_map[ikey] = itid
            imp_type_records.append((itid, k, d_sign, per, lbl))
        imp_records.append((idx, imp_type_map[ikey], *quad))

    # Build LAMMPS data text
    lines: list[str] = [
        "LAMMPS data file generated by GO Workbench (units real, atom_style full)",
        "",
        f"{len(structure.atoms)} atoms",
        f"{len(bond_records)} bonds",
        f"{len(angle_records)} angles",
        f"{len(dihe_records)} dihedrals",
        f"{len(imp_records)} impropers",
        "",
        f"{len(atom_type_records)} atom types",
        f"{len(bond_type_records)} bond types",
        f"{len(angle_type_records)} angle types",
        f"{len(dihe_type_records)} dihedral types",
        f"{len(imp_type_records)} improper types",
        "",
        f"{box['xlo']:.10f} {box['xhi']:.10f} xlo xhi",
        f"{box['ylo']:.10f} {box['yhi']:.10f} ylo yhi",
        f"{box['zlo']:.10f} {box['zhi']:.10f} zlo zhi",
    ]
    if triclinic:
        lines.append(f"{box['xy']:.10f} {box['xz']:.10f} {box['yz']:.10f} xy xz yz")
    lines.append("")

    lines.append("Masses")
    lines.append("")
    for tid, atype, mass, _, _, elem_lbl in atom_type_records:
        lines.append(f"{tid} {_fmt_f(mass)} # {atype} ({elem_lbl})")
    lines.append("")

    lines.append("Pair Coeffs # lj/cut/coul/*")
    lines.append("")
    for tid, atype, _, eps, sig, elem_lbl in atom_type_records:
        lines.append(f"{tid} {_fmt_f(eps)} {_fmt_f(sig)} # {atype} ({elem_lbl})")
    lines.append("")

    if bond_type_records:
        lines.append("Bond Coeffs # harmonic")
        lines.append("")
        for btid, k, req, lbl in bond_type_records:
            lines.append(f"{btid} {_fmt_f(k)} {_fmt_f(req)} # {lbl}")
        lines.append("")

    if angle_type_records:
        lines.append("Angle Coeffs # harmonic")
        lines.append("")
        for atid, k, theteq, lbl in angle_type_records:
            lines.append(f"{atid} {_fmt_f(k)} {_fmt_f(theteq)} # {lbl}")
        lines.append("")

    if dihe_type_records:
        lines.append("Dihedral Coeffs # fourier")
        lines.append("")
        for dtid, sig, lbl in dihe_type_records:
            terms_str = " ".join(
                f"{_fmt_f(k)} {per} {_fmt_f(ph)}" for k, per, ph in sig
            )
            lines.append(f"{dtid} {len(sig)} {terms_str} # {lbl}")
        lines.append("")

    if imp_type_records:
        lines.append("Improper Coeffs # cvff")
        lines.append("")
        for itid, k, d_sign, per, lbl in imp_type_records:
            lines.append(f"{itid} {_fmt_f(k)} {d_sign} {per} # {lbl}")
        lines.append("")

    lines.append("Atoms # full")
    lines.append("")
    for idx, atom in enumerate(structure.atoms):
        aid = idx + 1
        mol_id = int(atom.residue.idx) + 1
        tid = atom_type_indices[idx]
        q = float(atom.charge)
        x, y, z = coords[idx]
        lines.append(
            f"{aid} {mol_id} {tid} {q:.10f} {x:.10f} {y:.10f} {z:.10f} # {atom.residue.name}:{atom.name}"
        )
    lines.append("")

    if bond_records:
        lines.append("Bonds")
        lines.append("")
        for bid, btid, a1, a2 in bond_records:
            lines.append(f"{bid} {btid} {a1} {a2}")
        lines.append("")

    if angle_records:
        lines.append("Angles")
        lines.append("")
        for aid_n, atid, a1, a2, a3 in angle_records:
            lines.append(f"{aid_n} {atid} {a1} {a2} {a3}")
        lines.append("")

    if dihe_records:
        lines.append("Dihedrals")
        lines.append("")
        for did, dtid, a1, a2, a3, a4 in dihe_records:
            lines.append(f"{did} {dtid} {a1} {a2} {a3} {a4}")
        lines.append("")

    if imp_records:
        lines.append("Impropers")
        lines.append("")
        for iid, itid, a1, a2, a3, a4 in imp_records:
            lines.append(f"{iid} {itid} {a1} {a2} {a3} {a4}")
        lines.append("")

    data_text = "\n".join(lines)

    # Matching LAMMPS input template
    if is_periodic:
        boundary_line = "boundary        p p p"
        pair_line = f"pair_style      lj/cut/coul/long {cutoff_ang:.4f}"
        kspace_line = "kspace_style    pppm 1.0e-5"
        tail_line = "pair_modify     mix arithmetic tail yes"
    else:
        boundary_line = "boundary        f f f"
        pair_line = f"pair_style      lj/cut/coul/cut {cutoff_ang:.4f}"
        kspace_line = "# non-periodic finite solute"
        tail_line = "pair_modify     mix arithmetic tail no"

    in_lines = [
        "# Ready-to-use LAMMPS input template (AMBER/GAFF2 conventions)",
        "units           real",
        "atom_style      full",
        boundary_line,
        "",
        pair_line,
        tail_line,
        "special_bonds   amber",
        "bond_style      harmonic",
        "angle_style     harmonic",
        "dihedral_style  fourier",
        "improper_style  cvff",
        "",
        "read_data       __DATA_FILE__",
        kspace_line,
        "",
        "thermo_style    custom step pe ebond eangle edihed eimp evdwl ecoul elong etail",
        "run             0",
        "",
    ]
    in_text = "\n".join(in_lines)

    counts = {
        "atom_count": len(structure.atoms),
        "bond_count": len(bond_records),
        "angle_count": len(angle_records),
        "proper_dihedral_count": len(dihe_records),
        "proper_dihedral_terms_count": sum(len(v) for v in proper_groups.values()),
        "improper_dihedral_count": len(imp_records),
        "atom_types": len(atom_type_records),
        "bond_types": len(bond_type_records),
        "angle_types": len(angle_type_records),
        "dihedral_types": len(dihe_type_records),
        "improper_types": len(imp_type_records),
        "triclinic": triclinic,
        "box_A": box,
    }
    return data_text, in_text, counts


def _stage_and_write_lammps(
    output_prefix: str | Path,
    *,
    data_text: str,
    in_text: str,
    prov_payload: dict[str, Any],
    manifest_payload: dict[str, Any],
) -> dict[str, Any]:
    prefix = Path(output_prefix).resolve()
    paths = {
        "data": prefix.with_suffix(".data"),
        "in": prefix.with_suffix(".in"),
        "manifest": prefix.with_suffix(".manifest.json"),
        "provenance": prefix.with_suffix(".provenance.json"),
    }
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("LAMMPS export files already exist; choose a fresh output prefix")

    data_bytes = data_text.encode("utf-8")
    in_bytes = in_text.replace("__DATA_FILE__", paths["data"].name).encode("utf-8")
    data_sha = hashlib.sha256(data_bytes).hexdigest()

    prov_payload["data_sha256"] = data_sha
    prov_payload["outputs"] = {k: str(v) for k, v in paths.items() if k != "provenance"}
    manifest_payload.update({k: str(v) for k, v in paths.items()})
    manifest_payload["data_sha256"] = data_sha

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
        with tempfile.TemporaryDirectory(dir=prefix.parent, prefix=f".{prefix.name}.tmp.") as stage_dir_str:
            stage_dir = Path(stage_dir_str)
            staged = {k: stage_dir / fp.name for k, fp in paths.items()}
            staged["data"].write_bytes(data_bytes)
            staged["in"].write_bytes(in_bytes)
            _json(staged["provenance"], prov_payload)
            _json(staged["manifest"], manifest_payload)
            for key in ("data", "in", "provenance", "manifest"):
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
    return manifest_payload


def export_periodic_lammps(
    model: dict[str, Any],
    output_prefix: str | Path,
    *,
    charge_model: str = "auto",
    charges_json: str | Path | None = None,
    reference_dir: str | Path | None = None,
    current_reference_dir: str | Path | None = None,
    cutoff_ang: float | None = None,
    allow_finite: bool = True,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Export a periodic or finite GO model directly with ``gaff2_v3`` to LAMMPS."""
    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh output prefix")
    if charges_json is not None and charge_model != "custom":
        raise ValueError(f"charges_json requires charge_model='custom' (got {charge_model!r})")

    has_pbc_bonds = _periodic_edges(model) or bool(model.get("box", {}).get("pbc", [False, False, False])[0]) or bool(model.get("box", {}).get("pbc", [False, False, False])[1])
    if not allow_finite and not has_pbc_bonds:
        raise ValueError("allow_finite=False requires a periodic GO model")

    amber_parm, assigned, meta = build_v3_parmed_structure(
        model,
        charge_model=charge_model,
        custom_charges=charges_json,
        reference_dir=reference_dir,
        current_reference_dir=current_reference_dir,
        allow_finite=allow_finite,
    )
    is_periodic = any(bool(x) for x in meta["pbc"])
    if is_periodic:
        resolved_cutoff_ang = (
            _resolve_periodic_cutoff_nm(meta, None if cutoff_ang is None else cutoff_ang * 0.1) * 10.0
        )
    else:
        resolved_cutoff_ang = 120.0 if cutoff_ang is None else float(cutoff_ang)

    data_text, in_text, counts = _serialize_amber_parm_to_lammps(
        amber_parm, is_periodic=is_periodic, cutoff_ang=resolved_cutoff_ang
    )
    ff_meta = assigned.get("metadata", {}).get("force_field", {})
    prov = {
        "schema_version": 1,
        "format": "lammps_full",
        "direct_from_model": True,
        "periodic": is_periodic,
        "preset": ff_meta.get("preset", "gaff2_v3"),
        "charge_model": ff_meta.get("charge_model"),
        "frcmod_sha256": meta["frcmod_sha256"],
        "bundle_provenance_sha256": meta["bundle_provenance_sha256"],
        "cutoff_ang": resolved_cutoff_ang,
        "conventions": {
            "units": "real",
            "atom_style": "full",
            "pair_modify": "mix arithmetic",
            "special_bonds": "amber (lj=0,0,0.5; coul=0,0,5/6)",
            "dihedral_style": "fourier",
            "improper_style": "cvff (amber quartet order preserved)",
        },
        "counts": counts,
    }
    manifest = {
        "status": "exported",
        "parameterized": True,
        "direct_periodic": is_periodic,
        "periodic": is_periodic,
        "cutoff_ang": resolved_cutoff_ang,
        "preset": ff_meta.get("preset", "gaff2_v3"),
        "charge_model": ff_meta.get("charge_model"),
        "net_charge": float(sum(a["charge"] for a in assigned["atoms"])),
        **counts,
        "engine_sources": LAMMPS_SOURCES,
    }
    return _stage_and_write_lammps(
        output_prefix, data_text=data_text, in_text=in_text, prov_payload=prov, manifest_payload=manifest
    )


def export_lammps(
    topology: str | Path,
    coordinates: str | Path,
    output_prefix: str | Path,
    *,
    cutoff_ang: float = 120.0,
    overwrite: bool = False,
    model: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Export a finite AMBER prmtop/inpcrd topology into LAMMPS ``data`` + ``in`` files."""
    if overwrite:
        raise ValueError("overwrite is unsupported; choose a fresh output prefix")
    if model and (_periodic_edges(model) or model.get("metadata", {}).get("covalently_periodic")):
        raise ValueError(
            "Periodic LAMMPS export from an AMBER prmtop is unsupported via tleap; use export_periodic_lammps with model.json"
        )
    structure, top_sha, crd_sha = _load_amber_with_digests(topology, coordinates)
    if not structure.atoms:
        raise ValueError("Amber topology contains no atoms")
    upstream = _link_upstream_provenance(
        topology, coordinates, top_sha=top_sha, crd_sha=crd_sha
    )
    data_text, in_text, counts = _serialize_amber_parm_to_lammps(
        structure, is_periodic=False, cutoff_ang=cutoff_ang
    )
    prov = {
        "schema_version": 1,
        "format": "lammps_full",
        **upstream,
        "periodic": False,
        "cutoff_ang": float(cutoff_ang),
        "counts": counts,
    }
    manifest = {
        "status": "exported",
        "parameterized": True,
        "periodic": False,
        "cutoff_ang": float(cutoff_ang),
        "source_topology_sha256": upstream["source_topology_sha256"],
        "preset": upstream["preset"],
        "charge_model": upstream["charge_model"],
        "net_charge": float(sum(float(a.charge) for a in structure.atoms)),
        **counts,
        "engine_sources": LAMMPS_SOURCES,
    }
    return _stage_and_write_lammps(
        output_prefix, data_text=data_text, in_text=in_text, prov_payload=prov, manifest_payload=manifest
    )
