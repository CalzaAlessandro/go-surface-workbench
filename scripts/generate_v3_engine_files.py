#!/usr/bin/env python3
"""Generate standalone cross-engine force-field files for gaff2_v3 from GO_gaff2_v3.frcmod."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

WORKBENCH_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKBENCH_ROOT / "src"))

from go_builder.ff_v3 import V3_RESIDUE_NAMES  # noqa: E402
from go_builder.topology import _CUSTOM_LEAP_ATOM_TYPES  # noqa: E402

V3_DIR = WORKBENCH_ROOT / "resources" / "reference" / "gaff2_v3"
FRCMOD_PATH = V3_DIR / "GO_gaff2_v3.frcmod"
PROV_PATH = V3_DIR / "PROVENANCE.json"

KCAL_TO_KJ = 4.184


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_frcmod(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8").splitlines()
    section = None
    masses: list[tuple[str, float, float]] = []
    bonds: list[tuple[str, str, float, float]] = []
    angles: list[tuple[str, str, str, float, float]] = []
    dihes: list[tuple[tuple[str, str, str, str], int, float, float, float]] = []
    impropers: list[tuple[tuple[str, str, str, str], float, float, float]] = []
    nonbon: dict[str, tuple[float, float]] = {}

    headers = {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}
    for raw in lines[1:]:
        s = raw.strip()
        if not s:
            continue
        if s in headers:
            section = s
            continue
        if section == "MASS":
            p = s.split()
            masses.append((p[0], float(p[1]), float(p[2])))
        elif section == "BOND":
            head, rest = raw[:5].strip(), raw[5:].split()
            a, b = [x.strip() for x in head.split("-")]
            bonds.append((a, b, float(rest[0]), float(rest[1])))
        elif section == "ANGLE":
            head, rest = raw[:8].strip(), raw[8:].split()
            a, b, c = [x.strip() for x in head.split("-")]
            angles.append((a, b, c, float(rest[0]), float(rest[1])))
        elif section == "DIHE":
            head, rest = raw[:11].strip(), raw[11:].split()
            q = tuple(x.strip() for x in head.split("-"))
            assert len(q) == 4
            dihes.append((q, int(round(float(rest[0]))), float(rest[1]), float(rest[2]), float(rest[3])))
        elif section == "IMPROPER":
            head, rest = raw[:11].strip(), raw[11:].split()
            q = tuple(x.strip() for x in head.split("-"))
            assert len(q) == 4
            impropers.append((q, float(rest[0]), float(rest[1]), float(rest[2])))
        elif section == "NONBON":
            p = s.split()
            nonbon[p[0]] = (float(p[1]), float(p[2]))

    return {
        "masses": masses,
        "bonds": bonds,
        "angles": angles,
        "dihes": dihes,
        "impropers": impropers,
        "nonbon": nonbon,
    }


def build_leaprc() -> str:
    lines = [
        "# Standalone AMBER LEaP resource file for GO_gaff2_v3 (18 custom GO atom types)",
        "source leaprc.gaff2",
        "addAtomTypes {",
    ]
    for tname in sorted(k for k in _CUSTOM_LEAP_ATOM_TYPES if k.startswith("g")):
        elem, hyb = _CUSTOM_LEAP_ATOM_TYPES[tname]
        lines.append(f'  {{ "{tname}" "{elem}" "{hyb}" }}')
    lines.extend(["}", "loadAmberParams GO_gaff2_v3.frcmod"])
    for res in V3_RESIDUE_NAMES:
        lines.append(f"loadoff {res}_gaff2_v3.lib")
    return "\n".join(lines) + "\n"


def _elem_atnum(tname: str) -> tuple[str, int]:
    elem, _ = _CUSTOM_LEAP_ATOM_TYPES[tname]
    atnum = {"C": 6, "H": 1, "O": 8}[elem]
    return elem, atnum


def build_gromacs_ffnonbonded(data: dict) -> str:
    lines = [
        "; Standalone GROMACS [ atomtypes ] for GO_gaff2_v3 (explicit atomic numbers for gmx / gmx_MMPBSA / ParmEd ChRad)",
        "[ atomtypes ]",
        "; name    at.num      mass    charge  ptype          sigma        epsilon",
    ]
    for tname, mass, _ in data["masses"]:
        _, atnum = _elem_atnum(tname)
        rmin2, eps_kcal = data["nonbon"][tname]
        sigma_nm = 2.0 * rmin2 * (2.0 ** (-1.0 / 6.0)) * 0.1
        eps_kj = eps_kcal * KCAL_TO_KJ
        lines.append(
            f"{tname:<6s} {atnum:6d} {mass:10.4f}  0.000000      A {sigma_nm:16.10e} {eps_kj:16.10e}"
        )
    return "\n".join(lines) + "\n"


def build_gromacs_ffbonded(data: dict) -> str:
    lines = [
        "; Standalone GROMACS bonded parameters for GO_gaff2_v3 (converted verbatim from GO_gaff2_v3.frcmod)",
        "[ bondtypes ]",
        "; i      j  func         b0_nm            kb_kj_mol_nm2",
    ]
    for a, b, k_amber, req in data["bonds"]:
        lines.append(f"{a:<5s} {b:<5s}    1 {req * 0.1:14.8f} {2.0 * k_amber * 418.4:18.8f}")

    lines.extend(["", "[ angletypes ]", "; i      j      k  func      th0_deg          cth_kj_mol_rad2"])
    for a, b, c, k_amber, th0 in data["angles"]:
        lines.append(f"{a:<5s} {b:<5s} {c:<5s}    1 {th0:12.6f} {2.0 * k_amber * KCAL_TO_KJ:18.8f}")

    lines.extend(["", "[ dihedraltypes ]", "; Proper torsions (func 9, specific then wildcard)"])
    lines.append("; i      j      k      l  func    phase_deg          kd_kj_mol  pn")
    spec_d = [d for d in data["dihes"] if "X" not in d[0]]
    wild_d = [d for d in data["dihes"] if "X" in d[0]]
    for q, idivf, pk, phase, pn in (*spec_d, *wild_d):
        per = abs(int(round(pn)))
        kd_kj = (pk / idivf) * KCAL_TO_KJ
        lines.append(
            f"{q[0]:<5s} {q[1]:<5s} {q[2]:<5s} {q[3]:<5s}    9 {phase:12.6f} {kd_kj:18.10f} {per:3d}"
        )

    lines.extend(["", "[ dihedraltypes ]", "; Improper torsions (func 4, central atom at position 3)"])
    lines.append("; i      j      k      l  func    phase_deg          kd_kj_mol  pn")
    for q, pk, phase, pn in data["impropers"]:
        per = abs(int(round(pn)))
        kd_kj = pk * KCAL_TO_KJ
        lines.append(
            f"{q[0]:<5s} {q[1]:<5s} {q[2]:<5s} {q[3]:<5s}    4 {phase:12.6f} {kd_kj:18.10f} {per:3d}"
        )
    return "\n".join(lines) + "\n"


def build_openmm_xml(data: dict) -> str:
    def _om_cls(t: str) -> str:
        return "" if t == "X" else t

    lines = [
        "<ForceField>",
        "  <!-- Standalone OpenMM ForceField for GO_gaff2_v3 (converted verbatim from GO_gaff2_v3.frcmod) -->",
        "  <!-- Parity note: HarmonicBondForce, HarmonicAngleForce, NonbondedForce, and all Proper torsions match -->",
        "  <!-- tleap AmberPrmtopFile to < 1e-4 kcal/mol (0.000000 on proper torsions). Note that OpenMM -->",
        "  <!-- <ForceField> wildcard <Improper> ordering='amber' breaks (I, J, L) peripheral-slot ties via -->",
        "  <!-- residue/template/element rules, whereas AMBER tleap breaks ties by atom index; use -->",
        "  <!-- AmberPrmtopFile or 08_export_openmm.py / export_periodic_openmm (<System> XML, < 1e-6 kJ/mol) -->",
        "  <!-- whenever exact prmtop improper quartet identity on puckered sp2 sites is required. -->",
        "  <AtomTypes>",
    ]
    for tname, mass, _ in data["masses"]:
        elem, _ = _elem_atnum(tname)
        lines.append(f'    <Type name="gaff2_v3-{tname}" class="{tname}" element="{elem}" mass="{mass:.4f}"/>')
    lines.append("  </AtomTypes>")

    lines.append("  <HarmonicBondForce>")
    for a, b, k_amber, req in data["bonds"]:
        length_nm = req * 0.1
        k_omm = 2.0 * k_amber * 418.4
        lines.append(f'    <Bond class1="{a}" class2="{b}" length="{length_nm:.8f}" k="{k_omm:.8f}"/>')
    lines.append("  </HarmonicBondForce>")

    lines.append("  <HarmonicAngleForce>")
    for a, b, c, k_amber, th0 in data["angles"]:
        th_rad = math.radians(th0)
        k_omm = 2.0 * k_amber * KCAL_TO_KJ
        lines.append(f'    <Angle class1="{a}" class2="{b}" class3="{c}" angle="{th_rad:.10f}" k="{k_omm:.8f}"/>')
    lines.append("  </HarmonicAngleForce>")

    lines.append('  <PeriodicTorsionForce ordering="amber">')
    # Group proper multi-term torsions preserving specific-first order
    grouped_proper: dict[tuple[str, str, str, str], list[tuple[int, float, float]]] = defaultdict(list)
    proper_order: list[tuple[str, str, str, str]] = []
    spec_d = [d for d in data["dihes"] if "X" not in d[0]]
    wild_d = [d for d in data["dihes"] if "X" in d[0]]
    for q, idivf, pk, phase, pn in (*spec_d, *wild_d):
        if q not in grouped_proper:
            proper_order.append(q)
        per = abs(int(round(pn)))
        k_kj = (pk / idivf) * KCAL_TO_KJ
        grouped_proper[q].append((per, math.radians(phase), k_kj))

    for q in proper_order:
        attrs = [f'class{i}="{_om_cls(t)}"' for i, t in enumerate(q, 1)]
        for idx, (per, ph_rad, k_kj) in enumerate(grouped_proper[q], 1):
            attrs.extend([f'periodicity{idx}="{per}"', f'phase{idx}="{ph_rad:.10f}"', f'k{idx}="{k_kj:.10f}"'])
        lines.append(f'    <Proper {" ".join(attrs)}/>')

    # Amber IMPROPER c1-c2-c3(central)-c4 -> OpenMM class1=c3(central), class2=c1, class3=c2, class4=c4
    for q, pk, phase, pn in data["impropers"]:
        per = abs(int(round(pn)))
        k_kj = pk * KCAL_TO_KJ
        ph_rad = math.radians(phase)
        c1, c2, c3, c4 = _om_cls(q[2]), _om_cls(q[0]), _om_cls(q[1]), _om_cls(q[3])
        lines.append(
            f'    <Improper class1="{c1}" class2="{c2}" class3="{c3}" class4="{c4}" '
            f'periodicity1="{per}" phase1="{ph_rad:.10f}" k1="{k_kj:.10f}"/>'
        )
    lines.append("  </PeriodicTorsionForce>")

    lines.append('  <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">')
    lines.append('    <UseAttributeFromResidue name="charge"/>')
    for tname, _, _ in data["masses"]:
        rmin2, eps_kcal = data["nonbon"][tname]
        sigma_nm = 2.0 * rmin2 * (2.0 ** (-1.0 / 6.0)) * 0.1
        eps_kj = eps_kcal * KCAL_TO_KJ
        lines.append(f'    <Atom class="{tname}" sigma="{sigma_nm:.10f}" epsilon="{eps_kj:.10f}"/>')
    lines.append("  </NonbondedForce>")
    lines.append("</ForceField>")
    return "\n".join(lines) + "\n"


def build_lammps_inc(data: dict) -> str:
    lines = [
        "# Standalone LAMMPS reference settings & 18-type LJ table for GO_gaff2_v3",
        "# Intended for use alongside GO Workbench exported .data files (atom_style full).",
        "units           real",
        "atom_style      full",
        "boundary        p p p",
        "bond_style      harmonic",
        "angle_style     harmonic",
        "dihedral_style  fourier",
        "improper_style  cvff",
        "special_bonds   amber",
        "pair_modify     mix arithmetic",
        "# Pair Coeffs table (units real: epsilon [kcal/mol], sigma [Angstrom]):",
    ]
    for tname, mass, _ in data["masses"]:
        elem, _ = _elem_atnum(tname)
        rmin2, eps_kcal = data["nonbon"][tname]
        sig_ang = 2.0 * rmin2 * (2.0 ** (-1.0 / 6.0))
        lines.append(
            f"#   type={tname:<4s} elem={elem} mass={mass:7.3f} eps={eps_kcal:8.4f} sigma={sig_ang:10.6f} (Rmin/2={rmin2:7.4f})"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    data = parse_frcmod(FRCMOD_PATH)
    outputs = {
        "leaprc.go_gaff2_v3": build_leaprc(),
        "GO_gaff2_v3_ffnonbonded.itp": build_gromacs_ffnonbonded(data),
        "GO_gaff2_v3_ffbonded.itp": build_gromacs_ffbonded(data),
        "GO_gaff2_v3_openmm.xml": build_openmm_xml(data),
        "GO_gaff2_v3_lammps.inc": build_lammps_inc(data),
    }

    prov = json.loads(PROV_PATH.read_text(encoding="utf-8"))
    for fname, content in outputs.items():
        p = V3_DIR / fname
        raw = content.encode("utf-8")
        p.write_bytes(raw)
        prov["files"][fname] = {"sha256": _sha256_bytes(raw), "bytes": len(raw)}
        print(f"Wrote {fname} ({len(raw)} bytes, sha256={prov['files'][fname]['sha256'][:12]}...)")

    PROV_PATH.write_text(json.dumps(prov, indent=2) + "\n", encoding="utf-8")
    print(f"Updated {PROV_PATH.name} ({len(prov['files'])} files total)")


if __name__ == "__main__":
    main()
