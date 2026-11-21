"""Block 2 versioned AMBER GAFF2 v3 force field & dual-regime L0_v2 RESP charges.

Implements:
- Integrity-verified resource loading from ``resources/reference/gaff2_v3/``.
- Graph-only chemical perception covering Priority 1+2 GO functional groups
  (GGG, GRT, HYD, EPX, EPH, EPT, CAB, CBX, CBD), invariant to atom numbering,
  residue labels, and Kekulé bond-order representation.
- Explicit rejection of Priority 3 unsupported chemistry (phenols, planar ethers,
  lactones/esters, lactols, aldehydes) and unparameterized covalent junctions.
- Deterministic 4-color DSATUR graph coloring of epoxide units (gx, gy, gz, ge)
  eliminating the intra- vs inter-epoxide 58.5 deg / 115 deg angle ambiguity.
- Dual-regime L0_v2 charge calculator (fixed role charges when N_COO == 0;
  graph-distance sp2 sheet delocalization exp(-0.30 * d_graph) when N_COO > 0)
  plus auditable per-atom custom charge injection.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import deque
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

from .ff import FRAGMENTS as CURRENT_FRAGMENTS
from .ff import _assign_dedicated, load_fragment_library

DEFAULT_V3_REFERENCE = Path(__file__).resolve().parents[2] / "resources/reference/gaff2_v3"
DEFAULT_CURRENT_REFERENCE = Path(__file__).resolve().parents[2] / "resources/reference/current"

EPX_COLOR_TYPES = ("gx", "gy", "gz", "ge")
RES_COLORS_EPX = ("EPX", "EP2", "EP3", "EP4")
RES_COLORS_EPH = ("EPH", "EH2", "EH3", "EH4")
RES_COLORS_EPT = ("EPT", "ET2", "ET3", "ET4")

V3_RESIDUE_NAMES = (
    "GGG",
    "GRT",
    "HYD",
    "CAB",
    "CBX",
    "CBD",
    *RES_COLORS_EPX,
    *RES_COLORS_EPH,
    *RES_COLORS_EPT,
)

FORCE_FIELD_PRESETS: dict[str, dict[str, Any]] = {
    "gaff2_v3": {
        "name": "gaff2_v3",
        "preset_id": "gaff2_v3",
        "version": "3.0.0",
        "tier": "recommended",
        "default": True,
        "charge_model": "l0_v2",
        "frcmod_filename": "GO_gaff2_v3.frcmod",
        "bonded_model": "GAFF2 + GO_gaff2_v3.frcmod (18 lowercase GO-only types, 4-color DSATUR epoxides)",
        "status": "production_validated_block2",
        "description": (
            "Block 2 validated AMBER GAFF2 v3 (18 GO-only types, 4-color epoxides) "
            "+ L0_v2 dual-regime HF/6-31G* MK RESP charges."
        ),
    },
    "gaff2_current": {
        "name": "gaff2_current",
        "preset_id": "gaff2_current",
        "version": "1.0.0",
        "tier": "legacy",
        "default": False,
        "charge_model": "current_lib",
        "frcmod_filename": "GO_gaff2.frcmod",
        "bonded_model": "GAFF2 + GO_gaff2.frcmod (legacy 10-type fragment libraries)",
        "status": "legacy_deprecated",
        "description": (
            "Deprecated legacy 10-type GO_gaff2.frcmod + fragment library charges "
            "(known defects on COO-/COOH; retained only for backwards compatibility)."
        ),
    },
}

_PRESET_ALIASES = {
    "gaff2_v3": "gaff2_v3",
    "gaff2-v3": "gaff2_v3",
    "v3": "gaff2_v3",
    "recommended": "gaff2_v3",
    "resp-l0-v2": "gaff2_v3",
    "gaff2_current": "gaff2_current",
    "gaff2-current": "gaff2_current",
    "current": "gaff2_current",
    "resp-current": "gaff2_current",
    "legacy": "gaff2_current",
}

_ROLE_TO_CURRENT_LIB_KEY = {
    "scaffold_C": ("GGG", "CX"),
    "edge_C": ("GRT", "CH"),
    "edge_H": ("GRT", "HC"),
    "hydroxyl_C": ("HYD", "CL"),
    "hydroxyl_O": ("HYD", "OL"),
    "hydroxyl_H": ("HYD", "HL"),
    "epoxide_C": ("EPX", "CY"),
    "epoxide_O": ("EPX", "OE"),
    "EPH_C": ("EPH", "CZ"),
    "EPH_CH": ("EPH", "CK"),
    "EPH_H": ("EPH", "HK"),
    "EPH_O": ("EPH", "OE"),
    "EPT_C": ("EPT", "CP"),
    "EPT_H": ("EPT", "HP"),
    "EPT_O": ("EPT", "OE"),
    "carbonyl_C": ("CAB", "CN"),
    "carbonyl_O": ("CAB", "ON"),
    "COOH_anchor": ("CBX", "CY"),
    "COOH_C": ("CBX", "C4"),
    "COOH_O": ("CBX", "OJ"),
    "COOH_OH": ("CBX", "OK"),
    "COOH_H": ("CBX", "HK"),
    "COO_anchor": ("CBD", "CY"),
    "COO_C": ("CBD", "C4"),
    "COO_O": ("CBD", "OJ"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_force_field_preset(name: str | None) -> str:
    """Resolve a user preset name or alias into a canonical preset identifier."""
    if name is None:
        return "gaff2_v3"
    key = str(name).strip().lower()
    canonical = _PRESET_ALIASES.get(key)
    if canonical is None:
        choices = ", ".join(sorted(FORCE_FIELD_PRESETS))
        raise ValueError(f"Unknown force-field preset {name!r}; choose from {choices}")
    return canonical


def list_force_field_presets() -> list[dict[str, Any]]:
    """Return metadata records for packaged AMBER GO force-field presets."""
    return [dict(FORCE_FIELD_PRESETS[k]) for k in ("gaff2_v3", "gaff2_current")]


def verify_v3_bundle(reference_dir: str | Path | None = None) -> dict[str, Any]:
    """Verify presence and SHA-256 checksums of every file in resources/reference/gaff2_v3."""
    ref_dir = Path(reference_dir or DEFAULT_V3_REFERENCE).resolve()
    prov_path = ref_dir / "PROVENANCE.json"
    if not prov_path.is_file():
        raise ValueError(f"Missing gaff2_v3 PROVENANCE.json in {ref_dir}")
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    files = provenance.get("files")
    if not isinstance(files, dict):
        raise ValueError("Malformed PROVENANCE.json: missing files manifest")
    required_files = {"GO_gaff2_v3.frcmod", "L0_v2_charges.json"} | {
        f"{res}_gaff2_v3.lib" for res in V3_RESIDUE_NAMES
    }
    missing_required = sorted(required_files - set(files))
    if missing_required:
        raise ValueError(
            f"Incomplete gaff2_v3 bundle manifest; missing required entries: {missing_required}"
        )
    for fname, meta in files.items():
        target = ref_dir / fname
        if not target.is_file():
            raise ValueError(f"Missing gaff2_v3 bundle file: {target}")
        if _sha256(target) != meta["sha256"]:
            raise ValueError(f"gaff2_v3 SHA-256 mismatch for {fname}")
    return {
        "reference_dir": str(ref_dir),
        "provenance": provenance,
        "provenance_sha256": _sha256(prov_path),
        "frcmod_path": ref_dir / "GO_gaff2_v3.frcmod",
        "frcmod_sha256": files["GO_gaff2_v3.frcmod"]["sha256"],
        "l0_v2_path": ref_dir / "L0_v2_charges.json",
        "l0_v2_sha256": files["L0_v2_charges.json"]["sha256"],
    }


def _check_v3_preconditions(model: dict[str, Any]) -> None:
    metadata = model.get("metadata", {})
    if metadata.get("requires_parameterization") or any(
        a.get("metadata", {}).get("junction") for a in model.get("atoms", [])
    ):
        raise ValueError(
            "New covalent junction requires dedicated charge/parameterization template and frcmod"
        )
    atoms = model.get("atoms", [])
    if not atoms:
        raise ValueError("Model contains no atoms")
    bad_el = sorted({a.get("element") for a in atoms if a.get("element") not in ("C", "H", "O")})
    if bad_el:
        raise ValueError(
            f"Unsupported elements in gaff2_v3 GO domain (expected C/H/O only): {bad_el}"
        )


def _classify_single_epoxy(
    c1: int, c2: int, o_idx: int, atoms: list[dict[str, Any]], adj: list[list[int]]
) -> dict[int, str]:
    el = [a["element"] for a in atoms]
    sub: dict[int, str] = {}
    hs_by_c = {}
    for c in (c1, c2):
        if len(adj[c]) != 4:
            raise ValueError(f"Epoxide carbon {atoms[c]['id']} must have coordination 4")
        c_nbs = [j for j in adj[c] if el[j] == "C"]
        o_nbs = [j for j in adj[c] if el[j] == "O"]
        h_nbs = [j for j in adj[c] if el[j] == "H"]
        if o_nbs != [o_idx] or len(h_nbs) > 1 or len(c_nbs) + len(h_nbs) != 3:
            raise ValueError(
                f"Unsupported Priority-3 substitution at epoxide carbon {atoms[c]['id']}"
            )
        hs_by_c[c] = h_nbs
    n_h1, n_h2 = len(hs_by_c[c1]), len(hs_by_c[c2])
    if n_h1 == 0 and n_h2 == 0:
        sub[o_idx] = "epoxide_O"
        sub[c1] = "epoxide_C"
        sub[c2] = "epoxide_C"
    elif n_h1 + n_h2 == 1:
        sub[o_idx] = "EPH_O"
        if n_h1 == 1:
            sub[c1], sub[c2] = "EPH_CH", "EPH_C"
            sub[hs_by_c[c1][0]] = "EPH_H"
        else:
            sub[c2], sub[c1] = "EPH_CH", "EPH_C"
            sub[hs_by_c[c2][0]] = "EPH_H"
    elif n_h1 == 1 and n_h2 == 1:
        sub[o_idx] = "EPT_O"
        sub[c1], sub[c2] = "EPT_C", "EPT_C"
        sub[hs_by_c[c1][0]] = "EPT_H"
        sub[hs_by_c[c2][0]] = "EPT_H"
    else:
        raise ValueError("Unsupported >1 H substitution on epoxide pair")
    return sub


def _classify_hydroxyl(
    c: int, o_idx: int, h_idx: int, atoms: list[dict[str, Any]], adj: list[list[int]]
) -> dict[int, str]:
    el = [a["element"] for a in atoms]
    if len(adj[c]) == 4 and all(el[j] == "C" for j in adj[c] if j != o_idx):
        return {c: "hydroxyl_C", o_idx: "hydroxyl_O", h_idx: "hydroxyl_H"}
    if len(adj[c]) == 3 and sorted(el[j] for j in adj[c]) == ["C", "O", "O"]:
        return _classify_carboxyl(c, atoms, adj)
    if len(adj[c]) == 3:
        raise ValueError(
            f"Unsupported Priority-3 phenolic sp2-OH environment at carbon {atoms[c]['id']}"
        )
    raise ValueError(
        f"Unsupported Priority-3 hydroxyl environment (lactol/gem-diol/CHOH) at atom {atoms[c]['id']}"
    )


def _classify_carboxyl(c: int, atoms: list[dict[str, Any]], adj: list[list[int]]) -> dict[int, str]:
    el = [a["element"] for a in atoms]
    if len(adj[c]) != 3:
        raise ValueError(f"Carboxyl carbon {atoms[c]['id']} must be trigonal")
    c_nbs = [j for j in adj[c] if el[j] == "C"]
    o_nbs = [j for j in adj[c] if el[j] == "O"]
    if len(c_nbs) != 1 or len(o_nbs) != 2:
        raise ValueError(f"Unsupported carboxyl connectivity at carbon {atoms[c]['id']}")
    anchor = c_nbs[0]
    if len(adj[anchor]) != 3 or any(el[j] != "C" for j in adj[anchor]):
        raise ValueError(
            f"Carboxyl anchor {atoms[anchor]['id']} must be a bare sp2 ring carbon"
        )
    oh_oxy, c_o_oxy = [], []
    for oi in o_nbs:
        if len(adj[oi]) == 1:
            c_o_oxy.append(oi)
        elif len(adj[oi]) == 2 and any(el[j] == "H" for j in adj[oi]):
            oh_oxy.append(oi)
        else:
            raise ValueError(
                f"Unsupported Priority-3 ester/lactone/peroxide O at atom {atoms[oi]['id']}"
            )
    if len(oh_oxy) == 1 and len(c_o_oxy) == 1:
        oi = oh_oxy[0]
        hi = next(j for j in adj[oi] if el[j] == "H")
        return {
            anchor: "COOH_anchor",
            c: "COOH_C",
            c_o_oxy[0]: "COOH_O",
            oi: "COOH_OH",
            hi: "COOH_H",
        }
    if len(c_o_oxy) == 2:
        return {anchor: "COO_anchor", c: "COO_C", c_o_oxy[0]: "COO_O", c_o_oxy[1]: "COO_O"}
    raise ValueError(f"Unsupported dihydroxy/ester carboxyl motif at carbon {atoms[c]['id']}")


def classify_go_v3_roles(model: dict[str, Any]) -> list[str]:
    """Perceive all 25 Priority 1+2 GO chemical roles strictly from (element, bonds) graph."""
    _check_v3_preconditions(model)
    atoms = model["atoms"]
    n = len(atoms)
    id2i = {a["id"]: i for i, a in enumerate(atoms)}
    adj: list[list[int]] = [[] for _ in range(n)]
    for b in model["bonds"]:
        i, j = id2i[b["a"]], id2i[b["b"]]
        adj[i].append(j)
        adj[j].append(i)

    roles = [""] * n

    def assign(idx: int, role: str) -> None:
        if roles[idx] and roles[idx] != role:
            raise ValueError(
                f"Conflicting chemical environments at atom {atoms[idx]['id']}: "
                f"{roles[idx]} vs {role}"
            )
        roles[idx] = role

    for i, a in enumerate(atoms):
        if a["element"] != "O":
            continue
        c_nbs = [j for j in adj[i] if atoms[j]["element"] == "C"]
        h_nbs = [j for j in adj[i] if atoms[j]["element"] == "H"]
        if len(c_nbs) + len(h_nbs) != len(adj[i]):
            raise ValueError(f"Unsupported O-O or non-CH bond at oxygen {a['id']}")
        if len(c_nbs) == 2 and not h_nbs:
            c1, c2 = c_nbs
            if c2 not in adj[c1]:
                raise ValueError(
                    f"Unsupported Priority-3 non-epoxide ether/lactone oxygen at atom {a['id']}"
                )
            sub = _classify_single_epoxy(c1, c2, i, atoms, adj)
        elif len(c_nbs) == 1 and len(h_nbs) == 1:
            sub = _classify_hydroxyl(c_nbs[0], i, h_nbs[0], atoms, adj)
        elif len(c_nbs) == 1 and not h_nbs:
            c = c_nbs[0]
            el_nbs = sorted(atoms[j]["element"] for j in adj[c])
            if len(adj[c]) == 3 and el_nbs == ["C", "C", "O"]:
                sub = {c: "carbonyl_C", i: "carbonyl_O"}
            elif len(adj[c]) == 3 and el_nbs == ["C", "O", "O"]:
                sub = _classify_carboxyl(c, atoms, adj)
            else:
                raise ValueError(
                    f"Unsupported Priority-3 C=O environment (e.g. aldehyde) at atom {atoms[c]['id']}"
                )
        else:
            raise ValueError(f"Unsupported oxygen valence at atom {a['id']}")
        for idx, role in sub.items():
            assign(idx, role)

    for i, a in enumerate(atoms):
        if roles[i]:
            continue
        if a["element"] == "C":
            els = sorted(atoms[j]["element"] for j in adj[i])
            if els == ["C", "C", "H"]:
                hi = next(j for j in adj[i] if atoms[j]["element"] == "H")
                assign(i, "edge_C")
                assign(hi, "edge_H")
            elif els == ["C", "C", "C"]:
                assign(i, "scaffold_C")
            else:
                raise ValueError(
                    f"Unsupported carbon environment at atom {a['id']} (neighbors={els})"
                )
        elif a["element"] == "H":
            if len(adj[i]) != 1 or atoms[adj[i][0]]["element"] != "C":
                raise ValueError(f"Unsupported hydrogen valence at atom {a['id']}")

    if any(not r for r in roles):
        bad = [atoms[i]["id"] for i, r in enumerate(roles) if not r]
        raise ValueError(f"Unclassified atoms in GO v3 perception: {bad}")
    return roles


def color_epoxides_v3(
    model: dict[str, Any], roles: list[str]
) -> tuple[dict[int, int], list[dict[str, Any]]]:
    """Deterministically color epoxide units {0,1,2,3} via DSATUR so adjacent epoxides differ."""
    atoms = model["atoms"]
    id2i = {a["id"]: i for i, a in enumerate(atoms)}
    adj: list[set[int]] = [set() for _ in range(len(atoms))]
    for b in model["bonds"]:
        i, j = id2i[b["a"]], id2i[b["b"]]
        adj[i].add(j)
        adj[j].add(i)

    epx_units: list[dict[str, Any]] = []
    atom_to_epx: dict[int, int] = {}
    for o_idx, r in enumerate(roles):
        if r in ("epoxide_O", "EPH_O", "EPT_O"):
            cs = sorted(j for j in adj[o_idx] if atoms[j]["element"] == "C")
            if len(cs) != 2:
                raise ValueError(f"Expected 2 epoxide carbons at atom {atoms[o_idx]['id']}")
            uid = len(epx_units)
            epx_units.append({"uid": uid, "o_idx": o_idx, "c_idxs": cs, "role": r})
            for c_idx in cs:
                if c_idx in atom_to_epx:
                    raise ValueError(f"Atom {atoms[c_idx]['id']} belongs to multiple epoxides")
                atom_to_epx[c_idx] = uid

    if not epx_units:
        return {}, []

    eadj: dict[int, set[int]] = {u: set() for u in range(len(epx_units))}
    for u, info in enumerate(epx_units):
        for c_idx in info["c_idxs"]:
            for nb in adj[c_idx]:
                v = atom_to_epx.get(nb)
                if v is not None and v != u:
                    eadj[u].add(v)

    colors: dict[int, int] = {}
    uncolored = set(range(len(epx_units)))
    while uncolored:
        def key_fn(u: int) -> tuple[int, int, int, int]:
            sat = len({colors[v] for v in eadj[u] if v in colors})
            deg = len(eadj[u])
            return (-sat, -deg, epx_units[u]["c_idxs"][0], epx_units[u]["c_idxs"][1])

        u = min(uncolored, key=key_fn)
        used_c = {colors[v] for v in eadj[u] if v in colors}
        col = next((c for c in range(4) if c not in used_c), None)
        if col is None:
            raise ValueError("Epoxide adjacency graph requires more than 4 colors (exceeds v3 palette)")
        colors[u] = col
        uncolored.remove(u)

    for u, nbs in eadj.items():
        for v in nbs:
            if colors[u] == colors[v]:
                raise RuntimeError("Adjacent epoxides received identical DSATUR color")

    atom_color = {c_idx: colors[u] for c_idx, u in atom_to_epx.items()}
    return atom_color, epx_units


def decompose_go_v3(
    model: dict[str, Any], roles: list[str] | None = None
) -> tuple[list[str], list[tuple[int, str, list[tuple[int, str]]]], list[str], dict[int, int]]:
    """Decompose GO graph into clean 1-unit v3 residues & return (roles, res_list, v3_types, atom_color)."""
    roles = classify_go_v3_roles(model) if roles is None else roles
    atoms = model["atoms"]
    n = len(atoms)
    id2i = {a["id"]: idx for idx, a in enumerate(atoms)}
    adj: list[list[int]] = [[] for _ in range(n)]
    for b in model["bonds"]:
        i, j = id2i[b["a"]], id2i[b["b"]]
        adj[i].append(j)
        adj[j].append(i)

    atom_color, _ = color_epoxides_v3(model, roles)
    used: set[int] = set()
    res_list: list[tuple[int, str, list[tuple[int, str]]]] = []
    v3_types = [""] * n
    rnum = 1

    for i in range(n):
        if i in used:
            continue
        r = roles[i]
        if r == "edge_C":
            hs = [j for j in adj[i] if roles[j] == "edge_H"]
            res_list.append((rnum, "GRT", [(i, "CH"), (hs[0], "HC")]))
            v3_types[i], v3_types[hs[0]] = "gb", "gh"
            used.update((i, hs[0]))
            rnum += 1
        elif r == "hydroxyl_C":
            o = next(j for j in adj[i] if roles[j] == "hydroxyl_O")
            h = next(j for j in adj[o] if roles[j] == "hydroxyl_H")
            res_list.append((rnum, "HYD", [(i, "CL"), (o, "OL"), (h, "HL")]))
            v3_types[i], v3_types[o], v3_types[h] = "g6", "gm", "gw"
            used.update((i, o, h))
            rnum += 1
        elif r == "epoxide_O":
            cs = sorted(j for j in adj[i] if roles[j] == "epoxide_C")
            col = atom_color[cs[0]]
            ct = EPX_COLOR_TYPES[col]
            res_list.append((rnum, RES_COLORS_EPX[col], [(cs[0], "CY"), (cs[1], "CZ"), (i, "OE")]))
            v3_types[cs[0]] = v3_types[cs[1]] = ct
            v3_types[i] = "gp"
            used.update((cs[0], cs[1], i))
            rnum += 1
        elif r == "EPH_O":
            c_ch = next(j for j in adj[i] if roles[j] == "EPH_CH")
            c_c = next(j for j in adj[i] if roles[j] == "EPH_C")
            h_k = next(j for j in adj[c_ch] if roles[j] == "EPH_H")
            col = atom_color[c_c]
            ct = EPX_COLOR_TYPES[col]
            res_list.append(
                (rnum, RES_COLORS_EPH[col], [(c_c, "CZ"), (c_ch, "CK"), (h_k, "HK"), (i, "OE")])
            )
            v3_types[c_c] = v3_types[c_ch] = ct
            v3_types[h_k], v3_types[i] = "g1", "gp"
            used.update((c_c, c_ch, h_k, i))
            rnum += 1
        elif r == "EPT_O":
            cs = sorted(j for j in adj[i] if roles[j] == "EPT_C")
            h0 = next(j for j in adj[cs[0]] if roles[j] == "EPT_H")
            h1 = next(j for j in adj[cs[1]] if roles[j] == "EPT_H")
            col = atom_color[cs[0]]
            ct = EPX_COLOR_TYPES[col]
            res_list.append(
                (
                    rnum,
                    RES_COLORS_EPT[col],
                    [(cs[0], "CP"), (cs[1], "CK"), (h1, "HK"), (i, "OE"), (h0, "HP")],
                )
            )
            v3_types[cs[0]] = v3_types[cs[1]] = ct
            v3_types[h0] = v3_types[h1] = "g1"
            v3_types[i] = "gp"
            used.update((cs[0], cs[1], h0, h1, i))
            rnum += 1
        elif r == "COOH_anchor":
            c4 = next(j for j in adj[i] if roles[j] == "COOH_C")
            oj = next(j for j in adj[c4] if roles[j] == "COOH_O")
            ok = next(j for j in adj[c4] if roles[j] == "COOH_OH")
            hk = next(j for j in adj[ok] if roles[j] == "COOH_H")
            res_list.append((rnum, "CBX", [(i, "CY"), (c4, "C4"), (oj, "OJ"), (ok, "OK"), (hk, "HK")]))
            v3_types[i], v3_types[c4], v3_types[oj], v3_types[ok], v3_types[hk] = (
                "gb",
                "gc",
                "go",
                "gm",
                "gw",
            )
            used.update((i, c4, oj, ok, hk))
            rnum += 1
        elif r == "COO_anchor":
            c4 = next(j for j in adj[i] if roles[j] == "COO_C")
            os_ = sorted(j for j in adj[c4] if roles[j] == "COO_O")
            res_list.append((rnum, "CBD", [(i, "CY"), (c4, "C4"), (os_[0], "OJ"), (os_[1], "OK")]))
            v3_types[i], v3_types[c4] = "gb", "gd"
            v3_types[os_[0]] = v3_types[os_[1]] = "gu"
            used.update((i, c4, os_[0], os_[1]))
            rnum += 1
        elif r == "carbonyl_C":
            o = next(j for j in adj[i] if roles[j] == "carbonyl_O")
            res_list.append((rnum, "CAB", [(i, "CN"), (o, "ON")]))
            v3_types[i], v3_types[o] = "gk", "gq"
            used.update((i, o))
            rnum += 1

    for i in range(n):
        if i in used:
            continue
        if roles[i] == "scaffold_C":
            res_list.append((rnum, "GGG", [(i, "CX")]))
            v3_types[i] = "ga"
            used.add(i)
            rnum += 1

    if len(used) != n or any(not t for t in v3_types):
        missing = [atoms[i]["id"] for i in range(n) if i not in used]
        raise ValueError(f"Unassigned atoms in GO v3 residue decomposition: {missing}")

    return roles, res_list, v3_types, atom_color


def _expected_formal_charge(model: dict[str, Any], roles: list[str]) -> int:
    n_coo = sum(1 for r in roles if r == "COO_anchor")
    expected = -n_coo
    if "formal_charge" in model and model["formal_charge"] is not None:
        if int(model["formal_charge"]) != expected:
            raise ValueError(
                f"Model formal_charge={model['formal_charge']} disagrees with deprotonated carboxylate count (-{n_coo})"
            )
    elif not model.get("metadata", {}).get("formal_charges_unknown"):
        fc_vals = [int(a.get("formal_charge") or 0) for a in model["atoms"]]
        if any(v != 0 for v in fc_vals) and sum(fc_vals) != expected:
            raise ValueError(
                f"Per-atom formal_charge sum ({sum(fc_vals)}) disagrees with carboxylate count (-{n_coo})"
            )
    return expected


def compute_l0_v2_charges(
    model: dict[str, Any],
    roles: list[str] | None = None,
    *,
    reference_dir: str | Path | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Compute exact L0_v2 (L0_v2c_a030) dual-regime HF/6-31G* MK RESP charges."""
    info = verify_v3_bundle(reference_dir)
    l0_data = json.loads(Path(info["l0_v2_path"]).read_text(encoding="utf-8"))
    roles = classify_go_v3_roles(model) if roles is None else roles
    atoms = model["atoms"]
    n_atoms = len(atoms)
    id2i = {a["id"]: i for i, a in enumerate(atoms)}
    adj: list[list[int]] = [[] for _ in range(n_atoms)]
    for b in model["bonds"]:
        ia, ib = id2i[b["a"]], id2i[b["b"]]
        adj[ia].append(ib)
        adj[ib].append(ia)

    role_q = l0_data["role_charges_e"]
    acceptor_roles = set(l0_data["acceptor_roles"])
    exclude_roles = set(l0_data["exclude_bfs_roles"])
    alpha = float(l0_data["alpha"])
    d_max = int(l0_data["d_max"])
    delta_sheet = float(l0_data["delta_sheet_e"])

    q = np.array([float(role_q[r]) for r in roles], dtype=float)
    n_coo = sum(1 for r in roles if r == "COO_anchor")
    expected_fc = _expected_formal_charge(model, roles)

    if abs(delta_sheet) > 0.0 and n_coo > 0:
        m_mat = np.zeros((n_atoms, 4), dtype=float)
        for i, r in enumerate(roles):
            if r == "COO_C":
                m_mat[i, 0] = 1.0
            elif r == "COO_O":
                m_mat[i, 1] = 1.0
            elif r == "COO_anchor":
                m_mat[i, 2] = 1.0
                dist: dict[int, int] = {i: 0}
                dq = deque([i])
                while dq:
                    u = dq.popleft()
                    if dist[u] >= d_max:
                        continue
                    for v in adj[u]:
                        if v not in dist and roles[v] not in exclude_roles:
                            dist[v] = dist[u] + 1
                            dq.append(v)
                targets = [
                    (j, 1.0 if alpha == 0.0 else math.exp(-alpha * d))
                    for j, d in dist.items()
                    if 1 <= d <= d_max and roles[j] in acceptor_roles
                ]
                if not targets:
                    raise ValueError(
                        f"No eligible sp2 acceptor carbon within d_max={d_max} bonds of carboxylate anchor {atoms[i]['id']}"
                    )
                wsum = sum(w for _, w in targets)
                for j, w in targets:
                    m_mat[j, 3] += w / wsum
        for i, r in enumerate(roles):
            if r in ("COO_C", "COO_O", "COO_anchor"):
                q[i] = 0.0
        x_coo = np.array(
            [role_q["COO_C"], role_q["COO_O"], role_q["COO_anchor"], delta_sheet], dtype=float
        )
        q = q + m_mat @ x_coo
        regime = "anionic_sp2_delocalization"
    else:
        regime = "neutral_fixed_roles"

    residual = float(q.sum()) - float(expected_fc)
    if abs(residual) >= 1e-9:
        raise ValueError(
            f"L0_v2 charge conservation check failed: sum={float(q.sum())}, expected={expected_fc}"
        )

    charge_meta = {
        "charge_model": "L0_v2c_a030",
        "regime": regime,
        "n_coo": n_coo,
        "formal_charge_e": expected_fc,
        "net_charge_e": float(q.sum()),
        "charge_residual_e": residual,
        "alpha": alpha,
        "d_max": d_max,
        "delta_sheet_e": delta_sheet,
        "l0_v2_sha256": info["l0_v2_sha256"],
    }
    return q, charge_meta


def compute_current_lib_role_charges(
    model: dict[str, Any],
    roles: list[str] | None = None,
    *,
    reference_dir: str | Path | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Compute legacy fixed-library charges directly from resources/reference/current/*_gaff2.lib."""
    roles = classify_go_v3_roles(model) if roles is None else roles
    cur_dir = Path(reference_dir or DEFAULT_CURRENT_REFERENCE).resolve()
    libs = {
        name: load_fragment_library(cur_dir / f"{name}_gaff2.lib", name)
        for name in CURRENT_FRAGMENTS
    }
    lookup = {
        role: next(a["charge"] for a in libs[resname]["atoms"] if a["name"] == aname)
        for role, (resname, aname) in _ROLE_TO_CURRENT_LIB_KEY.items()
    }
    q = np.array([float(lookup[r]) for r in roles], dtype=float)
    expected_fc = _expected_formal_charge(model, roles)
    residual = float(q.sum()) - float(expected_fc)
    if abs(residual) >= 1e-6:
        raise ValueError(f"current_lib net charge {float(q.sum())} != expected {expected_fc}")
    return q, {
        "charge_model": "current_lib",
        "regime": "fixed_fragment_libraries",
        "n_coo": -expected_fc,
        "formal_charge_e": expected_fc,
        "net_charge_e": float(q.sum()),
        "charge_residual_e": residual,
        "source_dir": str(cur_dir),
    }


def _to_charge_float(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"Custom charge value must be a finite number, got {type(value).__name__}")
    val = float(value)
    if not math.isfinite(val):
        raise ValueError("Custom charge value must be finite")
    return val


def _to_atom_id(key: Any) -> int:
    if isinstance(key, bool):
        raise ValueError("Custom charge atom ID cannot be a boolean")
    try:
        return int(key)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Custom charge atom ID must be an integer, got {key!r}") from exc


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise ValueError(f"Duplicate JSON object key in custom charges: {k!r}")
        out[k] = v
    return out


def _parse_id_charge_map(raw_map: dict[Any, Any]) -> dict[int, float]:
    by_id: dict[int, float] = {}
    for k, v in raw_map.items():
        aid = _to_atom_id(k)
        if aid in by_id:
            raise ValueError(
                f"Duplicate atom ID {aid} after normalization in custom charges map (raw key {k!r})"
            )
        by_id[aid] = _to_charge_float(v)
    return by_id


def load_custom_charges(
    charges_source: str | Path | dict[str, Any] | list[float] | np.ndarray,
    model: dict[str, Any],
    roles: list[str] | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load and validate explicit per-atom partial charges with full provenance tracking."""
    atoms = model["atoms"]
    n = len(atoms)
    roles = classify_go_v3_roles(model) if roles is None else roles
    expected_fc = _expected_formal_charge(model, roles)

    source_sha256 = None
    source_desc = "inline_custom"
    external_prov: dict[str, Any] | None = None

    if isinstance(charges_source, (str, Path)):
        p = Path(charges_source).resolve()
        if not p.is_file():
            raise FileNotFoundError(f"Custom charges JSON not found: {p}")
        source_sha256 = _sha256(p)
        source_desc = str(p)
        payload = json.loads(
            p.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    elif isinstance(charges_source, (list, np.ndarray)):
        payload = [_to_charge_float(v) for v in charges_source]
    elif isinstance(charges_source, dict):
        payload = charges_source
    else:
        raise TypeError(f"Unsupported custom charges type: {type(charges_source).__name__}")

    if isinstance(payload, dict):
        external_prov = payload.get("provenance") if isinstance(payload.get("provenance"), dict) else None
        if "charges_e" in payload:
            raw_list = payload["charges_e"]
            if not isinstance(raw_list, list):
                raise ValueError("Custom 'charges_e' must be a JSON list of finite numbers")
            q = np.array([_to_charge_float(v) for v in raw_list], dtype=float)
        elif "charges_by_id" in payload:
            raw_map = payload["charges_by_id"]
            if not isinstance(raw_map, dict):
                raise ValueError(
                    "Custom 'charges_by_id' must be a JSON object mapping atom IDs to numbers"
                )
            by_id = _parse_id_charge_map(raw_map)
            if set(by_id) != {int(a["id"]) for a in atoms}:
                raise ValueError("Custom charges_by_id keys must match complete model atom IDs")
            q = np.array([by_id[int(a["id"])] for a in atoms], dtype=float)
        else:
            by_id = _parse_id_charge_map(payload)
            if set(by_id) != {int(a["id"]) for a in atoms}:
                raise ValueError("Custom charge dict keys must match model atom IDs or contain charges_e")
            q = np.array([by_id[int(a["id"])] for a in atoms], dtype=float)
    elif isinstance(payload, list):
        q = np.array([_to_charge_float(v) for v in payload], dtype=float)
    else:
        raise ValueError("Custom charges JSON must be a list or object")

    if q.shape != (n,) or not np.all(np.isfinite(q)):
        raise ValueError(f"Custom charges must have shape ({n},) with all finite floats")

    residual = float(q.sum()) - float(expected_fc)
    if abs(residual) > 1e-4:
        raise ValueError(
            f"Custom charge sum ({float(q.sum()):.6f} e) differs from formal charge ({expected_fc} e) beyond 1e-4 e"
        )

    return q, {
        "charge_model": "custom",
        "regime": "per_atom_override",
        "n_coo": -expected_fc,
        "formal_charge_e": expected_fc,
        "net_charge_e": float(q.sum()),
        "charge_residual_e": residual,
        "source": source_desc,
        "source_sha256": source_sha256,
        "external_provenance": external_prov,
    }


def assign_gaff2_v3(
    model: dict[str, Any],
    *,
    reference_dir: str | Path | None = None,
    charge_model: str = "auto",
    custom_charges: str | Path | dict[str, Any] | list[float] | np.ndarray | None = None,
    current_reference_dir: str | Path | None = None,
    parameter_templates: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a deep copy with GAFF2 v3 18-type assignments, clean v3 residues, and auditable charges."""
    assigned = deepcopy(model)
    metadata = assigned.setdefault("metadata", {})
    metadata["parameterized"] = False
    if custom_charges is not None and charge_model != "custom":
        raise ValueError(
            f"custom_charges/charges_json requires charge_model='custom' (got {charge_model!r})"
        )
    if parameter_templates is not None:
        if charge_model != "auto" or custom_charges is not None:
            raise ValueError(
                "parameter_templates supplies dedicated charges and cannot be combined with explicit charge_model/custom_charges overrides"
            )
        out = _assign_dedicated(assigned, parameter_templates)
        for atom in out["atoms"]:
            ameta = atom.setdefault("metadata", {})
            ameta["charge_source"] = "dedicated_parameter_templates"
            ameta.pop("v3_role", None)
            ameta.pop("epoxide_color", None)
        out.setdefault("metadata", {})["force_field"] = {
            "preset": "gaff2_v3",
            "version": FORCE_FIELD_PRESETS["gaff2_v3"]["version"],
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
        return out

    bundle = verify_v3_bundle(reference_dir)
    roles, res_list, v3_types, atom_color = decompose_go_v3(assigned)

    if charge_model in ("auto", "l0_v2", "L0_v2"):
        q_vec, charge_meta = compute_l0_v2_charges(
            assigned, roles, reference_dir=bundle["reference_dir"]
        )
        custom_combo = False
        q_tag = "gaff2_v3:L0_v2c_a030"
    elif charge_model in ("current", "current_lib"):
        q_vec, charge_meta = compute_current_lib_role_charges(
            assigned, roles, reference_dir=current_reference_dir
        )
        custom_combo = True
        q_tag = "custom_combo:gaff2_v3+current_lib"
    elif charge_model == "custom":
        if custom_charges is None:
            raise ValueError("custom_charges is required when charge_model='custom'")
        q_vec, charge_meta = load_custom_charges(custom_charges, assigned, roles)
        custom_combo = True
        q_tag = f"custom:{charge_meta['source']}"
    else:
        raise ValueError(
            f"Unsupported charge_model {charge_model!r} for gaff2_v3; choose l0_v2, current_lib, or custom"
        )

    for rnum, rname, pairs in res_list:
        for idx, aname in pairs:
            atom = assigned["atoms"][idx]
            atom["resid"] = rnum
            atom["resname"] = rname
            atom["name"] = aname
            atom["atom_type"] = v3_types[idx]
            atom["charge"] = float(q_vec[idx])
            ameta = atom.setdefault("metadata", {})
            ameta["v3_role"] = roles[idx]
            ameta["charge_source"] = f"{q_tag}:{roles[idx]}"
            if idx in atom_color:
                ameta["epoxide_color"] = int(atom_color[idx])

    metadata["force_field"] = {
        "preset": "gaff2_v3",
        "version": bundle["provenance"]["version"],
        "tier": "custom_combination" if custom_combo else bundle["provenance"]["tier"],
        "custom_combination": custom_combo,
        "name": "GAFF2 + GO_gaff2_v3.frcmod",
        "bonded_model": "GO_gaff2_v3 (18 lowercase GO-only types, 4-color DSATUR epoxides)",
        "charge_model": charge_meta["charge_model"],
        "charge_regime": charge_meta["regime"],
        "charge_info": charge_meta,
        "frcmod_sha256": bundle["frcmod_sha256"],
        "provenance_sha256": bundle["provenance_sha256"],
        "epoxide_colors_used": len(set(atom_color.values())),
    }
    return assigned
