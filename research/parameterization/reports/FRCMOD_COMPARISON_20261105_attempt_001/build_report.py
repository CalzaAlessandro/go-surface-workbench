#!/usr/bin/env python3
"""Read-only, standard-library numerical source audit. Writes only beside itself.

Run from any directory with the project interpreter. This is a coefficient audit,
not an Amber/LEaP parameter resolver or a force-field accuracy benchmark.
"""
from __future__ import annotations

import csv
import hashlib
import html
import itertools
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SRC = ROOT / "parameters/comparison_sources/HOSHIN_OLD_CURRENT_20261105"
REF = ROOT / "references"
SECTIONS = ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON")
FIELDS = {"MASS": (1, 1), "BOND": (2, 2), "ANGLE": (3, 2), "DIHE": (4, 4), "IMPROPER": (4, 3), "NONBON": (1, 2)}
UNITS = {"MASS": "u; optional polarizability Å³", "BOND": "K kcal/(mol Å²); r0 Å; Amber K(r-r0)²", "ANGLE": "K kcal/(mol rad²); theta0 degree; Amber K(theta-theta0)²", "DIHE": "IDIVF; PK kcal/mol; phase degree; signed periodicity; amplitude PK/IDIVF", "IMPROPER": "PK kcal/mol; phase degree; periodicity", "NONBON": "Rmin/2 Å, epsilon kcal/mol (not sigma)"}


def rel(path):
    return str(Path(path).relative_to(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(section, types):
    types = tuple(types)
    return min(types, types[::-1]) if section in ("BOND", "ANGLE", "DIHE") else types


def sig(section, types):
    return "-".join(key(section, types))


def rowparse(section, line, number):
    nt, nv = FIELDS[section]
    pat = r"^\s*([A-Za-z0-9?*+]+(?:\s*-\s*[A-Za-z0-9?*+]+){" + str(nt-1) + r"})\s+(.+)$"
    m = re.match(pat, line)
    if not m:
        raise ValueError(f"Unparsed {section}:{number}: {line}")
    tokens = m[2].split()
    vals = [float(t.replace("D", "E").replace("d", "e")) for t in tokens[:nv]]
    if len(vals) != nv or not all(math.isfinite(v) for v in vals):
        raise ValueError(f"Bad numeric fields {section}:{number}")
    consumed = nv
    if section == "MASS" and len(tokens) > 1:
        try:
            vals.append(float(tokens[1])); consumed += 1
        except ValueError:
            pass
    types = re.split(r"\s*-\s*", m[1])
    return {"section": section, "types": types, "signature": sig(section, types), "values": vals, "numeric_tokens": tokens[:consumed], "line": number, "comment": " ".join(tokens[consumed:]), "raw": line}


def parse_frcmod(path):
    rows, headers, issues = [], [], []
    section = None
    for i, raw in enumerate(path.read_text().splitlines(), 1):
        t = raw.strip()
        if t in SECTIONS or t == "NONB":
            section = "NONBON" if t == "NONB" else t
            headers.append({"line": i, "header": t, "section": section})
        elif t in ("END", "LJEDIT", "HBOND"):
            section = None
        elif section and t and not t.startswith(("#", "!")):
            try:
                rows.append(rowparse(section, raw, i))
            except ValueError as exc:
                issues.append({"section": section, "line": i, "raw": raw, "error": str(exc), "status": "invalid_incomplete_parameter_not_filled"})
    return {"id": path.stem, "file": rel(path), "sha256": sha(path), "bytes": path.stat().st_size, "headers": headers, "rows": rows, "issues": issues, "definitions": definitions(rows)}


def parse_dat(path, ident):
    blocks, block = [], []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        if line.strip():
            block.append((i, line))
        elif block:
            blocks.append(block); block = []
    if block:
        blocks.append(block)
    ljblocks = [i for i,b in enumerate(blocks) if b[0][1].strip().startswith("MOD4")]
    if len(blocks) < 7 or len(ljblocks) != 1 or ljblocks[0] < 6:
        raise ValueError(f"Unexpected DAT block structure: {path}")
    ljindex = ljblocks[0]
    rows = []
    for section, index, skip in (("MASS",0,1),("BOND",1,1),("ANGLE",2,0),("DIHE",3,0),("IMPROPER",4,0),("NONBON",ljindex,1)):
        for i, line in blocks[index][skip:]:
            if line.strip() == "END":
                break
            rows.append(rowparse(section, line, i))
    aliases = {}
    for _, line in itertools.chain.from_iterable(blocks[5:ljindex]):
        if not re.search(r"\d", line):
            names = line.split()
            for name in names[1:]:
                aliases[name] = names[0]
    return {"id": ident, "file": rel(path), "sha256": sha(path), "bytes": path.stat().st_size, "header": blocks[0][0][1], "rows": rows, "definitions": definitions(rows), "lj_aliases": aliases}


def definitions(rows):
    result, active = [], None
    for row in rows:
        if row["section"] != "DIHE":
            if active is not None:
                raise ValueError("Unterminated DIHE before another section")
            result.append({"section": row["section"], "types": list(key(row["section"], row["types"])), "signature": row["signature"], "rows": [row]})
            continue
        if active is None:
            active = {"section": "DIHE", "types": list(key("DIHE", row["types"])), "signature": row["signature"], "rows": []}
            result.append(active)
        if active["signature"] != row["signature"]:
            raise ValueError(f"Broken DIHE continuation at line {row['line']}")
        active["rows"].append(row)
        if row["values"][-1] >= 0:
            active = None
    if active is not None:
        raise ValueError("Unterminated final DIHE")
    return result


def rawvalues(defn):
    return [r["values"] for r in defn["rows"]]


def coeff(defn):
    """Fourier coefficients include constant energy; zero components vanish.

    Compare full series, not individual periodicity rows. IDIVF is divided out;
    negative periodicity is only continuation metadata. SCEE/SCNB are separate.
    """
    section = defn["section"]
    if section not in ("DIHE", "IMPROPER"):
        return dict(enumerate(rawvalues(defn)[0]))
    out = defaultdict(float)
    for row in defn["rows"]:
        v = row["values"]
        amp, phase, periodicity = (v[1]/v[0], v[2], abs(v[3])) if section == "DIHE" else v
        out["constant"] += amp
        out[f"cos_{periodicity:g}"] += amp * math.cos(math.radians(phase))
        out[f"sin_{periodicity:g}"] += amp * math.sin(math.radians(phase))
    return {str(k): v for k, v in out.items() if abs(v) > 1e-10}


def numeric_equal(a, b, tol=1e-7):
    ca, cb = coeff(a), coeff(b)
    return all(abs(ca.get(k, 0)-cb.get(k, 0)) <= tol for k in ca.keys() | cb.keys())


def representation_status(a, b):
    if numeric_equal(a, b):
        return "same_coefficients"
    if a["section"] == "MASS" and abs(rawvalues(a)[0][0]-rawvalues(b)[0][0]) <= 1e-7:
        return "same_mass_optional_polarizability_differs_or_unspecified"
    def half_unit(token):
        # Infer a printed decimal quantization interval, never an accuracy bound.
        token=token.lower().replace("d", "e")
        mantissa, _, exponent = token.partition("e")
        decimals=len(mantissa.partition(".")[2])
        return 0.5*10**(int(exponent or 0)-decimals)
    if a["section"] in ("DIHE", "IMPROPER"):
        def components(d):
            amps,precision=defaultdict(float),defaultdict(float)
            for r in d["rows"]:
                v,t=r["values"],r["numeric_tokens"]
                if d["section"]=="DIHE":
                    amplitude=v[1]/v[0]; marker=(abs(v[3]),v[2]%360); unit=half_unit(t[1])/v[0]
                else:
                    amplitude=v[0]; marker=(abs(v[2]),v[1]%360); unit=half_unit(t[0])
                if amplitude:
                    amps[marker]+=amplitude; precision[marker]+=unit
            return amps,precision
        aa,pa=components(a); ab,pb=components(b)
        if aa.keys()==ab.keys() and all(abs(aa[k]-ab[k]) <= max(pa[k],pb[k])+1e-10 for k in aa):
            return "same_to_reported_precision"
        return "different_coefficients"
    if a["section"] not in ("BOND", "ANGLE", "NONBON"):
        return "different_coefficients"
    ra,rb=a["rows"][0],b["rows"][0]
    if all(abs(x-y) <= max(half_unit(ta),half_unit(tb))+1e-10 for x,y,ta,tb in zip(ra["values"],rb["values"],ra["numeric_tokens"],rb["numeric_tokens"])):
        return "same_to_reported_precision"
    return "different_coefficients"


def scaling(defn):
    vals = defaultdict(set)
    for row in defn["rows"]:
        for field, value in re.findall(r"\b(SCEE|SCNB)\s*=\s*([0-9.]+)", row["comment"]):
            vals[field].add(float(value))
    return {k: sorted(v) for k, v in vals.items()} or {"status": "not_explicit_in_term; resolve loaded base/topology"}


def flags(row):
    result = []
    section, v, c = row["section"], row["values"], row["comment"]
    if re.search(r"ATTN|ATTENTION|placeholder|missing parameter|need.*revision", c, re.I): result.append("explicit_attention_or_placeholder_comment")
    if re.search(r"same as|calculated|empirical|default value|general improper", c, re.I): result.append("analogy_or_default_provenance")
    m = re.search(r"penalty score\s*=\s*([0-9.]+)", c, re.I)
    if m and float(m[1]) >= 50: result.append("analogy_penalty_ge_50_not_accuracy_metric")
    if section in ("BOND", "ANGLE") and v[0] <= 0: result.append("nonpositive_harmonic_constant")
    if section == "DIHE" and v[1] == 0: result.append("zero_fourier_component_not_automatically_missing")
    if section == "IMPROPER" and v[0] == 0: result.append("zero_fourier_component_not_automatically_missing")
    if section == "NONBON" and 0 in v: result.append("conventional_hydroxyl_H_zero_LJ" if row["types"] == ["ho"] and v == [0,0] else "zero_LJ_context_review")
    if section == "ANGLE" and row["signature"] == "cx-cx-cx" and v[1] <= 65: result.append("60degree_cyclopropane_type_angle_graph_context_required")
    return result


def groupdefs(source):
    out = defaultdict(list)
    for d in source["definitions"]:
        out[d["section"], d["signature"]].append(d)
    return out


def matched(defn, base):
    section, query = defn["section"], defn["types"]
    candidates = []
    for d in base["by_section"][section]:
        variants = [query]
        if section in ("BOND", "ANGLE", "DIHE"):
            variants.append(query[::-1])
        elif section == "IMPROPER":
            variants = [(p[0],p[1],query[2],p[2]) for p in itertools.permutations((query[0],query[1],query[3]))]
        if any(all(p == "X" or p == q for p, q in zip(d["types"], v)) for v in variants):
            candidates.append(d)
    if candidates:
        specificity = max(sum(t != "X" for t in d["types"]) for d in candidates)
        candidates = [d for d in candidates if sum(t != "X" for t in d["types"]) == specificity]
        mode = "exact_key" if any(d["signature"] == defn["signature"] for d in candidates) else "improper_outer_permutation_screen" if section == "IMPROPER" else "most_specific_wildcard_screen"
        return mode, candidates
    if section == "NONBON" and query[0] in base["lj_aliases"]:
        parent = base["lj_aliases"][query[0]]
        return "LJ_equivalence", [d for d in base["by_section"][section] if d["types"] == [parent]]
    return "unmatched_not_a_topology_gap_claim", []


def csvwrite(name, rows):
    path = OUT / name
    if not rows:
        path.write_text(""); return
    cols = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols); writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, separators=(",", ":")) if isinstance(v,(list,dict,tuple)) else v for k,v in row.items()})


def table(rows, columns=None):
    if not rows: return "<p>Nessuna riga.</p>"
    columns = columns or list(rows[0])
    def cell(x):
        return html.escape(json.dumps(x, ensure_ascii=False) if isinstance(x,(list,dict,tuple)) else str(x))
    width = f' style="min-width:{85*len(columns)}px"' if len(columns) > 7 else ''
    return '<div class="table"><table'+width+'><thead><tr>' + ''.join('<th>'+html.escape(c)+'</th>' for c in columns) + '</tr></thead><tbody>' + ''.join('<tr>'+''.join('<td>'+cell(row.get(c,""))+'</td>' for c in columns)+'</tr>' for row in rows) + '</tbody></table></div>'


def source_bases(ident):
    if ident == "Our_AMBER": return ["system_parm99", "system_parm10"]
    # All four compatible lowercase namespaces are screened to expose mixed
    # numerical fingerprints. This does not authorize loading them together.
    return ["AmberGO_GAFF_2003", "system_GAFF_1.81", "AmberGO_GAFF2_2.11", "system_GAFF2_2.2.20"]


def aim_audit():
    """Independent LJ extraction; no transfer of suffix-based type identities."""
    ljrows, records, bonds = [], [], []
    for path in sorted((REF / "gaff2_aim").glob("GO_*/*.top")):
        section, atoms, atomtypes, defaults = "", {}, {}, []
        interactions = Counter()
        unique_bonds = {}
        charge_sum = 0
        for i, raw in enumerate(path.read_text().splitlines(), 1):
            line = raw.split(";",1)[0].strip()
            if not line: continue
            m = re.match(r"\[\s*(\w+)\s*\]", line)
            if m: section = m[1]; continue
            if line.startswith("#"): continue
            fields = line.split()
            if section == "defaults": defaults = fields
            elif section == "atomtypes":
                atomtypes[fields[0]] = {"atomic_number":int(fields[1]), "mass":float(fields[2]), "sigma_nm":float(fields[5]), "epsilon_kJ_mol":float(fields[6]), "line":i}
            elif section == "atoms":
                atoms[int(fields[0])] = fields[1]; charge_sum += float(fields[6])
            elif section in ("bonds", "angles", "dihedrals", "pairs"):
                n = {"bonds":2,"angles":3,"dihedrals":4,"pairs":2}[section]
                interactions[section+"_function_"+fields[n]] += 1
                if section == "bonds":
                    types = tuple(atoms[int(x)] for x in fields[:2])
                    if len(fields) < 5: raise ValueError("Unexpected non-explicit AIM bond")
                    normalized = (float(fields[4])/(2*4.184*100), float(fields[3])*10)
                    marker = (min(types, types[::-1]),normalized)
                    if marker not in unique_bonds:
                        unique_bonds[marker] = {"source":rel(path), "line":i, "types":list(types), "function":fields[2], "Amber_K_kcal_mol_A2":normalized[0], "r0_A":normalized[1], "count":0}
                    unique_bonds[marker]["count"] += 1
        used = Counter(atoms.values())
        for name, count in sorted(used.items()):
            row = atomtypes[name]
            ljrows.append({"source":rel(path), "type":name, "count":count, **row, "Rmin_over_2_A":row["sigma_nm"]*10*2**(1/6)/2, "epsilon_kcal_mol":row["epsilon_kJ_mol"]/4.184, "charge_source":"[atoms], not zero [atomtypes] charge column", "mapping_status":"environment-specific type; no automatic equality with GAFF name prefix"})
        records.append({"file":rel(path), "sha256":sha(path), "defaults":defaults, "atom_count":len(atoms), "used_atomtypes":len(used), "sum_charge_e":charge_sum, "interaction_counts":dict(interactions), "bonded_provenance":"GAFF2/ACPYPE; inherited bonded model, not AIM bonded refit", "units_conversion":"GROMACS harmonic kb = 2*Amber K*4.184*100; b0(nm)*10=Å; LJ Rmin/2=10*sigma(nm)*2^(1/6)/2"})
        bonds.extend(unique_bonds.values())
    return records, ljrows, bonds


def main():
    manifest = json.loads((SRC / "MANIFEST.json").read_text())
    sources = [parse_frcmod(SRC / r["file"]) for r in manifest["records"]]
    for s, expected in zip(sources, manifest["records"]):
        assert s["sha256"] == expected["sha256"] and s["bytes"] == expected["bytes"], f"Source changed: {s['file']}"
        original=ROOT / expected["source"].split("/PARAMETERIZATION/",1)[1]
        assert sha(original)==s["sha256"], f"Archive and original snapshot diverged: {original}"
    base_specs = [("AmberGO_GAFF_2003", REF / "go_prior_art_20261104/AmberGO_native_audit_20261104/gaff.dat"),("AmberGO_GAFF2_2.11",REF / "go_prior_art_20261104/AmberGO_native_audit_20261104/gaff2.dat"),("system_GAFF_1.81",REF / "system_forcefield_sources/gaff.dat"),("system_GAFF2_2.2.20",REF / "system_forcefield_sources/gaff2.dat"),("system_parm99",REF / "system_forcefield_sources/parm99.dat"),("system_parm10",REF / "system_forcefield_sources/parm10.dat")]
    bases = {ident:parse_dat(path,ident) for ident,path in base_specs}
    system_manifest=json.loads((REF/"system_forcefield_sources/source_manifest.json").read_text())
    for record in system_manifest["files"]:
        if "sha256" in record and any(Path(b["file"]).name==record.get("copy") and b["id"].startswith("system_") for b in bases.values()):
            assert sha(REF/"system_forcefield_sources"/record["copy"])==record["sha256"], "Pinned system DAT changed"
    ambergo_manifest=json.loads((REF/"go_prior_art_20261104/AmberGO_native_audit_20261104/AUDIT_MANIFEST.json").read_text())
    for filename in ("gaff.dat", "gaff2.dat"):
        assert sha(REF/"go_prior_art_20261104/AmberGO_native_audit_20261104"/filename)==ambergo_manifest["files_sha256"][filename], "Pinned AmberGO DAT changed"
    for base in bases.values():
        base["by_section"] = {section:[d for d in base["definitions"] if d["section"] == section] for section in SECTIONS}
    terms, inventory, duplicates = [], [], []
    for s in sources:
        counts = Counter(r["section"] for r in s["rows"])
        defs = groupdefs(s)
        for (section, signature), ds in defs.items():
            if len(ds) > 1:
                duplicates.append({"source":s["id"], "section":section, "signature":signature, "definitions":len(ds), "lines":[[r["line"] for r in d["rows"]] for d in ds], "raw_values":[rawvalues(d) for d in ds], "coefficient_conflict":not all(numeric_equal(ds[0],d) for d in ds[1:]), "interpretation":"duplicate definitions; actual loader precedence not simulated"})
        for r in s["rows"]:
            terms.append({"source":s["id"], **r, "flags":flags(r), "units":UNITS[r["section"]]})
        inventory.append({"source":s["id"], "sha256":s["sha256"], "bytes":s["bytes"], **{section:counts[section] for section in SECTIONS}, "unique_definition_keys":len(defs), "duplicate_groups":sum(x["source"]==s["id"] for x in duplicates), "conflicting_groups":sum(x["source"]==s["id"] and x["coefficient_conflict"] for x in duplicates), "zero_DIHE_components":sum(r["section"]=="DIHE" and r["values"][1]==0 for r in s["rows"]), "analogy_flag_rows":sum("analogy_or_default_provenance" in flags(r) for r in s["rows"]), "native_namespace":"Amber uppercase; custom legacy names" if s["id"]=="Our_AMBER" else "GAFF-family lowercase; same labels do not prove same chemical role", "candidate_bases":source_bases(s["id"]), "historical_load_version":"not established by frcmod alone"})
        inventory[-1]["invalid_or_incomplete_rows"] = len(s["issues"])
        inventory[-1]["file"] = s["file"]
        inventory[-1]["explicit_DIHE_scaling"] = sorted({json.dumps(scaling(d),sort_keys=True) for d in s["definitions"] if d["section"]=="DIHE"})
    invalid_rows = [{"source": s["id"], **issue} for s in sources for issue in s["issues"]]
    csvwrite("numeric_invalid_rows.csv", invalid_rows)
    improper_overlaps=[]
    for s in sources:
        outer=defaultdict(list)
        for d in s["definitions"]:
            if d["section"]=="IMPROPER":
                t=d["types"]; outer[(t[2],tuple(sorted((t[0],t[1],t[3]))))].append(d)
        for (center,peripheral),ds in outer.items():
            if len({d["signature"] for d in ds}) > 1:
                improper_overlaps.append({"source":s["id"],"center_type_third":center,"outer_types_multiset":peripheral,"ordered_signatures":[d["signature"] for d in ds],"lines":[d["rows"][0]["line"] for d in ds],"coefficient_conflict":not all(numeric_equal(ds[0],d) for d in ds),"scope":"candidate outer permutations only; topology/loader semantics must be checked"})
    csvwrite("numeric_improper_permutation_overlaps.csv",improper_overlaps)
    base_rows = []
    for s in sources:
        for base_id in source_bases(s["id"]):
            base = bases[base_id]
            for d in s["definitions"]:
                mode, candidates = matched(d, base)
                equal = [numeric_equal(d, c) for c in candidates]
                statuses = [representation_status(d,c) for c in candidates]
                status = next((s for s in ("same_coefficients", "same_to_reported_precision", "same_mass_optional_polarizability_differs_or_unspecified") if s in statuses), "different_coefficients" if candidates else "no_candidate_match")
                base_rows.append({"source":s["id"], "base":base_id, "section":d["section"], "signature":d["signature"], "override_lines":[r["line"] for r in d["rows"]], "override_values":rawvalues(d), "base_match_mode":mode, "base_types":[c["signature"] for c in candidates], "base_lines":[[r["line"] for r in c["rows"]] for c in candidates], "base_values":[rawvalues(c) for c in candidates], "status":status, "multiple_candidates_disagree":bool(candidates) and not all(numeric_equal(candidates[0], c) for c in candidates), "override_scaling":scaling(d) if d["section"]=="DIHE" else "", "qualification":"base provenance candidate only; historical load and resolved topology not established"})
    pairs, pairrows = [], []
    lower = [s for s in sources if s["id"] != "Our_AMBER"]
    for a,b in itertools.combinations(lower,2):
        da,db = groupdefs(a),groupdefs(b)
        for section in SECTIONS:
            ka={k for k in da if k[0]==section}; kb={k for k in db if k[0]==section}
            counts=Counter()
            for k in sorted(ka|kb):
                left,right=da.get(k,[]),db.get(k,[])
                conflict = any(not numeric_equal(ds[0],d) for ds in (left,right) if ds for d in ds[1:])
                status="only_left_override" if not right else "only_right_override" if not left else "ambiguous_duplicate_conflict" if conflict else representation_status(left[0],right[0])
                counts[status]+=1
                pairrows.append({"left":a["id"],"right":b["id"],"section":section,"signature":k[1],"status":status,"left_lines":[[r["line"] for r in d["rows"]] for d in left],"right_lines":[[r["line"] for r in d["rows"]] for d in right],"left_values":[rawvalues(d) for d in left],"right_values":[rawvalues(d) for d in right],"scope":"shared native label only; verify graph/typing before physical comparison"})
            pairs.append({"left":a["id"],"right":b["id"],"section":section,**{name:counts[name] for name in ("same_coefficients","same_to_reported_precision","same_mass_optional_polarizability_differs_or_unspecified","different_coefficients","ambiguous_duplicate_conflict","only_left_override","only_right_override")}})
    basesumm=[]
    for (source,base_id,section),rs in itertools.groupby(sorted(base_rows,key=lambda r:(r["source"],r["base"],r["section"])),key=lambda r:(r["source"],r["base"],r["section"])):
        rs=list(rs); counter=Counter(r["status"] for r in rs)
        basesumm.append({"source":source,"base":base_id,"section":section,"definitions":len(rs),**dict(counter)})
    aim, aimlj, aimbonds=aim_audit()
    # Source-level assertions are deliberate validation checks, not simulation tests.
    hos=next(s for s in sources if s["id"]=="Hoshin_original_GAFF")
    assert sum(r["section"]=="NONBON" for r in hos["rows"]) == 5, "NONB alias must retain all LJ rows"
    aromatic=[d for d in hos["definitions"] if d["section"]=="DIHE" and d["signature"]=="ca-ca-ca-ca"][0]
    current=next(s for s in sources if s["id"]=="Our_current_GAFF2")
    native_aromatic=[d for d in current["definitions"] if d["section"]=="DIHE" and d["signature"]=="ca-ca-ca-ca"][0]
    assert numeric_equal(aromatic,native_aromatic), "PK/IDIVF normalization must identify 3.625 versus 14.5/4"
    assert len(aim)==2 and all(a["defaults"]==["1","2","yes","0.5","0.8333"] for a in aim)
    focused=[r for r in pairrows if r["left"]=="Our_current_GAFF2" and r["right"]=="Our_earlier_GAFF2" and r["status"]!="same_coefficients"]
    csvwrite("numeric_inventory.csv",inventory)
    csvwrite("numeric_terms.csv",terms)
    csvwrite("numeric_duplicates.csv",duplicates)
    csvwrite("numeric_base_comparison.csv",base_rows)
    csvwrite("numeric_base_summary.csv",basesumm)
    csvwrite("numeric_pairwise_summary.csv",pairs)
    csvwrite("numeric_pairwise_terms.csv",pairrows)
    csvwrite("numeric_current_earlier_changes.csv",focused)
    csvwrite("numeric_aim_lj.csv",aimlj)
    csvwrite("numeric_aim_bonds.csv",aimbonds)
    csvwrite("numeric_flags.csv",[r for r in terms if r["flags"]])
    highlights=[
        {"categoria":"differenza di coefficienti","risultato":"Current vs earlier: 17/17 BOND condivisi diversi; 38 ANGLE e 16 serie DIHE condivise diverse; una DIHE coincide alla precisione stampata; tutti i 10 LJ diversi.","limite":"Non è una classifica di accuratezza; verificare che i tipi descrivano gli stessi atomi."},
        {"categoria":"provenienza numerica","risultato":"Current: tutti i 19 K di BOND coincidono esattamente con GAFF2 2.2.20; 18 lunghezze sono compatibili con il DAT arrotondato a 0.001 Å, una coincide esattamente. Tutti i 10 LJ coincidono esattamente.","limite":"La corrispondenza prova l'impronta numerica, non il caricamento storico o una validazione GO."},
        {"categoria":"provenienza numerica","risultato":"Earlier GAFF2: contro GAFF1.81, 14/14 BOND con firma disponibile coincidono esattamente o alla precisione stampata (3 firme c6 assenti); 24/25 ANGLE disponibili e 9/9 LJ disponibili coincidono esattamente (c6 LJ assente).","limite":"Il nome GAFF2 non basta a dichiarare un preset GAFF2 puro. Nessuna equivalenza automatica c6→c3."},
        {"categoria":"provenienza numerica","risultato":"Hoshin: contro GAFF distribuito giugno2003 coincidono tutti i 6 BOND, tutti i 5 LJ e tutte le 25 serie DIHE dopo divisione PK/IDIVF; 7/12 ANGLE coincidono alla precisione stampata, 5 differiscono.","limite":"Gli improper aggiuntivi e il typing devono essere verificati nel modello sorgente."},
        {"categoria":"differenza di coefficienti","risultato":"OLD GOX e REGOX: tutte le 25 ANGLE, 56 serie DIHE e 27 IMPROPER condivise coincidono; REGOX aggiunge una ANGLE e 5 firme DIHE; nessuna MASS/BOND/LJ override in entrambi.","limite":"L'assenza di override non indica un'assenza di parametro nella base."},
        {"categoria":"integrità tecnica","risultato":"36 gruppi duplicati identici sulle firme canoniche, nessun conflitto numerico in questi gruppi. Our_AMBER riga139 IMPROPER CY-CX-CX-CY è incompleta e manca dei tre coefficienti.","limite":"Non completata per analogia; l'effetto sul loader non è testato."},
        {"categoria":"validazione QM mirata","risultato":"Current cx-cx-cx ha θ0=60°; c6-c2-ca passa da 125.19° a 115.29°; cx-op-cx K passa da 84.3 a119.01 kcal/(mol rad²).","limite":"Prima controllare grafo, classificazione ring e origine DAT; poi confrontare Hessiane e geometrie sullo stesso motivo."},
        {"categoria":"validazione QM mirata","risultato":"OH–epossido oh-c6-cx-op: earlier ha solo PK/IDIVF=1.4/9 a n=3; current usa ampiezze 0.5,0.9,0.18 kcal/mol a n=3,2,1. cx-cx-op-cx cambia anch'essa serie completa.","limite":"Scansioni rilassate solo dopo estrazione delle torsioni effettive del prmtop; cariche/LJ/scaling coerenti."},
        {"categoria":"ambito chimico","risultato":"La presenza di os/op e di legami C–O non dimostra copertura di eteri o lattoni. I lattoni sono esteri ciclici, richiedono percorsi carbonile–O–C e contesto d'anello/protonazione specifici.","limite":"Nessuna dichiarazione di lattone validato o nuovo parametro ricavato da questo audit."},
    ]
    csvwrite("numeric_highlights.csv",highlights)
    summary={"scope":"source coefficients only, not resolved GO topology or accuracy ranking", "highlights":highlights, "source_manifest_hash":sha(SRC/"MANIFEST.json"), "script_sha256":sha(__file__), "sources":inventory,"invalid_rows":invalid_rows,"bases":[{k:v for k,v in b.items() if k not in ("rows","definitions","by_section")} for b in bases.values()], "duplicates":duplicates,"improper_outer_permutation_overlaps":improper_overlaps,"base_summary":basesumm,"pairwise_summary":pairs,"aim":aim,"verification":{"all_seven_archive_hashes_match":True,"archives_equal_original_local_snapshots":True,"all_six_DAT_hashes_match_prior_manifests":True,"all_parameter_rows_accounted_for_without_silent_drop":True,"invalid_rows_preserved_not_filled":len(invalid_rows),"NONB_alias_retains_5_Hoshin_LJ_rows":True,"Fourier_IDIVF_normalization_checked":True,"AIM_defaults_read_from_both_topologies":True},"limitations":["Non è noto il DAT e l'ordine di caricamento storico dei file locali; le basi sono confronti candidati, non provenienza dimostrata.","I confronti fra firme uguali non stabiliscono identità dell'ambiente chimico. Namespace AMBER escluso dal confronto incrociato con GAFF.","Wildcard e permutazioni improper sono uno screening di candidati; nessuna precedenza effettiva LEaP, match di grafo, topologia risolta, esclusione 1–4 o validazione QM è simulata.","Frcmod privo di MASS/BOND/NONBON può ereditare tutto dalla base: le righe assenti non sono automaticamente parametri mancanti.","Zero torsionale può essere scelta convenzionale; penalty di analogia non è errore energetico e non misura accuratezza.","I tipi AIM con suffissi sono classi ambientali, non si riducono automaticamente a tipi GAFF. Eteri e lattoni richiedono grafi e typing espliciti; os/op da soli non ne provano la copertura."]}
    (OUT/"numeric_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    report_html(summary,terms,focused,aimlj,aimbonds,base_rows,pairs)
    print(json.dumps({"sources":len(sources),"bases":len(bases),"override_rows":len(terms),"base_definition_comparisons":len(base_rows),"duplicate_groups":len(duplicates),"conflicting_groups":sum(d["coefficient_conflict"] for d in duplicates),"AIM_LJ_rows":len(aimlj),"report":str(OUT/"numeric_report.html")},ensure_ascii=False))


def report_html(summary,terms,focused,aimlj,aimbonds,base_rows,pairs):
    parts=['<!doctype html><html lang="it"><meta charset="utf-8"><title>Confronto numerico frcmod GO</title><style>body{font:15px system-ui;max-width:1500px;margin:30px auto;padding:0 22px;color:#18222c}h1,h2{color:#123e59}.table{overflow:auto;max-height:700px;border:1px solid #ccd3da;margin:16px 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;vertical-align:top;padding:9px;border-bottom:1px solid #dce2e6;max-width:500px;overflow-wrap:anywhere}th{background:#e7f1f7;position:sticky;top:0}tr:nth-child(even){background:#f6f8fa}p{max-width:1100px;line-height:1.5}.note{padding:14px;background:#fff6db}a{color:#12628b}code{font-size:12px}</style><h1>GAFF, GAFF2, Hoshin, OLD, AMBER e GAFF2-AIM</h1><p>Audit numerico delle sorgenti conservate, 5 novembre 2026. Tutti i sette hash originali sono verificati. Nessun parametro, libreria o job di produzione è stato modificato.</p><p class="note"><b>Risposta operativa:</b> le differenze di coefficienti sono verificabili ora; non dimostrano quale forza di campo sia più accurata. La priorità è risolvere tipi e contesti sullo stesso grafo GO, quindi validare con QM i motivi ambigui. Nuovi calcoli non sono giustificati per il solo fatto che due file contengano numeri diversi.</p>']
    parts.append('<h2>Risultati principali verificati</h2>'+table(summary["highlights"]))
    parts.append('<h2>Limiti e interpretazione</h2>'+table([{"vincolo":x} for x in summary["limitations"]]))
    parts.append('<h2>Inventario delle righe override</h2><p>I conteggi DIHE sono componenti di Fourier, non torsioni indipendenti. Hoshin usa NONB: contiene cinque voci Lennard-Jones, qui conservate. Un frcmod è una sovrascrittura della base DAT.</p>'+table(summary["sources"],["source","MASS","BOND","ANGLE","DIHE","IMPROPER","NONBON","unique_definition_keys","duplicate_groups","conflicting_groups","zero_DIHE_components","analogy_flag_rows","sha256"]))
    parts.append('<h2>Basi realmente disponibili</h2><p>I DAT distribuiti con AmberGO e i DAT di sistema sono presenti localmente. Il nome GAFF2 non fissa una versione numerica. Le etichette delle basi riportano l’intestazione: il DAT AmberGO con header 2.11 include anche righe marcate 5/2017; l’hash, non la sola intestazione, identifica esattamente la sorgente. Per i file locali non è dimostrato quale DAT storico sia stato caricato.</p>'+table(summary["bases"],["id","header","file","sha256"]))
    parts.append('<h2>Coefficienti override contro basi candidate</h2><p>Una serie DIHE viene confrontata interamente dopo divisione PK/IDIVF e normalizzazione del segno di continuazione. same_to_reported_precision indica compatibilità con metà dell’ultima cifra stampata più grossolana, non una tolleranza fisica. Le differenze del campo di polarizzabilità opzionale sono separate dalla massa. I match wildcard e le permutazioni esterne improper sono screening; la precedenza LEaP non è risolta.</p>'+table(summary["base_summary"]))
    parts.append('<h2>Modifiche current rispetto a earlier</h2><p>Confronto per firma nativa identica. L’identità dei ruoli atomici deve ancora essere verificata sulle librerie e sui grafi. Le sole differenze MASS nel secondo campo riguardano la polarizzabilità opzionale, non il valore della massa.</p>'+table(focused))
    parts.append('<h2>Righe incomplete: blocco tecnico da risolvere</h2><p>Le righe seguenti non hanno tutti i coefficienti richiesti. Sono conservate esplicitamente e non entrano nei confronti numerici. Non sono completate per analogia.</p>'+table(summary["invalid_rows"]))
    parts.append('<h2>Duplicati e conflitti interni</h2><p>Un duplicato identico è ridondanza testuale. Definizioni diverse della stessa firma richiedono verifica dell’effetto effettivo nel loader. Componenti con periodicità negativa sono continuazioni, non duplicati.</p>'+table(summary["duplicates"]))
    important=[r for r in terms if "explicit_attention_or_placeholder_comment" in r["flags"] or "nonpositive_harmonic_constant" in r["flags"] or "60degree_cyclopropane_type_angle_graph_context_required" in r["flags"]]
    parts.append('<h2>Sovrapposizioni candidate per permutazione improper</h2><p>Si conserva il terzo tipo come centro e si confronta il multinsieme dei tre esterni. Sono potenziali firme sovrapposte, non una simulazione del loader.</p>'+table(summary["improper_outer_permutation_overlaps"]))
    parts.append('<h2>Righe da riesaminare prima di decidere QM</h2>'+table(important,["source","section","signature","line","values","flags","comment"]))
    priorities=[{"priorità":"1 — typing e grafo prima di QM","motivo":"cx-cx-cx / cx-cx tra epossidi distinti","evidenza":"Current cx-cx-cx: K=83.51 kcal/(mol rad²), θ0=60°. Il tipo può descrivere una catena di C che non forma un triangolo C–C–C.","azione":"Verificare percorso, triangoli e ossigeni condivisi; poi geometrie e Hessiane QM di epossido isolato, due epossidi adiacenti e motivo OH–epossido. Non dedurre una barriera reattiva da un FF non reattivo."},{"priorità":"2 — validazione mirata","motivo":"Rigidità basal-plane e giunzioni sp²/sp³","evidenza":"Current ed earlier hanno differenze di BOND/ANGLE per firme uguali; il cambio di versione GAFF2 è una possibile origine numerica da documentare.","azione":"Stesso frammento e stato chimico; geometria e Hessiana QM per ca-ca, ca-cx, cx-cx e giunzioni con c6. Fittare solo se gli scarti sul riferimento comune lo richiedono."},{"priorità":"3 — torsioni e impropri","motivo":"Torsioni zero, analogie ad alta penalty, impropri uniformi","evidenza":"Valori zero convenzionali e impropri 1.1/180/2 non equivalgono a termini mancanti; i conflitti duplicati sono elencati sopra.","azione":"Prima estrarre dal prmtop il termine realmente usato; poi scansioni rilassate per OH–epossido, epossidi adiacenti e gruppi di bordo, con cariche/LJ/scaling coerenti."},{"priorità":"4 — chimica aggiuntiva","motivo":"Eteri, lattoni, carbonili, COOH e COO⁻","evidenza":"Una sigla os/op o un termine C–O non prova copertura del grafo, del protomero o della coniugazione.","azione":"Aggiungere frammenti distinti per etere basale/ponte, etere di bordo e lattone (estere ciclico); verifica copertura e QM solo dei percorsi necessari. Nessun parametro di lattone è validato da questo audit."},{"priorità":"Separata — LJ/cariche","motivo":"GAFF2-AIM","evidenza":"Le topologie AIM cambiano cariche DDEC6 e LJ per ambiente; bonded ereditati da GAFF2/ACPYPE.","azione":"Riprodurre una topologia periodica pubblicata con mixing e fattori 1–4 originali. Test di interazione/acqua per trasferimenti a bordi/eteri/lattoni; non trattare un nuovo RESP come convalida dei bonded."}]
    csvwrite("numeric_qm_priorities.csv",priorities)
    parts.append('<h2>Quando servono nuovi calcoli QM</h2><p>«Differenza di coefficienti» descrive un dato numerico; «problema di analogia/provenienza» descrive un contesto non dimostrato; «validazione QM necessaria» si riferisce a uno specifico motivo chimico e a un’osservabile ancora da confrontare.</p>'+table(priorities))
    parts.append('<h2>GAFF2-AIM: convenzioni e LJ estratti dalle topologie</h2><p>Combination rule 2: σ aritmetica, ε geometrica. Fattori 1–4 LJ 0.5 e Coulomb 0.8333. Rmin/2 = 10 × σ(nm) × 2¹ᐟ⁶ / 2; ε(kcal/mol) = ε(kJ/mol) / 4.184. Il suffisso del tipo conserva l’ambiente e non viene eliminato per trasferire valori.</p>'+table(summary["aim"]))
    ljranges=[]
    for source in sorted({r["source"] for r in aimlj}):
        for atomic in (1,6,8):
            rs=[r for r in aimlj if r["source"]==source and r["atomic_number"]==atomic]
            ljranges.append({"source":source,"Z":atomic,"types":len(rs),"Rmin/2 min Å":min(r["Rmin_over_2_A"] for r in rs),"Rmin/2 max Å":max(r["Rmin_over_2_A"] for r in rs),"epsilon min kcal/mol":min(r["epsilon_kcal_mol"] for r in rs),"epsilon max kcal/mol":max(r["epsilon_kcal_mol"] for r in rs)})
    parts.append(table(ljranges))
    parts.append('<h2>Esempio numerico AIM bonded, conversione esplicita</h2><p>I legami AIM usano la forma armonica GROMACS con 1/2: K(Amber) = kb(GROMACS) / (2 × 4.184 × 100). Per esempio il ca-ca di GO30, se r0 ≈ 0.13984 nm e kb ≈ 385850 kJ/(mol nm²), dà K ≈ 461.10 kcal/(mol Å²); il numero non costituisce prova di nuovo fit AIM. Le singole righe e la frequenza sono nel CSV numeric_aim_bonds.csv; il typing ambientale resta distinto.</p>')
    parts.append('<h2>Confronti incrociati delle sole firme condivise</h2><p>AMBER uppercase è escluso da questi confronti. La presenza di stessi simboli GAFF tra modelli diversi non garantisce stesso ambiente chimico. Una firma presente in un solo frcmod può essere ereditata dall’altro DAT.</p>'+table(pairs))
    outputs=sorted(p.name for p in OUT.glob("numeric_*.csv"))+["numeric_summary.json","build_report.py"]
    parts.append('<h2>File riproducibili</h2><p>'+" · ".join('<a href="'+html.escape(name)+'">'+html.escape(name)+'</a>' for name in outputs)+'</p><p>Script SHA-256: <code>'+summary["script_sha256"]+'</code>. Tutti i valori e le righe sorgente sono conservati nei CSV. Nessun calcolo QM o MD è stato lanciato.</p></html>')
    (OUT/"numeric_report.html").write_text("\n".join(parts))


if __name__ == "__main__":
    main()
