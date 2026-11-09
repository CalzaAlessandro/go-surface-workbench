#!/usr/bin/env python3
"""Bounded, offline audit of already fitted molecular RESP charges.

Run from any directory with the existing workspace Python via uv run --no-project.
All outputs are confined to this script's directory. No fits or QM jobs are run.
"""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from collections import defaultdict, deque
from pathlib import Path
from statistics import mean, median

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
SOURCES = {
    "HF": ROOT / "PARAMETERIZATION/reports/CHARGE_COMPARISON_20261103_attempt_001/context_HF_16/context_analysis.json",
    "DFT": ROOT / "PARAMETERIZATION/reports/CHARGE_CONTROLS_UPDATE_20261103_attempt_001/context_DFT_19/context_analysis.json",
}
DOC = ROOT / "PARAMETERIZATION/docs/RESP_CONTEXT_ANALYSIS.md"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stats(values):
    if not values:
        return {k: None for k in ("n", "mean_e", "mean_abs_e", "median_abs_e", "rms_e", "min_e", "max_e", "max_abs_e", "n_abs_gt_0p05", "n_abs_gt_0p10")}
    return {"n": len(values), "mean_e": mean(values), "mean_abs_e": mean(map(abs, values)),
            "median_abs_e": median(map(abs, values)), "rms_e": math.sqrt(mean(x*x for x in values)),
            "min_e": min(values), "max_e": max(values), "max_abs_e": max(map(abs, values)),
            "n_abs_gt_0p05": sum(abs(x) > .05 for x in values),
            "n_abs_gt_0p10": sum(abs(x) > .10 for x in values)}


def csv_write(name, rows):
    assert rows, name
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()})


def dump(name, data):
    (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def distances(anchors, adjacency):
    result = {a: 0 for a in anchors}
    queue = deque(anchors)
    while queue:
        a = queue.popleft()
        for b in adjacency[a]:
            if b not in result:
                result[b] = result[a] + 1
                queue.append(b)
    return result


def signature(model):
    # Excludes workflow metadata, includes chemical identity and the exact geometry.
    return {"formal_charge": model["formal_charge"], "multiplicity": model["multiplicity"],
            "atoms": [{k: a.get(k) for k in ("id", "element", "name", "role", "roles", "formal_charge", "position", "symmetry_class", "chiral_tag", "cip_code")} for a in model["atoms"]],
            "bonds": sorted((min(b["a"], b["b"]), max(b["a"], b["b"]), b["order"]) for b in model["bonds"])}


def table(rows, columns):
    def fmt(v):
        if v is None:
            return "—"
        if isinstance(v, float):
            return f"{v:.4f}"
        if isinstance(v, (list, dict)):
            return json.dumps(v, ensure_ascii=False)
        return str(v)
    return '<div class="scroll"><table><thead><tr>' + ''.join('<th>'+html.escape(title)+'</th>' for _, title in columns) + '</tr></thead><tbody>' + ''.join('<tr>'+''.join('<td>'+html.escape(fmt(row.get(k)))+'</td>' for k, _ in columns)+'</tr>' for row in rows)+'</tbody></table></div>'


def main():
    manifest = []
    def record(path, kind, **kwargs):
        manifest.append({"path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path), "kind": kind, "sha256": sha(path), "bytes": path.stat().st_size, **kwargs})
    record(DOC, "schema_and_scope_documentation")
    record(Path(__file__).resolve(), "analysis_script")
    datasets = {method: json.loads(path.read_text()) for method, path in SOURCES.items()}
    for method, path in SOURCES.items():
        record(path, "context_analysis", method=method)
    context_rows, atom_rows, bg_rows, bg_summary, shell_rows, groups, pairs, neighborhoods = [], [], [], [], [], [], [], []
    originals, contexts = {}, {}
    source_gates = []
    for method, data in datasets.items():
        assert data["coverage"]["all_models_passed_charge_identity_and_total_checks"]
        assert data["context_count"] == data["model_count"] == len(data["models"])
        for m in data["models"]:
            mid, prov = m["model_id"], m["provenance"]
            key = (method, mid)
            assert key not in contexts
            contexts[key] = m
            modelpath, chargepath, evidencepath = (Path(prov[k]) for k in ("model_path", "charges_path", "fit_evidence_path"))
            original, chargefile, evidence = (json.loads(p.read_text()) for p in (modelpath, chargepath, evidencepath))
            originals[key] = original
            for path, kind in ((modelpath, "model"), (chargepath, "charges"), (evidencepath, "fit_evidence")):
                record(path, kind, method=method, model_id=mid)
            atoms = {a["atom_id"]: a for a in m["per_atom_charges"]}
            assert len(atoms) == len(original["atoms"]) == len(chargefile["atoms"])
            assert [(a["atom_id"], a["element"], a["charge_e"]) for a in m["per_atom_charges"]] == [(a["atom_id"], a["element"], a["charge_e"]) for a in chargefile["atoms"]]
            assert [(a["atom_id"], a["element"], a["position_angstrom"]) for a in m["per_atom_charges"]] == [(a["id"], a["element"], a["position"]) for a in original["atoms"]]
            assert m["charge_check"]["passed"] and m["coverage"]["atom_identity_order_elements_match"]
            observed = sum(a["charge_e"] for a in atoms.values())
            assert abs(observed - original["formal_charge"]) <= m["charge_check"]["rounding_tolerance_e"]
            assert abs(observed - m["charge_check"]["observed_charge_sum_e"]) < 1e-12
            embedded = evidence.get("artifacts_sha256", {})
            for p in (modelpath, chargepath):
                expected = embedded.get(p.name)
                source_gates.append({"method": method, "model_id": mid, "artifact": p.name, "embedded_sha256": expected, "actual_sha256": sha(p), "match": expected == sha(p) if expected else None})
            scaf = set(m["scaffold_carbon_ids"])
            bg = m["background_carbon_distance_charges"]
            assert {a["atom_id"] for a in bg} == set(m["background_carbon_atom_ids"])
            anchors = set().union(*(set(g["scaffold_carbon_anchor_ids"]) for g in m["functional_groups"])) if m["functional_groups"] else set()
            prefix = {"method": method, "model_id": mid, "scaffold_C": len(scaf), "context_id": m["context_id"]}
            context_rows.append({**prefix, "n_atoms": len(atoms), "n_groups": len(m["functional_groups"]), "functionalized_anchor_C": len(anchors), "anchor_loading_fraction": len(anchors)/len(scaf), "formal_charge_e": original["formal_charge"], "charge_sum_e": observed, "charge_error_e": observed-original["formal_charge"], "physical_conformers": prov["physical_conformers"], "orientations": prov["orientations"], "is_diagnostic": prov["is_diagnostic"], "reported_wavefunction_stability_validated": prov["reported_wavefunction_stability_validated"], "fit_status": prov["fit_status"], "fit_status_is_acceptance": prov["fit_status_is_acceptance"], "final_parameter_validation": evidence.get("final_parameter_validation"), "qm_method": prov["qm_method"], "scope": prov["scope"], "model_scope": prov["model_scope"], "protocol_scope": prov["protocol_scope"], "upstream_qm_gate": prov["upstream_qm_gate"]})
            for a in atoms.values():
                atom_rows.append({**prefix, **a})
            adjacency = {a: set() for a in scaf}
            neighbors = defaultdict(set)
            for bond in original["bonds"]:
                a, b = bond["a"], bond["b"]
                neighbors[a].add(b); neighbors[b].add(a)
                if a in scaf and b in scaf:
                    adjacency[a].add(b); adjacency[b].add(a)
            gd = {g["group_id"]: distances(g["scaffold_carbon_anchor_ids"], adjacency) for g in m["functional_groups"]}
            for a in bg:
                assert a["charge_e"] == atoms[a["atom_id"]]["charge_e"]
                expected_distance = min((ds[a["atom_id"]] for ds in gd.values()), default=None)
                assert expected_distance == a["nearest_functional_group_graph_distance"]
                bg_rows.append({**prefix, **a})
            for kind in ("interior_scaffold", "hydrogen_terminated_edge"):
                selected = [a for a in bg if a["background_kind"] == kind]
                bg_summary.append({**prefix, "background_kind": kind, **stats([a["charge_e"] for a in selected])})
                bins = defaultdict(list)
                for a in selected:
                    bins[a["nearest_functional_group_graph_distance"]].append(a["charge_e"])
                for distance, values in sorted(bins.items(), key=lambda x: -1 if x[0] is None else x[0]):
                    shell_rows.append({**prefix, "background_kind": kind, "distance_CC_bonds": distance, **stats(values)})
            for g in m["functional_groups"]:
                ids, anchorids = set(g["atom_ids"]), set(g["scaffold_carbon_anchor_ids"])
                q = sum(atoms[a]["charge_e"] for a in ids)
                qa = sum(atoms[a]["charge_e"] for a in anchorids)
                assert abs(q-g["net_charge_e"]) < 1e-12 and abs(qa-g["scaffold_anchor_charge_e_sum"]) < 1e-12
                assert anchorids <= ids, "This report requires anchor-included blocks"
                groups.append({**prefix, "group_id": g["group_id"], "block_type": g["block_type"], "n_groups_in_model": len(m["functional_groups"]), "atom_ids_with_anchor": sorted(ids), "anchor_atom_ids": sorted(anchorids), "block_charge_e": q, "target_formal_charge_e": g["target_formal_charge_e"], "block_minus_formal_e": q-g["target_formal_charge_e"], "anchor_charge_e": qa, "group_without_anchor_charge_e": q-qa, "remaining_molecule_charge_e": observed-q, "max_abs_atomic_charge_in_block_e": max(abs(atoms[a]["charge_e"]) for a in ids), "atomic_charges_e": {a: atoms[a]["charge_e"] for a in sorted(ids)}})
                for radius in range(4):
                    included_scaf = {a for a, dist in gd[g["group_id"]].items() if dist <= radius}
                    region = ids | included_scaf
                    # Include each added scaffold carbon's actual bonded H cap.
                    region |= {b for a in included_scaf for b in neighbors[a] if atoms[b]["element"] == "H"}
                    overlap = [other["group_id"] for other in m["functional_groups"] if other is not g and region.intersection(other["atom_ids"])]
                    charge = sum(atoms[a]["charge_e"] for a in region)
                    neighborhoods.append({**prefix, "group_id": g["group_id"], "block_type": g["block_type"], "n_groups_in_model": len(m["functional_groups"]), "radius_CC_bonds": radius, "charge_e": charge, "minus_group_formal_e": charge-g["target_formal_charge_e"], "remaining_molecule_charge_e": observed-charge, "n_atoms": len(region), "atom_ids": sorted(region), "other_groups_intersected": overlap})
            for p in m["functional_group_pairs"]:
                expected = min(gd[p["group_a_id"]][a] for a in next(g["scaffold_carbon_anchor_ids"] for g in m["functional_groups"] if g["group_id"] == p["group_b_id"]))
                assert p["graph_distance"] == expected
                pairs.append({**prefix, **p})

    paired_models, paired_atoms, paired_groups, mapping = [], [], [], []
    hf = {m["model_id"]: m for m in datasets["HF"]["models"]}
    dft = {m["model_id"]: m for m in datasets["DFT"]["models"]}
    for mid in sorted(hf.keys() & dft.keys()):
        hm, dm = hf[mid], dft[mid]
        identical = signature(originals["HF", mid]) == signature(originals["DFT", mid])
        assert identical, f"Unmatched chemistry or geometry: {mid}"
        assert hm["scaffold_carbon_ids"] == dm["scaffold_carbon_ids"]
        assert [(g["group_id"], g["atom_ids"]) for g in hm["functional_groups"]] == [(g["group_id"], g["atom_ids"]) for g in dm["functional_groups"]]
        mapping.append({"model_id": mid, "chemical_signature_exact_match": identical, "coordinates_exact_match": True, "atom_order_ids_elements_roles_formal_charges_bonds_exact_match": True, "raw_model_sha256_identical": sha(Path(hm["provenance"]["model_path"])) == sha(Path(dm["provenance"]["model_path"]))})
        hbg = {a["atom_id"]: a["background_kind"] for a in hm["background_carbon_distance_charges"]}
        ds = []
        for ha, da in zip(hm["per_atom_charges"], dm["per_atom_charges"], strict=True):
            assert (ha["atom_id"], ha["element"], ha["role"]) == (da["atom_id"], da["element"], da["role"])
            delta = da["charge_e"]-ha["charge_e"]
            ds.append(delta)
            paired_atoms.append({"model_id": mid, "scaffold_C": len(hm["scaffold_carbon_ids"]), "atom_id": ha["atom_id"], "element": ha["element"], "role": ha["role"], "background_kind": hbg.get(ha["atom_id"]), "HF_charge_e": ha["charge_e"], "DFT_charge_e": da["charge_e"], "DFT_minus_HF_e": delta})
        paired_models.append({"model_id": mid, "scaffold_C": len(hm["scaffold_carbon_ids"]), **stats(ds)})
        for hg, dg in zip(hm["functional_groups"], dm["functional_groups"], strict=True):
            paired_groups.append({"model_id": mid, "scaffold_C": len(hm["scaffold_carbon_ids"]), "group_id": hg["group_id"], "block_type": hg["block_type"], "HF_block_charge_e": hg["net_charge_e"], "DFT_block_charge_e": dg["net_charge_e"], "DFT_minus_HF_e": dg["net_charge_e"]-hg["net_charge_e"], "HF_anchor_charge_e": hg["scaffold_anchor_charge_e_sum"], "DFT_anchor_charge_e": dg["scaffold_anchor_charge_e_sum"]})

    group_summary = []
    for method, nscaf, block in sorted({(g["method"],g["scaffold_C"],g["block_type"]) for g in groups}):
        selected = [g for g in groups if (g["method"], g["scaffold_C"], g["block_type"]) == (method,nscaf,block)]
        group_summary.append({"method": method, "scaffold_C": nscaf, "block_type": block, "n_contexts": len({g["model_id"] for g in selected}), "n_blocks": len(selected), "block_min_e": min(g["block_charge_e"] for g in selected), "block_max_e": max(g["block_charge_e"] for g in selected), "mean_abs_deviation_from_formal_e": mean(abs(g["block_minus_formal_e"]) for g in selected), "max_abs_deviation_from_formal_e": max(abs(g["block_minus_formal_e"]) for g in selected)})
    inventory = []
    for method in SOURCES:
        for nscaf in sorted({r["scaffold_C"] for r in context_rows if r["method"] == method}):
            selected = [r for r in context_rows if r["method"] == method and r["scaffold_C"] == nscaf]
            inventory.append({"method": method, "scaffold_C": nscaf, "n_contexts": len(selected), "n_pristine": sum(r["n_groups"] == 0 for r in selected), "n_functionalized": sum(r["n_groups"] > 0 for r in selected), "anchor_loading_min": min(r["anchor_loading_fraction"] for r in selected), "anchor_loading_max": max(r["anchor_loading_fraction"] for r in selected)})
    nonmonotone = []
    for c in context_rows:
        bins = sorted([r for r in shell_rows if r["method"] == c["method"] and r["model_id"] == c["model_id"] and r["background_kind"] == "interior_scaffold" and r["distance_CC_bonds"] is not None], key=lambda r:r["distance_CC_bonds"])
        if len(bins) >= 2:
            nonmonotone.append({"method": c["method"], "model_id": c["model_id"], "n_distance_bins": len(bins), "mean_abs_monotone_nonincreasing": all(b["mean_abs_e"] <= a["mean_abs_e"] for a,b in zip(bins,bins[1:])), "farthest_distance": bins[-1]["distance_CC_bonds"], "farthest_n_atoms": bins[-1]["n"], "farthest_mean_abs_e": bins[-1]["mean_abs_e"], "farthest_max_abs_e": bins[-1]["max_abs_e"]})
    summary = {"schema_version": 1, "analysis": "descriptive_molecular_RESP_assumption_audit", "generated_by": str(Path(__file__).relative_to(ROOT)), "inputs": {k:str(v.relative_to(ROOT)) for k,v in SOURCES.items()}, "inventory": inventory, "context_count": len(context_rows), "matched_HF_DFT_models": len(mapping), "DFT_only_models": sorted(dft.keys()-hf.keys()), "group_summary": group_summary, "distance_monotonicity_descriptive": nonmonotone, "background_summary_by_context": bg_summary, "paired_model_differences": paired_models, "COO_blocks": [g for g in groups if g["block_type"] == "COO"], "all_source_charge_identity_total_checks_revalidated": True, "all_group_sums_and_graph_distances_recomputed": True, "source_artifact_hash_gates": source_gates, "scope": {"new_QM_or_fit": False, "physical_fit_is_unit": True, "rotations_are_independent_samples": False, "HF_ground_state_or_stability_claimed": False, "MD_validation_claimed": False, "transferable_distance_rule_inferred": False, "human_scientific_approval": "pending"}}
    context_summary = []
    for c in context_rows:
        row = dict(c)
        for b in bg_summary:
            if (b["method"], b["model_id"]) == (c["method"], c["model_id"]):
                label = "interior" if b["background_kind"] == "interior_scaffold" else "H_edge"
                row.update({f"{label}_{k}": v for k,v in b.items() if k in stats([])})
        context_summary.append(row)
    outputs = {"summary.csv": context_summary, "contexts.csv": context_rows, "per_atom_charges.csv": atom_rows, "background_carbon_charges.csv": bg_rows, "background_summary_by_context.csv": bg_summary, "background_by_distance.csv": shell_rows, "functional_group_blocks.csv": groups, "functional_group_summary.csv": group_summary, "functional_group_pairs.csv": pairs, "local_neighborhood_charges.csv": neighborhoods, "paired_HF_DFT_models.csv": paired_models, "paired_HF_DFT_atoms.csv": paired_atoms, "paired_HF_DFT_groups.csv": paired_groups, "pairing_validation.csv": mapping, "source_artifact_hash_gates.csv": source_gates, "distance_monotonicity.csv": nonmonotone}
    for name, rows in outputs.items():
        csv_write(name, rows)
    dump("summary.json", summary)
    dump("source_manifest.json", {"algorithm": "SHA-256", "sources": manifest, "memory_used_for_workflow_only": {"file": "MEMORY.md", "lines": "96-110", "rollout_id": "01a07be4-18ac-7111-bcc4-ab20caec3b9d"}})
    # Tables deliberately retain context sizes and denominators rather than pool atoms.
    bgcols = [("method","Metodo"),("model_id","Modello"),("n","n C"),("mean_e","media q / e"),("mean_abs_e","media |q| / e"),("min_e","q min / e"),("max_e","q max / e"),("max_abs_e","max |q| / e"),("n_abs_gt_0p10","n |q| > 0,10 e")]
    interior = [r for r in bg_summary if r["background_kind"] == "interior_scaffold"]
    pristine = [r for r in interior if not contexts[r["method"], r["model_id"]]["functional_groups"]]
    coo = [g for g in groups if g["block_type"] in ("COOH","COO")]
    worst_interior = max((r for r in bg_rows if r["background_kind"]=="interior_scaffold"), key=lambda r: abs(r["charge_e"]))
    nondecaying = sum(not r["mean_abs_monotone_nonincreasing"] for r in nonmonotone)
    max_delta = max(paired_models, key=lambda r:r["max_abs_e"])
    primary = [
        {"question": "I C interni hanno q = 0?", "result": f"No nei fit disponibili: max |q| = {abs(worst_interior['charge_e']):.6f} e ({worst_interior['method']}, {worst_interior['model_id']}, C{worst_interior['atom_id']}, d = {worst_interior['nearest_functional_group_graph_distance']}). Anche i pristini hanno cariche interne non nulle."},
        {"question": "|q| cala sempre con d?", "result": f"No: la media |q| per guscio non è monotona non crescente in {nondecaying}/{len(nonmonotone)} contesti con almeno due distanze interne. Conteggio descrittivo su campioni dipendenti."},
        {"question": "COO completo porta −1 e?", "result": "HF: carbossile −0,766142 + anchor +0,086049 = −0,680093 e. DFT: carbossile −0,760300 + anchor +0,114968 = −0,645332 e. Il restante −0,319908/−0,354667 e è sul resto del modello."},
        {"question": "Quanto cambia il fit HF/DFT?", "result": f"Nei 16 modelli appaiati la media atomica |Δq| varia da {min(r['mean_abs_e'] for r in paired_models):.6f} a {max(r['mean_abs_e'] for r in paired_models):.6f} e per modello; il massimo |Δq| atomico è {max_delta['max_abs_e']:.6f} e ({max_delta['model_id']})."},
    ]
    summary["primary_findings"] = primary
    dump("summary.json", summary)
    sections = ['<h1>Audit delle assunzioni RESP molecolari</h1><p class="lead">Analisi descrittiva dei fit completati disponibili al 5 novembre 2026. Le cariche sono parametri dell’ESP di modelli PAH finiti con H al bordo. Tutti i valori di carica sono in e. Revisione scientifica umana: da completare.</p>',
        table(primary,[("question","Assunzione"),("result","Risultato quantitativo")]),
        '<h2>Campione e limiti del confronto</h2><p>35 contesti di fit: 16 HF e 19 DFT. Ogni contesto usa un conformero fisico e tre orientazioni rigide nello stesso fit; le orientazioni non sono repliche indipendenti. I 16 modelli comuni hanno identità chimica, ordine atomico, connettività, ruoli, classi di simmetria e coordinate esattamente coincidenti nei file locali. Gli hash completi dei model.json differiscono per metadati di workflow; la verifica di appaiamento usa il contenuto chimico. I tre C54 hanno solo DFT.</p>',
        table(inventory,[("method","Metodo"),("scaffold_C","C scaffold"),("n_contexts","n fit"),("n_pristine","pristini"),("n_functionalized","funzionalizzati"),("anchor_loading_min","frazione anchor min"),("anchor_loading_max","frazione anchor max")]),
        '<p>La frazione anchor è il numero di carboni scaffold di ancoraggio distinti diviso per i carboni scaffold. È una descrizione del caricamento, non una densità superficiale sperimentale. C16, C24 e C54 hanno dimensione, bordo, numero di casi e caricamenti diversi: non sono repliche intercambiabili.</p>',
        '<h2>1. Carboni interni: l’ipotesi q = 0</h2><p>Questa sezione comprende solo carboni scaffold non appartenenti ai blocchi funzionali e classificati interior_scaffold. I carboni con H al bordo sono riportati separatamente. Già i riferimenti pristini presentano cariche interne non nulle; lo scostamento nei sistemi funzionalizzati non può essere attribuito automaticamente a un trasferimento indotto dal gruppo.</p>',
        '<h3>Riferimenti pristini</h3>', table(pristine,bgcols),
        '<details><summary>Tutti i contesti: carboni interni</summary>'+table(interior,bgcols)+'</details>',
        '<h3>Distanza nel grafo dello scaffold</h3><p>d è il minimo numero di legami C–C dello scaffold da uno degli anchor del gruppo più vicino. Non è una distanza euclidea e non include O o carboni carbossilici esociclici. Le medie per guscio sono descrittive entro il singolo contesto; gli atomi dello stesso fit non sono osservazioni indipendenti. Nessun test di significatività o legge di decadimento è stimato. d assente indica un modello pristino.</p>',
        table(nonmonotone,[("method","Metodo"),("model_id","Modello"),("n_distance_bins","n distanze"),("mean_abs_monotone_nonincreasing","media |q| non crescente"),("farthest_distance","d massimo"),("farthest_n_atoms","n C a d massimo"),("farthest_mean_abs_e","media |q| lontano / e"),("farthest_max_abs_e","max |q| lontano / e")]),
        '<details><summary>Tutti i gusci: solo carboni interni</summary>'+table([r for r in shell_rows if r["background_kind"]=="interior_scaffold"],[("method","Metodo"),("model_id","Modello"),("distance_CC_bonds","d / legami C–C"),("n","n C"),("mean_e","media q / e"),("mean_abs_e","media |q| / e"),("max_abs_e","max |q| / e")])+'</details>',
        '<h2>2. Neutralità e localizzazione dei blocchi</h2><p>OH comprende C–O–H; EPX comprende i due C e O; CAB è il carbonile C=O; COOH e COO comprendono il carbossile e il C scaffold di ancoraggio. EPH/EPT sono controlli epoxide al bordo e rimangono separati da EPX basale. La carica formale bersaglio è 0 per OH, EPX, CAB e COOH e −1 per COO. Non è un vincolo imposto separatamente ai blocchi dal fit molecolare. Le righe aggregano intervalli fra blocchi dello stesso metodo e della stessa dimensione, mantenendo n contesti e n blocchi.</p>',
        table(group_summary,[("method","Metodo"),("scaffold_C","C scaffold"),("block_type","Blocco"),("n_contexts","n fit"),("n_blocks","n blocchi"),("block_min_e","q min / e"),("block_max_e","q max / e"),("mean_abs_deviation_from_formal_e","media |q − formale| / e"),("max_abs_deviation_from_formal_e","max |q − formale| / e")]),
        '<h3>Carbossile da solo e anchor</h3><p>Per COO il riferimento molecolare è −1 e. La carica residua sul resto della molecola è la somma dei parametri RESP esterni al blocco; non identifica un trasferimento elettronico osservabile o univoco.</p>',
        table(coo,[("method","Metodo"),("model_id","Modello"),("block_type","Tipo"),("group_without_anchor_charge_e","carbossile / e"),("anchor_charge_e","anchor / e"),("block_charge_e","blocco completo / e"),("target_formal_charge_e","formale / e"),("block_minus_formal_e","blocco − formale / e"),("remaining_molecule_charge_e","resto molecola / e")]),
        '<h3>Somme di vicinato per modelli con un solo gruppo</h3><p>Si sommano il blocco completo e i C scaffold entro r legami C–C dagli anchor, includendo gli H direttamente legati ai C aggiunti. Non si aggiungono altri O automaticamente. La lista degli atomi è nel CSV. r = 0 ricostruisce il blocco in questi modelli. È una misura operativa della distribuzione dei parametri, non una definizione unica della localizzazione elettronica. I modelli multigruppo sono disponibili nel CSV con eventuali intersezioni fra blocchi esplicite.</p>',
        table([r for r in neighborhoods if r["n_groups_in_model"]==1],[("method","Metodo"),("model_id","Modello"),("block_type","Blocco"),("radius_CC_bonds","r / legami C–C"),("n_atoms","n atomi"),("charge_e","q vicinato / e"),("minus_group_formal_e","q − formale / e"),("remaining_molecule_charge_e","resto molecola / e")]),
        '<details><summary>Tutti i blocchi e le cariche atomiche interne</summary>'+table(groups,[("method","Metodo"),("model_id","Modello"),("group_id","Gruppo"),("block_charge_e","q blocco / e"),("anchor_charge_e","q anchor / e"),("group_without_anchor_charge_e","q senza anchor / e"),("block_minus_formal_e","q − formale / e"),("atomic_charges_e","atomo: q / e")])+'</details>',
        '<h2>3. Confronto HF/DFT appaiato</h2><p>Δq = q(DFT) − q(HF), per lo stesso modello e lo stesso atomo sulla stessa geometria DFT. HF usa HF/6-31G*; DFT usa wB97XD/6-311++G(d,p); entrambi usano ESP MK e RESP iperbolico Amber. La tabella descrive la sensibilità al metodo del fit entro ogni modello, senza appaiare C54 a C24 o C16. La media firmata vicina a zero è dovuta anche al vincolo di carica molecolare e non implica accordo atomico.</p>',
        table(paired_models,[("model_id","Modello"),("scaffold_C","C scaffold"),("n","n atomi"),("mean_abs_e","media |Δq| / e"),("rms_e","RMS Δq / e"),("max_abs_e","max |Δq| / e")]),
        '<details><summary>Differenze appaiate dei blocchi funzionali</summary>'+table(paired_groups,[("model_id","Modello"),("group_id","Gruppo"),("HF_block_charge_e","HF blocco / e"),("DFT_block_charge_e","DFT blocco / e"),("DFT_minus_HF_e","DFT − HF / e")])+'</details>',
        '<h2>4. Interpretazione scientifica consentita</h2><ul><li>I dati consentono di verificare se i parametri dei carboni interni sono nulli, se i blocchi rispettano una neutralità locale scelta e quanto variano con metodo e contesto.</li><li>Un piccolo totale di gruppo può nascondere compensazioni fra cariche atomiche grandi. Una grande carica RESP su un C interno non dimostra da sola un trasferimento fisico di elettroni: la ripartizione atomica deriva dal modello di fit dell’ESP.</li><li>I carboni poco esposti possono essere scarsamente identificati dall’ESP esterno. Questo audit non calcola il condizionamento del fit, intervalli di identificabilità o una densità elettronica partizionata; i valori grandi non sono validati come osservabili atomiche.</li><li>Le tabelle a distanza e le coppie di gruppi non separano in modo causale distanza, faccia, geometria, caricamento, bordo e dimensione. Non giustificano una regola trasferibile o un raggio universale oltre cui porre q = 0.</li><li>I fit HF mantengono il caveat della stabilità del determinante RHF non validata; non si afferma lo stato fondamentale accettato. La stabilità DFT è riportata dai metadati sorgente e non ricalcolata. Il completamento RESP non implica approvazione del gate QM a monte.</li><li>Nessun nuovo calcolo QM, fit, assegnamento di force field, test MD, proprietà bulk o trasferibilità modulare è stato eseguito. Non è una parametrizzazione pronta per produzione.</li></ul>',
        '<h2>5. Controlli e riproducibilità</h2><p>Sono stati ricalcolati da file locali i totali molecolari, le somme dei blocchi, gli anchor e le distanze del grafo; identità, ordine, cariche atomiche e coordinate dei riepiloghi sono stati confrontati con gli input. Gli hash SHA-256 dei sorgenti sono in source_manifest.json. Il confronto con gli hash incorporati nei FIT_COMPLETE.json è separato in source_artifact_hash_gates.csv: eventuali differenze vengono conservate, senza promuovere i metadati a validazione del metodo.</p>',
        '<p>Unità statistica: un fit molecolare. Soglie 0,05 e 0,10 e nei CSV: descrittive, non criteri fisici di accettazione. Nessuna fonte scientifica esterna viene usata per sostituire i dati locali. Procedura redazionale assistita dalla skill locale scientific-writing v2.3; responsabilità e approvazione scientifica restano umane.</p>',
        '<pre>uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python VS_POLL_HIGH_COVERAGE/.venv/bin/python PARAMETERIZATION/reports/RESP_ASSUMPTIONS_20261105_attempt_001/build_report.py</pre>',
        '<h3>File del rapporto</h3><ul>'+''.join('<li><a href="'+name+'">'+name+'</a></li>' for name in ["summary.json","source_manifest.json",*outputs])+ '</ul>',
        '<details><summary>Gate e scope di ogni contesto</summary>'+table(context_rows,[("method","Metodo"),("model_id","Modello"),("is_diagnostic","diagnostico"),("reported_wavefunction_stability_validated","stabilità riportata"),("fit_status_is_acceptance","fit = accettazione"),("final_parameter_validation","parametri finali validati"),("scope","scope"),("upstream_qm_gate","gate QM")])+'</details>',
        '<details><summary>Carboni al bordo con H: tabella separata</summary>'+table([r for r in bg_summary if r["background_kind"]=="hydrogen_terminated_edge"],bgcols)+'</details>']
    from build_extra_loading import build as build_extra
    appendix, extra_files = build_extra(ROOT, OUT, stats, table, csv_write, dump, sha)
    summary["additional_HF_loading"] = {"n_contexts": 2, "included_in_original_35": False, "paired_HF_DFT": False, "L0_baseline_holdout": False, "raw_molecular_RESP_unprojected": True, "artifacts": extra_files}
    dump("summary.json", summary)
    sections.append(appendix)
    page = '<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Audit assunzioni RESP — 2026-11-05</title><style>body{font:16px/1.5 system-ui,sans-serif;margin:2rem auto;max-width:1500px;padding:0 1.2rem;color:#17212b;background:#fbfcfd}h1{font-size:2rem}h2{margin-top:2.5rem;border-bottom:2px solid #b5c6d5;padding-bottom:.3rem}h3{margin-top:1.6rem}.lead{font-size:1.1rem;max-width:1050px}.scroll{overflow:auto;max-height:650px;border:1px solid #c4d0da;margin:1rem 0}table{border-collapse:collapse;width:100%;font-size:.86rem;background:white}th,td{padding:.5rem .65rem;border-bottom:1px solid #dce3e9;text-align:right;white-space:nowrap}th{background:#e5edf4;position:sticky;top:0}td:nth-child(2),th:nth-child(2){text-align:left}body>.scroll:first-of-type td{white-space:normal;text-align:left;min-width:220px}tbody tr:nth-child(even){background:#f3f6f9}details{margin:1.2rem 0}summary{cursor:pointer;font-weight:600}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#e9eff4;padding:1rem}a{color:#164c82}li{margin:.35rem 0}</style><body>'+''.join(sections)+'</body></html>'
    (OUT / "REPORT.html").write_text(page, encoding="utf-8")
    product_files = [OUT / name for name in ["REPORT.html", "summary.json", "source_manifest.json", *outputs, *extra_files]]
    dump("output_manifest.json", {"outputs": [{"file": p.name, "sha256": sha(p), "bytes": p.stat().st_size} for p in product_files]})
    print(json.dumps({"contexts": len(context_rows), "matched_models": len(mapping), "interior_atom_rows": sum(a["background_kind"]=="interior_scaffold" for a in bg_rows), "group_blocks": len(groups), "hash_gate_mismatches": sum(g["match"] is False for g in source_gates), "output": str(OUT)},ensure_ascii=False))


if __name__ == "__main__":
    main()
