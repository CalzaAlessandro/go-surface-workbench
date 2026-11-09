"""Separate descriptive appendix for two newly recovered HF loading fits.

Uses the existing local context parser and the workspace's existing NumPy.
Never projects raw charges onto the L0 library or changes fit inputs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def build(root, out, stats, table, csv_write, dump, sha):
    sys.dont_write_bytecode = True
    code = root / "PARAMETERIZATION/src"
    sys.path.insert(0, str(code))
    from go_param.context import analyze_dataset

    names = ("C24_EPX3_loading_s20261308", "C24_OH2_loading_s20261304")
    archive = root / "PARAMETERIZATION/calculations/expanded_resp_completed_20261105_snapshot_001/cases"
    specs = []
    paths = []
    for name in names:
        fit = archive / name / "HF_RESP_attempt_001/fits/attempt_001"
        spec = {"model_path": str(fit / "model.json"), "charges_path": str(fit / "charges.json"), "fit_evidence_path": str(fit / "FIT_COMPLETE.json")}
        assert all(Path(p).is_file() for p in spec.values())
        specs.append(spec)
        paths.extend(Path(p) for p in spec.values())
    context = analyze_dataset(specs)
    assert context["model_count"] == context["context_count"] == 2
    assert context["coverage"]["all_models_passed_charge_identity_and_total_checks"]
    rows, blocks, shells = [], [], []
    gates = []
    for m in context["models"]:
        atoms = {a["atom_id"]: a for a in m["per_atom_charges"]}
        prefix = {"method": "HF_extra_loading", "model_id": m["model_id"], "scaffold_C": len(m["scaffold_carbon_ids"])}
        anchors = set().union(*(set(g["scaffold_carbon_anchor_ids"]) for g in m["functional_groups"]))
        evidence = json.loads(Path(m["provenance"]["fit_evidence_path"]).read_text())
        for label in ("model", "charges"):
            p = Path(m["provenance"][f"{label}_path"])
            expected = evidence["artifacts_sha256"][p.name]
            assert sha(p) == expected
            gates.append({"model_id": m["model_id"], "artifact": p.name, "embedded_sha256": expected, "actual_sha256": sha(p), "match": True})
        for kind in ("interior_scaffold", "hydrogen_terminated_edge"):
            background = [a for a in m["background_carbon_distance_charges"] if a["background_kind"] == kind]
            rows.append({**prefix, "n_groups": len(m["functional_groups"]), "anchor_loading_fraction": len(anchors)/len(m["scaffold_carbon_ids"]), "background_kind": kind, **stats([a["charge_e"] for a in background]), "physical_conformers": m["provenance"]["physical_conformers"], "orientations": m["provenance"]["orientations"], "scope": m["provenance"]["scope"], "reported_wavefunction_stability_validated": m["provenance"]["reported_wavefunction_stability_validated"]})
            for distance in sorted({a["nearest_functional_group_graph_distance"] for a in background}):
                shells.append({**prefix, "background_kind": kind, "distance_CC_bonds": distance, **stats([a["charge_e"] for a in background if a["nearest_functional_group_graph_distance"] == distance])})
        for g in m["functional_groups"]:
            q = sum(atoms[a]["charge_e"] for a in g["atom_ids"])
            assert abs(q-g["net_charge_e"]) < 1e-12
            blocks.append({**prefix, **g})
    files = ["extra_loading_context_analysis.json", "extra_loading_summary.csv", "extra_loading_blocks.csv", "extra_loading_background_by_distance.csv", "extra_loading_source_manifest.json"]
    dump(files[0], context)
    csv_write(files[1], rows)
    csv_write(files[2], blocks)
    csv_write(files[3], shells)
    paths += [Path(__file__).resolve(), *(code / "go_param" / f for f in ("context.py", "modular.py", "resp.py", "protocols.py"))]
    dump(files[4], {"separate_from_original_35_contexts": True, "paired_HF_DFT": False, "L0_baseline_holdout": False, "raw_molecular_RESP_unprojected": True, "source_artifact_hash_gates": gates, "sources": [{"path": str(p.relative_to(root)), "sha256": sha(p), "bytes": p.stat().st_size} for p in paths]})
    markup = '<h2>Appendice: due nuovi fit HF di caricamento</h2><p>Questi due fit sono stati archiviati separatamente il 5 novembre 2026. Sono una raccolta aggiuntiva, esclusa dai 35 contesti e dai 16 confronti HF/DFT delle sezioni principali. Non hanno qui una controparte DFT e non costituiscono un holdout della baseline L0. Sono cariche molecolari RESP grezze: non è imposto q = 0 ai carboni interni e non è applicata la neutralizzazione dei blocchi della libreria L0. Il caveat sul determinante RHF resta invariato. Ogni caso usa un conformero e tre orientazioni nello stesso fit.</p>'
    markup += table([r for r in rows if r["background_kind"] == "interior_scaffold"],[("model_id","Modello"),("n_groups","n gruppi"),("anchor_loading_fraction","frazione anchor"),("n","n C interni"),("mean_abs_e","media |q| / e"),("max_abs_e","max |q| / e")])
    markup += table(blocks,[("model_id","Modello"),("group_id","Gruppo"),("net_charge_e","q blocco / e"),("target_formal_charge_e","formale / e"),("scaffold_anchor_charge_e_sum","q anchor / e")])
    markup += '<p>I quattro hash model/charges corrispondono ai manifest dei fit. I CSV separano anche i C al bordo con H e i gusci di distanza. Questi due modelli aggiungono contesti di caricamento; non isolano un effetto causale del numero di gruppi.</p><ul>'+''.join('<li><a href="'+name+'">'+name+'</a></li>' for name in files)+'</ul>'
    return markup, files
