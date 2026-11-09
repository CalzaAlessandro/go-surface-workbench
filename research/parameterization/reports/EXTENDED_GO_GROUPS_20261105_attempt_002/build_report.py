#!/usr/bin/env python3
"""Rebuild this report and 2D graph plate from the exported candidate catalog.

Run with the same existing workspace interpreter; never submits QM or changes
the candidates. To regenerate candidates use the module CLI in a NEW directory.
"""
from __future__ import annotations

import hashlib
import html
import json
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    catalog = json.loads((OUT / "catalog.json").read_text())
    rows, molecules, legends = [], [], []
    for entry in catalog["models"]:
        record = json.loads((OUT / entry["model_path"]).read_text())
        mol = Chem.MolFromSmiles(entry["smiles"])
        rdDepictor.Compute2DCoords(mol)
        molecules.append(mol)
        legends.append(entry["name"]+" | "+entry["motif"])
        rows.append({**entry, "electron_count": record["electron_count"], "atom_count": len(record["atoms"]), "graph_gate": record["metadata"]["motif_graph_gate"]["passed"]})
    svg = Draw.MolsToGridImage(molecules, molsPerRow=3, subImgSize=(350, 240), legends=legends, useSVG=True)
    (OUT / "motifs.svg").write_text(svg)
    readme = """# Nuovi gruppi GO: riferimenti di ricerca non parametrizzati

Sono stati preparati sette riferimenti C/H/O neutri, senza radicali nel grafo e con elettroni pari. Il singoletto chiuso è un’ipotesi iniziale: nessun QM, verifica dello stato fondamentale, fit RESP o test MD è stato eseguito.

| Riferimento | Motivo | Formula | Fonte dello SMILES |
|---|---|---|---|
"""
    for row in rows:
        readme += f"| `{row['model_id']}` | {row['motif']}, anello {row['ring_size'] or 'non richiesto'} | {row['formula']} | [{row['name']}]({row['source']}) |\n"
    readme += """
Il dibenzofurano è un piccolo riferimento aromatico fuso per un ambiente etereo di bordo/difetto; il difeniletere è un controllo di ponte diarile. Ftalide, cumarina e 1,4-naftochinone sono riferimenti aromatici fusi; fenolo è un controllo aromatico monociclico. Oxolan-2-ol è soltanto un controllo **alifatico** della connettività lattolica. Nessuno è un sostituto equivalente di un foglio GO completo o una prova della presenza/abbondanza dei gruppi nei nostri campioni XPS.

Le fonti ChEBI sono state aperte il 5 novembre 2026. Per oxolan-2-ol lo SMILES e InChIKey sono stati recuperati dalla pagina J-GLOBAL/JST indicizzata; l’accesso diretto ha restituito HTTP 503. Questo limite di accesso è conservato nel manifest delle fonti. Le identità grafiche dei sette input sono state controllate anche mediante formula e InChIKey RDKit, prima di assegnare un singolo stereoisomero geometrico dove lo SMILES non specifica la stereochimica.

Il gate separa etere non epossidico, estere ciclico a cinque/sei termini, emiacetale ciclico senza carbonile, fenolo neutro e chinone para coniugato. Nei modelli ogni atomo ha ID stabile, elemento, carica formale intera, ruolo `extended_…`, equivalenza del grafo e coordinate. **Non ci sono cariche parziali.** I nuovi ruoli sono `UNPARAMETERIZED_RESEARCH_REFERENCE`; `modular.partition_model` e il builder non sono dichiarati compatibili.

Questa seconda esportazione ripristina la rappresentazione aromatica RDKit dopo MMFF94s, prima di assegnare ruoli e scrivere JSON/SDF. Il test di round trip della cumarina verifica che non cambino gli ordini di legame fra i due formati. La prima esportazione rimane storica e non va usata come pannello corrente.

I sette `model.json`/`initial.sdf`/`initial.xyz` contengono geometrie ETKDGv3 + MMFF94s convergenti con seme 20261105. Sono esclusivamente geometrie iniziali. `geometry_opt/reference.gjf` prepara Opt/Freq/MK e Stable a wB97XD/6-311++G(d,p). `HF_RESP_template/reference.gjf` prepara un **template** HF/6-31G* MK/Stable con coordinate iniziali: va rigenerato usando una geometria DFT accettata prima di qualunque uso scientifico. Nessuna orientazione, sottomissione o carica è stata generata.

HierGO comprende ETH per l’etere non epossidico, ma non un set esplicito di lattoni o COO−. I lattoni e il lattolo qui sono estensioni di ricerca ispirate alla letteratura, non blocchi già stabiliti di HierGO. Le motivazioni GO e i limiti di attribuzione sono nella revisione locale `docs/MISSING_GROUPS_LITERATURE.md`; l’identità dei piccoli riferimenti non prova l’identità dei siti nel materiale.

Per arrivare a parametri trasferibili occorrono QM accettato, conservazione della connettività e dell’ordine atomico, fit e validazione su contesti esclusi dalla calibrazione, definizione di blocchi/anchor senza sovrapposizioni e bonded/LJ/scaling coerenti. Non si impone automaticamente la neutralità di sottoinsiemi a cariche RESP molecolari future e non si trasferiscono parametri EPX agli eteri o COOH ai lattoni.

La trasferibilità da circa **4 nm a circa 10 nm**, necessaria per proteine/DNA, resta un requisito da verificare. Il metodo finale dovrà conservare le stesse classi di grafo e gli stessi vincoli di carica molecolare/gruppo tra dimensioni. Servono verifiche tenute fuori dalla calibrazione su fogli più grandi, variando bordo/interno, ossidazione bassa e alta a dimensione fissa, dimensioni diverse a composizione comparabile, gruppi vicini e protonazione, con controlli elettrostatici/MD pertinenti. **Non si richiede nuovo QM dell’intero foglio a ogni dimensione produttiva.** Un nuovo ambiente non supportato richiede riferimenti mirati e una nuova verifica della copertura.

I test reali del grafo sono stati scritti e visti fallire prima dell’implementazione. Verifica completa: 265 test passati, di cui 20 nell’estensione. Le cariche, gli input scientifici precedenti e il codice di assegnamento esistente restano fuori dall’estensione.

Per rigenerare il pannello, dalla directory PARAMETERIZATION, usare una nuova directory:

```sh
PYTHONPATH=src uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python -m go_param.extended_fragments --output-dir reports/EXTENDED_GO_GROUPS_20261105_attempt_003 --cpus 8 --mem-gb 6
```

Per ricostruire solo questo README, la tavola 2D e il rapporto:

```sh
uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python reports/EXTENDED_GO_GROUPS_20261105_attempt_002/build_report.py
```

Vedere anche `catalog.json`, `source_manifest.json`, `artifact_manifest.json` e `docs/EXTENDED_GO_GROUPS_20261105.md`.
"""
    (OUT / "README.md").write_text(readme)
    sources = [{"model_id": row["model_id"], "url": row["source"], "input_smiles": row["smiles"], "formula": row["formula"], "input_inchikey": row["inchikey"], "access_date": "2026-11-05", "source_access": "indexed J-GLOBAL content retrieved; direct request HTTP 503" if row["motif"] == "lactol" else "ChEBI primary record opened", "identity_scope": "molecular reference identity, not presence in GO"} for row in rows]
    code = [PROJECT / "src/go_param/extended_fragments.py", PROJECT / "tests/test_extended_fragments.py", PROJECT / "docs/EXTENDED_GO_GROUPS_20261105.md", PROJECT / "src/go_param/qm.py", PROJECT / "src/go_param/reference_qm.py", PROJECT / "src/go_param/resp.py", PROJECT / "src/go_param/protocols.py", PROJECT / "docs/MISSING_GROUPS_LITERATURE.md"]
    manifest = {"sources": sources, "local_files": [{"path": str(p.relative_to(PROJECT)), "sha256": digest(p)} for p in code], "source_identity_human_verification": "pending", "no_claim_of_GO_population": True}
    (OUT / "source_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n")
    cells = []
    for row in rows:
        vals = [row["model_id"], row["motif"], row["ring_size"] or "—", row["formula"], row["atom_count"], row["electron_count"], row["reference_scope"]]
        cells.append("<tr>"+"".join("<td>"+html.escape(str(v))+"</td>" for v in vals)+"</tr>")
    paragraphs = []
    for block in readme.split("\n\n"):
        if block.startswith("#") or block.startswith("|"):
            continue
        if block.startswith("```sh"):
            paragraphs.append("<pre><code>"+html.escape(block.removeprefix("```sh\n").removesuffix("\n```"))+"</code></pre>")
            continue
        text = html.escape(block)
        text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
        text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
        text = re.sub(r"\[([^]]+)\]\((https://[^)]+)\)", r'<a href="\2">\1</a>', text)
        paragraphs.append("<p>"+text+"</p>")
    source_links = '<ul>'+''.join('<li><a href="'+html.escape(row['source'])+'">'+html.escape(row['name'])+'</a>: <code>'+html.escape(row['smiles'])+'</code></li>' for row in rows)+'</ul>'
    page = '<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Nuovi gruppi GO — riferimenti</title><style>body{font:16px/1.55 system-ui,sans-serif;max-width:1300px;margin:2rem auto;padding:0 1rem;color:#17212b}table{border-collapse:collapse;font-size:.84rem;width:100%}td,th{border-bottom:1px solid #c9d4de;text-align:left;padding:.5rem}th{background:#e7eff5}.scroll{overflow:auto}img{width:100%;max-width:1050px}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f5f8;padding:1rem}a{color:#155d91}</style><body><h1>Sette nuovi riferimenti GO non parametrizzati</h1><p>Grafi sanitizzati e geometrie iniziali. Nessun risultato QM, carica RESP o parametro produttivo. La trasferibilità da circa 4 nm a circa 10 nm resta da validare.</p><figure><img src="motifs.svg" alt="Sette grafi di riferimento distinti: due eteri, due lattoni, un lattolo alifatico, fenolo e naftochinone"><figcaption>Connettività 2D dai sette SMILES documentati; non rappresenta geometrie QM né densità dei gruppi nei fogli.</figcaption></figure><div class="scroll"><table><thead><tr>'+''.join('<th>'+x+'</th>' for x in ['Modello','Motivo','Anello','Formula','Atomi con H','Elettroni','Scope'])+'</tr></thead><tbody>'+''.join(cells)+'</tbody></table></div><h2>Protocollo, limiti e riproducibilità</h2>'+''.join(paragraphs)+'<h2>SMILES documentati</h2>'+source_links+'<p><a href="catalog.json">Catalogo</a> · <a href="source_manifest.json">Fonti</a> · <a href="artifact_manifest.json">Hash degli artefatti</a> · <a href="README.md">README</a></p></body></html>'
    (OUT / "REPORT.html").write_text(page)
    files = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "artifact_manifest.json" and "__pycache__" not in p.parts)
    (OUT / "artifact_manifest.json").write_text(json.dumps({"algorithm": "SHA-256", "files": [{"path": str(p.relative_to(OUT)), "sha256": digest(p), "bytes": p.stat().st_size} for p in files]}, indent=2)+"\n")
    print(json.dumps({"models": len(rows), "artifact_files": len(files), "new_QM_or_fit": False}))


if __name__ == "__main__":
    main()
