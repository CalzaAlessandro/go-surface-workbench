"""Reproduce a matched L0/L1 research comparison from completed native fits."""
from __future__ import annotations

import csv
import hashlib
import html
import json
import pathlib
import statistics

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / 'calculations/environment_resp_results_20261105_snapshot_001/results'
hashes = {}


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p, csv_rows=False):
    hashes[str(p.relative_to(ROOT))] = digest(p)
    if csv_rows:
        with p.open() as stream:
            return list(csv.DictReader(stream))
    return json.loads(p.read_text())


def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                      allow_nan=False) + '\n')


def table(name, rows):
    with (HERE / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


complete = read(DATA / 'COMPLETE.json')
assert complete['status'] == 'RESEARCH_RESP_ENVIRONMENT_FITS_COMPLETE'
for rel, expected in complete['artifacts_sha256'].items():
    assert digest(DATA / rel) == expected, rel
coverage = read(DATA / 'class_coverage.csv', True)
scores = read(DATA / 'context_scores.csv', True)
old = read(ROOT / 'calculations/complementary_results_20261103_snapshot_001/l0_hf_train_C16/shared_library.json')['shared_role_charges_e']
new = read(DATA / 'C16_only_L0/FIT_COMPLETE.json')['shared_role_charges_e']
reproduction = max(abs(old[r] - new[r]) for r in old)
assert reproduction == 0.
matched = []; support = []; class_rows = []
for mode in ('C16_only', 'C16_plus_loading'):
    mids = sorted(r['model_id'] for r in coverage
                  if r['candidate'] == mode + '_L1' and r['split'] == 'held_out'
                  and r['original_panel'] == 'True' and r['fully_supported'] == 'True')
    assert len(mids) == (2 if mode == 'C16_only' else 4)
    for kind in ('L0', 'L1'):
        candidate = mode + '_' + kind
        covered = [r for r in coverage if r['candidate'] == candidate
                   and r['original_panel'] == 'True' and r['split'] == 'held_out'
                   and r['fully_supported'] == 'True']
        support.append({'training': mode, 'method': kind,
                        'supported_original_C24_models': len(covered),
                        'original_C24_models': 9})
        for ref in ('HF', 'DFT'):
            rows = [r for r in scores if r['candidate'] == candidate
                    and r['model_id'] in mids and r['reference'] == ref]
            assert len(rows) == len(mids)
            matched.append({'training': mode, 'method': kind, 'QM_reference': ref,
                            'N_matched_models': len(mids), 'model_ids': '|'.join(mids),
                            'median_ESP_RRMS': statistics.median(float(r['mean_ESP_RRMS']) for r in rows)})
        problem = read(DATA / candidate / 'problem.json')
        fitted = read(DATA / candidate / 'FIT_COMPLETE.json')
        for role, value in fitted['shared_role_charges_e'].items():
            models = sorted({m['model_id'] for m in problem['models']
                             if any(a['role'] == role for a in m['atoms'])})
            count = sum(a['role'] == role for m in problem['models'] for a in m['atoms'])
            class_rows.append({'candidate': candidate, 'class_id': role,
                               'charge_e': value, 'physical_model_support': len(models),
                               'physical_atom_occurrences': count,
                               'model_ids': '|'.join(models)})
table('matched_scores.csv', matched); table('coverage.csv', support)
table('fitted_class_charges.csv', class_rows)
assumptions = read(ROOT / 'reports/RESP_ASSUMPTIONS_20261105_attempt_001/summary.json')
save('SUMMARY.json', {'job_id': '1469243', 'native_Amber_RESP_fits': complete['fits'],
                    'raw_grids_verified': complete['raw_grids_verified'],
                    'previous_L0_reproduction_max_diff_e': reproduction,
                    'coverage': support, 'matched_scores': matched,
                    'L1_outperforms_L0_on_this_matched_panel': False,
                    'reference_or_production_validated': False,
                    'unknown_classes_not_filled': True,
                    'comparison_is_retrospective_not_blind': True,
                    'new_loading_cases_are_not_original_C24_test_cases': True,
                    'size_and_oxidation_acceptance': '../../docs/SCALABILITY_ACCEPTANCE_20261105.md'})


def esc(s):
    return html.escape(str(s))


score_rows = ''.join(f'<tr><td>{esc(r["training"])}</td><td>{r["method"]}</td>'
                     f'<td>{r["N_matched_models"]}</td><td>{r["QM_reference"]}</td>'
                     f'<td>{100*r["median_ESP_RRMS"]:.1f}%</td></tr>' for r in matched)
support_rows = ''.join(f'<tr><td>{esc(r["training"])}</td><td>{r["method"]}</td>'
                       f'<td>{r["supported_original_C24_models"]}/9</td></tr>' for r in support)
fit_rows = ''.join(f'<tr><td>{esc(r["candidate"])}</td><td>{r["training_models"]}</td>'
                   f'<td>{r["fitted_classes"]}</td><td>{r["free_dimensions"]}</td></tr>' for r in complete['fits'])
findings = ''.join(f'<tr><td>{esc(r["question"])}</td><td>{esc(r["result"])}</td></tr>'
                   for r in assumptions['primary_findings'])
chart_rows = [r for r in matched if r['training'] == 'C16_plus_loading' and r['QM_reference'] == 'DFT']
chart = ''.join(f'<div class="barrow"><span>{r["method"]}</span><span class="track">'
                f'<span class="bar" style="width:{100*r["median_ESP_RRMS"]:.2f}%"></span></span>'
                f'<strong>{100*r["median_ESP_RRMS"]:.1f}%</strong></div>' for r in chart_rows)
style = '''<style>:root{--foreground:#202831;--muted:#596575;--border:#d6dce4;--chart-1:#356f9b}*{box-sizing:border-box}html,body{background:transparent}body{font:15px/1.5 system-ui;color:var(--foreground);margin:0}main{max-width:1050px;margin:auto;padding:22px}h1{font-size:25px}h2{font-size:19px;margin-top:28px}table{border-collapse:collapse;width:100%;font-size:14px;margin:12px 0}th,td{border-bottom:1px solid var(--border);padding:8px;text-align:left;vertical-align:top}a{color:var(--chart-1)}.muted,small{color:var(--muted)}.barrow{display:grid;grid-template-columns:50px 1fr 60px;gap:10px;align-items:center;margin:10px 0}.track{height:14px;background:#e9edf2}.bar{display:block;height:14px;background:var(--chart-1)}code{overflow-wrap:anywhere}@media(max-width:600px){main{padding:12px}table{font-size:12px}th,td{padding:7px 4px}}</style>'''
report = '<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">' + style + '''<main>
<h1>GO: primo fit RESP per intorni e audit delle ipotesi</h1>
<p>5 novembre 2026 · quattro fit Amber RESP reali, job 1469243 · nessuna modifica delle librerie di produzione.</p>
<p><strong>Risultato:</strong> questa prima versione L1, dipendente dai primi vicini, non migliora L0 sui pochi contesti comuni di verifica. La copertura è ancora limitata. Non è un confronto del preset GAFF2–AIM completo, né una bocciatura di tutte le possibili descrizioni dell'intorno.</p>
<h2>Cosa abbiamo fittato</h2>
<table><thead><tr><th>Candidato</th><th>Modelli training</th><th>Variabili di carica</th><th>Gradi liberi</th></tr></thead><tbody>''' + fit_rows + '''</tbody></table>
<p>L0 usa ruoli condivisi. L1 distingue l'intorno dei carboni del gruppo attraverso i ruoli dei primi vicini; O/H ereditano il contesto del gruppo. Numerazione, nome del modello, coordinate, orientazione e rappresentazione di risonanza non entrano nella chiave. Gli O del COO restano equivalenti. I test verificano queste invarianti.</p>
<p>Per isolare il contributo delle classi, entrambe le varianti mantengono carbonio basale a zero, blocchi neutri o COO −1 e bordi C–H condivisi. Questa scelta controllata non dimostra che tali vincoli siano fisicamente esatti. Le cariche sono state ottenute dall'eseguibile esterno Amber RESP già verificato; il nuovo L0 riproduce il precedente fit con differenza massima 0.000000 e.</p>
<h2>Copertura delle nove strutture C24 originali escluse dal fit</h2>
<table><thead><tr><th>Training</th><th>Metodo</th><th>Strutture interamente supportate</th></tr></thead><tbody>''' + support_rows + '''</tbody></table>
<p>C16_only usa sette frammenti piccoli. C16_plus_loading aggiunge due nuovi C24 di loading, senza usare i nove C24 originali nel fit. Questo secondo caso verifica strutture escluse dal training, ma non è più un test puro del trasferimento C16→C24. Le classi non viste non ricevono una carica di ripiego.</p>
<h2>Confronto su sottoinsiemi identici</h2>
<table><thead><tr><th>Training</th><th>Metodo</th><th>N comuni</th><th>Potenziale QM</th><th>Mediana RRMS</th></tr></thead><tbody>''' + score_rows + '''</tbody></table>
<p>Ogni valore è la mediana dopo media sulle tre orientazioni di ogni struttura. Non confrontare direttamente le due serie di training: i sottoinsiemi comuni sono diversi, rispettivamente due e quattro strutture. I risultati sono esplorativi e retrospettivi, non un test cieco indipendente di selezione del metodo.</p>
<p>Confronto DFT, training ampliato, stesse quattro strutture; scala 0–100% RRMS:</p>''' + chart + '''
<h2>Cosa dicono le RESP molecolari sui vincoli</h2>
<table><thead><tr><th>Domanda</th><th>Risultato descrittivo</th></tr></thead><tbody>''' + findings + '''</tbody></table>
<p>Le cariche atomiche RESP non sono osservabili uniche. Carboni poco esposti possono essere scarsamente identificati dal potenziale; valori non nulli nei frammenti finiti non bastano a invalidare una libreria con fondo zero. Servono potenziali, momenti, controlli dei vincoli e trasferimento a dimensioni maggiori. La non monotonia osservata non esclude un effetto della distanza; esclude l'assunzione che una sola funzione decrescente sia già dimostrata.</p>
<h2>Dimensione, ossidazione e gruppi mancanti</h2>
<p>La validazione finale deve incrociare fogli di circa 4–10 nm con ossidazione bassa/alta, gruppi isolati/raggruppati e contesti di bordo. I confronti a numero fisso di gruppi su frammenti più grandi cambiano anche la copertura; non dimostrano da soli trasferibilità nella dimensione. <a href="../../docs/SCALABILITY_ACCEPTANCE_20261105.md">Criteri di accettazione e limiti attuali</a>.</p>
<p>Il prototipo SQE usa matrici dense: costo e risposta alla dimensione devono essere verificati prima di proporlo per grandi fogli. Per L0/L1 serve copertura delle nuove classi e una validazione di assegnazione separata dai limiti dell'input RESP. Nessun metodo è già dichiarato validato sui 10 nm ad alta ossidazione.</p>
<p>Gli eteri non epossidici sono documentati nei sorgenti HierGO. Il lattone non ha una parametrizzazione dedicata verificata in quei sorgenti: lattone (estere ciclico) e lattolo (emiacetale ciclico) entrano come riferimenti distinti. <a href="../EXTENDED_GO_GROUPS_20261105_attempt_002/catalog.json">Catalogo dei nuovi riferimenti</a>; nessuna nuova carica o parametro finale è assegnato a questi gruppi.</p>
<h2>Confronto frcmod e calcoli ancora necessari</h2>
<p>La verifica dei coefficienti e delle loro basi è separata dal confronto delle cariche. Un frcmod è un insieme di override e deve essere risolto con base DAT, tipi atomici e ordine di caricamento. Le coincidenze con GAFF/GAFF2 documentano provenienza, non validazione fisica.</p>
<ul><li><a href="../FRCMOD_COMPARISON_20261105_attempt_001/numeric_report.html">Confronto numerico GAFF/GAFF2/Hoshin/OLD/AIM</a>.</li>
<li><a href="../FRCMOD_QM_EVIDENCE_20261105_attempt_001/REPORT.html">Evidenza e priorità dei controlli QM</a>: riutilizzare gli Hessiani validi, verificare prima le topologie, poi introdurre controlli mirati su veri contesti mancanti. Nessun refit globale o job QM aggiuntivo per i frcmod è stato avviato qui.</li></ul>
<h2>Dati e stato della ricerca</h2>
<p>HF/6-31G* è ancora il ramo convenzionale con caveat di stabilità RHF; il controllo DFT usa ωB97XD/6-311++G(d,p). Il riferimento indipendente, i fogli interi e il pannello ampliato restano in validazione. Nessun fit completato viene promosso automaticamente a cariche di produzione.</p>
<ul><li><a href="../RESP_ASSUMPTIONS_20261105_attempt_001/REPORT.html">Audit completo delle ipotesi RESP</a>.</li>
<li><a href="matched_scores.csv">Punteggi appaiati</a> · <a href="coverage.csv">copertura</a> · <a href="fitted_class_charges.csv">cariche e supporto per classe</a>.</li>
<li><a href="SUMMARY.json">Sintesi e ambito</a> · <a href="INPUT_MANIFEST.json">hash degli input</a> · <a href="../../docs/ENVIRONMENT_RESP_PILOT_20261105.md">disegno e registro di esecuzione</a>.</li></ul>
</main></html>'''
(HERE / 'REPORT.html').write_text(report)
compact = '''<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>html,body{background:transparent;margin:0}body{font:14px/1.5 var(--font-sans,system-ui);color:var(--foreground)}h2{font-size:18px;margin:5px 0}p{margin:8px 0}.barrow{display:grid;grid-template-columns:115px 1fr 58px;gap:10px;align-items:center;margin:13px 0}.track{height:14px;background:var(--border)}.bar{display:block;height:14px;background:var(--chart-1)}strong{font-variant-numeric:tabular-nums;text-align:right}small{font-size:12px}</style><h2>Primo confronto RESP: blocchi fissi e intorni</h2><p>Stesse 4 strutture C24 escluse dal fit · errore relativo del potenziale DFT · minore è meglio</p>''' + ''.join(
    f'<div class="barrow"><span>{"L0 · blocchi fissi" if r["method"] == "L0" else "L1 · primi vicini"}</span><span class="track"><span class="bar" style="width:{100*r["median_ESP_RRMS"]:.2f}%"></span></span><strong>{100*r["median_ESP_RRMS"]:.1f}%</strong></div>'
    for r in chart_rows) + '''<small>Scala comune 0–100% RRMS. Training: 7 C16 + 2 nuovi controlli di loading; 3 orientazioni per struttura. L1 supporta soltanto 4 dei 9 C24 originali, L0 ne supporta 7: il confronto usa esclusivamente i 4 comuni. Questo primo L1 non migliora L0 e non rappresenta il preset GAFF2–AIM completo. Trasferibilità 4–10 nm e bassa/alta ossidazione ancora da validare.</small></html>'''
(HERE / 'comparison.html').write_text(compact)
save('INPUT_MANIFEST.json', hashes)
save('REPORT_MANIFEST.json', {p.name: digest(p) for p in HERE.iterdir()
                             if p.is_file() and p.name != 'REPORT_MANIFEST.json'})
print(json.dumps({'fit_count': 4, 'previous_L0_max_diff_e': reproduction,
                  'matched_scores': matched}, indent=2))
