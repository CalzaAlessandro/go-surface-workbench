"""Join frozen completed controls on identical contexts; no fitting or new QM."""
from __future__ import annotations

import csv
import datetime
import hashlib
import html
import json
import pathlib
import statistics

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SNAP = ROOT / 'calculations/BCC_NAGL_20261105_snapshot_001'
NAGL = SNAP / 'calculations/nagl_benchmark_20261105_attempt_001/results'
BCC = SNAP / 'calculations/bcc_family_benchmark_20261105_attempt_001/results'
OLD = ROOT / 'reports/DIRECT_QM_ESP_LITERATURE_20261104_attempt_001'
LIT = ROOT / 'reports/AMBERGO_QM_BENCHMARK_20261104_attempt_001'
INPUTS: dict[str, str] = {}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path, *, csv_rows=False):
    INPUTS[str(path.relative_to(ROOT))] = digest(path)
    if csv_rows:
        with path.open() as stream:
            return list(csv.DictReader(stream))
    return json.loads(path.read_text())


def save(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                                      allow_nan=False) + '\n')


def write_csv(name, rows):
    with (HERE / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


done_n = read(NAGL / 'COMPLETE.json')
done_b = read(BCC / 'COMPLETE.json')
assert done_n['technical_successes'] == 16 and not done_n['failures']
assert len(done_b['records']) == 24 and len(done_b['failures']) == 8
assert all(r['native_BCC_all_bonds_verified'] for r in done_b['records'])
values = {}


def add(mid, ref, method, value):
    key = mid, ref, method
    value = float(value)
    if key in values:
        assert abs(values[key] - value) < 1e-12, key
    values[key] = value


for directory in (OLD, LIT):
    for r in read(directory / 'context_scores.csv', csv_rows=True):
        add(r['model_id'], r['QM_reference'], r['candidate'], r['mean_ESP_RRMS'])
for r in read(NAGL / 'QM_context_scores.csv', csv_rows=True):
    add(r['model_id'], r['QM_reference'], 'NAGL_AM1BCC_ELF10', r['mean_ESP_RRMS'])
for r in read(BCC / 'QM_context_scores.csv', csv_rows=True):
    method = 'ABCG2_integer_charge_control' if r['method'] == 'abcg2' else 'BCC_integer_charge_control'
    add(r['model_id'], r['QM_reference'], method, r['mean_ESP_RRMS'])

seven = sorted({mid for mid, ref, method in values
                if mid.startswith('C24_') and method == 'NAGL_AM1BCC_ELF10'})
five = sorted({mid for mid, ref, method in values
               if mid in seven and method == 'ABCG2_integer_charge_control'})
assert len(seven) == 7 and len(five) == 5
assert all((mid, ref, method) in values for mid in five for ref in ('HF', 'DFT')
           for method in ('BCC_integer_charge_control', 'ABCG2_integer_charge_control'))
methods = [
    ('L0_HF_RESP_C16', 'Libreria L0–RESP, training C16'),
    ('SQE_HF_RESP_C16', 'SQE–RESP, training C16'),
    ('current_library', 'LIB originale'),
    ('BCC_integer_charge_control', 'AM1-BCC classico · controllo Q intero'),
    ('ABCG2_integer_charge_control', 'ABCG2 · controllo Q intero'),
    ('NAGL_AM1BCC_ELF10', 'NAGL 1.0 · surrogato ELF10'),
    ('Espaloma_BCC_ELF10', 'EspalomaCharge · modello provato'),
    ('literature_hiergo_opls', 'HierGO / OPLS'),
    ('Stauffer_ESP_current_caps', 'Stauffer ESP + bordi della LIB'),
    ('Stauffer_SESP_current_caps', 'Stauffer SESP + bordi della LIB'),
    ('AmberGO_current_caps', 'AmberGO + bordi della LIB'),
]
aggregate = []
for subset, mids in [('matched_C24_5', five), ('C24_7', seven)]:
    for method, label in methods:
        if not all((mid, ref, method) in values for mid in mids for ref in ('HF', 'DFT')):
            continue
        aggregate.append({'subset': subset, 'N_contexts': len(mids), 'method': method,
                          'label': label,
                          **{ref + '_median_ESP_RRMS': statistics.median(values[mid, ref, method]
                               for mid in mids) for ref in ('HF', 'DFT')}})
write_csv('aggregate_scores.csv', aggregate)
write_csv('context_scores.csv', [{'model_id': mid, 'QM_reference': ref, 'method': method,
                                 'mean_ESP_RRMS': value}
                                for (mid, ref, method), value in sorted(values.items())])

availability = [
    {'method': 'Amber AM1-BCC', 'status': '12/16 coperti; 5 C24 confrontabili', 'meaning': 'Native correzioni BCC classiche sulle popolazioni AM1 archiviate.'},
    {'method': 'ABCG2 / AM1-BCC-GAFF2', 'status': '12/16 coperti; stessi 5 C24', 'meaning': 'Native correzioni ABCG2 sullo stesso AM1. Tabella identificata da hash; non assunta equivalente a ogni revisione pubblicata.'},
    {'method': 'OpenEye AM1-BCC Sym', 'status': 'Non eseguito: licenza assente', 'meaning': 'Variante nativa Quacpac con simmetrizzazione.'},
    {'method': 'OpenEye AM1-BCC NoSym', 'status': 'Non eseguito: licenza assente', 'meaning': 'Variante nativa senza simmetrizzazione.'},
    {'method': 'OpenEye AM1-BCC SymSPt', 'status': 'Non eseguito: licenza assente', 'meaning': 'Variante nativa single point con simmetrizzazione.'},
    {'method': 'OpenEye AM1-BCC NoSymSPt', 'status': 'Non eseguito: licenza assente', 'meaning': 'Variante nativa single point senza simmetrizzazione.'},
    {'method': 'OpenEye AM1-BCC ELF10', 'status': 'Non eseguito: licenza assente', 'meaning': 'Selezione e media su conformeri; non ricostruito mediante un surrogato ML.'},
    {'method': 'NAGL openff-gnn-am1bcc-1.0.0', 'status': '16/16 controlli tecnici superati', 'meaning': 'Predittore ML di cariche tipo ELF10. Nessun fit ai nostri frammenti.'},
    {'method': 'EspalomaCharge', 'status': 'Controllo già eseguito, riportato', 'meaning': 'Altro surrogato AM1-BCC ELF10; distinto da NAGL e da ABCG2.'},
    {'method': 'Stauffer modified AM1-BCC', 'status': 'Solo valori di motivi pubblicati', 'meaning': 'Procedura YASARA; non eseguita come metodo completo sui nostri frammenti.'},
    {'method': 'AM1-BCC2 / BCC-GAFF', 'status': 'Nome da disambiguare', 'meaning': 'Nessuna riga numerica duplicata senza una distinta implementazione e tabella. -at gaff2 non seleziona ABCG2.'},
]
save('METHOD_AVAILABILITY.json', availability)
save('SUMMARY.json', {
    'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'NAGL_job_id': '1469240', 'BCC_ABCG2_job_id': '1469241',
    'NAGL_technical_successes': 16, 'native_BCC_covered': 12, 'native_ABCG2_covered': 12,
    'matched_subsets': {'C24_7': seven, 'matched_C24_5': five},
    'aggregate_scores': aggregate,
    'OpenEye_license_absent_confirmed_by_user': True,
    'whole_sheet_NAGL_predictions': False,
    'whole_sheet_exclusions': done_n['whole_sheet_exclusions'],
    'production_force_field_validated': False,
    'production_libraries_modified': False,
    'new_QM_or_charge_fit': False,
    'HF_caveat': 'Conventional HF/6-31G* branch; RHF determinant stability caveat retained.',
    'DFT_reference': 'wB97XD/6-311++G(d,p), reference validation still running.',
    'BCC_rounding_policy': 'Raw charges preserved. Integer-control values use explicitly documented minimum-L2 correction for three-decimal AM1 population rounding, not missing-parameter imputation.',
    'scope': 'Same archived AM1 populations, different native corrections; not a new independent full geometry-optimized workflow for each method.',
})


def esc(value):
    return html.escape(str(value))


def score_table(subset):
    rows = [r for r in aggregate if r['subset'] == subset]
    return '<table><thead><tr><th>Cariche</th><th>HF · RRMS</th><th>DFT · RRMS</th></tr></thead><tbody>' + ''.join(
        f'<tr><td>{esc(r["label"])}</td><td>{100*r["HF_median_ESP_RRMS"]:.1f}%</td><td>{100*r["DFT_median_ESP_RRMS"]:.1f}%</td></tr>'
        for r in rows) + '</tbody></table>'


style = '''<style>:root{color-scheme:light;--foreground:#202831;--border:#d6dce4;--muted:#596575;--chart-1:#356f9b;--chart-2:#8b5aa3;--radius:8px}*{box-sizing:border-box}body{font:15px/1.5 system-ui;color:var(--foreground);margin:0}main{max-width:1050px;margin:auto;padding:24px}h1{font-size:25px}h2{font-size:20px;margin-top:30px}p{max-width:90ch}table{border-collapse:collapse;width:100%;margin:14px 0;font-size:14px}th,td{border-bottom:1px solid var(--border);padding:9px 8px;vertical-align:top;text-align:left}td:nth-child(n+2){font-variant-numeric:tabular-nums}small,.muted{color:var(--muted)}a{color:var(--chart-1)}.bars{margin:16px 0;max-width:700px}.barrow{display:grid;grid-template-columns:minmax(140px,1fr) 2fr 58px;align-items:center;gap:10px;margin:9px 0}.track{height:14px;background:#e9edf2}.bar{height:14px;background:var(--chart-1)}code{overflow-wrap:anywhere} @media(max-width:600px){main{padding:12px}table{font-size:12px}th,td{padding:7px 4px}.barrow{grid-template-columns:125px 1fr 47px}}</style>'''
chart_rows = [r for r in aggregate if r['subset'] == 'matched_C24_5' and r['method'] in {
    'L0_HF_RESP_C16', 'SQE_HF_RESP_C16', 'current_library', 'BCC_integer_charge_control',
    'ABCG2_integer_charge_control', 'NAGL_AM1BCC_ELF10', 'Espaloma_BCC_ELF10'}]
chart = '<div class="bars" role="img" aria-label="Errore relativo nel potenziale DFT sui medesimi cinque frammenti C24; valori più bassi indicano un accordo migliore.">' + ''.join(
    f'<div class="barrow"><span>{esc(r["label"])}</span><span class="track"><span class="bar" style="display:block;width:{min(100,100*r["DFT_median_ESP_RRMS"]/2):.2f}%"></span></span><strong>{100*r["DFT_median_ESP_RRMS"]:.1f}%</strong></div>'
    for r in chart_rows) + '</div><small>Scala comune: 0–200% RRMS; nessuna barra tagliata.</small>'
assert all(r['DFT_median_ESP_RRMS'] <= 2 for r in chart_rows)
method_table = '<table><thead><tr><th>Metodo</th><th>Disponibilità</th><th>Significato</th></tr></thead><tbody>' + ''.join(
    f'<tr><td>{esc(r["method"])}</td><td>{esc(r["status"])}</td><td>{esc(r["meaning"])}</td></tr>' for r in availability) + '</tbody></table>'
report = '''<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">''' + style + '''<main>
<h1>Cariche GO: confronto esteso BCC, ABCG2 e NAGL</h1>
<p>5 novembre 2026 · risultati completati sui frammenti, nessuna sostituzione delle librerie di produzione.</p>
<p><strong>Risultato:</strong> NAGL supera i controlli tecnici sui 16 frammenti. BCC e ABCG2 hanno copertura verificata su 12/16 ciascuno. Su questi frammenti, ABCG2 riproduce il potenziale QM meno bene del BCC classico. Questo confronto non misura l'accuratezza delle energie di solvatazione, bersaglio della parametrizzazione ABCG2.</p>
<h2>Stesse strutture per tutti: cinque C24</h2>
<p>Errore relativo quadratico medio del potenziale elettrostatico, espresso in percentuale. Prima si media sulle tre orientazioni rigide di ciascun frammento, poi si prende la mediana dei cinque frammenti. Valori inferiori indicano accordo migliore; non sono errori percentuali delle singole cariche.</p>
''' + score_table('matched_C24_5') + chart + '''
<h2>Set più ampio: sette C24</h2>
<p>BCC e ABCG2 sono esclusi da questa tabella perché non coprono tutti i sette contesti. Non confrontare direttamente una mediana a cinque casi con una a sette. L0 e SQE sono stati calibrati su C16, non sui C24 mostrati; NAGL ed Espaloma sono usati senza riaddestramento.</p>
''' + score_table('C24_7') + '''
<h2>Cosa è stato eseguito</h2>''' + method_table + '''
<p><strong>Copertura:</strong> entrambi i metodi nativi sono stati esclusi per pyrene_dione, C24_single_basal_EPX_control, C24_EPH_edge_control e C24_OH_OH_near_trans. Manca la correzione per il legame tra tipi BCC 12 e 16 di ordine 2. L'assenza viene rifiutata esplicitamente: non viene sostituita con zero.</p>
<p><strong>Protocollo BCC:</strong> sono riutilizzate le stesse popolazioni AM1 archiviate; cambia la tabella e la tipizzazione nativa delle correzioni. Le cariche native grezze sono conservate. Per il punteggio si usa un controllo separato con totale intero, correggendo soltanto il residuo da arrotondamento delle popolazioni AM1 a tre decimali. La correzione è dichiarata per ogni atomo; non è una nuova parametrizzazione.</p>
<p><strong>Protocollo NAGL:</strong> modello openff-gnn-am1bcc-1.0.0, openff-nagl 0.6.1, PyTorch 2.7.1. Controllati mappa atomica, connettività, carica totale, ripetibilità, permutazione degli atomi e rappresentazione di risonanza del carbossilato. Nessuna rinormalizzazione successiva delle cariche. Passare questi controlli non prova trasferibilità al GO.</p>
<p>I fogli interi GOX935 e REGOX980 restano esclusi dall'assegnazione NAGL: il grafo di valenza e lo stato elettronico non sono ancora stabiliti. I conteggi di siti di carbonio irrisolti sono diagnostiche del grafo, non molteplicità fisiche calcolate. Non sono stati aggiunti idrogeni o modificati i fogli per forzare il modello.</p>
<h2>Come leggere il risultato scientifico</h2>
<p>HF/6-31G* è qui il ramo convenzionale già calcolato, con la nota instabilità del determinante RHF da conservare. Il secondo potenziale è ωB97XD/6-311++G(d,p). La scelta del riferimento QM è ancora in validazione. Nessuna delle due colonne identifica una carica atomica unica o dimostra da sola accuratezza in MD.</p>
<p>L0 e SQE rimangono promettenti per l'elettrostatica di questo pannello. NAGL migliora sensibilmente rispetto al modello Espaloma provato, ma non supera la LIB originale sulle mediane di questo pannello. ABCG2 è un controllo importante del diverso obiettivo di fit: proprietà in soluzione e riproduzione del potenziale del gas richiedono verifiche distinte.</p>
<h2>RESP oppure ABCG2 per le small molecules?</h2>
<p>RESP è previsto sia con GAFF sia con GAFF2. ABCG2 è la variante BCC ottimizzata per GAFF2. Nel confronto pubblicato del 2025 su 642 soluti neutri FreeSolv, il RMSE dell'energia libera di idratazione è 0.99 kcal/mol per GAFF2/ABCG2, 1.48 per GAFF2/RESP e 1.71 per GAFF2/AM1-BCC. FreeSolv è anche un bersaglio di ottimizzazione: questi numeri non sono interamente una validazione esterna. Il confronto non stabilisce una superiorità per GO, ioni o GBSA.</p>
<p>Per il prossimo database pollutants è annotato GAFF2 + ABCG2 come candidato, con GAFF2 + RESP come confronto. Nessun cambiamento retroattivo: <a href="../../docs/POLLUTANTS_CHARGE_METHOD_NEXT_PHASE_20261105.md">nota operativa per la prossima fase</a>.</p>
<h2>Confronto con Stauffer e parametri conservati</h2>
<p>Stauffer ESP e SESP non sono RESP. SESP scala di 1.5 le cariche ESP periodiche per il contesto del modello d'acqua. Il nostro test aggiunge bordi finiti della LIB, assenti nella parametrizzazione periodica: è un controllo ibrido dichiarato. L'accordo del modello non scalato con un potenziale del gas non dimostra che sia migliore in acqua.</p>
<p>I sette frcmod originali Hoshin, AmberGO, OLD e attuali sono conservati byte per byte con manifest SHA256 in <a href="../../parameters/comparison_sources/HOSHIN_OLD_CURRENT_20261105/README.md">HOSHIN_OLD_CURRENT_20261105</a>. Il confronto e la scelta dei parametri bonded/LJ rimangono successivi alla fase sulle cariche.</p>
<h2>Fonti e dati verificabili</h2>
<ul><li><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7728379/">He et al. 2020: sviluppo AM1-BCC-GAFF2</a>.</li>
<li><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11948320/">He et al. 2025: estensione ABCG2 e confronti con RESP</a>.</li>
<li><a href="https://github.com/openforcefield/openff-nagl-models/blob/main/docs/models/openff-gnn-am1bcc-1.0.0/index.md">OpenFF: modello NAGL 1.0 e bersaglio ELF10</a>.</li>
<li><a href="https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonConstants/OECharges.html">OpenEye: varianti native AM1-BCC</a>.</li>
<li><a href="https://orubel.mcmaster.ca/publications/2014/Stauffer_JCP_141_044705.pdf">Stauffer et al. 2014</a>.</li>
<li><a href="SUMMARY.json">Sintesi, ambito e limitazioni</a> · <a href="aggregate_scores.csv">mediane CSV</a> · <a href="context_scores.csv">singoli contesti CSV</a> · <a href="INPUT_MANIFEST.json">hash degli input</a> · <a href="METHOD_AVAILABILITY.json">disponibilità dei metodi</a>.</li></ul>
</main></html>'''
(HERE / 'REPORT.html').write_text(report)
small_molecule_rows = [
    ('GAFF2 + RESP', 1.48),
    ('GAFF2 + AM1-BCC classico', 1.71),
    ('GAFF2 + ABCG2', 0.99),
]
small_chart = '''<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<style>*{box-sizing:border-box}html,body{background:transparent;margin:0}body{font:14px/1.5 var(--font-sans,system-ui);color:var(--foreground)}h2{font-size:18px;margin:6px 0}p{margin:8px 0;color:var(--muted-foreground,var(--foreground))}.row{display:grid;grid-template-columns:minmax(145px,1fr) 1.3fr 48px;gap:10px;align-items:center;margin:12px 0}.track{height:14px;background:var(--border)}.fill{height:14px;background:var(--chart-1);display:block}strong{font-variant-numeric:tabular-nums;text-align:right}small{font-size:12px}a{color:var(--primary)}@media(max-width:430px){.row{grid-template-columns:135px 1fr 42px;gap:6px}}</style>
<h2>Small molecules: confronto pubblicato in acqua</h2>
<p>Errore RMSE dell'energia libera di idratazione · kcal/mol · minore è meglio</p>
<div role="img" aria-label="GAFF2 RESP: 1.48; GAFF2 AM1-BCC: 1.71; GAFF2 ABCG2: 0.99 kcal/mol. Scala comune da zero a due.">''' + ''.join(
    f'<div class="row"><span>{esc(label)}</span><span class="track"><span class="fill" style="width:{100*value/2:.1f}%"></span></span><strong>{value:.2f}</strong></div>'
    for label, value in small_molecule_rows) + '''</div>
<small>642 soluti neutri FreeSolv, GAFF2 e acqua TIP3P. Il dataset è anche usato per ottimizzare ABCG2: non è una verifica interamente indipendente. Nessuna conclusione automatica per ioni, GO o GBSA. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11948320/">He et al., JCTC 2025, tabella 2</a>.</small></html>'''
(HERE / 'small_molecules_comparison.html').write_text(small_chart)
save('INPUT_MANIFEST.json', INPUTS)
save('REPORT_MANIFEST.json', {p.name: digest(p) for p in HERE.iterdir()
                             if p.is_file() and p.name != 'REPORT_MANIFEST.json'})
print(json.dumps({'matched_C24_5': five, 'aggregate_scores': aggregate}, indent=2))
