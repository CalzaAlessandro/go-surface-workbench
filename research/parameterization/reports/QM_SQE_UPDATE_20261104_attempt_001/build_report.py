"""Rebuild the bounded SQE comparison from frozen, matched reference artifacts."""
from __future__ import annotations

import base64
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import statistics
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CALC = ROOT/'calculations'
FIT = CALC/'sqe_resp_fit_results_20261104_snapshot_001/results'
OLD = CALC/'complementary_results_20261103_snapshot_001'
os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='go-sqe-mpl-'))
os.environ.setdefault('XDG_CACHE_HOME', tempfile.mkdtemp(prefix='go-sqe-font-'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def table(path, data):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def main():
    complete = read(FIT/'CALIBRATION_COMPLETE.json')
    assert complete['status'] == 'SQE_RESP_CALIBRATION_COMPLETE'
    for name, checksum in complete['artifacts_SHA256'].items():
        assert sha(FIT/name) == checksum, name
    source = read(CALC/'sqe_resp_fit_20261104_attempt_001/inputs.json')
    assert sha(CALC/'sqe_resp_fit_20261104_attempt_001/inputs.json') == complete['input_SHA256']
    # The L0 training grids must be exactly the HF grids declared by the SQE inputs.
    for entry in read(OLD/'l0_hf_train_C16/inputs.json'):
        name = Path(entry['source_path']).name
        assert entry['sha256'] == source['HF_fit_sources'][entry['model_id']]['artifacts_sha256'][name]
    scores = rows(FIT/'per_context_charge_scores.csv')
    held = sorted(r['model_id'] for r in scores if r['candidate'] == 'default_C16_training'
                  and r['split'] == 'held_out_size_transfer' and r['status'] == 'available')
    assert len(held) == 7
    excluded = [dict(model_id=r['model_id'], reason=r['reason']) for r in scores
                if r['candidate'] == 'default_C16_training' and r['status'] != 'available']
    new = rows(FIT/'HF_ESP_scores.csv')
    old = rows(OLD/'cross_method_ESP_scores.csv')
    lookup = {(r['model_id'], r['orientation'], r['candidate']): r for r in new}
    selected = []
    for row in old:
        if row['model_id'] not in held or row['reference_ESP'] != 'HF/6-31G*':
            continue
        corresponding = {'atomic_hf_on_dft': 'molecular_HF_RESP', 'current_full_block_library': 'current_library'}
        if row['candidate'] in corresponding and row['status'] == 'matched_charge_control':
            other = lookup[row['model_id'], row['orientation'], corresponding[row['candidate']]]
            assert abs(float(row['relative_rmse'])-float(other['HF_ESP_RRMS'])) < 1e-10
        if row['candidate'] == 'l0_hf_train_C16':
            assert row['status'] == 'matched_charge_control' and row['C16_train_C24_holdout'] == 'True'
            selected.append({'model_id': row['model_id'], 'orientation': int(row['orientation']),
                             'candidate': row['candidate'], 'HF_ESP_RRMS': float(row['relative_rmse'])})
    assert len(selected) == 21
    names = ['current_library', 'l0_hf_train_C16', 'default_C16_training', 'molecular_HF_RESP']
    for row in new:
        if row['model_id'] in held:
            selected.append({'model_id': row['model_id'], 'orientation': int(row['orientation']),
                             'candidate': row['candidate'], 'HF_ESP_RRMS': float(row['HF_ESP_RRMS'])})
    per_context = []
    for mid in held:
        for candidate in names+['default_all16_training']:
            values = [r['HF_ESP_RRMS'] for r in selected if r['model_id'] == mid and r['candidate'] == candidate]
            assert len(values) == 3
            per_context.append({'model_id': mid, 'candidate': candidate, 'orientations': 3,
                                'mean_HF_ESP_RRMS': statistics.mean(values)})
    medians = {candidate: statistics.median(r['mean_HF_ESP_RRMS'] for r in per_context
                                          if r['candidate'] == candidate)
               for candidate in names+['default_all16_training']}
    summary = {
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'physical_C24_contexts': held, 'excluded_from_C16_transfer_test': excluded,
        'aggregation': 'mean over three rigid HF MK orientations within each molecular context; median over seven contexts',
        'unit_of_observation': 'molecular context, not ESP rotation',
        'median_HF_ESP_RRMS': medians,
        'all16_fit_C24_evaluation_is_in_sample': True,
        'molecular_HF_RESP_is_fitted_in_sample_reference': True,
        'COOH_COO_independent_transfer_validation': False,
        'production_parameter_validation': False,
        'fit_source_files_SHA256_verified': len(complete['artifacts_SHA256']),
    }
    save(HERE/'SUMMARY.json', summary)
    table(HERE/'matched_orientation_scores.csv', selected)
    table(HERE/'matched_context_scores.csv', per_context)
    for name in ['CALIBRATION_COMPLETE.json', 'per_context_charge_scores.csv', 'default_per_atom_charges.csv',
                 'default_C16_training.fit.json', 'default_all16_training.fit.json']:
        shutil.copy2(FIT/name, HERE/name)
    status = read(HERE/'SCHEDULER_STATUS.json')
    expanded = read(CALC/'expanded_resp_panel_20261104_attempt_001/PREPARATION_COMPLETE.json')
    submission = read(CALC/'expanded_resp_panel_20261104_attempt_001/SUBMISSION.json')
    save(HERE/'EXPANDED_PREPARATION.json', expanded)
    save(HERE/'EXPANDED_SUBMISSION.json', submission)

    labels = ['Libreria\nattuale', 'L0 HF–RESP\ntraining C16', 'SQE–RESP\ntraining C16', 'RESP molecolare\nriferimento fittato']
    colors = ['#666666', '#0072B2', '#D55E00', '#009E73']
    markers = ['o', 's', '^', 'D']
    with plt.rc_context({'font.size': 11, 'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'font.family': 'DejaVu Sans'}):
        fig, ax = plt.subplots(figsize=(7.1, 4.4), layout='constrained')
        maximum = 0
        for index, candidate in enumerate(names):
            values = [r['mean_HF_ESP_RRMS'] for r in per_context if r['candidate'] == candidate]
            maximum = max(maximum, max(values))
            ax.scatter(index+np.linspace(-.12, .12, len(values)), values, color=colors[index],
                       marker=markers[index], s=37, alpha=.8, zorder=3)
            median = medians[candidate]
            ax.plot([index-.25, index+.25], [median, median], color=colors[index], linewidth=3)
            ax.annotate(f'{median:.3f}', (index+.27, median), va='center', fontsize=10)
        ax.set(xticks=range(4), xticklabels=labels, ylabel='HF ESP RRMS (adimensionale)',
               ylim=(0, maximum*1.1), xlim=(-.45, 3.7))
        ax.set_title('Stessi 7 contesti C24 · valori inferiori = residuo minore', fontsize=12, pad=13)
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
        fig.savefig(HERE/'HF_ESP_transfer.png', dpi=300)
        fig.savefig(HERE/'HF_ESP_transfer.svg')
        plt.close(fig)
    image = 'data:image/png;base64,'+base64.b64encode((HERE/'HF_ESP_transfer.png').read_bytes()).decode()
    caption = ('Ogni simbolo è un contesto C24, dopo media sulle tre rotazioni MK; le linee sono mediane, '
               'senza intervalli di confidenza. C24 è escluso dal training L0/SQE; le RESP molecolari '
               'sono invece fittate su quei riferimenti. EPH/EPT esclusi perché assenti nel training C16.')
    css = '''<style>
    :root{font-family:var(--font-sans,system-ui,sans-serif);color:var(--foreground,#202832)}
    body{margin:0}main{max-width:1000px;margin:auto}h1{font-size:1.7rem}h2{font-size:1.2rem;margin-top:1.8rem}
    p,li{line-height:1.55}img,svg{display:block;width:100%;height:auto}.note{color:var(--muted-foreground,#53616d);font-size:.91rem}
    a{color:var(--primary,#005a8c)}.table{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:.92rem}
    td,th{text-align:left;padding:.6rem;border-bottom:1px solid var(--border,#d4dce3)}th{font-weight:650}
    .callout{padding:16px;border-radius:var(--radius,8px);background:var(--secondary,#eef2f6)}
    pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.85rem}code{overflow-wrap:anywhere}
    </style>'''
    svg = (HERE/'HF_ESP_transfer.svg').read_text()
    svg = svg[svg.index('<svg'):]
    visual = f'<!doctype html><html lang="it"><head><meta charset="utf-8">{css}</head><body><main><figure style="margin:0" aria-label="Dispersione e mediana del residuo HF ESP per quattro modelli sugli stessi sette frammenti C24">{svg}</figure><p class="note">{caption}</p></main></body></html>'
    (HERE/'comparison.html').write_text(visual)
    context_table = []
    for mid in held:
        values = {r['candidate']: r['mean_HF_ESP_RRMS'] for r in per_context if r['model_id'] == mid}
        context_table.append('<tr><td>'+html.escape(mid)+'</td>'+''.join(f'<td>{values[c]:.4f}</td>' for c in names)+'</tr>')
    cases_text = ''.join('<li>'+html.escape(r['case'])+': '+
                         (f'iterazione SCF {r["last_SCF_iteration"]}; non ancora convergente' if r['last_SCF_iteration'] is not None and not r['SCF_convergence_marker_present']
                          else 'vedere marcatori e log nella snapshot' if r['log_exists'] else 'non ancora avviato')+'</li>'
                         for r in status['whole_sheet'])
    counts = status['expanded_RESP']['counts']
    document = f'''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">{css}<title>GO · HF–RESP, SQE e QM intero · 4 novembre 2026</title></head><body><main>
    <h1>GO: primo fit SQE, nuovi riferimenti RESP e QM intero</h1>
    <p class="note">Snapshot del cluster: {html.escape(status['verified_at_utc'])}. Risultati, input pronti e job in corso sono distinti. Nessuna nuova libreria finale selezionata.</p>
    <div class="callout"><strong>Il primo SQE migliora la libreria attuale, ma L0 ricalibrato resta migliore sul potenziale HF dei sette contesti C24.</strong> I nuovi riferimenti serviranno a verificare se una correzione dipendente dall’intorno giustifica la maggiore complessità, soprattutto per carbonio basale e carbossilato.</div>
    <h2>Confronto sullo stesso riferimento HF</h2><img src="{image}" alt="Confronto dei residui del potenziale HF, con sette punti per metodo e mediana."><p class="note">{caption}</p>
    <p>Mediane RRMS: libreria attuale {medians['current_library']:.6f}; L0 HF–RESP {medians['l0_hf_train_C16']:.6f}; SQE {medians['default_C16_training']:.6f}; riferimento molecolare fittato {medians['molecular_HF_RESP']:.6f}. Il test non include carbossili o carbossilati indipendenti dal training.</p>
    <p>L0 significa una carica condivisa per ruolo, con i vincoli dei blocchi modulari. SQE usa invece la libreria attuale come q0 e fitta elettronegatività di ruolo alle RESP molecolari: il grafo e la geometria determinano la risposta statica. Carica totale e componenti disconnesse sono conservate; le cariche restano fisse in MD. Non è QM né un calcolo di charge embedding.</p>
    <p>La mediana RMS delle differenze atomiche dalle RESP molecolari passa da {complete['median_C24_baseline_charge_RMS_difference_e']:.6f} e a {complete['median_C24_charge_RMS_difference_e']:.6f} e. È un obiettivo diverso dal residuo del potenziale. Infatti il fit SQE su tutti i 16 contesti dà RRMS C24 {medians['default_all16_training']:.6f}, peggiore del fit C16: questo dato, <strong>in sample</strong>, non è nascosto né usato per scegliere il default.</p>
    <div class="table"><table><thead><tr><th>Contesto</th><th>Attuale</th><th>L0 C16</th><th>SQE C16</th><th>RESP fittata</th></tr></thead><tbody>{''.join(context_table)}</tbody></table></div>
    <h2>Protocollo e limiti del fit</h2>
    <p>Riferimento: HF/6-31G* MK e Amber RESP esterno su geometrie ωB97X-D/6-311++G(d,p) accettate. L’instabilità esterna RHF non è cancellata e non si dichiara una funzione d’onda RHF stabile. Training: sette contesti C16. Test: sette C24 supportati; EPH/EPT esclusi, non riempiti artificialmente. Default predefinito: hardness 4, bond hardness 10, width scale 1, ridge 0.0001; solo elettronegatività di ruolo fittate. Le 11 configurazioni di sensibilità sono esplorative.</p>
    <p>Il solver è verificato contro tre risultati nativi ChargeFW e con test di conservazione, permutazione, orientazione dei legami, rotazioni e risonanza COO. Queste verifiche numeriche non provano trasferibilità fisica né qualità di adsorbimento. I coefficienti pubblicati NPA non sono stati usati come coefficienti RESP.</p>
    <h2>QM sui fogli corrected originali</h2><p>Array <code>1468880</code>: PBE/ma-def2-SVP, RI-J/def2/J, UKS e CPCM acqua, single point sulle coordinate originali. GOX: 1007 atomi, Q = −16 e; REGOX: 913 atomi, Q = −6 e. Doppietto/quartetto per GOX, singoletto/tripletto per REGOX; le prove non esauriscono la ricerca dello stato fondamentale.</p>
    <ul>{cases_text}</ul>
    <p>Il precedente array <code>1468875</code> è andato fuori memoria prima di SCF. Il recupero cambia soltanto risorse: 32 rank MPI, 1024 GB per nodo, maxcore 24000 MB/rank, QoS CPU autorizzata, due casi simultanei e cinque giorni per caso. Nessuna carica finale è ancora derivata da questi fogli. Hirshfeld sarà un diagnostico della densità, non una RESP.</p>
    <h2>Espansione dei frammenti</h2><p>Array <code>1468993</code>: 55 casi, 27 C24 e 28 C54, quattro simultanei e cinque giorni per caso. OH fino a 6; EPX fino a 4; miscele OH/EPX, due carbonili, COOH e COO singoli/doppi, protonazione mista, combinazioni vicine/lontane con OH/EPX e controlli terminali EPH/EPT.</p>
    <p>Vicino: distanza minima di grafo ≤ 2; lontano: ≥ 4 legami C–C tra ancoraggi. Le coppie selezionate e tutte le distanze sono registrate nei modelli. 37 casi sono neutri, 14 hanno Q = −1 e, 4 hanno Q = −2 e. Tutte le geometrie iniziali hanno conteggi esatti, matching π completo e minimizzazione MMFF vincolata convergente; solo DFT potrà validare il minimo e la stabilità.</p>
    <p>Una richiesta C24 con due EPX lontani non è stata ottenuta nel numero limitato di tentativi; resta nel registro delle esclusioni. Nessun gruppo è stato eliminato per farla passare. Alla snapshot: {counts['reference_logs_started']} riferimenti Gaussian avviati, {counts['HF_RESP_complete']} nuove RESP complete, {counts['HF_RESP_failed_or_rejected']} fallimenti/rifiuti della fase HF–RESP. Il percorso è Opt/Freq/Stable DFT → gate → tre MK HF → RESP Amber.</p>
    <h2>Preset e prossimo criterio di scelta</h2><p>GAFF2–AIM e HierGO/OPLS rimangono bundle indivisibili di cariche, bonded, LJ, mixing e convenzioni 1–4. Le loro tabelle non sono convertite in cariche GAFF2 generiche. GAFF2–AIM non copre automaticamente bordi e carbossilati; nella tabella HierGO archiviata manca COO⁻. La mappatura generica completa e la conversione Amber restano da validare. Il catalogo descrive questa politica, non implementa ancora un selettore universale.</p>
    <p>Manterrei L0 come baseline semplice. Con i nuovi riferimenti confronterei A: libreria fissa; B: classi o incrementi conservativi per intorno di grafo; C: SQE, verificando anche un obiettivo diretto sul potenziale con regolarizzazione RESP. Occorre escludere interi contesti dal fit, inclusi carbossilati, distanze e dimensioni; non scegliere le impostazioni sul test. Il confronto con LJ/adsorbimento viene dopo la scelta provvisoria delle cariche.</p>
    <h2>APBS: utile dopo l’assegnazione</h2><p>APBS parte da cariche e raggi assegnati e calcola elettrostatica e risposta di solvente/ioni in un modello continuo. Può mostrare la sensibilità dei modelli al sale e alla solvatazione, ma non fornisce da solo un riferimento indipendente per derivare o validare quelle stesse cariche. Per il nostro obiettivo resta primario HF–RESP; il QM intero controlla effetti elettronici e di scala. Non fittiamo le cariche Amber a un potenziale PB già schermato dal solvente. <a href="https://apbs.readthedocs.io/en/latest/using/index.html">Workflow APBS</a>; <a href="https://apbs.readthedocs.io/en/latest/background.html">modello e limiti</a>.</p>
    <h2>Dati e provenienza</h2><p><a href="SUMMARY.json">Riepilogo numerico</a> · <a href="matched_context_scores.csv">valori per contesto</a> · <a href="default_per_atom_charges.csv">cariche atomiche SQE</a> · <a href="SCHEDULER_STATUS.json">snapshot del cluster</a> · <a href="EXPANDED_SUBMISSION.json">input e lancio dei 55 casi</a> · <a href="REPORT_MANIFEST.json">checksum</a>.</p>
    <p class="note">Gli output scientifici originali, le topologie, le traiettorie e i job di altri progetti restano invariati. Nessun nuovo frcmod o parametro LJ è stato fittato. <a href="https://doi.org/10.1186/s13321-021-00528-w">Riferimento SQE</a>.</p></main></body></html>'''
    (HERE/'REPORT.html').write_text(document)
    files = [p for p in HERE.iterdir() if p.is_file() and p.name != 'REPORT_MANIFEST.json']
    save(HERE/'REPORT_MANIFEST.json', {
        'generated_at_utc': summary['generated_at_utc'],
        'files_SHA256': {p.name: sha(p) for p in sorted(files)},
        'original_source_SHA256': {
            str(FIT/'CALIBRATION_COMPLETE.json'): sha(FIT/'CALIBRATION_COMPLETE.json'),
            str(FIT/'HF_ESP_scores.csv'): sha(FIT/'HF_ESP_scores.csv'),
            str(OLD/'cross_method_ESP_scores.csv'): sha(OLD/'cross_method_ESP_scores.csv'),
            str(OLD/'l0_hf_train_C16/inputs.json'): sha(OLD/'l0_hf_train_C16/inputs.json'),
        },
        'figure': {'width_inches': 7.1, 'height_inches': 4.4, 'PNG_dpi': 300,
                   'matplotlib': matplotlib.__version__, 'numpy': np.__version__,
                   'points': '7 molecular contexts per method; deterministic horizontal displacement',
                   'lines': 'median; no uncertainty interval or smoothing'},
        'production_parameter_validation': False,
    })
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
