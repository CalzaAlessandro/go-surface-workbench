"""Render verified Espaloma results; charge differences are not FF accuracy."""
import argparse
import base64
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import statistics
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
source, out = args.input.resolve(), args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
root = source/'results'
summary = json.loads((root/'BENCHMARK_COMPLETE.json').read_text())
assert summary['status'] == 'BENCHMARK_COMPLETE' and summary['successful_contexts'] == 16
for name, sha in summary['artifacts_SHA256'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest() == sha
for path in root.iterdir():
    if path.is_file(): shutil.copy2(path, out/path.name)
shutil.copy2(source/'requirements.release.frozen.txt', out/'requirements.release.frozen.txt')
probes = json.loads((root/'CONTROL_PROBES.json').read_text())
cases = [json.loads(p.read_text())['summary'] for p in root.glob('*.json')
         if p.name not in ['BENCHMARK_COMPLETE.json','CONTROL_PROBES.json']]
matched = [r for r in cases if r['native_BCC_primary_available']]
assert len(matched) == 12
keys = ['native_AM1_BCC__current_library','espaloma__native_AM1_BCC','espaloma__current_library']
names = ['BCC native ↔ libreria','Espaloma ↔ BCC native','Espaloma ↔ libreria']
values = [summary['aggregate_common_subset'][k]['median_context_charge_RMS_difference_e'] for k in keys]
with tempfile.TemporaryDirectory(prefix='go_espa_report_') as cache:
    os.environ['MPLCONFIGDIR'] = cache; os.environ['XDG_CACHE_HOME'] = cache
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with plt.rc_context({'font.size':10,'svg.fonttype':'none'}):
        fig, ax = plt.subplots(figsize=(7.6,3.1), layout='constrained')
        ax.barh(names[::-1], values[::-1], color='#0072B2', height=.55)
        ax.set_xlim(0, max(values)*1.2)
        for i, value in enumerate(values[::-1]): ax.text(value+.002, i, f'{value:.3f} e', va='center')
        ax.set_xlabel('Scarto RMS atomico: mediana sui medesimi 12 contesti (e)')
        ax.spines[['top','right']].set_visible(False)
        for suffix in ['png','svg']:fig.savefig(out/('charge_difference.'+suffix),dpi=300,facecolor='white')
        plt.close(fig)
pair_rows=[]
for case in sorted(matched,key=lambda r:r['model_id']):
    pair_rows.append({'model_id':case['model_id'],**{
        k:case['comparisons'][k]['charge_RMS_difference_e'] for k in keys}})
with (out/'matched_context_RMS_differences.csv').open('x',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(pair_rows[0]));writer.writeheader();writer.writerows(pair_rows)
role_rows=list(csv.DictReader((root/'per_context_role_comparison.csv').open()))
selected=[r for r in role_rows if (r['model_id']=='pyrene_COO' and r['role'] in ['COO_C','COO_O','COO_anchor'])
          or (r['model_id']=='pyrene_diol_epoxide' and r['role'] in ['epoxide_O','hydroxyl_O'])]
def number(value):return f'{float(value):+.4f}'.replace('-', '−').replace('.',',')
table_rows=''.join('<tr><td>'+html.escape(r['model_id'])+'</td><td>'+html.escape(r['role'])+'</td>'
                   +''.join('<td>'+number(r[k+'_mean_charge_e'])+'</td>' for k in
                            ['current_library','native_AM1_BCC','espaloma'])+'</tr>' for r in selected)
img='data:image/png;base64,'+base64.b64encode((out/'charge_difference.png').read_bytes()).decode()
max_repeat=max(r['repeat_max_difference_e'] for r in cases)
max_permutation=max(r['permutation_max_difference_e'] for r in cases)
max_q=max(abs(r['Espaloma_charge_sum_e']-r['formal_charge']) for r in cases)
page='''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GO · Controllo EspalomaCharge · 3 novembre 2026</title><style>
:root{--foreground:#182532;--muted-foreground:#526272;--secondary:#f1f4f7;--border:#ccd4dc;--ring:#2866a8;--font-sans:system-ui,sans-serif}body{margin:0;font:15px var(--font-sans);line-height:1.55;color:var(--foreground)}main{max-width:1000px;margin:auto}h2{font-size:1.15rem;margin:1.3em 0 .5em}p{margin:.7em 0}.note{color:var(--muted-foreground);font-size:.85rem}img{display:block;width:100%;max-width:850px;height:auto}a{color:var(--ring)}table{border-collapse:collapse;width:100%;font-size:.85rem;font-variant-numeric:tabular-nums}th,td{padding:7px 9px;border-bottom:1px solid var(--border);text-align:right}th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}thead{background:var(--secondary)}.scroll{overflow:auto}code{font-size:.85em;overflow-wrap:anywhere}details{margin:14px 0}summary{font-weight:600;cursor:pointer}li{margin:8px 0}</style></head><body><main>
<p class="note">3 novembre 2026 · controllo di assegnazione delle cariche · non selezione del force field</p>
<p><strong>EspalomaCharge è stato eseguito su 16 grafi GO finiti.</strong> Il confronto con BCC native usa i medesimi 12 casi con parametri BCC completi; i quattro casi esclusi restano espliciti. HF–RESP resta il riferimento principale del progetto.</p>
<h2>La vicinanza BCC–libreria non implica equivalenza con Espaloma</h2>
<img src="__IMAGE__" alt="Scarto RMS delle cariche: mediana su dodici contesti comuni. BCC e libreria 0,039 e, Espaloma e BCC 0,091 e, Espaloma e libreria 0,092 e.">
<p>La libreria è più vicina alle BCC native di quanto Espaloma lo sia su questo pannello. È una differenza tra vettori di carica, non una misura dell’accuratezza fisica. Espaloma approssima un protocollo AM1-BCC ELF10; il controllo nativo qui proviene da Amber/SQM. Differenze di protocollo e trasferimento di dominio vanno distinti: questo test non stabilisce da solo la causa dello scarto.</p>
<p class="note">Per ogni molecola si calcola lo scarto RMS tra atomi corrispondenti, poi la mediana fra contesti. I contesti hanno lo stesso peso. I casi BCC non coperti non sono imputati con zero.</p>
<h2>Esempi atomici reali dei risultati</h2>
<div class="scroll"><table><thead><tr><th>Contesto</th><th>Ruolo</th><th>Libreria (e)</th><th>BCC native (e)</th><th>Espaloma (e)</th></tr></thead><tbody>__ROLE_TABLE__</tbody></table></div>
<p class="note">Quando un ruolo compare più volte nella molecola, è mostrata la media degli atomi di quel ruolo in quel contesto. I valori individuali sono nel CSV.</p>
<h2>Controllo del totale: adattatore necessario nella release provata</h2>
<p>La release v0.0.8 calcola un valore <code>total_charge</code> nell’API, ma non lo inoltra al modello; il featurizer memorizza <code>formal_charge</code>, mentre il readout usa <code>q_ref</code>. Il probe acetato con l’API non adattata restituisce un totale vicino a zero. Nel benchmark il campo <code>q_ref</code> riceve le cariche formali già presenti nel grafo: il readout pubblicato ne usa la somma. I pesi non sono rifittati e il totale è rispettato senza un riscalamento a posteriori.</p>
<p>Il modello è stato caricato in modalità inferenza; sono registrate ripetizione della predizione e permutazione dell’ordine atomico. Questa è una prova del modello pubblicato con adattatore di ingresso dichiarato, non una certificazione dell’API standard per GO carico.</p>
<ul><li>Massimo scarto sul totale formale nei 16 casi: __Q_ERROR__ e.</li><li>Massima differenza ripetendo l’inferenza: __REPEAT__ e.</li><li>Massima differenza dopo permutazione e ripristino degli IDs: __PERMUTATION__ e.</li></ul>
<h2>Limiti del test</h2>
<p>Il dataset contiene frammenti finiti, non grandi fogli o proteine. Il test conferma fattibilità tecnica e conservazione del totale con l’adattatore; non valida adsorbimento, solvente o parametri LJ. Otto contesti hanno metadati R/S ereditati che non coincidono con quelli ricavati dalle coordinate ottimizzate. Sono registrati senza modificare le sorgenti. I descrittori di questa versione non includono coordinate o stereochimica; il benchmark non costituisce una validazione delle etichette stereochimiche.</p>
<h2>Accoppiamento con GAFF e Lennard-Jones</h2>
<p>Cariche e LJ influenzano insieme idratazione, contatti e adsorbimento. Anche i bonded possono dipendere dal modello elettrostatico con cui sono stati sviluppati. Studiare prima le cariche è utile per isolare le ipotesi; la validazione finale deve mantenere identificata la coppia cariche–LJ e verificare i bonded pertinenti.</p>
<p>In <a href="https://doi.org/10.1021/acs.jpcb.2c06037">GAFF2–AIM</a> i bonded provengono dalla parametrizzazione GAFF2 dei frammenti. Dalla densità DFT ottengono DDEC6 per le cariche e volumi atomici per i LJ tramite la procedura descritta nel SI. Non descrivono un nuovo fit bonded AIM. La derivazione comune di cariche e LJ non equivale a un fit simultaneo di tutti i parametri contro esperimenti.</p>
<p>Vedi la <a href="../../docs/CHARGES_BONDED_LJ_COUPLING.md">nota sul loro accoppiamento</a>. Questo controllo Espaloma non modifica alcun parametro GAFF, LJ o frcmod.</p>
<h2>Cosa richiede la ricalibrazione SQE</h2>
<ol><li>Definire poche famiglie atomiche e di legame che coprano basale, ossidazione, bordo e protonazione.</li><li>Usare i riferimenti molecolari HF–RESP già calcolati; ampliare soprattutto COOH/COO⁻ e i contesti misti. Una sola molecola carbossilato non basta a validare la trasferibilità.</li><li>Selezionare una baseline di carica e fissare inizialmente parte dei parametri. Stimare pochi coefficienti di risposta con regolarizzazione, prima di liberare tutte le famiglie di parametri.</li><li>Confrontare cariche e potenziale prodotto su contesti interi esclusi dal fit, verificando totale, stabilità numerica e dipendenza da dimensione.</li></ol>
<p>Un set SQE adattato alle nostre RESP sarebbe una nuova calibrazione GO. I set originali del lavoro <a href="https://link.springer.com/article/10.1186/s13321-021-00528-w">SQE+qp</a> usano riferimenti NPA. Il software e l’architettura si possono riutilizzare; i coefficienti non diventano adatti al nostro target semplicemente sostituendo le cariche di partenza.</p>
<details><summary>Provenienza e dati</summary>
<p>Job Slurm 1463920, completato in 10 s. Release v0.0.8, commit <code>b32649c8df2029fbd25dd18716fe838e0c62d6ae</code>; i metadati di versione del pacchetto da archivio riportano 0+unknown, perciò il commit e l’hash dei pesi sono l’identificazione riproducibile.</p>
<p>Il branch di sviluppo inizialmente esaminato puntava a un checkpoint v0.1.0 non pubblicato. Quel percorso non ha prodotto predizioni; ambiente e log sono preservati separatamente nel tentativo.</p>
<p>Archivio: <code>espaloma_results_20261103_snapshot_001.tar.gz</code><br>SHA-256: <code>4e71ebb063a17e137ed7763bb0034ca585fc197467d60ca3986674e8fe0c524e</code></p>
<ul><li><a href="BENCHMARK_COMPLETE.json">Esito e hash</a>; <a href="CONTROL_PROBES.json">probe dell’API e dell’adattatore</a>.</li><li><a href="per_atom_charge_comparison.csv">Cariche per atomo</a>; <a href="per_context_role_comparison.csv">cariche per ruolo e contesto</a>; <a href="matched_context_RMS_differences.csv">scarti sui 12 contesti comuni</a>.</li><li><a href="requirements.release.frozen.txt">Dipendenze</a>; <a href="charge_difference.svg">figura vettoriale</a>.</li></ul>
<p>Riferimenti: <a href="https://pubs.acs.org/doi/10.1021/acs.jpca.4c01287">EspalomaCharge</a>; <a href="https://github.com/choderalab/espaloma-charge/tree/v0.0.8">sorgente della release</a>; <a href="https://github.com/choderalab/espaloma-charge/releases/tag/v0.0.8">pesi pubblicati</a>. Originali e installazioni scientifiche precedenti preservati.</p>
</details></main></body></html>'''
for key,val in {'__IMAGE__':img,'__ROLE_TABLE__':table_rows,'__Q_ERROR__':f'{max_q:.3g}',
                '__REPEAT__':f'{max_repeat:.3g}','__PERMUTATION__':f'{max_permutation:.3g}'}.items():page=page.replace(key,val)
(out/'REPORT.html').write_text(page)
(out/'REPORT_MANIFEST.json').write_text(json.dumps({'HTML_SHA256':hashlib.sha256(page.encode()).hexdigest(),
  'source_archive_SHA256':'4e71ebb063a17e137ed7763bb0034ca585fc197467d60ca3986674e8fe0c524e',
  'aggregation':'median of per-context atom-matched RMS differences; same 12 native-BCC-covered contexts',
  'charge_similarity_is_physical_validation':False,'final_parameter_validation':False},indent=2)+'\n')
print(json.dumps({'report':str(out/'REPORT.html'),'successful_contexts':len(cases),'matched_contexts':len(matched),
                  'max_charge_total_error_e':max_q,'repeat_difference_e':max_repeat,
                  'permutation_difference_e':max_permutation},indent=2))
