"""Summarize verified, explicitly diagnostic whole-sheet charge predictions."""
from collections import defaultdict
import base64
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CALC=HERE.parent
OUT=ROOT/'reports/CORRECTED_SHEETS_CHARGE_CONTROL_20261103_attempt_001'
OUT.mkdir(exist_ok=False)
FIRST=CALC/'corrected_sheets_charge_results_20261103_snapshot_001'
SECOND=CALC/'corrected_sheets_resonance_results_20261103_snapshot_001'

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def table(p,rows):
    with p.open('x',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

expected={'corrected_sheets_charge_results_20261103_snapshot_001.tar.gz':'0327903353e74a287ba689aab5299f7b3ad68c291d892908fd5116c76a393ed9',
          'corrected_sheets_resonance_results_20261103_snapshot_001.tar.gz':'e59ddf633b1e865b60eb8b74e92c3aab10cacb2b327fe7db98121b40f03dbc92'}
for name,checksum in expected.items():assert sha(CALC/name)==checksum
inputs=read(FIRST/'inputs.json');atoms=[];roles=[];records=[];verified=0
for source in [FIRST,SECOND]:
    for folder in source.glob('results_*'):
        d=read(folder/'BENCHMARK_COMPLETE.json')
        assert d['status']=='DIAGNOSTIC_PREDICTIONS_COMPLETE' and not d['failures']
        assert all(sha(folder/k)==v for k,v in d['artifacts_SHA256'].items())
        verified+=len(d['artifacts_SHA256']);records.extend(d['records'])
        atoms.extend(list(csv.DictReader((folder/'per_atom_charges.csv').open())))
        roles.extend(list(csv.DictReader((folder/'per_context_role_charges.csv').open())))
assert len(records)==18
table(OUT/'per_atom_all_methods.csv',atoms);table(OUT/'per_role_all_methods.csv',roles)
table(OUT/'prediction_timings_and_checks.csv',records)
save(OUT/'GRAPH_QC.json',[m['graph_QC'] for m in inputs['models']])
shutil.copy2(HERE/'INDEPENDENT_GRAPH_AND_RESONANCE_CHECKS.json',OUT/'INDEPENDENT_GRAPH_AND_RESONANCE_CHECKS.json')
shutil.copy2(SECOND/'SCHEDULER_RESONANCE_STATUS.json',OUT/'SCHEDULER_STATUS.json')
for m in inputs['models']:
    target=OUT/'source_and_diagnostic_graphs'/m['name'];target.mkdir(parents=True)
    for p in (FIRST/'sources'/m['name']).iterdir():shutil.copy2(p,target/p.name)

lookup={(r['model'],r['method'],r['variant'],r['atom_id']):float(r['charge_e']) for r in atoms}
methods=['espaloma','sqeq0','sqeqp'];labels=['Espaloma','SQE+q0','SQE+qp']
numbers=[];sensitivity=[];widget=[]
for model in inputs['models']:
    name=model['name'];n=len(model['atoms'])
    for method in methods:
        a=np.array([lookup[name,method,'A',str(i+1)] for i in range(n)])
        for alt,label in [('B','Alternative maximum-valence assignment'),('C_COO_resonance_swap','COO resonance-only swap')]:
            other=np.array([lookup[name,method,alt,str(i+1)] for i in range(n)])
            diff=a-other
            sensitivity.append({'model':name,'method':method,'test':label,'RMS_difference_e':float(np.sqrt(np.mean(diff**2))),
                                'max_atomic_difference_e':float(np.max(abs(diff)))})
        basal=np.array([a[i] for i,at in enumerate(model['atoms']) if at['role']=='scaffold_C'])
        res=defaultdict(float)
        for at,value in zip(model['atoms'],a):
            if at['resname']=='CBD':res[at['resid']]+=float(value)
        numbers.append({'model':name,'method':method,'variant':'A','basal_atoms':len(basal),
                        'basal_mean_e':float(basal.mean()),'basal_min_e':float(basal.min()),'basal_max_e':float(basal.max()),
                        'basal_fraction_abs_q_above_0_1e':float(np.mean(abs(basal)>.1)),
                        'mean_CBD_residue_charge_e':float(np.mean(list(res.values()))),
                        'minimum_CBD_residue_charge_e':min(res.values()),'maximum_CBD_residue_charge_e':max(res.values())})
    groups=defaultdict(list)
    for i,at in enumerate(model['atoms']):groups[at['role']].append(i)
    for variant in ['A','B','C_COO_resonance_swap']:
        for role,indices in groups.items():
            row={'model':name,'variant':variant,'role':role,'count':len(indices),
                 'source':float(np.mean([model['atoms'][i]['source_partial_charge_e'] for i in indices]))}
            for method in methods:
                values=np.array([lookup[name,method,variant,str(i+1)] for i in indices])
                row[method]=float(values.mean());row[method+'_q05']=float(np.quantile(values,.05));row[method+'_q95']=float(np.quantile(values,.95))
            widget.append(row)
table(OUT/'basal_and_CBD_summary.csv',numbers);table(OUT/'representation_sensitivity.csv',sensitivity)
table(OUT/'mean_charge_table.csv',widget)
with tempfile.TemporaryDirectory(prefix='go_sheets_plot_') as cache:
    os.environ['MPLCONFIGDIR']=cache;os.environ['XDG_CACHE_HOME']=cache
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=['#0072B2','#009E73','#D55E00']
    with plt.rc_context({'font.size':10,'svg.fonttype':'none'}):
        fig,axes=plt.subplots(1,2,figsize=(9,3.6),sharex=True,sharey=True,layout='constrained')
        for ax,model in zip(axes,inputs['models']):
            ids=[at['id'] for at in model['atoms'] if at['role']=='scaffold_C']
            for method,label,color in zip(methods,labels,colors):
                values=[lookup[model['name'],method,'A',str(i)] for i in ids]
                ax.hist(values,bins=np.arange(-.6,.601,.03),density=True,histtype='step',linewidth=1.8,label=label,color=color)
            ax.axvline(0,color='#555555',ls='--',lw=1,label='Libreria: 0 e')
            ax.set_title(model['name']+f' · {len(ids)} carboni GGG')
            ax.set_xlabel('Carica atomica (e)');ax.set_xlim(-.6,.6);ax.spines[['top','right']].set_visible(False)
        axes[0].set_ylabel('Densità di probabilità (1/e)')
        axes[0].legend(frameon=False,fontsize=8)
        for suffix in ['png','svg']:fig.savefig(OUT/('basal_distributions.'+suffix),dpi=300,facecolor='white')
        plt.close(fig)
        fig,ax=plt.subplots(figsize=(8,3.4),layout='constrained')
        for j,(method,label,color) in enumerate(zip(methods,labels,colors)):
            values=[r['max_atomic_difference_e'] for name in ['GOX935','REGOX980'] for r in sensitivity
                    if r['model']==name and r['method']==method and r['test']=='COO resonance-only swap']
            y=np.arange(2)+(j-1)*.2;ax.barh(y,values,height=.18,color=color,label=label)
            for yy,v in zip(y,values):ax.text(v+.01,yy,f'{v:.3f} e',va='center',fontsize=9)
        ax.set_yticks([0,1],['GOX935','REGOX980']);ax.invert_yaxis();ax.set_xlim(0,.76)
        ax.set_xlabel('Massima variazione atomica scambiando la sola risonanza COO⁻ (e)')
        ax.legend(frameon=False,ncols=3,loc='lower right',fontsize=9);ax.spines[['top','right']].set_visible(False)
        for suffix in ['png','svg']:fig.savefig(OUT/('COO_resonance_sensitivity.'+suffix),dpi=300,facecolor='white')
        plt.close(fig)

def image(name):return 'data:image/png;base64,'+base64.b64encode((OUT/name).read_bytes()).decode()
css='''<style>:root{--foreground:#192936;--muted-foreground:#536270;--secondary:#f0f4f7;--border:#ccd5dc;--ring:#2866a8;--font-sans:system-ui,sans-serif}body{margin:0;font:15px var(--font-sans);line-height:1.55;color:var(--foreground)}main{max-width:1000px;margin:auto}h2{font-size:1.15rem;margin:1.3em 0 .4em}img{display:block;width:100%;height:auto}a{color:var(--ring)}.note{font-size:.85rem;color:var(--muted-foreground)}.scroll{overflow:auto}table{border-collapse:collapse;font-size:.85rem;font-variant-numeric:tabular-nums;width:100%}th,td{padding:7px;border-bottom:1px solid var(--border);text-align:right}th:first-child,td:first-child{text-align:left}select{font:inherit;margin:5px 14px 5px 0;padding:5px;color:var(--foreground);background:var(--secondary);border:1px solid var(--border);border-radius:4px}code{overflow-wrap:anywhere}</style>'''
page='''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GOX e REGOX corrected · prove diagnostiche di carica</title>'''+css+'''</head><body><main>
<h2>GOX e REGOX corrected: prove completate, cariche ancora diagnostiche</h2>
<p>Espaloma, SQE+q0 e SQE+qp sono stati eseguiti sui due fogli originali corrected del cluster, estratti dalle topologie dei complessi con le maschere <code>:1-494</code> e <code>:1-628</code>. Atomi, geometrie, connettività e cariche di produzione non sono stati modificati. Tutte le 18 predizioni tecniche sono terminate e conservano la carica netta.</p>
<div class="scroll"><table><tr><th>Struttura</th><th>Atomi</th><th>Formula</th><th>Carica</th><th>Residui CBD</th><th>Elettroni</th></tr><tr><td>GOX935 corrected</td><td>1007</td><td>C₆₁₃H₁₄₁O₂₅₃</td><td>−16 e</td><td>16</td><td>5859</td></tr><tr><td>REGOX980 corrected</td><td>913</td><td>C₆₄₇H₁₄₆O₁₂₀</td><td>−6 e</td><td>6</td><td>4994</td></tr></table></div>
<h2>Prima delle cariche: il grafo elettronico richiede verifica</h2>
<p>Amber memorizza i legami, ma non un grafo di Lewis completo con tutti gli ordini di legame. Ho mantenuto i gruppi della libreria, un COO⁻ per ogni CBD e carboni formalmente neutri, cercando il massimo numero possibile di doppi legami C=C compatibili con le valenze.</p>
<p>Restano almeno <strong>63 siti carboniosi non completati in GOX e 28 in REGOX</strong> sotto queste ipotesi. Il risultato è verificato sia con un’ottimizzazione intera sia con un algoritmo indipendente di matching bipartito. GOX ha inoltre un numero dispari di elettroni: con formula e carica attuali non può essere un sistema interamente a guscio chiuso.</p>
<p>Questi conteggi non determinano il numero fisico di spin spaiati, la molteplicità o il numero di idrogeni da aggiungere. Una diversa descrizione elettronica può includere delocalizzazione, separazione di carica o stati a guscio aperto. Per le predizioni ho usato due completamenti massimali A e B, rendendo espliciti i siti residui nel grafo diagnostico. Nessuna passivazione o correzione chimica è stata applicata agli originali.</p>
<h2>Cariche per gruppo: medie e intervalli nel foglio</h2>
<p>Selezionare il foglio e la rappresentazione. A/B cambiano la collocazione dei doppi legami e dei siti carboniosi residui; C conserva A e scambia soltanto i due ossigeni nella rappresentazione di risonanza di ogni COO⁻. Le celle riportano la media; il tooltip mostra il 5°–95° percentile. Questi valori non sono nuove cariche universali di LIB.</p>
<label>Foglio <select id="model" onchange="update()"><option>GOX935</option><option>REGOX980</option></select></label><label>Rappresentazione <select id="variant" onchange="update()"><option value="A">A · completamento primario</option><option value="B">B · completamento alternativo</option><option value="C_COO_resonance_swap">C · sola risonanza COO⁻</option></select></label>
<div class="scroll"><table><thead><tr><th>Ruolo</th><th>N</th><th>Libreria attuale</th><th>Espaloma</th><th>SQE+q0</th><th>SQE+qp</th></tr></thead><tbody id="charges"></tbody></table></div>
<h2>Carboni basali: una media quasi zero può nascondere distribuzioni ampie</h2>
<img src="'''+image('basal_distributions.png')+'''" alt="Distribuzioni delle cariche sui carboni GGG nei due fogli, per i tre metodi diagnostici.">
<p>Con Espaloma il carbonio GGG medio vale −0,004 e in GOX e −0,040 e in REGOX, ma circa il 56% e il 55% dei siti supera |q| = 0,1 e. Gli intervalli A sono rispettivamente −0,414…+0,358 e e −0,380…+0,539 e. SQE produce distribuzioni più strette. Poiché il grafo elettronico non è validato, questi dati non bastano a confutare l’approssimazione della tua LIB con carbonio basale zero: mostrano dove confrontare riferimenti QM ben definiti.</p>
<p>La carica media dell’intero residuo CBD, incluso il carbonio di ancoraggio, è circa −0,717/−0,430 e con Espaloma, −1,027/−1,044 e con SQE+q0 e −0,106/−0,110 e con SQE+qp (GOX/REGOX). Il resto della carica è distribuito sul foglio. La conservazione di −16 e o −6 e da sola non garantisce una localizzazione fisicamente corretta.</p>
<h2>Controllo cruciale: la sola risonanza del carbossilato</h2>
<img src="'''+image('COO_resonance_sensitivity.png')+'''" alt="Variazione massima delle cariche quando si scambiano gli ordini di legame e la carica formale tra ossigeni COO, mantenendo immutata la struttura fisica.">
<p>Scambiando soltanto O=C–O⁻ con ⁻O–C=O, con gli stessi atomi, coordinate e carica totale, SQE+q0 cambia fino a circa <strong>0,65 e</strong> su un ossigeno; SQE+qp fino a 0,24 e; Espaloma fino a 0,043 e. La variante C non sposta i siti carboniosi residui: separa quindi questo effetto dall’ambiguità A/B del resto del foglio. SQE+q0 assegna ai due ossigeni anche valori vicini a −1,23 e e −0,58 e, nonostante una media di circa −0,90 e.</p>
<p>Questo controllo rende Espaloma più interessante degli attuali set SQE per la robustezza alla notazione del COO⁻, senza provarne l’accuratezza fisica. Per la futura calibrazione occorre trattare simmetria e risonanza esplicitamente: tipizzazione coerente, eventuale media di rappresentazioni validate e confronto con RESP, includendo test su molecole intere escluse dal fit.</p>
<h2>Tempi, metodi e verifiche</h2>
<p>Le singole inferenze primarie hanno richiesto circa 0,08–0,46 s per Espaloma e 0,10–0,18 s per SQE, su CPU. Il picco cumulativo del processo è circa 290 MiB per Espaloma e 75 MiB per SQE. Questi tempi escludono installazione, preparazione e caricamento dei pesi; non dimostrano scalabilità lineare a fogli arbitrariamente grandi.</p>
<p>Espaloma v0.0.8 usa il modello pubblicato invariato e il totale formale esplicito nel campo q_ref; massimo errore sul totale in questi test circa 5,1 × 10⁻⁶ e. SQE usa <a href="https://pypi.org/project/chargefw/0.1.7/">ChargeFW 0.1.7</a> nativo, set CCD_gen di <a href="https://doi.org/10.1186/s13321-021-00528-w">Schindler 2021</a>, calibrati su NPA, non sulle nostre RESP. Matching rigoroso, soluzione completa e nessun parametro rifittato. Totali SQE conservati entro circa 1,5 × 10⁻¹² e. Permutazione degli atomi e, per SQE, rototraslazione delle coordinate passano i controlli numerici. La riproducibilità numerica non risolve l’ambiguità chimica del grafo.</p>
<p>Job array 1464904 e 1464946: tutti e quattro i task completati, uscita 0. Le ultime analisi C54 sono anch’esse concluse: <strong>16 HF–RESP e 19 DFT–RESP</strong> disponibili sui frammenti. Nessun QM/RESP nuovo è stato eseguito sui fogli completi. Il controllo scheduler salvato qui riporta zero job attivi o in coda.</p>
<h2>Un QM diretto sull’intero foglio è possibile?</h2>
<p>La dimensione, da sola, non lo rende impossibile. <a href="https://www.cp2k.org/quickstep">CP2K/Quickstep</a> è progettato anche per sistemi dell’ordine di mille atomi. Il Supporting Information del <a href="https://doi.org/10.1021/acs.jpcb.2c06037">modello GAFF2–AIM</a> descrive single-point DFT su GO di oltre 10000 atomi con PBE, base DZVP-MOLOPT-SR-GTH e metodi a scaling ridotto; dichiara però anche un limite di 25 cicli SCF e una convergenza energetica approssimata. Quel protocollo produce DDEC6 e volumi atomici, non RESP HF/6-31G*.</p>
<p>Per il nostro obiettivo partirei dalla verifica dello stato elettronico e da un single-point DFT ben controllato su un foglio selezionato, come validazione della distribuzione di carica e del comportamento a lunga distanza. Il costo di un HF/6-31G* completo, della geometria e delle griglie multi-orientazione RESP richiede un benchmark separato prima di stimare se bastino i cinque giorni della coda. Le RESP sui frammenti restano il riferimento coerente con Amber per addestrare un assegnatore riproducibile e scalabile. Per un QM sui fogli con carica −16/−6 occorre inoltre definire ambiente, controioni o embedding e condizioni al contorno: un calcolo isolato in fase gas e una superficie in soluzione rispondono a problemi diversi. Non è stata stimata la stabilità dei polianioni isolati.</p>
<h2>File e provenienza</h2>
<ul><li><a href="per_atom_all_methods.csv">Tutte le cariche atomiche, varianti A/B/C</a>; <a href="per_role_all_methods.csv">distribuzioni per ruolo</a>; <a href="representation_sensitivity.csv">sensibilità alla rappresentazione</a>; <a href="basal_and_CBD_summary.csv">carboni basali e carica dei residui CBD</a>.</li>
<li><a href="GRAPH_QC.json">Audit chimico e hash delle topologie originali</a>; <a href="INDEPENDENT_GRAPH_AND_RESONANCE_CHECKS.json">verifica indipendente dei matching e del controllo COO</a>; <a href="SCHEDULER_STATUS.json">stato Slurm verificato</a>.</li>
<li><a href="../CHARGE_CONTROLS_UPDATE_20261103_attempt_001/REPORT.html">Confronto precedente sui frammenti e ultimo C54 recuperato</a>; <a href="../ESPALOMA_CONTROL_20261103_attempt_001/REPORT.html">metodi e controlli Espaloma</a>.</li>
<li>I file originali parm7/rst7 e le rappresentazioni SDF diagnostiche sono conservati in <code>source_and_diagnostic_graphs</code>. I risultati non sostituiscono le cariche di produzione.</li></ul>
<script>const data='''+json.dumps(widget,allow_nan=False)+''';function number(v){return (v>=0?'+':'−')+Math.abs(v).toFixed(3).replace('.',',')}function update(){const model=document.getElementById('model').value,variant=document.getElementById('variant').value;document.getElementById('charges').innerHTML=data.filter(r=>r.model===model&&r.variant===variant).map(r=>'<tr><td>'+r.role+'</td><td>'+r.count+'</td><td>'+number(r.source)+'</td>'+['espaloma','sqeq0','sqeqp'].map(k=>'<td title="5°–95° percentile: '+number(r[k+'_q05'])+' / '+number(r[k+'_q95'])+' e">'+number(r[k])+'</td>').join('')+'</tr>').join('')}update();</script>
</main></body></html>'''
(OUT/'REPORT.html').write_text(page)
comparison='<!doctype html><html lang="it"><head><meta charset="utf-8">'+css+'</head><body><main><img src="'+image('basal_distributions.png')+'" alt="Distribuzioni diagnostiche delle cariche sui carboni basali GGG."><p class="note">Varianti A: stesso foglio, tre metodi. La linea a 0 e indica la tua LIB. Il grafo elettronico dei fogli richiede verifica: queste distribuzioni sono diagnostiche.</p></main></body></html>'
(OUT/'comparison.html').write_text(comparison)
save(OUT/'REPORT_MANIFEST.json',{'report_SHA256':sha(OUT/'REPORT.html'),'source_archives':expected,
    'verified_backend_artifact_files':verified,'diagnostic_predictions':len(records),
    'production_charge_model_eligible':False,'new_QM_or_RESP_calculation':False})
shutil.copy2(Path(__file__),OUT/'build_report.py')
print(json.dumps({'report':str(OUT/'REPORT.html'),'verified_files':verified,'predictions':len(records)}))
