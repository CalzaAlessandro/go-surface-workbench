"""Report verified SQE/Espaloma controls and the final recovered C54 reference."""
import base64
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CALC = HERE.parent
OUT = ROOT/'reports/CHARGE_CONTROLS_UPDATE_20261103_attempt_001'
OUT.mkdir(exist_ok=False)

def read(p):
    return json.loads(p.read_text())

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def save(p, d):
    p.write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def table(p, rows):
    with p.open('x', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

sqe = CALC/'sqe_results_20261103_snapshot_001'
summary = read(sqe/'results/BENCHMARK_COMPLETE.json')
assert summary['status'] == 'BENCHMARK_COMPLETE' and not summary['failures']
for name, checksum in summary['artifacts_SHA256'].items():
    assert sha(sqe/'results'/name) == checksum
for name, checksum in read(sqe/'INPUT_MANIFEST.json').items():
    assert sha(sqe/name) == checksum
for name in ['BENCHMARK_COMPLETE.json', 'per_atom_charges.csv', 'per_context_role_charges.csv']:
    shutil.copy2(sqe/'results'/name, OUT/name)
shutil.copytree(sqe/'results/parameters', OUT/'parameters')
shutil.copy2(sqe/'requirements.frozen.txt', OUT/'requirements.frozen.txt')

from go_param.context import analyze_dataset, write_outputs
primary = CALC/'complementary_results_20261103_snapshot_001/atomic_inputs'
entries = [{'model_path': str(p/'model.json'), 'charges_path': str(p/'charges.json'),
            'fit_evidence_path': str(p/'FIT_COMPLETE.json')} for p in sorted(primary.glob('*/DFT'))]
extra = sorted((CALC/'c54_recovered_results_snapshot_20261103_attempt_001/outputs').glob('*/fits/attempt_001'))
extra.append(CALC/'c54_single_results_20261103_snapshot_001/fits/attempt_001')
for p in extra:
    marker = read(p/'FIT_COMPLETE.json')
    for name, checksum in marker['artifacts_sha256'].items():
        # Earlier compact snapshots include only their declared subset of files.
        if (p/name).exists():
            assert sha(p/name) == checksum
    entries.append({'model_path': str(p/'model.json'), 'charges_path': str(p/'charges.json'),
                    'fit_evidence_path': str(p/'FIT_COMPLETE.json')})
assert len(entries) == 19
context = analyze_dataset(entries)
assert context['coverage']['all_models_passed_charge_identity_and_total_checks']
write_outputs(context, OUT/'context_DFT_19')
save(OUT/'context_DFT_19/input_spec.json', entries)
final = extra[-1]
for name in ['FIT_COMPLETE.json', 'charges.csv']:
    shutil.copy2(final/name, OUT/('C54_final_'+name))

inputs = read(HERE/'inputs.json')
methods = [('current_library', 'Libreria attuale'), ('espaloma', 'Espaloma'),
           ('sqeq0', 'SQE+q0'), ('sqeqp', 'SQE+qp'), ('molecular_HF_RESP', 'HF–RESP molecolare')]
q = {}
for case in inputs['cases']:
    mid = case['model']['model_id']
    for key, values in case['references'].items():
        if values is None:
            continue
        by_role = {}
        for atom, value in zip(case['model']['atoms'], values):
            by_role.setdefault(atom['role'], []).append(value)
        for role, values in by_role.items():
            q[mid, role, key] = sum(values)/len(values)
for row in csv.DictReader((sqe/'results/per_context_role_charges.csv').open()):
    q[row['model_id'], row['role'], row['method']] = float(row['mean_charge_e'])

selected = [('pyrene_diol_cis', 'hydroxyl_O', 'OH · ossigeno'),
            ('pyrene_diol_cis', 'hydroxyl_H', 'OH · idrogeno'),
            ('pyrene_diol_epoxide', 'epoxide_O', 'Epossido · ossigeno, misto OH/EPX'),
            ('C24_single_basal_EPX_control', 'epoxide_O', 'Epossido · ossigeno, C24 singolo'),
            ('pyrene_dione', 'carbonyl_O', 'Carbonile · ossigeno'),
            ('pyrene_COOH', 'COOH_O', 'COOH · ossigeno carbonilico'),
            ('pyrene_COO', 'COO_C', 'COO⁻ · carbonio'),
            ('pyrene_COO', 'COO_O', 'COO⁻ · ossigeno'),
            ('pyrene_COO', 'COO_anchor', 'COO⁻ · carbonio di ancoraggio')]
selected_rows = [{'model_id': mid, 'role': role, **{key: q[mid,role,key] for key,_ in methods}}
                 for mid,role,_ in selected]
table(OUT/'selected_role_charges.csv', selected_rows)
oxygen = [(mid, role, label) for mid,role,label in selected if role.endswith('_O')]
with tempfile.TemporaryDirectory(prefix='go_sqe_plot_') as cache:
    os.environ['MPLCONFIGDIR'] = cache
    os.environ['XDG_CACHE_HOME'] = cache
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = ['#666666','#0072B2','#009E73','#D55E00','#CC79A7']
    with plt.rc_context({'font.size':9, 'svg.fonttype':'none'}):
        fig, ax = plt.subplots(figsize=(9,4.2), layout='constrained')
        for j,(key,label) in enumerate(methods):
            ys = [i + (j-2)*.12 for i in range(len(oxygen))]
            ax.scatter([q[mid,role,key] for mid,role,_ in oxygen], ys, label=label,
                       color=colors[j], marker=['o','s','^','D','x'][j], s=35)
        ax.set_yticks(range(len(oxygen)), [label for _,_,label in oxygen])
        ax.invert_yaxis();ax.set_xlabel('Carica media del ruolo nel contesto indicato (e)')
        ax.set_xlim(-1.02,-.02);ax.grid(axis='x',alpha=.2)
        ax.spines[['top','right']].set_visible(False)
        ax.legend(loc='upper center', bbox_to_anchor=(.4,1.18),ncols=3,frameon=False,fontsize=8)
        for suffix in ['png','svg']:fig.savefig(OUT/f'oxygen_charges.{suffix}',dpi=300,facecolor='white')
        plt.close(fig)

def number(v):
    return f'{v:+.3f}'.replace('-', '−').replace('.', ',')

rows = ''.join('<tr><td>'+html.escape(label)+'</td>'+''.join(
    '<td>'+number(q[mid,role,key])+'</td>' for key,_ in methods)+'</tr>' for mid,role,label in selected)
headers = ''.join('<th>'+label+'</th>' for _,label in methods)
img = 'data:image/png;base64,'+base64.b64encode((OUT/'oxygen_charges.png').read_bytes()).decode()
css = '''<style>:root{--foreground:#192936;--muted-foreground:#536270;--border:#ccd5dc;--ring:#2866a8;--font-sans:system-ui,sans-serif}body{margin:0;font:15px var(--font-sans);line-height:1.55;color:var(--foreground)}main{max-width:1060px;margin:auto}h2{font-size:1.15rem;margin:1.3em 0 .4em}img{display:block;width:100%;height:auto}a{color:var(--ring)}.note{font-size:.85rem;color:var(--muted-foreground)}.scroll{overflow-x:auto}table{border-collapse:collapse;font-size:.85rem;font-variant-numeric:tabular-nums;width:100%}th,td{padding:7px;border-bottom:1px solid var(--border);text-align:right}th:first-child,td:first-child{text-align:left}code{overflow-wrap:anywhere}</style>'''
page = '''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GO · Espaloma, SQE e RESP · aggiornamento 3 novembre 2026</title>'''+css+'''</head><body><main>
<h2>Espaloma e SQE: risultati effettivi sui frammenti GO</h2>
<p>Espaloma resta un candidato interessante. Il modello pubblicato restituisce cariche OH e COO⁻ plausibili come controllo, ma sposta molto le cariche epossidiche rispetto alla libreria attuale. Questo test valuta le cariche; non seleziona una nuova coppia cariche–LJ per dinamica.</p>
<p>SQE+q0 e SQE+qp hanno completato <strong>19/19</strong> frammenti. SQE standard ha completato i 18 neutri ed escluso correttamente il carbossilato −1. Nessun parametro mancante è stato sostituito con zero. Job Slurm <strong>1464485</strong>, completato in 5 s, uscita 0.</p>
<img src="'''+img+'''" alt="Confronto delle cariche di ossigeno nei contesti indicati, con cinque metodi distinti.">
<p class="note">Ogni riga è un contesto specifico, non una nuova carica universale. Sono medie per ruolo nel singolo frammento; i valori atomici e gli intervalli sono nei CSV. La somiglianza tra cariche non misura l’accuratezza fisica del force field.</p>
<div class="scroll"><table><thead><tr><th>Ruolo / contesto</th>'''+headers+'''</tr></thead><tbody>'''+rows+'''</tbody></table></div>
<h2>Interpretazione e decisione di lavoro</h2>
<p>Per l’epossido nel frammento misto, l’ossigeno passa da −0,370 e della libreria a −0,138 e con Espaloma e −0,607 e con SQE+qp; HF–RESP molecolare dà −0,419 e. I controlli non convergono a una sola distribuzione. Il confronto con HF–RESP riguarda l’intero frammento, con vincoli diversi da una libreria di blocchi neutri.</p>
<p>Nel carbossilato, SQE+q0 dà mediamente −0,908 e per ossigeno e SQE+qp −0,472 e. Il totale del frammento resta −1 e in entrambi: cambia quanto la carica è localizzata sul COO⁻ o distribuita sullo scheletro. La carica del sottoinsieme COO⁻ non è stata forzata a −1 nei riferimenti molecolari. È proprio uno dei casi da ampliare prima di scegliere un modello trasferibile.</p>
<p>Terrei Espaloma come controllo rapido orientato alle AM1-BCC. Per un futuro modello di redistribuzione calibrato su Amber, confronterei SQE+q0 e SQE+qp sulle nostre HF–RESP, con validazione su famiglie di frammenti escluse dal fit. Non adotterei direttamente i coefficienti SQE pubblicati come preset RESP per GO.</p>
<h2>Riproducibilità: dipendenza dalla forma di Kekulé</h2>
<p>La libreria nativa richiede legami con ordine intero. I grafi aromatici validati sono stati kekulizzati con RDKit senza cambiare atomi, idrogeni, geometria, connettività o carica formale. Una seconda rappresentazione equivalente, ottenuta invertendo l’ordine degli atomi prima della kekulizzazione, cambia alcuni ordini singolo/doppio.</p>
<p>Questa variazione cambia una carica atomica fino a <strong>0,0961 e per SQE+qp</strong> nel carbossilato, 0,0451 e per SQE+q0 e 0,0580 e per SQE standard. Anche i fogli pristini mostrano una dipendenza più piccola. La permutazione degli atomi mantenendo gli stessi ordini di legame e la rotazione/traslazione rigida passano invece i controlli numerici. Serve quindi definire e validare una gestione coerente dell’aromaticità/risonanza prima di farne un assegnatore generale.</p>
<p>Conservazione del totale: massimo errore 2,46 × 10⁻¹⁴ e per SQE+qp; nessuna correzione manuale applicata alle predizioni. Le cariche resterebbero fisse durante la MD. La redistribuzione uniforme iniziale dei q0 tabulati è parte della formulazione SQE+qp pubblicata; non è una correzione RESP né un uso della tua LIB come baseline.</p>
<h2>Metodi, copertura e limiti</h2>
<p>Implementazione C++ nativa <a href="https://pypi.org/project/chargefw/0.1.7/">ChargeFW 0.1.7</a>, distribuita dal manutentore di ChargeFW2; Python 3.11 in ambiente isolato. Set espliciti <code>SQE/SQEq0/SQEqp_Schindler2021_CCD_gen</code>, matching rigoroso, soluzione completa, nessuna selezione automatica o rifit. I coefficienti pubblicati sono calibrati su <strong>B3LYP/6-311G NPA</strong>, secondo i metadati del set e il <a href="https://doi.org/10.1186/s13321-021-00528-w">lavoro originale</a>. Non sono cariche RESP.</p>
<p>Espaloma usa i pesi pubblicati v0.0.8 e il campo q_ref del readout per fornire il totale formale. L’API non adattata di quella release dava circa 0 per acetato invece di −1; il controllo con l’adattatore conserva −1. Pesi invariati; 16/16 frammenti, nessun nuovo fit. Il modello è orientato alle AM1-BCC ELF10, che non coincidono necessariamente con la procedura nativa Amber/SQM del nostro controllo.</p>
<p>Le 12 comparazioni con AM1-BCC mantengono le esclusioni dei quattro casi privi di parametri BCC. I metadati stereochimici ereditati discordanti sono conservati: non è dichiarata validata la stereochimica né la trasferibilità a grandi superfici. SQE usa una soluzione densa, con memoria quadratica e costo della soluzione cubico nel numero di legami: i tempi sui frammenti non dimostrano scalabilità lineare.</p>
<h2>Le ultime analisi in coda sono finite</h2>
<p>Controllo scheduler del <strong>3 novembre 2026, 21:47 UTC</strong>: nessun job dell’utente in coda o in esecuzione. La sintesi 1452975, i recuperi C54 1453077 e 1453111, Espaloma 1463920 e SQE 1464485 risultano completati. Due task originali 1406363 restano marcati FAILED a causa del precedente postprocessing; i recuperi separati hanno prodotto i risultati verificati.</p>
<p>L’ultimo <code>C54_single_basal_EPX_control</code> è recuperato e verificato: 21 hash dei file, tre griglie complete da 121147, 121157 e 121397 punti, nessun nuovo Gaussian. Il pannello contiene ora <strong>19 contesti DFT–RESP e 16 HF–RESP</strong>; non sono 35 strutture chimiche indipendenti. HF/6-31G* resta il riferimento per le cariche Amber. DFT–RESP usa ωB97X-D/6-311++G(d,p). I riferimenti RHF mantengono l’avvertenza sulla stabilità elettronica già documentata.</p>
<p>Il report precedente a 18 contesti resta uno snapshot storico. La nuova <a href="context_DFT_19/context_analysis.json">analisi DFT dei 19 contesti</a> include cariche atomiche, distanze tra gruppi e carboni basali; <a href="C54_final_charges.csv">cariche dell’ultimo C54</a> e <a href="C54_final_FIT_COMPLETE.json">protocollo e metriche RESP</a>.</p>
<h2>Dati e provenienza</h2>
<ul><li><a href="per_atom_charges.csv">Cariche SQE per atomo, incluse le due rappresentazioni Kekulé</a>; <a href="per_context_role_charges.csv">medie e intervalli per ruolo/contesto</a>.</li>
<li><a href="selected_role_charges.csv">Valori dei cinque metodi riportati qui</a>; <a href="BENCHMARK_COMPLETE.json">esito, controlli numerici e hash SQE</a>; <a href="requirements.frozen.txt">dipendenze bloccate</a>.</li>
<li><a href="../ESPALOMA_CONTROL_20261103_attempt_001/REPORT.html">Report e dati Espaloma completi</a>; <a href="../CHARGE_COMPARISON_20261103_attempt_001/RESP_AMBER_REPORT.html">report RESP precedente</a>.</li>
<li><a href="oxygen_charges.svg">Figura vettoriale</a>. Originali scientifici e librerie di produzione preservati. Nessuna modifica a cariche, LJ, bonded o frcmod di produzione.</li></ul>
</main></body></html>'''
(OUT/'REPORT.html').write_text(page)
visual = '<!doctype html><html lang="it"><head><meta charset="utf-8">'+css+'</head><body><main><img src="'+img+'" alt="Cariche di ossigeno dei cinque metodi nei frammenti GO indicati."><p class="note">Medie per ruolo e contesto (e). Le predizioni SQE usano parametri NPA pubblicati, senza rifit. Espaloma usa il modello pubblicato con totale formale esplicito. La somiglianza delle cariche non valida il force field.</p></main></body></html>'
(OUT/'comparison.html').write_text(visual)
save(OUT/'REPORT_MANIFEST.json', {'report_SHA256':sha(OUT/'REPORT.html'),
     'SQE_archive_SHA256':sha(CALC/'sqe_results_20261103_snapshot_001.tar.gz'),
     'last_C54_archive_SHA256':sha(CALC/'c54_single_results_20261103_snapshot_001.tar.gz'),
     'DFT_contexts':19,'HF_contexts':16,'SQE_native_result_files_verified':len(summary['artifacts_SHA256']),
     'final_parameter_validation':False, 'figure_aggregation':'within-context means for a chemical role; not universal parameters'})
print(json.dumps({'report':str(OUT/'REPORT.html'),'DFT_contexts':context['context_count'],
                  'verified_SQE_result_files':len(summary['artifacts_SHA256'])}))
