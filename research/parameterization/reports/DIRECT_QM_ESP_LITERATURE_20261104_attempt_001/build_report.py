"""Report independent electronic-potential controls, mixing and literature scope."""
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import shutil
import statistics

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CALC=ROOT/'calculations'
EXPERIMENT=CALC/'direct_qm_esp_literature_20261104_attempt_001'
RESULTS=CALC/'direct_qm_esp_literature_results_20261104_snapshot_001/results'
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def save(path,data):path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def main():
    complete=read(RESULTS/'COMPLETE.json')
    assert complete['status']=='DIRECT_QM_ESP_CONTROL_COMPLETE' and complete['raw_QM_grids_verified']==54
    for name,s in complete['artifacts_SHA256'].items():assert sha(RESULTS/name)==s
    assert sha(EXPERIMENT/'inputs.json')==complete['input_SHA256']
    for name in ['COMPLETE.json','orientation_scores.csv','context_scores.csv','aggregate_scores.csv','QM_SOURCE_EVIDENCE.json','BCC_CHARGE_POLICIES.json']:
        shutil.copy2(RESULTS/name,HERE/name)
    for name in ['LITERATURE_COVERAGE.json','GAFF2_AIM_NATIVE_PARAMETERS.csv','GAFF2_AIM_NATIVE_PROVENANCE.json']:
        shutil.copy2(EXPERIMENT/name,HERE/name)
    context=rows(RESULTS/'context_scores.csv');aggregate=rows(RESULTS/'aggregate_scores.csv')
    inputs=read(EXPERIMENT/'inputs.json')
    held=[c['model']['model_id'] for c in inputs['cases'] if c['scope']=='held_out_C24' and c['BCC_source']]
    assert len(held)==5
    names=[('molecular_HF_RESP','RESP molecolare HF'),('BCC_integer_charge_control','AM1-BCC · controllo Q intero'),
           ('RESP_BCC_50_50_diagnostic','Miscela RESP/BCC 50:50'),('SQE_HF_RESP_C16','SQE–RESP · training C16'),
           ('L0_HF_RESP_C16','L0 RESP · training C16'),('Espaloma_BCC_ELF10','EspalomaCharge pubblicato'),
           ('literature_hiergo_opls','Tabella HierGO/OPLS'),('current_library','Libreria attuale'),
           ('molecular_DFT_RESP','RESP DFT · fit sul riferimento')]
    by_context={(r['model_id'],r['QM_reference'],r['candidate']):float(r['mean_ESP_RRMS']) for r in context}
    mix_cases=[]
    for mid in held:
        mix_cases.append({'model_id':mid,'HF_RESP':by_context[mid,'DFT','molecular_HF_RESP'],
                         'BCC':by_context[mid,'DFT','BCC_integer_charge_control'],
                         'mixture':by_context[mid,'DFT','RESP_BCC_50_50_diagnostic']})
    summary={'generated_at_utc':datetime.now(timezone.utc).isoformat(),'source_job_id':'1469066',
             'QM_grids_verified':54,'physical_contexts':9,'matched_C24_for_mixing':held,
             'mixing_diagnostic_DFT':mix_cases,
             'mixture_worse_than_HF_RESP_C24':sum(r['mixture']>r['HF_RESP'] for r in mix_cases),
             'mixture_better_than_BCC_C24':sum(r['mixture']<r['BCC'] for r in mix_cases),
             'aggregate_scores':aggregate,
             'QM_reference':'raw HF/6-31G* and wB97XD/6-311++G(d,p) electrostatic potentials, not fitted charge vectors',
             'independence':'DFT is independent of the HF-RESP fitting target, but remains an approximate electronic method; all candidates share the evaluated geometry',
             'GAFF2_AIM_finite_fragment_ESP_score_available':False,
             'GAFF2_AIM_scope':'native periodic charge/LJ inventory; finite edges and COO lack a validated published mapping',
             'new_QM_or_charge_parameter_fit':False,'production_validation':False}
    save(HERE/'SUMMARY.json',summary)
    native=rows(EXPERIMENT/'GAFF2_AIM_NATIVE_PARAMETERS.csv')
    def family(r):
        t=r['atom_type'];z=int(r['atomic_number'])
        if z==1:return 'Idrogeni (tipi H nativi)'
        if z==6:return 'Carboni, tipi ca*' if t.startswith('ca') else 'Carboni, tipi c3*' if t.startswith('c3') else 'Carboni, altri tipi'
        return 'Ossigeni, tipi oh*' if t.startswith('oh') else 'Ossigeni, tipi os*' if t.startswith('os') else 'Ossigeni, altri tipi'
    def decimal(v,digits=3):return f'{v:.{digits}f}'.replace('.',',')
    def percent(v):return decimal(100*v,1)+'%'
    css='<style>:root{font-family:var(--font-sans,system-ui,sans-serif);color:var(--foreground,#202832)}body{margin:0;font-size:16px}main{max-width:1100px;margin:auto}p,li{line-height:1.6}h1{font-size:1.6rem}h2{font-size:1.2rem;margin-top:1.6rem}h3{font-size:1.05rem}a{color:var(--primary,#005a8c)}.note{font-size:.94rem;color:var(--muted-foreground,#53616d)}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:1rem;margin:.7rem 0 1rem}td,th{padding:.7rem .6rem;text-align:left;vertical-align:top;border-bottom:1px solid var(--border,#d4dce3)}th{font-weight:650}.num{white-space:nowrap;font-variant-numeric:tabular-nums}code{overflow-wrap:anywhere}summary{cursor:pointer;padding:.6rem 0}</style>'
    caption='Errore relativo RMS del potenziale, espresso in percentuale: più basso significa maggiore fedeltà al DFT adottato. Ogni molecola pesa una volta: media delle tre rotazioni MK, poi mediana fra gli stessi cinque C24. Non è un errore sulle energie di adsorbimento. C24 è escluso dal training SQE/L0.'
    scopes={'molecular_HF_RESP':'Fit individuale sul potenziale HF; qui valutato su un diverso potenziale DFT',
        'BCC_integer_charge_control':'Antechamber; correzione del totale per gli arrotondamenti documentata',
        'RESP_BCC_50_50_diagnostic':'Media diagnostica a peso fisso; usa il fit individuale RESP',
        'SQE_HF_RESP_C16':'Parametri condivisi, fittati sui C16',
        'L0_HF_RESP_C16':'Cariche condivise per ruolo, fittate sui C16',
        'Espaloma_BCC_ELF10':'Pesi pubblicati con target AM1-BCC ELF10; nessun fit GO',
        'literature_hiergo_opls':'Tabella di cariche pubblicata; non intero force field',
        'current_library':'Blocchi della libreria corrente',
        'molecular_DFT_RESP':'Fit individuale sullo stesso riferimento DFT: controllo in sample'}
    comparisons=[]
    for name,label in names:
        d={r['QM_reference']:float(r['median_mean_ESP_RRMS']) for r in aggregate if r['subset']=='five_BCC_covered_C24' and r['candidate']==name}
        values=[by_context[mid,'DFT',name] for mid in held]
        comparisons.append(f'<tr><td>{html.escape(label)}</td><td class="num">{percent(d["DFT"])}</td><td class="num">{percent(min(values))}–{percent(max(values))}</td><td>{scopes[name]}</td></tr>')
    benchmark_table='<div class="scroll"><table><thead><tr><th>Candidato</th><th>Errore mediano DFT</th><th>Intervallo nei 5 casi</th><th>Come è ottenuto</th></tr></thead><tbody>'+''.join(comparisons)+'</tbody></table></div>'
    native_tables=[];family_rows=[]
    for source in sorted({r['source_topology'] for r in native}):
        rs=[r for r in native if r['source_topology']==source];body=[]
        for fam in sorted({family(r) for r in rs}):
            group=[r for r in rs if family(r)==fam];qs=[float(r['charge_e']) for r in group]
            n=sum(int(r['count']) for r in group);mean=sum(float(r['charge_e'])*int(r['count']) for r in group)/n
            family_rows.append({'source_topology':source,'family':fam,'native_types':len(group),'atoms':n,'mean_charge_e_atom_weighted':mean,'min_type_charge_e':min(qs),'max_type_charge_e':max(qs)})
            body.append(f'<tr><td>{fam}</td><td>{len(group)}</td><td>{n}</td><td class="num">{decimal(mean)}</td><td class="num">{decimal(min(qs))}–{decimal(max(qs))}</td></tr>')
        native_tables.append('<h3>Foglio periodico C/O ≈ '+decimal(float(rs[0]['C_over_O']),1)+'</h3><div class="scroll"><table><thead><tr><th>Famiglia nativa</th><th>Tipi</th><th>Atomi</th><th>q media (e)</th><th>q minima–massima (e)</th></tr></thead><tbody>'+''.join(body)+'</tbody></table></div>')
    with (HERE/'GAFF2_AIM_FAMILY_SUMMARY.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(family_rows[0]));w.writeheader();w.writerows(family_rows)
    inventory_table=''.join(native_tables)
    inventory_caption='Le medie pesano i tipi per il numero di atomi nei rispettivi fogli nativi. Gli intervalli descrivono valori dei tipi, non incertezza statistica. I fogli sono diversi dai frammenti QM del confronto e non costituiscono una graduatoria di accuratezza.'
    (HERE/'comparison.html').write_text('<!doctype html><html lang="it"><head><meta charset="utf-8">'+css+'</head><body><main><h2>Confronto sul potenziale DFT diretto</h2>'+benchmark_table+'<p class="note">'+caption+'</p><h2>GAFF2–AIM: cariche native per famiglia</h2>'+inventory_table+'<p class="note">'+inventory_caption+'</p></main></body></html>')
    coverage=read(EXPERIMENT/'LITERATURE_COVERAGE.json')
    coverage_rows=[]
    for candidate in sorted({r['candidate'] for r in coverage}):
        cs=[r for r in coverage if r['candidate']==candidate]
        scored=sum(r['scoring_status'].startswith('complete_') for r in cs)
        explanation=('Inventario nativo cariche + LJ; mapping sui frammenti finiti non definito' if candidate=='pinto_gaff2_aim' else
                     'Controllo di carica completo nei contesti coperti; COO− non fornito' if candidate=='hiergo_opls' else
                     'Dati parziali, caps ibridi o mancata mappatura: escluso dal confronto primario')
        coverage_rows.append(f'<tr><td>{candidate}</td><td>{scored}/9</td><td>{explanation}</td></tr>')
    text=f'''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">{css}<title>GO · QM diretto, RESP/BCC e letteratura</title></head><body><main>
    <h1>Cariche GO: confronto con il potenziale QM diretto</h1>
    <p class="note">Postanalisi completata, job 1469066. Nessun nuovo QM o fit di carica. Nessuna libreria di produzione modificata.</p>
    <h2>Quale riferimento è scientificamente difendibile?</h2>
    <p><strong>Raccomandazione: usare proprietà elettroniche QM dirette come riferimento primario, verificando il livello DFT; poi validare separatamente il modello a cariche fisse in acqua.</strong> RESP è un metodo di fit e può utilizzare un potenziale DFT. Il livello HF/6-31G* è il riferimento convenzionale del ramo Amber, non una verità universale per GO.</p>
    <div class="scroll"><table><thead><tr><th>Domanda</th><th>Riferimento appropriato</th><th>Limite da rispettare</th></tr></thead><tbody>
    <tr><td>Riproduciamo l’elettrostatica del GO?</td><td>Potenziale e dipolo QM diretti, con metodo/base/stato controllati</td><td>DFT è approssimato; serve una verifica del riferimento stesso</td></tr>
    <tr><td>Il modello si trasferisce a un nuovo foglio?</td><td>Frammenti, disposizioni e dimensioni esclusi dal fit</td><td>Il buon fit di una singola molecola non dimostra trasferibilità</td></tr>
    <tr><td>Funziona come modello a cariche fisse in acqua?</td><td>Polarizzazione media controllata, LJ/acqua coerenti, interazioni e proprietà indipendenti</td><td>Il migliore potenziale nel vuoto non garantisce la migliore MD in acqua</td></tr>
    </tbody></table></div>
    <p>Il controllo attuale ωB97X-D/6-311++G(d,p) è un riferimento di lavoro: non abbiamo ancora dimostrato che sia il migliore per GO. Il gate successivo è verificare sensibilità a funzionale/base, stato elettronico e dimensione. <a href="https://arxiv.org/abs/2309.12184">Studio sulla dipendenza dal metodo QM</a>.</p>
    <p>Per derivare cariche fisse efficaci in acqua, <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8204736/">RESP2</a> è un precedente di uso di riferimenti gas/acqua insieme alla validazione LJ. È una direzione da testare per GO, non un protocollo già adottato qui. Dopo l’assegnazione le cariche restano costanti in MD e GBSA. <a href="../../docs/DFT_REFERENCE_FIXED_CHARGES_20261104.md">Protocollo e verifiche mancanti</a>.</p>
    <p><strong>La miscela RESP/BCC 50:50 non migliora le RESP molecolari in questo controllo DFT.</strong> Nei cinque C24 con BCC disponibili è peggiore delle RESP e migliore delle BCC in tutti e cinque i casi. È un risultato su questo pannello, non una legge generale sui modelli di carica.</p>
    <h2>Il riferimento è il potenziale elettronico</h2><p>Valutiamo direttamente le griglie Gaussian HF/6-31G* e ωB97X-D/6-311++G(d,p), non il potenziale ricostruito dalle loro cariche RESP. Sono verificati 54 file ESP e i corrispondenti log: hash, metodo, convergenza, totale/spin, identità atomica e stessa geometria entro la precisione delle coordinate.</p>
    <p>La verifica dei log conferma il metodo e la molteplicità dichiarati, non lo stato fondamentale. I riferimenti HF sono controlli RHF convenzionali; questa postanalisi non risolve le instabilità esterne già documentate né stabilisce una nuova stabilità elettronica.</p>
    <p>Il DFT è un controllo elettronico diverso dall’obiettivo HF–RESP, ma resta approssimato. Non è una verità esatta né una validazione automatica in acqua. Le cariche atomiche non sono univocamente determinate dal QM: potenziale, densità e momenti sono riferimenti più diretti dei singoli numeri di una partizione.</p>
    <h2>RESP, BCC, miscela e altri candidati</h2>{benchmark_table}<p class="note">{caption}</p>
    <p>RESP molecolare indica il fit sul singolo sistema. L0 e SQE condividono parametri derivati dai C16; questo confronto non dimostra che un fit individuale sia trasferibile come libreria. Sui sette C24 supportati, la mediana DFT di SQE è circa 0.223 e quella L0 circa 0.231: valori vicini che non selezionano da soli un vincitore.</p>
    <p>Espaloma usa pesi pubblicati con target AM1-BCC ELF10; SQE usa un fit alle nostre RESP. La prova valuta questi candidati specifici, con target di addestramento diversi, e non stabilisce una superiorità universale di un’architettura. I grandi residui di Espaloma in queste geometrie suggeriscono cautela nel trasferimento al GO.</p>
    <h2>Miscelare i target è una nuova scelta di modello</h2><p>Due calibrazioni separate SQE–RESP e SQE–BCC sono coerenti e confrontabili. Un fit con entrambi i target, o una media delle cariche, produce un terzo compromesso. La media conserva il totale se entrambi gli ingressi hanno lo stesso Q; il peso deve comunque essere motivato e validato su dati indipendenti. Qui il 50:50 è fissato a priori e solo diagnostico; non è stato scelto il peso migliore dopo aver visto il test.</p>
    <p>AM1-BCC è parametrizzato per emulare il potenziale HF/6-31G*, quindi il suo accordo con RESP non è una conferma indipendente di quel riferimento. Non si pretende che ogni carica BCC coincida con un particolare fit RESP. <a href="https://doi.org/10.1002/jcc.10128">Jakalian et al.</a></p>
    <h2>Cariche di letteratura: copertura esplicita</h2><div class="scroll"><table><thead><tr><th>Fonte</th><th>Contesti con score</th><th>Perimetro</th></tr></thead><tbody>{''.join(coverage_rows)}</tbody></table></div>
    <p>HierGO/OPLS è valutato mediante la tabella pubblicata, senza eseguire qui il suo intero force field. GOPY e Stauffer restano inventari/parziali nei casi non coperti: nessuna carica mancante è sostituita con una carica della nostra libreria. L’OLD appartiene ai fogli QM/MM originali, senza un mapping univoco su questi frammenti.</p>
    <h2>GAFF2–AIM sui propri sistemi nativi</h2>{inventory_table}<p class="note">{inventory_caption}</p>
    <p>Sono inclusi i 132 tipi delle topologie pubblicate C/O ≈ 1.8 e 3.0: 15597 e 12819 atomi. Per ogni tipo conserviamo carica, σ, ε e abbondanza. Non sono cariche attribuite ex novo ai nostri frammenti. Il modello è derivato da DDEC6/AIM-LJ e descrive i contesti basali periodici pubblicati; finite edges, carbonili e carbossilati richiedono una mappatura/estensione non inventata qui.</p>
    <p><strong>Non è ancora disponibile uno score ESP di GAFF2–AIM confrontabile su questi frammenti finiti.</strong> Per ottenerlo serve un sistema con chimica e condizioni al contorno compatibili, con un riferimento QM adeguato. Un confronto periodico richiede anche coerenza del riferimento del potenziale e delle condizioni elettrostatiche. Riempire i bordi con cariche RESP produrrebbe un modello ibrido diverso dal pubblicato.</p>
    <p>In MD GAFF2–AIM rimane il bundle pubblicato di cariche, LJ, bonded, mixing e convenzioni 1–4. L’inventario di cariche non ne autorizza l’abbinamento a LJ arbitrari. <a href="https://doi.org/10.1021/acs.jpcb.2c06037">Pinto et al.</a> · <a href="https://github.com/arvpinto/GO_sheets">Dati originali</a> · <a href="GAFF2_AIM_NATIVE_PARAMETERS.csv">Tabella nativa cariche/LJ</a>.</p>
    <h2>Quanto possiamo concludere</h2><p>Il controllo non indica che RESP sia intrinsecamente errato per GO. Mostra che la sua rappresentazione molecolare conserva un buon accordo con questo secondo potenziale QM nei contesti studiati. Il quadro resta limitato da scelta del funzionale/base, stabilità/spin e frammenti finiti. Nei carbossili i due esempi C16 sono diagnostici di training e non entrano nella mediana C24; nel COO− il termine di monopolo può inoltre rendere piccolo il RRMS relativo senza garantire ogni dettaglio locale.</p>
    <p>Il gas-phase QM testa la fedeltà elettronica. Le cariche fisse per acqua sono parametri efficaci: la nota sovrapolarizzazione HF/6-31G* può contribuire al comportamento in fase condensata. Un modello con migliore ESP nel vuoto deve ancora essere verificato con LJ e solvente coerenti, interazioni con sonde e osservabili indipendenti. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10867846/">Discussione sul riferimento QM e sulla polarizzazione</a>.</p>
    <h2>Dati e riproducibilità</h2><p><a href="SUMMARY.json">Riepilogo</a> · <a href="context_scores.csv">Score per contesto</a> · <a href="orientation_scores.csv">Griglie/orientazioni</a> · <a href="QM_SOURCE_EVIDENCE.json">54 sorgenti QM verificate</a> · <a href="LITERATURE_COVERAGE.json">Copertura letteratura</a> · <a href="BCC_CHARGE_POLICIES.json">Correzione BCC dichiarata</a> · <a href="REPORT_MANIFEST.json">Checksum</a>.</p>
    <p class="note">La prima postanalisi 1469054 si è fermata per un percorso BCC incompleto; il recupero 1469060 e l’estensione 1469066 usano la snapshot verificata. Gli input e i log precedenti restano conservati. Nessun calcolo elettronico è stato ripetuto per questa analisi.</p></main></body></html>'''
    (HERE/'REPORT.html').write_text(text)
    save(HERE/'REPORT_MANIFEST.json',{'generated_at_utc':summary['generated_at_utc'],
        'files_SHA256':{p.name:sha(p) for p in HERE.iterdir() if p.is_file() and p.name!='REPORT_MANIFEST.json'},
        'source_complete_SHA256':sha(RESULTS/'COMPLETE.json'),'input_SHA256':sha(EXPERIMENT/'inputs.json'),
        'literature_source_SHA256':{name:sha(ROOT/'references'/name) for name in ['hiergo_17_charge_rows.json','stauffer2014_table_I_GO_charges.json','gaff2_aim_source_parameters.json']},
        'presentation':'tables_only; previous figure artifacts retained but not embedded',
        'mixture_selection':'fixed 50:50 diagnostic, not fitted or selected on the test','production_validation':False})
    print(json.dumps({'verified_QM_grids':54,'mixture_worse_than_RESP_on_C24':summary['mixture_worse_than_HF_RESP_C24'],
                      'GAFF2_AIM_inventory_types':132,'GAFF2_AIM_comparable_finite_ESP_available':False},indent=2))


if __name__=='__main__':main()
