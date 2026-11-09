"""Build a self-contained, table-based research status report from artifacts."""
import datetime,hashlib,html,json,pathlib
HERE=pathlib.Path(__file__).resolve().parent
ROOT=pathlib.Path('/Users/alessandrocalza/Desktop/CURRENT_PROJECTS/GO_BUILDER/PARAMETERIZATION')
ESC=html.escape

def table(headers,rows):
 return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+ESC(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def link(label,url):return '<a href="'+ESC(url,quote=True)+'">'+ESC(label)+'</a>'

probe=json.loads((ROOT/'calculations/qm_reference_runtime_results_20261104_snapshot_001/VALIDATION.json').read_text())
integration=ROOT/'calculations/qm_reference_integration_results_20261104_snapshot_001/outputs/runtime_water'
r=json.loads((integration/'COMPARISON_STATUS.json').read_text())
pairs=json.loads((integration/'PAIRWISE_QM_PROPERTIES.json').read_text())
assert all(x['accepted'] and x['esp_extracted'] for x in r['methods']) and len(pairs)==3 and all(x['matched_grid'] for x in pairs)
refs=ROOT/'references/go_prior_art_20261104'
csiro=json.loads((refs/'external_structures/csiro_neutral_candidates/CANDIDATES.json').read_text())
mace=json.loads((refs/'external_structures/STRUCTURE_CATALOG.json').read_text())
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
summary={'built_utc':now,'scientific_scope':'Independent QM validation; conventional Amber charge branch preserved','reference_selected':False,'production_charges_changed':False,'runtime_probe_job':'1469071','runtime_probe_passed':True,'end_to_end_job':'1469072','end_to_end_passed':True,'new_GO_array':'1469073','new_anchor_array':'1469076','GO_models':8,'GO_property_levels':48,'anchor_models':6,'anchor_property_levels':54,'finite_field_anchors':2,'local_tests_passed':238,'external_candidates_not_yet_curated':{'CSIRO_finite':len(csiro['files']),'GO_MACE_periodic_frames':len(mace['records'][0]['frames']),'GO_MACE_large_sheet_atoms':mace['records'][1]['frames'][0]['atoms']}}
(HERE/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
parts=['''<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GO — riferimento QM e compatibilità AMBER</title><style>
body{font:17px/1.55 system-ui,sans-serif;color:#152b37;background:#f3f6f8;margin:0}main{max-width:1200px;margin:auto;padding:32px 26px 60px}h1{font-size:32px;line-height:1.2}h2{font-size:23px;margin-top:34px}p{max-width:1000px}.note{background:#e3eef5;border-left:5px solid #176b8c;padding:16px 20px}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;background:white;margin:16px 0}th,td{padding:12px 14px;text-align:left;vertical-align:top;border-bottom:1px solid #d4e0e7}th{background:#123f53;color:white}td:first-child{font-weight:600}a{color:#075d88}small{color:#526573}code{overflow-wrap:anywhere}.status{font-weight:700;color:#744600}@media(max-width:600px){main{padding:18px 12px}body{font-size:15px}th,td{padding:9px}h1{font-size:26px}}@media print{body{background:white;font-size:10pt}main{padding:0}th{color:black;background:#e9eef1}tr{break-inside:avoid}a{color:black}}</style><main><h1>Cariche del GO: riferimento QM e coerenza AMBER</h1>''',
'<p><small>Report del 4 novembre 2026. Generato: '+ESC(now)+'. Stati dei job: verifica di questa sessione; non sono un monitor live.</small></p>',
'''<p class="note"><b>Decisione attuale:</b> HF/6-31G*–RESP resta la baseline convenzionale per AMBER/GAFF. Il confronto DFT/correlato serve a misurare errori e limiti del riferimento elettronico. Non è stato selezionato un nuovo preset né modificata una carica di produzione. Le cariche finali saranno fisse durante MD e GBSA, con totale intero determinato dalla chimica.</p>
<h2>Precedenti specifici e database</h2>''']
parts.append(table(['Fonte','Cosa fornisce','Uso per il nostro progetto / limite'],[
[link('AmberGO — Data 2026','https://www.mdpi.com/2306-5729/11/1/18'), 'Workflow e 15 grandi modelli GO compatibili con AMBER; tutorial e file nei repository degli autori.', 'Precedente diretto per costruzione/typing/parametri. Il tutorial assegna valori rappresentativi fissi: OH C/O/H = +0.16/−0.57/+0.41 e; EPX C/O = +0.16/−0.32 e. Questi valori non costituiscono un nuovo benchmark HF–RESP per tutti gli ambienti chimici. '+link('Tutorial','https://github.com/yingling-group/AmberGO/blob/main/Building_AMBER_Compatible_Graphene_Oxide_Tutorial.pdf')],
[link('CSIRO — Neutral GO Data Set','https://data.csiro.au/collection/csiro:41659'), '20.396 nanoflakes finiti, C150–C1014, rilassati con SCC-DFTB; file XYZ, diverse forme e terminazioni.', 'Sorgente di geometrie e diversità. DFTB/Mulliken non equivale al target HF–RESP. Neutralità elettronica e stati non magnetici del dataset non stabiliscono da soli protonazione o spin fondamentale del modello scelto.'],
[link('GO-MACE-23','https://zenodo.org/records/14066557'), 'Geometrie GO periodiche, dati CASTEP/PBE di energia e forze e strutture dopo annealing.', 'Utile per difetti, distribuzioni e deformazioni. I metadati riportano popn_calculate=false: i dati pubblicati non sono un dataset RESP pronto. Per celle periodiche occorre gestire esplicitamente i bordi prima di calcoli molecolari.'],
[link('al-Badri et al., Carbon 2021','https://discovery.ucl.ac.uk/id/eprint/10164884/'), 'QM ONETEP/PBE su GO da 979 atomi; DDEC e LJ ricavati dai volumi atomici, su base bonded OPLS.', 'Dimostra l’importanza di derivare e validare insieme cariche e LJ, anche con acqua/ioni. È un confronto scientifico pertinente, non un set di sole cariche da trasferire automaticamente su GAFF2.'],
[link('GAFF2-AIM','https://doi.org/10.1021/acs.jpcb.2c06037'), 'Preset pubblicato con DDEC6/AIM-LJ e base bonded GAFF2, già censito nel progetto.', 'Resta un preset completo distinto. Copertura prevalentemente basale OH/EPX periodica: bordi finiti, COO⁻ e giunzioni di graft richiedono verifica dedicata.']]))
parts+=['<h2>Strutture esterne effettivamente conservate</h2><p>I file sono candidati per la selezione del dataset; nessuno è già stato usato come nuova etichetta RESP o come carica di produzione.</p>']
rows=[]
for f in csiro['files']:
 composition=' '.join(f'{k}{v}' for k,v in f['composition'].items())
 rows.append([ESC(f['filename']),str(f['atoms'])+' atomi; '+ESC(composition),'Finito; SCC-DFTB; legami, passivazione, carica chimica e spin ancora da curare.'])
rows += [['GO-MACE iter-0-test.xyz','75 frame; 154–280 atomi','Periodico; non interpretare frame correlati come test indipendenti.'],['GO-MACE 900K-optimized.xyz','1 foglio; 15.009 atomi','Copertura strutturale; dimensione e periodicità lo rendono diverso dai controlli molecolari piccoli.']]
parts.append(table(['File / pool','Dimensioni verificate','Stato'],rows))
parts+=['''<p>Il campione CSIRO è stato scelto per dimensione del file solo per verificare l’accesso e la fattibilità: non è ancora un campione rappresentativo. Il catalogo completo degli XYZ è locale; non è stato scaricato l’intero database. Revisioni, URL stabili e SHA-256 sono salvati nei cataloghi.</p>
<h2>Come li useremo per calibrare le cariche</h2><ol>
<li>Curare grafo, unità, gruppi, bordi/H, protonazione, carica totale e molteplicità; verificare duplicati e celle periodiche.</li>
<li>Selezionare per dimensione, ossidazione, gruppi, clustering e curvatura. Conservare frammenti controllati per separare effetti di distanza e di composizione.</li>
<li>Separare training e test per foglio/genitore e strutture simili prima del fit. Riservare famiglie esterne mai usate per scegliere il modello.</li>
<li>Ricalcolare il riferimento elettrostatico con un protocollo uniforme. Non mescolare cariche DFTB, DDEC, OPLS e RESP come se fossero la stessa etichetta.</li>
<li>Confrontare libreria fissa, correzioni di contesto e SQE sullo stesso target; imporre il totale Q chimico, senza ridistribuzioni arbitrarie sul carbonio basale.</li>
</ol>
<h2>Validazione QM avviata</h2>''']
parts.append(table(['Campagna','Disegno','Stato verificato'],[
['Controlli GO — 1469073','8 geometrie esistenti fisse × 6 livelli: HF/6-31G*, ωB97X-D/6-311++G(d,p), ωB97X-D/aug-cc-pVTZ, PBE0/aug-cc-pVTZ, PW6B95/aug-cc-pVTZ e /aug-cc-pVDZ.','Primi due task in esecuzione; altri in coda. Niente nuove geometrie GO casuali.'],
['Molecole di controllo — 1469076','Benzene, fenolo, benzene oxide, p-benzochinone, acido acetico e acetato. Ottimizzazione/frequenze/stabilità prima del confronto; livelli HF/DFT, MP2 con due basi e CCSD.','Primi due task in esecuzione; altri in coda. Benzene oxide e acetato prevedono anche CCSD(T) a campi ±0.0001 e ±0.0002 a.u. nelle tre direzioni.'],
['Fogli interi — 1468880','GOX/REGOX originali; controlli di spin già avviati con ORCA/PBE e CPCM acqua.','Due task attivi; ultima lettura: nessuna convergenza SCF e nessuna terminazione normale. Questi risultati non sono ancora riferimenti RESP.'],
['RESP estese — 1468993','55 nuovi contesti già in coda: carico, prossimità e combinazioni di gruppi.','Quattro task attivi, rimanenti in attesa. Nessun completamento scientifico nuovo acquisito in questo report.']]))
parts.append('''<p>Ogni nuovo task reale dispone di 16 CPU, 64 GB e limite di cinque giorni; massimo due task per ciascuna nuova campagna. Geometria, Q, spin e orientazione sono comuni nel confronto di metodo/base. I gate distinguono proprietà complete, stabilità del determinante e convergenza dai soli exit code.</p>
<p>Le prime RHF→UHF instabilità sono registrate, non corrette cambiando arbitrariamente spin o idrogeni. Gli output RHF convenzionali restano diagnosi AMBER consultabili, ma non vengono presentati come stati elettronici stabili. La calibrazione della libreria non viene aggiornata automaticamente da questi confronti.</p>
<h2>Controlli preparatori già superati</h2>''')
parts.append(table(['Controllo','Evidenza','Cosa dimostra'],[
['Gaussian runtime — 1469071','21 casi positivi accettati; controllo negativo Density=SCF respinto.','Sintassi, selezione delle densità MP2/CC e separazione dal successivo link HF di stabilità. Non accuratezza del GO.'],
['Campo finito CCSD su acqua','Dipolo analitico 1.7242 D; estrapolato 1.724234 D. Errore ≈ 0.000034 D; differenza tra passi ≈ 0.0000064 D.','Verifica di unità e segno del campo. Per CCSD(T) si usano derivate di energia, senza spacciare la densità CCSD per CCSD(T).'],
['Worker completo — 1469072','3 metodi accettati; espgen nativo; 17.719 punti coincidenti e 3 confronti di coppia.','Pipeline eseguibile e confronto su griglie/coordinate effettivamente comuni.'],
['Test locali','238 superati, compresi log Gaussian reali e controlli negativi.','Controlli software; distinti dalla validazione scientifica del futuro modello di carica.']]))
parts+=['''<h2>Il vincolo AMBER resta decisivo</h2><p>Seguire RESP non significa che ogni fit a una qualunque densità QM sia automaticamente compatibile con GAFF. Il target convenzionale HF/6-31G* e AM1-BCC sono documentati dal tutorial ufficiale AMBER. La compatibilità del motore con un file di cariche non dimostra il bilanciamento fisico con proteine, ligandi, acqua, ioni e termini LJ.</p>
<p>La scelta ha quindi due passaggi: validare il riferimento elettronico su proprietà QM (potenziale, dipoli, base, stato elettronico), poi validare il modello a cariche fisse con il force field previsto. Un eventuale target DFT/RESP alternativo resta sperimentale finché non supera entrambi. L’accordo fra funzionali da solo non basta; i controlli correlati sono limitati alle molecole e richiedono valutazione del carattere multiriferimento. Nessun T1 mancante viene interpretato come T1 basso.</p>
<p class="status">Ancora da decidere: riferimento QM alternativo, calibrazione definitiva SQE/contesto e scelta delle cariche di produzione. I job avviati non costituiscono già questa decisione.</p>
<p>Fonti metodologiche: '''+link('Tutorial ufficiale AMBER / GAFF','https://ambermd.org/tutorials/basic/tutorial4b/')+'; '+link('RESP2 e accoppiamento cariche–LJ','https://pmc.ncbi.nlm.nih.gov/articles/PMC8204736/')+'''.</p><p>Protocollo e log: <code>PARAMETERIZATION/docs/QM_REFERENCE_VALIDATION_PLAN_20261104.md</code>; <code>PARAMETERIZATION/calculations/qm_reference_*</code>. Cataloghi esterni: <code>PARAMETERIZATION/references/go_prior_art_20261104/</code>.</p></main></html>''']
(HERE/'REPORT.html').write_text('\n'.join(parts))
files={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in HERE.iterdir() if f.is_file() and f.name!='REPORT_MANIFEST.json'}
(HERE/'REPORT_MANIFEST.json').write_text(json.dumps({'built_utc':now,'files_sha256':files},indent=2)+'\n')
print('HTML report written:',HERE/'REPORT.html')
