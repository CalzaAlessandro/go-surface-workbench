import pathlib,csv,json,hashlib,html,datetime
P=pathlib.Path(__file__).resolve().parent;R=P.parents[1]
source=R/'calculations/ambergo_qm_results_20261104_snapshot_001'
raw=list(csv.DictReader((source/'summary.csv').open()))
rows={r['candidate']:{k:(v if k=='candidate' else int(v) if k=='contexts' else float(v)) for k,v in r.items()} for r in raw}
assert all(r['contexts']==7 for r in rows.values())
# Reproduce the earlier validated benchmark exactly for unchanged candidates.
old=list(csv.DictReader((R/'reports/DIRECT_QM_ESP_LITERATURE_20261104_attempt_001/aggregate_scores.csv').open()))
for r in old:
 if r['subset']=='all_seven_C24':
  assert abs(rows[r['candidate']][r['QM_reference']+'_potential_median_RRMS']-float(r['median_mean_ESP_RRMS']))<1e-12
amber=[rows[k] for k in ['AmberGO_OPLS_caps','AmberGO_L0_C16_caps','AmberGO_current_caps']]
main=[('Libreria RESP L0, training C16',rows['L0_HF_RESP_C16']),('SQE–RESP, training C16',rows['SQE_HF_RESP_C16']),('AmberGO + bordi ibridi (3 scelte)',None),('HierGO / OPLS',rows['literature_hiergo_opls']),('Tua LIB originale',rows['current_library']),('EspalomaCharge, modello provato',rows['Espaloma_BCC_ELF10'])]
fields=['HF_potential_median_RRMS','DFT_potential_median_RRMS','HF_RESP_atomic_median_RMSE_e']
summary=[]
for label,row in main:
 entry={'model':label}
 for key in fields:entry[key]=[min(r[key] for r in amber),max(r[key] for r in amber)] if row is None else row[key]
 summary.append(entry)
(P/'SUMMARY.json').write_text(json.dumps({'job_id':'1469100','completed':True,'no_new_QM':True,'no_new_charge_fit':True,'contexts':7,'verified_QM_grids':42,'training_transfer':'C16→C24 for L0/SQE; individual molecular RESP references fitted in-sample','QM_DFT_method':'wB97XD/6-311++G(d,p)','HF_reference':'HF/6-31G* conventional branch with known RHF stability caveat','metric':'mean of 3 rigid orientations per context, then median over 7 contexts','AmberGO_limit':'Hybrid finite H-edge caps not provided by AmberGO; intervals across 3 cap choices, not statistical uncertainty','GAFF2_AIM':'Not numerically ranked: no validated native assignment to these finite fragments','production_validation':False,'summary':summary,'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in source.iterdir() if f.is_file()}},indent=2)+'\n')
for f in ['summary.csv','context_scores.csv','charge_scores.csv','QM_EVIDENCE.json']: (P/f).write_bytes((source/f).read_bytes())
def esc(x):return html.escape(str(x))
def table(headers,rows):return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(s)+'</th>' for s in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def val(x,pct=False):
 def one(v):return f'{100*v:.1f}%' if pct else f'{v:.3f}'
 return ' – '.join(one(v) for v in x) if isinstance(x,list) else one(x)
parts=['''<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AmberGO: confronto QM e origine dei parametri</title><style>body{font:15px/1.5 var(--font-sans,system-ui);color:var(--foreground,#24313b)}h2{font-size:19px}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid var(--border,#ccd3d8);padding:8px;text-align:right;vertical-align:top}th:first-child,td:first-child{text-align:left}td{white-space:nowrap}td:first-child{white-space:normal}p{margin:12px 0}svg{width:100%;height:auto;display:block}svg text{fill:var(--foreground,#24313b);font:12px system-ui}a{color:var(--primary,#146b92)}.text td{text-align:left;white-space:normal}</style><body><h2>Confronto sui sette frammenti C24</h2><p>Postanalisi completata: job 1469100, 42 griglie QM esistenti verificate. Nessun nuovo QM o fit. Valori più bassi indicano un migliore accordo con il riferimento specificato. Le cariche atomiche RESP e il potenziale esterno sono due metriche distinte.</p>''']
parts.append(table(['Assegnazione','Errore potenziale HF','Errore potenziale DFT','Scarto da HF–RESP, e'],[[esc(r['model']),val(r[fields[0]],True),val(r[fields[1]],True),val(r[fields[2]])] for r in summary]))
parts.append('''<p>DFT = ωB97X-D/6-311++G(d,p), controllo attuale ancora in fase di verifica metodologica. Media sulle tre orientazioni rigide di ogni contesto, poi mediana sui sette contesti; non 21 campioni indipendenti. La branca HF convenzionale conserva il limite noto delle instabilità RHF nei PAH.</p><p><b>AmberGO è un confronto diagnostico ibrido:</b> i suoi gruppi basali sono quelli originali, ma le cariche terminali C/H sono ±0.115, ±0.119646 oppure ±0.150 e. Gli intervalli descrivono queste tre scelte, non un’incertezza statistica. Non è un test del force field AmberGO completo.</p>''')
svg=['<svg viewBox="0 0 720 310" role="img" aria-label="Errore relativo del potenziale rispetto al DFT attuale, con intervallo di sensibilità ai bordi per AmberGO"><text x="190" y="18">Errore sul potenziale DFT, % — scala da zero</text>']
left=190;right=655;x=lambda v:left+v*100/150*(right-left)
for tick in [0,25,50,75,100,125,150]:
 xx=x(tick/100);svg.append(f'<line x1="{xx}" x2="{xx}" y1="30" y2="270" stroke="var(--border,#ccd3d8)"/><text x="{xx}" y="294" text-anchor="middle">{tick}</text>')
short=['RESP L0','SQE–RESP','AmberGO + bordi','HierGO / OPLS','Tua LIB originale','EspalomaCharge']
for i,r in enumerate(summary):
 y=48+i*40;v=r['DFT_potential_median_RRMS'];svg.append(f'<text x="2" y="{y+4}">{short[i]}</text>')
 if isinstance(v,list):
  lo,hi=v;svg.append(f'<line x1="{x(lo)}" x2="{x(hi)}" y1="{y}" y2="{y}" stroke="var(--chart-2,#197c67)" stroke-width="7"/><text x="{x(hi)+8}" y="{y+4}">{100*lo:.1f}–{100*hi:.1f}</text>')
 else:svg.append(f'<line x1="{left}" x2="{x(v)}" y1="{y}" y2="{y}" stroke="var(--chart-1,#286d9b)" stroke-width="7"/><text x="{x(v)+8}" y="{y+4}">{100*v:.1f}</text>')
svg.append('</svg>');parts.append(''.join(svg))
parts.append('''<h2>Cosa possiamo concludere</h2><p>SQE è il più vicino ai singoli valori HF–RESP tra gli assegnatori trasferibili valutati. L0 riproduce meglio il potenziale HF, mentre SQE è leggermente migliore sul DFT attuale. AmberGO è promettente per OH/epossidi, ma la sensibilità ai bordi impedisce una classifica generale del preset. La distanza fra vettori di carica non è identica alla qualità del potenziale: cariche interne differenti possono compensarsi.</p><p>I RESP molecolari sono controlli fittati sul singolo sistema, quindi non sono un test indipendente di trasferibilità della libreria. Il modello molecular_HF_RESP ha RRMS mediano HF 7.55%; molecular_DFT_RESP ha RRMS mediano DFT 10.98%. Non sono cariche già disponibili per un foglio arbitrario senza un nuovo calcolo/fit.</p><p>GAFF2-AIM e OLD non hanno una mappatura validata sugli stessi frammenti, quindi non ricevono una posizione numerica. GOPY grezzo ha problemi di bilancio dei blocchi ed è escluso dal confronto principale. AM1-BCC nativo ha copertura parametrica soltanto su cinque dei sette contesti.</p><h2>AmberGO: provenienza distinta per ogni componente</h2>''')
parts.append(table(['Componente','Origine verificata'],[
['Geometria / distribuzione dei gruppi','HierGO, con modifiche di costruzione e pulizia DS/VMD. La citazione di HierGO riguarda questo ruolo; le cariche non coincidono con la sua tabella OPLS.'],
['Cariche','L’articolo richiama Stauffer et al. 2014: ESP da DFT periodica PBE, con scaling SESP per la fase condensata. I valori MOL2 AmberGO differiscono dalla tabella originale: ad esempio O epossidico −0.32 contro −0.36 e SESP. La derivazione specifica di questa variazione non è stabilita dal protocollo consultato.'],
['Parametri GAFF','Base GAFF e GO_gaff_Hoshin.frcmod riutilizzato dal gruppo, attribuito a Hoshin Kim. Il tutorial riporta il riuso, non una nuova derivazione QM sistematica di tutti i termini.'],
['File GAFF2 aggiuntivo','GO_gaff2_Merve.frcmod è incluso, ma i MOL2 finali seguono il typing GAFF semplificato. Molti termini del file GAFF2 sono assegnati per analogia; la presenza del file non prova una nuova parametrizzazione interamente ottimizzata.']]))
parts.append('''<h2>ML dedicato alle cariche</h2><p>Esistono EspalomaCharge e OpenFF NAGL per cariche di tipo AM1-BCC. Espaloma è già stato testato qui, con prestazioni insufficienti su questo sottoinsieme GO: non significa che fallisca in generale sulle molecole del suo dominio. NAGL non è ancora stato testato sui nostri fogli. Un predittore addestrato a RESP è possibile, ma richiede dati RESP e verifiche di trasferibilità; non si ottiene dal normale GO-MACE senza cambiare modello e target.</p><p>GO-MACE viene lasciato fuori dal percorso di assegnazione delle cariche. Nessun parametro di produzione è stato modificato. I nuovi confronti tra livelli QM restano in corso.</p><p>Fonti: <a href="https://www.mdpi.com/2306-5729/11/1/18">AmberGO, articolo</a>; <a href="https://github.com/yingling-group/AmberGO">file originali</a>; <a href="https://orubel.mcmaster.ca/publications/2014/Stauffer_JCP_141_044705.pdf">Stauffer 2014</a>; <a href="https://pubs.acs.org/doi/10.1021/acs.jpca.4c01287">EspalomaCharge</a>; <a href="https://docs.openforcefield.org/projects/nagl-models/en/latest/">OpenFF NAGL</a>. Dati completi in summary.csv, dettagli per contesto e hash dei QM nei file adiacenti.</p></body></html>''')
text='\n'.join(parts).replace('<div class="scroll"><table><thead><tr><th>Componente','<div class="scroll"><table class="text"><thead><tr><th>Componente')
(P/'REPORT.html').write_text(text)
(P/'comparison.html').write_text(text.split('<h2>Cosa possiamo concludere</h2>')[0]+'<p>Confronto elettrostatico a sole cariche: non include LJ, dinamica o energie di adsorbimento. GAFF2-AIM non è classificato perché manca una mappatura validata sugli stessi frammenti.</p></body></html>')
(P/'MANIFEST.json').write_text(json.dumps({'files_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in P.iterdir() if f.is_file() and f.name!='MANIFEST.json'}},indent=2)+'\n')
print('Seven-context benchmark reproduced, AmberGO cap sensitivity quantified; no final force-field selection.')
