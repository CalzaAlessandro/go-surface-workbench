import pathlib,json,csv,hashlib,html,collections,datetime
from go_param.compare import _tables
P=pathlib.Path(__file__).resolve().parent;R=P.parents[1]
roles=['scaffold_C','hydroxyl_C','hydroxyl_O','hydroxyl_H','epoxide_C','epoxide_O']
labels=['C senza legame diretto con O','C dell’OH','O dell’OH','H dell’OH','C dell’epossido','O dell’epossido']
original=_tables('current_gaff2')[0];hier=_tables('hiergo_opls')[0]
for file in ['HYD_gaff2.lib','EPX_gaff2.lib','GGG_gaff2.lib']:
 a=R/'references/local_snapshot/current'/file;b=R.parent/'VS_POLL_WORKFLOW/resources/LIB'/file
 assert a.read_bytes()==b.read_bytes(),'live library differs from audited source'
fit={r['role']:float(r['charge_e']) for r in csv.DictReader((R/'reports/CHARGE_COMPARISON_20261103_attempt_001/L0_shared_role_charges.csv').open()) if r['fit']=='l0_hf_all16'}
aim=collections.defaultdict(list);source_files=[]
for f in [R/'references/gaff2_aim/GO_18/GO_sheet_1_8.top',R/'references/gaff2_aim/GO_30/GO_sheet_3_0.top']:
 atoms={};adj=collections.defaultdict(set);sec='';source_files.append(f)
 for line in f.read_text().splitlines():
  line=line.split(';')[0].strip()
  if not line:continue
  if line.startswith('['):sec=line.strip('[] ').strip();continue
  s=line.split()
  if sec=='atoms':atoms[int(s[0])]={'type':s[1],'q':float(s[6])}
  elif sec=='bonds':a,b=map(int,s[:2]);adj[a].add(b);adj[b].add(a)
 for i,a in atoms.items():
  typ=a['type'];role=None
  if typ.startswith('oh'):role='hydroxyl_O'
  elif typ.startswith('os'):role='epoxide_O'
  elif typ.lower().startswith('h'):role='hydroxyl_H'
  elif typ.startswith('c'):
   ns=[atoms[j]['type'] for j in adj[i]]
   role='hydroxyl_C' if any(t.startswith('oh') for t in ns) else 'epoxide_C' if any(t.startswith('os') for t in ns) else 'scaffold_C'
  if role:aim[role].append(a['q'])
amber=dict(zip(roles,[0.,.16,-.57,.41,.16,-.32]))
models=[('AmberGO',amber),('HierGO / OPLS',{r:hier[r][0] for r in roles}),('Tua LIB originale',{r:original[r][0] for r in roles}),('RESP L0 preliminare',{r:fit[r] for r in roles})]
rows=[]
for role,label in zip(roles,labels):
 row={'role':role,'label':label,**{name:data[role] for name,data in models},'GAFF2_AIM_min':min(aim[role]),'GAFF2_AIM_max':max(aim[role])};rows.append(row)
sources=source_files+[R/'references/hiergo_17_charge_rows.json',R/'references/library_charges.json',R/'references/go_prior_art_20261104/AmberGO_native_audit_20261104/FINAL_MOL2_CHARGE_AUDIT.json',R/'reports/CHARGE_COMPARISON_20261103_attempt_001/L0_shared_role_charges.csv']
(P/'DATA.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'units':'elementary charge e','meaning':'Native assignment comparison, not an accuracy or adsorption-energy ranking','AIM_scope':'Ranges of native atom types present in two original periodic GO_18 and GO_30 topologies; graph-based role assignment; not confidence intervals','RESP_L0':'Shared-role HF-RESP fit to 16 finite contexts with imposed constraints; preliminary, not production validated','rows':rows,'source_sha256':{str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}},indent=2)+'\n')
fmt=lambda q:('+' if q>0 else '−' if q<0 else '')+f'{abs(q):.3f}'
parts=['''<!doctype html><html lang="it"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Cariche GO e ruolo di GO-MACE</title><style>body{font:15px/1.5 var(--font-sans,system-ui);color:var(--foreground,#24313b)}h2{font-size:19px}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid var(--border,#ccd3d8);padding:8px;text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left;white-space:normal}svg{width:100%;height:auto;display:block}.usage{table-layout:fixed}.usage td,.usage th{text-align:left;white-space:normal;vertical-align:top}.usage th:first-child{width:22%}svg text{fill:var(--foreground,#24313b);font:12px system-ui}a{color:var(--primary,#146b92)}p{margin:12px 0}.muted{color:var(--muted-foreground,#68747c)}</style><body><h2>Cariche atomiche dei gruppi comuni</h2><p>Valori in e. La tua LIB originale e il nostro fit RESP preliminare sono riportati separatamente. Gli intervalli GAFF2-AIM descrivono i due fogli nativi esaminati: non sono barre d’errore.</p><div class="scroll"><table><thead><tr><th>Ambiente</th>''']
parts += [''.join('<th>'+html.escape(name)+'</th>' for name,data in models)+'<th>GAFF2-AIM, min…max</th></tr></thead><tbody>']
for r in rows:
 parts.append('<tr><td>'+html.escape(r['label'])+'</td>'+''.join('<td>'+fmt(r[name])+'</td>' for name,data in models)+'<td>'+fmt(r['GAFF2_AIM_min'])+' … '+fmt(r['GAFF2_AIM_max'])+'</td></tr>')
parts.append('</tbody></table></div><h2>Confronto dei soli ossigeni</h2><p class="muted">Scala comune da −0,75 a 0 e. Una carica più negativa non implica, da sola, maggiore accuratezza o adsorbimento.</p>')
# Directly labelled point/range chart, no averaging or inferred uncertainty.
svg=['<svg viewBox="0 0 720 280" role="img" aria-label="Cariche degli ossigeni OH ed epossido: punti per modelli fissi, intervalli dei tipi nativi per GAFF2-AIM">']
for panel,(role,title) in enumerate([('hydroxyl_O','O dell’OH'),('epoxide_O','O dell’epossido')]):
 left=145+panel*300;right=left+235;y0=43
 x=lambda q:left+(q+.75)/.75*(right-left)
 svg.append(f'<text x="{(left+right)/2}" y="18" text-anchor="middle">{title}</text>')
 for tick in [-.75,-.5,-.25,0]:
  xx=x(tick);svg.append(f'<line x1="{xx}" y1="28" x2="{xx}" y2="232" stroke="var(--border,#ccd3d8)"/><text x="{xx}" y="251" text-anchor="middle">{tick:g}</text>')
 for idx,(name,data) in enumerate(models):
  y=y0+idx*40
  if panel==0:svg.append(f'<text x="2" y="{y+4}">{html.escape(name)}</text>')
  xx=x(data[role]);svg.append(f'<circle cx="{xx}" cy="{y}" r="4" fill="var(--chart-{idx+1},#28678c)"/><text x="{xx+9}" y="{y+4}">{fmt(data[role])}</text>')
 y=y0+4*40;lo=min(aim[role]);hi=max(aim[role])
 if panel==0:svg.append(f'<text x="2" y="{y+4}">GAFF2-AIM</text>')
 svg.append(f'<line x1="{x(lo)}" x2="{x(hi)}" y1="{y}" y2="{y}" stroke="var(--chart-5,#386e5b)" stroke-width="6"/><text x="{right}" y="{y+18}" text-anchor="end">{fmt(lo)}…{fmt(hi)}</text>')
svg.append('</svg>');parts.append(''.join(svg))
parts.append('''<p>Questa è una comparazione descrittiva di assegnazioni ottenute con schemi e parametrizzazioni diversi. Il fit RESP L0 è vincolato e preliminare; GAFF2-AIM distribuisce carica anche sul carbonio senza O direttamente legato. Per una classifica servono confronti sullo stesso sistema e con i rispettivi parametri LJ/bonded.</p><h2>GO-MACE: tre usi diversi</h2><table class="usage"><thead><tr><th>Obiettivo</th><th>Uso possibile</th></tr></thead><tbody><tr><td>Predire cariche RESP</td><td>Il GO-MACE pubblicato non fornisce questa uscita. Servono potenziali/densità/dipoli QM o etichette RESP; energie e forze non separano univocamente cariche, LJ e termini bonded.</td></tr><tr><td>Ricavare parametri GAFF</td><td>Possibile nuovo fit di energie/forze/curvature dopo aver fissato cariche e LJ, entro il dominio del potenziale e con verifica QM indipendente. Non è una conversione diretta del file del modello in frcmod.</td></tr><tr><td>Calibrare SQE per AMBER</td><td>GO-MACE può fornire geometrie da campionare. Un sottoinsieme curato viene ricalcolato con il riferimento elettrostatico concordato; SQE viene calibrato contro quel riferimento, mantenendo il totale Q chimico.</td></tr></tbody></table><p>La superficie finale conserva cariche fisse durante MD/GBSA. Il normale GO-MACE, studiato per GO e trasformazioni in vuoto, non è già validato per la distribuzione elettrostatica di COO⁻, fogli a diversa carica o interazioni GO–proteina in soluzione.</p><p>Fonti: <a href="https://github.com/yingling-group/AmberGO">AmberGO</a>; <a href="https://github.com/IFM-molecular-simulation-group/HierGO/blob/master/Oplsatom2">HierGO OPLS</a>; <a href="https://github.com/arvpinto/GO_sheets">GAFF2-AIM</a>; <a href="https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202410088">GO-MACE</a>. Dati e hash nel file DATA.json accanto a questo report.</p></body></html>''')
(P/'REPORT.html').write_text('\n'.join(parts))
print('Verified live original libraries; wrote charge comparison and documented native AIM ranges.')
