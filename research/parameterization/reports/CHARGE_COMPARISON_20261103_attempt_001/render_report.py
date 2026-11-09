"""Create descriptive figures/tables from verified, completed charge results."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import tempfile

parser=argparse.ArgumentParser()
parser.add_argument('--input',type=Path,required=True)
parser.add_argument('--c54-input',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
source=args.input.resolve();c54=args.c54_input.resolve();out=args.output.resolve()
out.mkdir(parents=True,exist_ok=False)
cache=tempfile.TemporaryDirectory(prefix='go_charge_mpl_')
os.environ['MPLCONFIGDIR']=cache.name
os.environ.setdefault('XDG_CACHE_HOME',cache.name)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from go_param.context import analyze_dataset,write_outputs


def read(path):return json.loads(Path(path).read_text())
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def table(path,rows):
    with Path(path).open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


raw=list(csv.DictReader((source/'cross_method_ESP_scores.csv').open(newline='')))
assert len(raw)==1728
groups={}
for row in raw:groups.setdefault((row['model_id'],row['reference_ESP'],row['candidate']),[]).append(row)
averages=[];lookup={}
for key,rows in groups.items():
    assert len(rows)==3 and {r['orientation'] for r in rows}=={'1','2','3'}
    valid=[float(r['relative_rmse']) for r in rows if r['status']=='matched_charge_control']
    assert len(valid) in (0,3)
    item=dict(zip(['model_id','reference_ESP','candidate'],key))
    item.update(orientations=3,mean_RRMS=statistics.mean(valid) if valid else None,
                min_RRMS=min(valid) if valid else None,max_RRMS=max(valid) if valid else None,
                availability='available' if valid else 'unavailable',
                reason=';'.join(sorted({r['reason'] for r in rows if r['reason']})))
    averages.append(item);lookup[key]=item['mean_RRMS']
table(out/'RRMS_by_physical_context.csv',averages)
hf='HF/6-31G*';dft='wB97XD/6-311++G(d,p)'
models=list(dict.fromkeys(r['model_id'] for r in raw))
assert len(models)==16
keys=['atomic_hf_on_dft','atomic_dft','l0_hf_all16','l0_hf_train_C16','current_full_block_library','am1_bcc_constrained_rounding']
labels=['RESP HF\natomico','RESP DFT\natomico','L0 HF\n16 casi','L0 HF\nC16','Libreria\nattuale','AM1-BCC*']
short=['C16 · pristine','C16 · 2OH cis','C16 · 2OH trans','C16 · 2OH + EPX','C16 · carbonili','C16 · COOH','C16 · COO−',
       'C24 · pristine','C24 · EPX singolo','C24 · EPH bordo','C24 · EPT bordo','C24 · 2OH vicini',
       'C24 · 2OH lontani','C24 · 2EPX vicini','C24 · 2OH + EPX vicino','C24 · 2OH + EPX lontano']
finite=[lookup[m,g,k] for m in models for g in [hf,dft] for k in keys if lookup[m,g,k] is not None]
vmax=math.ceil(max(finite)*10)/10
for grid,name in [(hf,'HF_ESP'),(dft,'DFT_ESP')]:
    matrix=np.array([[lookup[m,grid,k] if lookup[m,grid,k] is not None else np.nan for k in keys] for m in models])
    with plt.rc_context({'font.size':8.5,'svg.fonttype':'none','pdf.fonttype':42}):
        fig,ax=plt.subplots(figsize=(7.6,6.9),layout='constrained')
        cmap=matplotlib.colormaps['Blues'].with_extremes(bad='#eeeeee')
        img=ax.imshow(matrix,vmin=0,vmax=vmax,cmap=cmap,aspect='auto')
        ax.set_xticks(range(len(keys)),labels);ax.xaxis.tick_top()
        ax.tick_params(axis='x',length=0,pad=8)
        ax.set_yticks(range(len(models)),short);ax.tick_params(axis='y',length=0)
        for i in range(len(models)):
            for j in range(len(keys)):
                v=matrix[i,j]
                ax.text(j,i,'N/D' if np.isnan(v) else f'{v:.3f}',ha='center',va='center',
                        color='white' if np.isfinite(v) and v>vmax*.58 else '#111111',fontsize=8.2)
        ax.axhline(6.5,color='#666666',linewidth=1)
        ax.set_xticks(np.arange(-.5,len(keys),1),minor=True)
        ax.set_yticks(np.arange(-.5,len(models),1),minor=True)
        ax.grid(which='minor',color='white',linewidth=.7);ax.tick_params(which='minor',length=0)
        fig.colorbar(img,ax=ax,shrink=.78,label='Errore relativo ESP (RRMS)')
        fig.suptitle('Riferimento ESP: '+('HF/6-31G*' if grid==hf else 'ωB97X-D/6-311++G(d,p)'),fontsize=11)
        fig.savefig(out/(name+'.png'),dpi=300,facecolor='white')
        fig.savefig(out/(name+'.svg'),facecolor='white')
        plt.close(fig)
common_holdout=[m for m in models if m.startswith('C24_') and lookup[m,hf,'l0_hf_train_C16'] is not None]
assert len(common_holdout)==7
heldout={k:statistics.median(lookup[m,hf,k] for m in common_holdout) for k in ['atomic_hf_on_dft','l0_hf_all16','l0_hf_train_C16','current_full_block_library']}
summary={'HF_contexts':16,'DFT_contexts_before_C54':16,'AM1_BCC_parameter_covered_contexts':12,
         'AM1_BCC_exclusions':read(source/'AM1_BCC_EXCLUSIONS.json'),
         'HF_ESP_all16_median_RRMS':{k:statistics.median(lookup[m,hf,k] for m in models) for k in ['atomic_hf_on_dft','atomic_dft','l0_hf_all16','current_full_block_library']},
         'C24_held_out_models':common_holdout,'C24_held_out_HF_ESP_median_RRMS':heldout,
         'three_orientations_are_independent_chemical_samples':False,'final_parameter_validation':False}
entries_hf=[];entries_dft=[];differences=[]
for mid in models:
    folders={method:source/'atomic_inputs'/mid/method for method in ['HF','DFT']}
    for method,entries in [('HF',entries_hf),('DFT',entries_dft)]:
        p=folders[method];entries.append({'model_path':str(p/'model.json'),'charges_path':str(p/'charges.json'),'fit_evidence_path':str(p/'FIT_COMPLETE.json')})
    model=read(folders['HF']/'model.json');qh=read(folders['HF']/'charges.json')['atoms'];qd=read(folders['DFT']/'charges.json')['atoms']
    role_diffs={}
    for atom,h,d in zip(model['atoms'],qh,qd):
        assert atom['id']==h['atom_id']==d['atom_id']
        role_diffs.setdefault(atom['role'],[]).append(h['charge_e']-d['charge_e'])
    for role,values in role_diffs.items():
        differences.append({'model_id':mid,'role':role,'atom_count':len(values),
                            'mean_HF_minus_DFT_charge_e':statistics.mean(values),'min_difference_e':min(values),'max_difference_e':max(values)})
for p in sorted((c54/'outputs').glob('*/fits/attempt_001')):
    entries_dft.append({'model_path':str(p/'model.json'),'charges_path':str(p/'charges.json'),'fit_evidence_path':str(p/'FIT_COMPLETE.json')})
assert len(entries_dft)==18
write_outputs(analyze_dataset(entries_hf),out/'context_HF_16')
write_outputs(analyze_dataset(entries_dft),out/'context_DFT_18')
table(out/'matched_HF_DFT_role_differences.csv',differences)
summary['DFT_contexts_including_recovered_C54']=18
summary['C54_ESP_input_IO_regression']=read(c54/'resp_header_io/IO_REGRESSION_GATE.json')
save(out/'SUMMARY_RESULTS.json',summary)
lib_rows=[]
for name,record in read(source/'SYNTHESIS_COMPLETE.json')['L0_fits'].items():
    for role,q in record['shared_roles'].items():
        lib_rows.append({'fit':name,'training_contexts':record['models'],'role':role,'charge_e':q,'production_validated':False})
table(out/'L0_shared_role_charges.csv',lib_rows)
caption='Media delle tre orientazioni rigide di ciascuna geometria. N/D = copertura assente; non zero. I fit atomici sono valutati anche sulle griglie usate per il loro fit: questo non classifica l’accuratezza del force field.'
html='''<!doctype html><html lang="it"><meta charset="utf-8"><style>
body{font-family:var(--font-sans);color:var(--foreground);margin:0}button{font:inherit;color:var(--foreground);background:var(--secondary);border:1px solid var(--border);border-radius:var(--radius);padding:7px 12px;cursor:pointer}button[aria-pressed=true]{outline:2px solid var(--ring)}.tabs{display:flex;gap:12px;flex-wrap:wrap;margin:8px 0}.figure{width:100%;max-width:760px}.figure img{display:block;width:100%;height:auto}.note{font-size:12px;line-height:1.45;color:var(--muted-foreground);margin:8px 0}
</style><div class="tabs"><button aria-pressed="true" id="bh" onclick="show('HF')">Griglia HF/6-31G*</button><button aria-pressed="false" id="bd" onclick="show('DFT')">Griglia DFT</button></div>
<div class="figure"><img id="hf" src="'''+str(out/'HF_ESP.png')+'''" alt="Errore ESP per sedici strutture e sei modelli di carica; ogni cella riporta il valore numerico e i casi non coperti sono N/D."><img id="dft" style="display:none" src="'''+str(out/'DFT_ESP.png')+'''" alt="Stesso confronto sulla griglia DFT, con identica scala dei colori."></div><p class="note">'''+caption+'''</p><p class="note">L0 HF C16: fit su 7 pireni; le righe C24 sono casi esclusi dall’addestramento. Due gruppi terminali non erano presenti nel training. AM1-BCC*: lieve correzione di arrotondamento della carica totale, distinta dai valori nativi conservati; 12 casi con parametri BCC completi.</p><script>function show(which){document.getElementById('hf').style.display=which==='HF'?'block':'none';document.getElementById('dft').style.display=which==='DFT'?'block':'none';document.getElementById('bh').setAttribute('aria-pressed',which==='HF');document.getElementById('bd').setAttribute('aria-pressed',which==='DFT')}</script></html>'''
(out/'comparison.html').write_text(html)
save(out/'FIGURE_MANIFEST.json',{'source':str(source/'cross_method_ESP_scores.csv'),
    'source_sha256':hashlib.sha256((source/'cross_method_ESP_scores.csv').read_bytes()).hexdigest(),
    'aggregation':'arithmetic mean of 3 rigid orientations per physical geometry, not replicates',
    'missing_data':'unavailable shown as N/D, not zero','normalization':{'RRMS_scale_min':0,'common_RRMS_scale_max':vmax},
    'formats':['PNG 300 dpi','SVG','HTML with selectable reference grid'],'journal_specific_validation':False})
print(json.dumps(summary,indent=2,ensure_ascii=False))
