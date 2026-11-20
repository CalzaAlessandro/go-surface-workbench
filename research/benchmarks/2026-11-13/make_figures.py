"""Rebuild README figures exclusively from the distributed benchmark tables."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

HERE=Path(__file__).resolve().parent
WORKBENCH=HERE.parents[2]
DATA=HERE/'data'
OUT=WORKBENCH/'docs/assets/benchmarks/2026-11-13'
LABELS={
 'GAFF_AM1BCC_closed':'GAFF 1.81 + AM1-BCC', 'GAFF_HF_RESP_closed':'GAFF 1.81 + HF-RESP',
 'GAFF2_AM1BCC_closed':'GAFF2 2.2.20 + AM1-BCC','GAFF2_HF_RESP_closed':'GAFF2 2.2.20 + HF-RESP',
 'OpenFF_2.0_AM1BCC_closed':'OpenFF 2.0 + AM1-BCC','Sage_2.3_AshGC':'Sage 2.3 + AshGC',
 'Grappa_1.4_OpenFF2_AM1BCC_closed':'Grappa 1.4 + AM1-BCC',
 'Espaloma_0.3.2_OpenFF2_AM1BCC_closed':'Espaloma 0.3.2 + AM1-BCC',
 'Espaloma_0.3.2_OpenFF2_EspalomaCharge':'Espaloma 0.3.2 + EspalomaCharge',
 'QForce_Hessian_only_GAFF2_LJ_AM1BCC':'Q-Force (Hessian-only, in-sample)*'}
ORDER=list(LABELS)
COLOURS=['#666666','#666666','#0072B2','#0072B2','#CC79A7','#E69F00','#009E73','#D55E00','#D55E00','#56B4E9']


def read(name):
    with (DATA/name).open() as handle:return list(csv.DictReader(handle))


def export(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=300,facecolor='white',metadata={'Software':'Matplotlib '+matplotlib.__version__})
    fig.savefig(OUT/(name+'.svg'),facecolor='white',metadata={'Date':None})
    plt.close(fig)


def mm_figure(rows,key,label,title,filename,scale=1):
    fig,ax=plt.subplots(figsize=(9.8,6.3));fig.subplots_adjust(left=.41,right=.97,bottom=.18,top=.83)
    names=sorted({r['model_id'] for r in rows});medians={}
    for method in ORDER:
        subset={r['model_id']:r for r in rows if r['method']==method};assert set(subset)==set(names)
        medians[method]=np.median([float(subset[name][key])*scale for name in names])
    ordered=sorted(ORDER,key=medians.get)
    for y,method in enumerate(ordered):
        value=medians[method];colour=COLOURS[ORDER.index(method)]
        ax.barh(y,value,height=.65,color=colour,zorder=3)
        ax.text(value+max(medians.values())*.025,y,f'{value:.2f}',va='center',fontsize=11,weight='bold')
    ax.set_yticks(range(len(ordered)),[LABELS[x] for x in ordered]);ax.invert_yaxis();ax.set_xlim(0,max(medians.values())*1.16)
    ax.set_xlabel(label);ax.grid(axis='x',alpha=.22);ax.set_axisbelow(True)
    fig.suptitle(title,x=.04,y=.97,ha='left',weight='bold')
    fig.text(.04,.90,'Median error across 6 finite fragments · Lower is better',fontsize=11,color='#333333')
    fig.text(.04,.025,'Common DFT geometries. Individual cases: PER_CASE_RESULTS.md.\n*Q-Force is fitted in-sample to each fragment\'s QM Hessian (no torsion scans; GAFF2 LJ/AM1-BCC), not a transferable preset.',fontsize=8.8,color='#555555')
    export(fig,filename)


def charge_figure():
    rows=read('charge_per_case.csv');names=['C54_EPX2_far_s20261407','C54_EPX4_loading_s20261409','C54_OH2_loading_s20261404','C54_OH4_loading_s20261405','C54_OH6_loading_s20261406','C54_COO1_loading_s20261417']
    methods=['L0_strict40','Local_BCC_ESP_strict40']
    labels=['L0 baseline','Local BCC fitted to ESP']
    colours=['#777777','#0072B2']
    fig,ax=plt.subplots(figsize=(9.5,5));fig.subplots_adjust(left=.09,right=.97,bottom=.19,top=.75)
    for j,(method,label,colour) in enumerate(zip(methods,labels,colours)):
        values=[100*float(next(r['median_ESP_RRMS'] for r in rows if r['model_id']==name and r['method']==method)) for name in names]
        bars=ax.bar(np.arange(6)+(j-.5)*.34,values,width=.32,color=colour,label=label,zorder=3)
        ax.bar_label(bars,fmt='%.1f',padding=3,fontsize=10)
    ax.axvline(4.5,color='#999999',ls=':',lw=1);ax.set_xticks(range(6),['EPX2','EPX4','OH2','OH4','OH6','COO⁻'])
    ax.set_ylabel('ESP relative RMSE (%)');ax.set_ylim(0,28);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    ax.set_title('Local charge correction: 5 of 6 C54 cases improve',loc='left',weight='bold',pad=55)
    ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),ncol=2,frameon=False,fontsize=9)
    ax.text(3,26,'OH4 worsens',ha='center',fontsize=10,color='#9B2323',weight='bold')
    fig.text(.09,.045,'Lower error is better. Median of 3 orientations per molecule. OH4: 21.44% → 22.50%.\n40 C16/C24 training molecules; these 6 C54 molecules excluded from fitting and tuning.',fontsize=9,color='#555555')
    export(fig,'charge_transfer_C54')


def pattern_figure():
    examples=json.loads((DATA/'pattern_examples.json').read_text());rows=read('pattern_per_case.csv')
    fig=plt.figure(figsize=(10,5.4));grid=fig.add_gridspec(2,3,height_ratios=[2,1.1],hspace=.35)
    for i,example in enumerate(examples):
        ax=fig.add_subplot(grid[0,i]);atoms={a['id']:a for a in example['atoms']};carbon={i:a for i,a in atoms.items() if a['element']=='C'}
        edges=[[carbon[b['a']]['position'][:2],carbon[b['b']]['position'][:2]] for b in example['bonds'] if b['a'] in carbon and b['b'] in carbon]
        ax.add_collection(LineCollection(edges,colors='#AAAAAA',linewidths=.6))
        xy=np.array([a['position'][:2] for a in carbon.values()]);ax.scatter(xy[:,0],xy[:,1],s=5,color='#555555')
        for group,colour,marker in [('OH','#0072B2','o'),('epoxide','#D55E00','^')]:
            oxygen=[a['position'][:2] for a in atoms.values() if a['element']=='O' and any(a['id'] in g.get('atom_ids',[]) for g in example['groups'] if g['type']==group)]
            # Native group records expose atom IDs through atom metadata rather than a shared group atom_ids field.
            if not oxygen:
                group_ids={g['id'] for g in example['groups'] if g['type']==group}
                oxygen=[a['position'][:2] for a in atoms.values() if a['element']=='O' and a.get('metadata',{}).get('group_id') in group_ids]
            if oxygen:
                oxy=np.array(oxygen);ax.scatter(oxy[:,0],oxy[:,1],s=45,color=colour,marker=marker,label=group,zorder=3)
        ax.autoscale();ax.set_aspect('equal');ax.set_title(example['mode'].capitalize());ax.set_xlabel('x (Å)')
        if i==0:ax.set_ylabel('y (Å)')
    ax=fig.add_subplot(grid[1,:])
    for i,mode in enumerate(['random','island','dispersed']):
        values=np.array([float(r['mean_anchor_pair_graph_distance']) for r in rows if r['mode']==mode]);assert len(values)==16
        value=float(np.mean(values));ax.bar(i,value,width=.55,color=['#666666','#0072B2','#D55E00'][i],zorder=3)
        ax.text(i,value+.2,f'{value:.2f}',ha='center',fontsize=11,weight='bold')
    ax.set_xticks(range(3),['Random','Island','Dispersed']);ax.set_ylabel('Group-pair distance\n(C-scaffold bonds)');ax.set_ylim(0,11);ax.grid(axis='y',alpha=.2)
    handles,labels=fig.axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.935),ncol=2,frameon=False,fontsize=10)
    fig.suptitle('Controlled oxidation patterns at identical composition',weight='bold',y=.99)
    fig.text(.08,.012,'C130H36O8; H omitted in top views. Bars = mean group-pair distance averaged over 16 seeds.\nShorter distance means more clustering. Starting structures, not equilibrium ensembles.',fontsize=9,color='#555555')
    fig.subplots_adjust(left=.09,right=.98,bottom=.15,top=.85)
    export(fig,'oxidation_patterns')


def charge_overview():
    rows=read('charge_per_case.csv')
    specs=[
        ('Local_BCC_ESP_strict40','Local BCC fitted to ESP',6,'#0072B2'),
        ('SQE_ESP_strict40','SQE fitted to ESP',6,'#0072B2'),
        ('L0_strict40','L0 baseline',6,'#0072B2'),
        ('HierGO_OPLS_COO_HF_ESP_extension','HierGO OPLS + COO⁻ HF-ESP ext.',6,'#009E73'),
        ('hiergo_opls','HierGO OPLS (5/6 neutral only)*',5,'#56B4E9'),
        ('current_gaff2','Current GAFF2 block table',6,'#E69F00'),
        ('AshGC','AshGC 1.0 (AM1-BCC target)',6,'#CC79A7'),
        ('EspalomaCharge','EspalomaCharge (ELF10 target)',6,'#CC79A7'),
        ('Atomic_RESP_in_sample_control','Molecular RESP: in-sample control',6,'#888888'),
    ]
    fig,ax=plt.subplots(figsize=(9.8,6.0));fig.subplots_adjust(left=.41,right=.97,bottom=.21,top=.81)
    for y,(method,label,expected_n,colour) in enumerate(specs):
        subset=sorted([r for r in rows if r['method']==method],key=lambda r:r['model_id']);assert len(subset)==expected_n
        values=[100*float(r['median_ESP_RRMS']) for r in subset]
        med=np.median(values)
        ax.barh(y,med,height=.65,color=colour,zorder=3)
        ax.text(med+3,y,f'{med:.1f}%',va='center',fontsize=11,weight='bold')
    ax.set_yticks(range(len(specs)),[s[1] for s in specs]);ax.invert_yaxis();ax.set_xlim(0,250)
    ax.set_xlabel('Median ESP relative RMSE (%)');ax.grid(axis='x',alpha=.2)
    fig.suptitle('C54 charge comparisons on the same HF ESP target',x=.04,y=.98,ha='left',weight='bold')
    fig.text(.04,.905,'Matched-charge C54 evaluations · Median of per-molecule orientation medians · Lower is better',fontsize=10,color='#333333')
    fig.text(.04,.025,'*HierGO OPLS lacks COO⁻ (covers 5 neutral C54 only). Unclosed GOPY raw (-0.45 to -1.48 e) and partial Stauffer/GAFF2-AIM\nare excluded from matched-charge scoring (see charge_coverage.csv). Different targets test HF ESP, not condensed-phase accuracy.',fontsize=8.6,color='#555555')
    export(fig,'charge_backends_C54')


def per_case_tables():
    lines=['# Individual benchmark cases','','All declared cases are retained. The overview bars show medians for QM/charge comparisons and means for pattern separation. No confidence intervals are inferred.','']
    names=['pyrene_pristine','pyrene_diol_cis','pyrene_diol_trans','pyrene_diol_epoxide','pyrene_COOH','pyrene_COO']
    rows=read('mm_per_case.csv')
    for title,key,scale in [('Force RMSE, kcal/(mol·Å)','force_component_RMSE_kcal_mol_A',1),('Cartesian-Hessian error, %','hessian_relative_Frobenius_error',100)]:
        lines.extend(['## '+title,'','| Combination | Pyrene | Diol cis | Diol trans | Diol + EPX | COOH | COO⁻ |','|---|---:|---:|---:|---:|---:|---:|'])
        for method in ORDER:
            values=[float(next(r[key] for r in rows if r['model_id']==name and r['method']==method))*scale for name in names]
            lines.append('| '+LABELS[method]+' | '+' | '.join(f'{v:.3f}' for v in values)+' |')
        lines.extend(['','`*` Q-Force (Hessian-only, in-sample) is fitted directly to each molecule\'s QM Hessian without torsion scans, using GAFF2 LJ/AM1-BCC.',''])
    lines.extend(['## C54 ESP relative RMSE, %','','Median of three orientations within each molecule; molecular RESP was fitted on the same grids. HierGO OPLS covers the five neutral C54 molecules only (`COO⁻` absent); `HierGO_OPLS_COO_HF_ESP_extension` adds a 2-parameter `COO⁻` block fitted on 10 non-C54 `COO⁻` fragments.','',
                  '| Assignment | EPX2 | EPX4 | OH2 | OH4 | OH6 | COO⁻ |','|---|---:|---:|---:|---:|---:|---:|'])
    rows=read('charge_per_case.csv')
    ids=['C54_EPX2_far_s20261407','C54_EPX4_loading_s20261409','C54_OH2_loading_s20261404','C54_OH4_loading_s20261405','C54_OH6_loading_s20261406','C54_COO1_loading_s20261417']
    for method in ['L0_strict40','Local_BCC_ESP_strict40','SQE_ESP_strict40','HierGO_OPLS_COO_HF_ESP_extension','hiergo_opls','current_gaff2','AshGC','EspalomaCharge','Atomic_RESP_in_sample_control']:
        cells=[]
        for name in ids:
            match=next((r for r in rows if r['model_id']==name and r['method']==method),None)
            cells.append(f"{100*float(match['median_ESP_RRMS']):.2f}" if match else 'not covered')
        lines.append('| '+method+' | '+' | '.join(cells)+' |')
    lines.extend(['','## Published charge-table coverage and charge-conservation audit on C54','',
                  '| Table / extension | Covered C54 molecules | Matched formal charge | Missing roles / scope limitation |',
                  '|---|---:|---:|---|'])
    cov=read('charge_coverage.csv')
    pub_order=['current_gaff2','hiergo_opls','HierGO_OPLS_COO_HF_ESP_extension','gopy_raw','stauffer_esp','stauffer_sesp','pinto_gaff2_aim']
    for method in pub_order:
        sub=[r for r in cov if r['method']==method]
        complete_n=sum(r['complete']=='True' for r in sub)
        matched_n=sum(r['charge_matches_formal']=='True' for r in sub)
        missing=sorted({role for r in sub if r['unsupported_roles'] for role in r['unsupported_roles'].split('|')})
        if method=='gopy_raw':
            errs=[float(r['charge_sum_error_e']) for r in sub if r['charge_sum_error_e']]
            note=f"Nonzero neutral motif sums ({min(errs):+.3f} to {max(errs):+.3f} e); missing: {', '.join(missing)}"
        elif method=='pinto_gaff2_aim':
            note='Periodic DDEC6 + AIM LJ (`GO_18`/`GO_30`); no automatic finite-fragment transfer'
        elif missing:
            note='Missing: '+', '.join(missing)
        elif method=='HierGO_OPLS_COO_HF_ESP_extension':
            note='Charge-only extension; OPLS-AA lacks aromatic `CA-C_3-O2` / `CA-CA-C_3-O2` bonded terms'
        else:
            note='Complete C54 role coverage'
        lines.append(f'| {method} | {complete_n}/6 | {matched_n}/6 | {note} |')
    lines.extend(['','## Pattern separation across 16 seeds','','Distance is the shortest C-scaffold path between group anchor sets, averaged over group pairs.','',
                  '| Pattern | Mean across seeds | Minimum seed value | Maximum seed value |','|---|---:|---:|---:|'])
    rows=read('pattern_per_case.csv')
    for mode in ['random','island','dispersed']:
        v=[float(r['mean_anchor_pair_graph_distance']) for r in rows if r['mode']==mode]
        lines.append(f'| {mode} | {np.mean(v):.3f} | {min(v):.3f} | {max(v):.3f} |')
    (HERE/'PER_CASE_RESULTS.md').write_text('\n'.join(lines)+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'go-20261113'}):
        rows=read('mm_per_case.csv')
        mm_figure(rows,'force_component_RMSE_kcal_mol_A','Force component RMSE [kcal/(mol·Å)]','MM forces versus the common QM reference','QM_force_comparison')
        mm_figure(rows,'hessian_relative_Frobenius_error','Relative Cartesian-Hessian error (%)','Local curvature versus the common QM reference','QM_hessian_comparison',100)
        charge_figure();charge_overview();pattern_figure();per_case_tables()
    manifest={str(p.relative_to(WORKBENCH)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.suffix in ['.png','.svg']}
    (HERE/'FIGURE_MANIFEST.json').write_text(json.dumps({'data_files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.is_file()},'figures':manifest,'matplotlib':matplotlib.__version__},indent=2)+'\n')
    print('Exported',len(manifest),'figure files to',OUT)


if __name__=='__main__':main()
