"""Describe existing RESP results; this script performs no new charge fitting."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--report', type=Path, required=True)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--replace-generated', action='store_true', help='Replace only this script’s derived outputs')
args = parser.parse_args()
report, source = args.report.resolve(), args.source.resolve()
inputs = [report/'context_HF_16/per_atom_charges.csv',
          report/'context_HF_16/background_carbon_distance_charges.csv',
          source/'l0_hf_all16/current_library.json',
          source/'l0_hf_all16/shared_library.json']
atomic = list(csv.DictReader(inputs[0].open()))
coo = [r for r in atomic if r['model_id'] == 'pyrene_COO' and r['role'].startswith('COO_')]
assert len(coo) == 4
roles = ['COO_anchor', 'COO_C', 'COO_O']
qa = {role: statistics.mean(float(r['charge_e']) for r in coo if r['role'] == role) for role in roles}
assert abs(qa['COO_C'] - .807216) < 1e-9
current = json.loads(inputs[2].read_text())['shared_role_charges_e']
l0 = json.loads(inputs[3].read_text())['shared_role_charges_e']
candidates = [('Libreria attuale', current), ('RESP HF atomico', qa), ('RESP HF L0', l0)]
rows = []
for label, values in candidates:
    region = values['COO_anchor'] + values['COO_C'] + 2*values['COO_O']
    for role in roles:
        rows.append({'candidate': label, 'model_or_fit': 'pyrene_COO' if label == 'RESP HF atomico' else 'shared_library',
                     'role': role, 'charge_e': values[role], 'multiplicity_in_region': 2 if role == 'COO_O' else 1,
                     'four_atom_region_charge_e': region, 'whole_pyrene_COO_target_e': -1,
                     'partial_region_charge_is_formal_charge': False})
basal = list(csv.DictReader(inputs[1].open()))
groups = {}
for row in basal:
    if row['role'] == 'scaffold_C':
        groups.setdefault((row['model_id'], row['nearest_functional_group_graph_distance']), []).append(float(row['charge_e']))
background = []
for (model, distance), values in groups.items():
    background.append({'model_id': model, 'nearest_group_graph_distance': distance,
                       'atom_count': len(values), 'mean_charge_e': statistics.mean(values),
                       'mean_absolute_charge_e': statistics.mean(abs(v) for v in values),
                       'min_charge_e': min(values), 'max_charge_e': max(values),
                       'scope': 'finite-fragment molecular HF-RESP; not physical charge-transfer measurement'})

def write_table(name, data):
    with (report/name).open('w' if args.replace_generated else 'x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)

write_table('COO_charge_diagnostic.csv', rows)
write_table('basal_C_summary_by_context.csv', background)
with tempfile.TemporaryDirectory(prefix='go_strategy_fig_') as cache:
    os.environ['MPLCONFIGDIR'] = cache
    os.environ['XDG_CACHE_HOME'] = cache
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    with plt.rc_context({'font.size': 10, 'svg.fonttype': 'none'}):
        fig, ax = plt.subplots(figsize=(8.2, 5.1), layout='constrained')
        fig.get_layout_engine().set(rect=(0, .14, 1, .86))
        for j, ((label, values), color) in enumerate(zip(candidates, ['#747474', '#0072B2', '#D55E00'])):
            x = np.arange(3)+(j-1)*.25
            y = [values[role] for role in roles]
            ax.bar(x, y, width=.23, color=color, label=label)
            for xi, v in zip(x, y):
                ax.text(xi, v+(.04 if v >= 0 else -.05), f'{v:+.3f}', ha='center',
                        va='bottom' if v >= 0 else 'top', fontsize=8.5)
        ax.axhline(0, color='#333333', linewidth=.8)
        ax.set_xticks(range(3), ['C di ancoraggio', 'C carbossilico', 'Ciascun O'])
        ax.set_ylim(-1.3, 2.18)
        ax.set_ylabel('Carica parziale (e)')
        ax.set_title('COO⁻: libreria attuale, RESP molecolare e RESP L0')
        ax.legend(loc='upper left', frameon=False, fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        fig.text(.5, .02, 'Il totale del modello è −1 e; il sottoinsieme COO può avere una somma parziale diversa.\nUn solo contesto fisico COO⁻: questa figura non stabilisce il vincolo responsabile della differenza.',
                 ha='center', fontsize=8.5)
        for extension in ['png', 'svg']:
            target = report/f'COO_charge_diagnostic.{extension}'
            if target.exists() and not args.replace_generated: raise FileExistsError(target)
            fig.savefig(target, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
manifest = {'scope': 'descriptive existing results; no new RESP fit or QM',
            'sources': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
            'rows_COO': len(rows), 'basal_context_distance_rows': len(background),
            'basal_atom_count': sum(r['atom_count'] for r in background),
            'error_bars': 'none; atoms and rigid rotations are not independent replicate fits',
            'final_parameter_validation': False}
with (report/'CHARGE_STRATEGY_DIAGNOSTIC_MANIFEST.json').open('w' if args.replace_generated else 'x') as handle:
    json.dump(manifest, handle, indent=2)
print(json.dumps({k: v for k, v in manifest.items() if k != 'sources'}, indent=2))
