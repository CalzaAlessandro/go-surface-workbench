"""Add coverage and matched-subset summaries to the descriptive charge report.

Run once after render_report.py. Scores are descriptive, not FF rankings.
"""
import argparse
import csv
import json
import statistics
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--report', type=Path, required=True)
args = parser.parse_args()
root = args.report
rows = list(csv.DictReader((root/'RRMS_by_physical_context.csv').open()))
assert len(rows) == 16 * 2 * 18
lookup = {(r['model_id'], r['reference_ESP'], r['candidate']): r for r in rows}
models = list(dict.fromkeys(r['model_id'] for r in rows))
grids = list(dict.fromkeys(r['reference_ESP'] for r in rows))
candidates = list(dict.fromkeys(r['candidate'] for r in rows))
coverage, scores = [], []
for candidate in candidates:
    subset = [m for m in models if lookup[m, grids[0], candidate]['availability'] == 'available']
    assert all(lookup[m, grids[1], candidate]['availability'] == 'available' for m in subset)
    coverage.append({'candidate': candidate, 'covered_contexts': len(subset),
                     'total_contexts': len(models), 'covered_model_ids': json.dumps(subset),
                     'missing_means_zero': False})
    if not subset:
        continue
    for grid in grids:
        comparisons = list(dict.fromkeys([candidate, 'atomic_hf_on_dft', 'l0_hf_all16',
                                          'current_full_block_library']))
        for evaluated in comparisons:
            values = [float(lookup[m, grid, evaluated]['mean_RRMS']) for m in subset]
            scores.append({'subset_defined_by': candidate, 'evaluated_candidate': evaluated,
                           'reference_ESP': grid, 'physical_contexts': len(values),
                           'median_RRMS': statistics.median(values),
                           'min_RRMS': min(values), 'max_RRMS': max(values),
                           'in_sample_or_held_out': 'see candidate training declaration; not independent FF validation'})
for name, data in [('CANDIDATE_COVERAGE.csv', coverage), ('MATCHED_SUBSET_SCORES.csv', scores)]:
    with (root/name).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
print(f'Wrote {len(coverage)} coverage rows and {len(scores)} matched-subset summaries.')
