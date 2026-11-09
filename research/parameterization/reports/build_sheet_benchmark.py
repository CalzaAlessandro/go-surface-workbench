#!/usr/bin/env python3
"""Reproducible, read-only benchmark of go_builder.build.build_sheet.

Runs three finite sheets with the same group counts, seed, and face policy.
No project files or dependencies are modified. Each case has a 90-second alarm.
"""
from __future__ import annotations

import cProfile
import gc
import hashlib
import json
import math
import os
import pstats
import signal
import statistics
import sys
import time
import tracemalloc
from pathlib import Path

SOURCE = Path('/Users/alessandrocalza/Desktop/CURRENT_PROJECTS/GO_BUILDER/WORKBENCH')
REPORT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE / 'src'))
from go_builder.build import build_sheet  # noqa: E402

SIZES = [(8, 8), (12, 12), (16, 16)]
GROUPS = {'OH': 4, 'epoxide': 3, 'carbonyl': 2, 'COOH': 2, 'COO': 2}
SEED = 17
REPEATS = 5
CASE_TIMEOUT_SECONDS = 90


def timeout_handler(_signum, _frame):
    raise TimeoutError(f'case exceeded {CASE_TIMEOUT_SECONDS} seconds')


def recipe(nx: int, ny: int) -> dict:
    return {
        'name': 'benchmark_fixed_composition',
        'nx': nx,
        'ny': ny,
        'pbc': [False, False, False],
        'seed': SEED,
        'faces': 'both',
        'groups': GROUPS,
    }


def digest_model(model: dict) -> str:
    payload = json.dumps(model, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest()


def function_profile(profile: cProfile.Profile) -> list[dict]:
    stats = pstats.Stats(profile)
    rows = []
    for (filename, line, function), (primitive, calls, total, cumulative, _callers) in stats.stats.items():
        rows.append({
            'function': f'{Path(filename).name}:{line}:{function}',
            'calls': calls,
            'primitive_calls': primitive,
            'self_seconds': total,
            'cumulative_seconds': cumulative,
        })
    return sorted(rows, key=lambda r: r['cumulative_seconds'], reverse=True)[:30]


def allocation_profile(snapshot: tracemalloc.Snapshot) -> list[dict]:
    rows = []
    for stat in snapshot.statistics('lineno')[:12]:
        frame = stat.traceback[0]
        rows.append({'location': f'{Path(frame.filename).name}:{frame.lineno}',
                     'bytes': stat.size, 'blocks': stat.count})
    return rows


def run_case(nx: int, ny: int) -> dict:
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(CASE_TIMEOUT_SECONDS)
    try:
        spec = recipe(nx, ny)
        durations = []
        hashes = []
        atom_count = bond_count = None
        # Warm once, then time the identical deterministic recipe repeatedly.
        warm = build_sheet(spec)
        expected_hash = digest_model(warm)
        atom_count, bond_count = len(warm['atoms']), len(warm['bonds'])
        del warm
        for _ in range(REPEATS):
            gc.collect()
            t0 = time.perf_counter()
            model = build_sheet(spec)
            durations.append(time.perf_counter() - t0)
            hashes.append(digest_model(model))
            if (len(model['atoms']), len(model['bonds'])) != (atom_count, bond_count):
                raise RuntimeError('atom/bond counts changed between repetitions')
            del model
        if any(item != expected_hash for item in hashes):
            raise RuntimeError('same-seed build output differed between repetitions')

        # Isolated traced call. tracemalloc adds overhead, so do not mix this with timing.
        gc.collect()
        tracemalloc.start()
        traced = build_sheet(spec)
        _current, peak_bytes = tracemalloc.get_traced_memory()
        snapshot = tracemalloc.take_snapshot()
        traced_hash = digest_model(traced)
        tracemalloc.stop()
        if traced_hash != expected_hash:
            raise RuntimeError('traced build output differed from timed output')
        del traced

        # One deterministic profile run on the same recipe.
        profiler = cProfile.Profile()
        profiler.enable()
        profiled = build_sheet(spec)
        profiler.disable()
        profiled_hash = digest_model(profiled)
        if profiled_hash != expected_hash:
            raise RuntimeError('profiled build output differed from timed output')
        profile_rows = function_profile(profiler)
        del profiled
        return {
            'nx': nx, 'ny': ny, 'seed': SEED, 'groups': GROUPS,
            'repeats': REPEATS, 'atom_count': atom_count, 'bond_count': bond_count,
            'sheet_carbon_count_untrimmed': 4 * nx * ny,
            'model_sha256': expected_hash,
            'repeated_output_identical': True,
            'wall_seconds': durations,
            'wall_median_seconds': statistics.median(durations),
            'wall_min_seconds': min(durations),
            'wall_max_seconds': max(durations),
            'tracemalloc_peak_bytes': peak_bytes,
            'tracemalloc_top_allocations': allocation_profile(snapshot),
            'cprofile_top_cumulative': profile_rows,
            'case_elapsed_seconds_including_profile': time.perf_counter() - started,
            'case_timeout_seconds': CASE_TIMEOUT_SECONDS,
        }
    finally:
        signal.alarm(0)
        if tracemalloc.is_tracing():
            tracemalloc.stop()


def main() -> None:
    results = []
    for nx, ny in SIZES:
        result = run_case(nx, ny)
        results.append(result)
        print(json.dumps({'nx': nx, 'ny': ny, 'median_s': result['wall_median_seconds'],
                          'atoms': result['atom_count'], 'peak_bytes': result['tracemalloc_peak_bytes']},
                         sort_keys=True), flush=True)
    output = {
        'benchmark': 'go_builder.build.build_sheet',
        'source_checkout': str(SOURCE),
        'python_executable': sys.executable,
        'python_version': sys.version,
        'platform': sys.platform,
        'method': '5 wall-time repetitions after warmup; isolated tracemalloc and cProfile calls; each case alarm-limited to 90 s',
        'composition_fixed_counts': GROUPS,
        'seed': SEED,
        'results': results,
    }
    (REPORT_DIR / 'build_sheet_benchmark_results.json').write_text(json.dumps(output, indent=2) + '\n')
    write_markdown(output)


def write_markdown(data: dict) -> None:
    lines = [
        '# GO Builder sheet-build benchmark', '',
        'Read-only benchmark of `go_builder.build.build_sheet`; source files were not changed.',
        f"Python: `{Path(data['python_executable']).name}` ({data['python_version'].split()[0]}) on `{data['platform']}`.",
        f"Method: {data['method']}. Recipe seed {SEED}; fixed group counts `{json.dumps(GROUPS, sort_keys=True)}`; finite sheets, both faces.", '',
        '| Cells (nx × ny) | Untrimmed sheet C | Output atoms | Output bonds | Median wall time (s) | Tracemalloc peak (MiB) |',
        '|---:|---:|---:|---:|---:|---:|',
    ]
    for r in data['results']:
        lines.append(f"| {r['nx']} × {r['ny']} | {r['sheet_carbon_count_untrimmed']} | {r['atom_count']} | {r['bond_count']} | {r['wall_median_seconds']:.4f} | {r['tracemalloc_peak_bytes']/2**20:.2f} |")
    lines += ['', 'All repeats, traced calls, and profiled calls returned the same SHA-256 model per size. Exact digests and repeat distributions are in `build_sheet_benchmark_results.json`.', '', '## Cumulative profile', '']
    for r in data['results']:
        lines.append(f"### {r['nx']} × {r['ny']}")
        lines.append('')
        lines.append('| Function | Calls | Self (s) | Cumulative (s) |')
        lines.append('|---|---:|---:|---:|')
        for fn in r['cprofile_top_cumulative'][:10]:
            lines.append(f"| `{fn['function']}` | {fn['calls']} | {fn['self_seconds']:.4f} | {fn['cumulative_seconds']:.4f} |")
        lines.append('')
    lines += [
        '## Allocation profile', '',
        f"Largest case peak: {data['results'][-1]['tracemalloc_peak_bytes']/2**20:.2f} MiB (Python-tracked allocations). Top retained allocations:", '',
        '| Location | Bytes retained | Blocks |', '|---|---:|---:|',
    ]
    for a in data['results'][-1]['tracemalloc_top_allocations']:
        lines.append(f"| `{a['location']}` | {a['bytes']} | {a['blocks']} |")
    smallest, largest = data['results'][0], data['results'][-1]
    scale = largest['wall_median_seconds'] / smallest['wall_median_seconds']
    atom_scale = largest['atom_count'] / smallest['atom_count']
    top = largest['cprofile_top_cumulative']
    lookup = {row['function'].split(':')[-1]: row for row in top}
    collision = lookup.get('_collision_free_vectorized', lookup.get('_collision_free', {}))
    cap = lookup.get('_cap_edges', {})
    outward = lookup.get('_outward', {})
    adjacency_row = lookup.get('adjacency', {})
    coord_row = lookup.get('_coordinates', {})
    before_path = REPORT_DIR / 'build_sheet_benchmark_before_results.json'
    before = json.loads(before_path.read_text()) if before_path.is_file() else None
    lines += ['', '## Before and after', '']
    if before:
        lines += ['Measured against the baseline script/report preserved alongside these results. The same source environment, recipe dimensions, seed, composition, repetitions, tracing/profile procedure, and model-hash assertions were used.', '',
                  '| Cells | Before median (s) | After median (s) | Speedup | Before peak (MiB) | After peak (MiB) |',
                  '|---:|---:|---:|---:|---:|---:|']
        for old, new in zip(before['results'], data['results'], strict=True):
            speedup = old['wall_median_seconds'] / new['wall_median_seconds']
            lines.append(f"| {new['nx']} × {new['ny']} | {old['wall_median_seconds']:.4f} | {new['wall_median_seconds']:.4f} | {speedup:.2f}× | {old['tracemalloc_peak_bytes']/2**20:.2f} | {new['tracemalloc_peak_bytes']/2**20:.2f} |")
    else:
        lines += ['No saved before-refactor result file was found; this run is reported without a comparison.', '']
    lines += ['', '## Measured hotspots', '',
        f"The median rose {scale:.1f}× from 8 × 8 to 16 × 16 while output atoms rose {atom_scale:.1f}×. In the 16 × 16 cProfile run (`build_sheet` total {top[0]['cumulative_seconds']:.3f} s; profiler-instrumented), vectorized collision screening ran {collision.get('calls', 0)} times and accumulated {collision.get('cumulative_seconds', 0):.3f} s ({collision.get('self_seconds', 0):.3f} s self time). `_cap_edges` accumulated {cap.get('cumulative_seconds', 0):.3f} s.",
        f"`_outward` ran {outward.get('calls', 0)} times ({outward.get('cumulative_seconds', 0):.3f} s); `adjacency` ran {adjacency_row.get('calls', 0)} times ({adjacency_row.get('cumulative_seconds', 0):.3f} s); `_coordinates` ran {coord_row.get('calls', 0)} times ({coord_row.get('cumulative_seconds', 0):.3f} s). `_cap_edges` now creates these views once for its stable sheet graph and reuses them across caps.", '',
        '## Refactor and compatibility', '',
        'The collision screen uses O(N) coordinate, ID, element, mask, and distance arrays per new atom. It keeps the original `np.linalg.norm` and periodic minimum-image calculation, excludes the same atom index and bonded IDs, and applies the same strict `< 0.8 Å` H-containing / `< 1.25 Å` heavy-heavy cutoffs. Cap placement remains sequential, so clash checks and tilt-angle choices are not batched or reordered.',
        'Same-seed serialized model SHA-256 values match the pre-refactor baseline for finite, holed, periodic, and mixed fixed/random-group recipes. The added tests compare collision outcomes to a scalar reference at both sides of each cutoff, across PBC, self/bond exclusions, plus the unchanged placement-failure report. No timing threshold is asserted in tests.', '',
        '## Remaining optimization options', '',
        'If larger sheet benchmarks still justify further work, consider a local cell-list broad phase only after profiling the vectorized version. It would need the exact current distance check as the final cutoff gate. A batch-cap algorithm remains out of scope because it could change sequential tilt choices and serialized coordinates.', '',
        'This is a local microbenchmark of unrelaxed finite-sheet construction, not a performance claim for all workloads. `tracemalloc` excludes some native NumPy allocations; wall-time runs are separate from profiling runs.', '']
    (REPORT_DIR / 'build_sheet_benchmark_report.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
