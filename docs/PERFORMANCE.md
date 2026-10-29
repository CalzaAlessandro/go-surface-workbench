# Graph construction performance

`build_sheet` now vectorizes the per-atom collision screen and reuses adjacency,
coordinate, atom-ID, and cell views while capping finite-sheet edges. The
collision screen still checks each new atom against the current model in the
original order. Edge hydrogens are placed and tilted sequentially; their checks
are not batched.

## Measured local result

The benchmark used Python 3.14.10 from the existing GO_BUILDER environment,
finite sheets, seed 17, both faces, and fixed group counts OH 4, epoxide 3,
carbonyl 2, COOH 2, and COO 2. It took five wall-time measurements per size
after warm-up, plus separate tracemalloc and cProfile runs. Each size had a
90-second timeout.

| Cells | Output atoms | Before median | After median | Speedup | Python-tracked peak before | After |
|---:|---:|---:|---:|---:|---:|---:|
| 8 × 8 | 303 | 0.0492 s | 0.0242 s | 2.04× | 1.02 MiB | 1.06 MiB |
| 12 × 12 | 639 | 0.1812 s | 0.0628 s | 2.88× | 2.84 MiB | 2.92 MiB |
| 16 × 16 | 1103 | 0.4667 s | 0.1352 s | 3.45× | 4.40 MiB | 4.56 MiB |

The pre- and post-change timed, traced, and profiled calls produced identical
model SHA-256 values for each benchmark recipe. Dedicated tests also pin
baseline output hashes for finite, holed, periodic, and mixed fixed/random-group
recipes, compare the vector collision outcomes with the previous scalar loop
around both strict cutoffs, and check the unchanged placement-failure report.
No timing threshold is asserted in tests.

The vectorized screen retains O(N) temporary arrays per new atom. The memory
increase above is about 3–4% in Python-tracked peak allocation for these sizes;
tracemalloc does not include every native NumPy allocation. This is a local
microbenchmark of unrelaxed finite-sheet construction, not a guarantee for all
recipes or larger workloads.

## Reproducing the benchmark

The original benchmark and raw before/after results are kept in the local
`Desktop/GO_PARAMETERIZATION/reports/` directory:

- `build_sheet_benchmark.py` — benchmark driver
- `build_sheet_benchmark_before_report.md` and
  `build_sheet_benchmark_before_results.json` — fresh pre-refactor baseline
- `build_sheet_benchmark_report.md` and
  `build_sheet_benchmark_results.json` — post-refactor report and raw results
- `build_sheet_benchmark_baseline_report.md` and
  `build_sheet_benchmark_baseline_results.json` — earlier baseline retained
  from the initial snapshot

Run the driver with the existing project Python environment. The script points
to the local `go-surface-workbench` checkout and writes its report/results next
to itself. The full project test suite passed before and after this refactor;
the new compatibility suite passed 14 tests.

## Remaining limits

Collision screening remains proportional to the number of candidates times the
current atom count, and each candidate still needs a minimum-image distance
array. A cell-list broad phase could help larger sheets, but would require a
separate equivalence check with the current exact distance/cutoff test as the
final gate. Batching edge capping was avoided because it could change the
sequential clash decisions and tilt angles.

A portable driver is versioned as `scripts/benchmark_build.py`; run it with
`uv run python scripts/benchmark_build.py`. Outputs go to ignored
`results/benchmarks/`. Copy-time numerical before/after evidence is versioned
under `docs/benchmarks/`; timings are machine-specific.
