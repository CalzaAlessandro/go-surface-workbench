# Graph construction performance

`build_sheet` now vectorizes the per-atom collision screen and reuses adjacency,
coordinate, atom-ID, and cell views while capping finite-sheet edges. The
collision screen still checks each new atom against the current model in the
original order. Edge hydrogens are placed and tilted sequentially; their checks
are not batched.

## Measured local result (historical pre-Blocco-1 vectorization baseline)

The initial vectorization benchmark used Python 3.14.10 from the existing GO_BUILDER environment,
finite sheets, seed 17, both faces, and fixed group counts OH 4, epoxide 3,
carbonyl 2, COOH 2, and COO 2. It took five wall-time measurements per size
after warm-up, plus separate tracemalloc and cProfile runs. Each size had a
90-second timeout.

| Cells | Output atoms | Before median | After median | Speedup | Python-tracked peak before | After |
|---:|---:|---:|---:|---:|---:|---:|
| 8 × 8 | 303 | 0.0492 s | 0.0242 s | 2.04× | 1.02 MiB | 1.06 MiB |
| 12 × 12 | 639 | 0.1812 s | 0.0628 s | 2.88× | 2.84 MiB | 2.92 MiB |
| 16 × 16 | 1103 | 0.4667 s | 0.1352 s | 3.45× | 4.40 MiB | 4.56 MiB |

During that historical refactor, pre- and post-change timed, traced, and profiled
calls produced identical model SHA-256 values for each benchmark recipe. Dedicated
tests in `tests/test_build_performance_equivalence.py` pin regression hashes for
finite, holed, periodic, and mixed fixed/random-group recipes and compare the
vector collision outcomes with the scalar reference around both cutoffs. (Three of
the four pinned hashes were intentionally updated during Blocco 1 for improved
constructed edge-cap/functional-group geometry and orthonormal post-puckering
frames, as documented inline in `tests/test_build_performance_equivalence.py`;
`tests/test_passo1e_audit.py` separately guards the dense `24×24` `O/C=0.35`,
`pH=7.4` build at `dt < 3.0 s`.)

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

## Blocco 1 $O(N)$ scalability (2026-11-20)

In Blocco 1, both `build_sheet` and `check_geometry_integrity` were upgraded with:
- an incremental **PBC-aware 3D spatial voxel grid** (`_SpatialGrid` in `build.py` and `_build_spatial_cells` in `model.py`) replacing $O(N^2)$ distance screening;
- **zero-recursion iterative bipartite matching** (`_iterative_kuhn_matching` for `<= 1000` left nodes preserving exact RNG equivalence, and `_hopcroft_karp_matching` above `1000` nodes);
- **fast structural cloning** (`_clone_model_fast`) inside the decoration loop and steric edge capacity estimation in `composition.py`.

- **translation-stable convex-hull scaffold-diameter calculation** (`_compute_scaffold_diameter_nm` in `model.py`): exact $O(N \log N)$ via 2D convex hull + rotating calipers for planar sheets, and 3D convex-hull vertex reduction with centered pairwise differences for puckered sheets.

Measured on Apple Silicon (`CPython 3.13.15`, single thread, uninstrumented wall time; `0` clashes / `Gate: PASS` on all cases; reproducible via `uv run python scripts/benchmark_build.py --blocco1`):

| Benchmark case | Output atoms | Groups | Build (s) | Geometry gate (s) | Chemistry summary (s) | Memory | Realized O/C |
|---|---:|---:|---:|---:|---:|---:|---:|
| Finite `3×3` (`O/C=0.35`, pH 7.4) | 58 | 9 | 0.017 | 0.001 | 0.000 | 0.2 MiB Python peak | 0.367 |
| Finite `6×6` (`O/C=0.35`, pH 7.4) | 231 | 40 | 0.059 | 0.011 | 0.001 | 0.9 MiB Python peak | 0.348 |
| Hexagonal zigzag (`d=4.0 nm`, pH 7.4) | 638 | 119 | 0.112 | 0.020 | 0.002 | 3.2 MiB Python peak | 0.350 |
| Finite `12×12` + 5% mixed holes + 2 islands | 899 | 161 | 0.246 | 0.027 | 0.002 | 4.1 MiB Python peak | 0.348 |
| Finite `12×12` (`O/C=0.35`, pH 7.4) | 913 | 170 | 0.159 | 0.030 | 0.005 | 4.3 MiB Python peak | 0.349 |
| Periodic `24×24` (`O/C=0.30`) | 3,456 | 691 | 0.764 | 0.172 | 0.021 | 17.3 MiB Python peak | 0.300 |
| Finite `24×24` (`O/C=0.35`, pH 7.4) | 3,666 | 729 | 0.671 | 0.116 | 0.025 | 17.5 MiB Python peak | 0.350 |
| Finite `36×36` (`O/C=0.35`, pH 7.4) | 8,294 | 1,697 | 1.595 | 0.258 | 0.068 | 39.3 MiB Python peak | 0.350 |
| Periodic `48×48` (`O/C=0.30`) | 13,824 | 2,765 | 3.293 | 0.748 | 0.215 | 69.3 MiB Python peak | 0.300 |
| Periodic `64×64` (`O/C=0.30`) | 24,576 | 4,915 | 6.304 | 1.331 | 0.394 | ~203.4 MiB total RSS | 0.300 |

Pinned Blocco 1 regression hashes in `tests/test_build_performance_equivalence.py` pass with zero mismatches.

A portable benchmark driver is versioned as `scripts/benchmark_build.py`:
- `uv run python scripts/benchmark_build.py` reproduces the 3-size fixed-count microbenchmark (`8×8`, `12×12`, `16×16`).
- `uv run python scripts/benchmark_build.py --blocco1` reproduces the 10-case Blocco 1 scaling table above (`58` to `24,576` atoms).

Outputs go to ignored `results/benchmarks/`. Copy-time numerical before/after evidence is versioned under `docs/benchmarks/`; timings are machine-specific.
