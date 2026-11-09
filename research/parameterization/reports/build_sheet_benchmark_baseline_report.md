# GO Builder sheet-build benchmark

Read-only benchmark of `go_builder.build.build_sheet`; source files were not changed.
Python: `python` (3.14.10) on `darwin`.
Method: 5 wall-time repetitions after warmup; isolated tracemalloc and cProfile calls; each case alarm-limited to 90 s. Recipe seed 17; fixed group counts `{"COO": 2, "COOH": 2, "OH": 4, "carbonyl": 2, "epoxide": 3}`; finite sheets, both faces.

| Cells (nx × ny) | Untrimmed sheet C | Output atoms | Output bonds | Median wall time (s) | Tracemalloc peak (MiB) |
|---:|---:|---:|---:|---:|---:|
| 8 × 8 | 256 | 303 | 403 | 0.0497 | 1.02 |
| 12 × 12 | 576 | 639 | 883 | 0.1808 | 2.84 |
| 16 × 16 | 1024 | 1103 | 1555 | 0.4617 | 4.40 |

All repeats, traced calls, and profiled calls returned the same SHA-256 model per size. Exact digests and repeat distributions are in `build_sheet_benchmark_results.json`.

## Cumulative profile

### 8 × 8

| Function | Calls | Self (s) | Cumulative (s) |
|---|---:|---:|---:|
| `build.py:547:build_sheet` | 1 | 0.0002 | 0.0758 |
| `build.py:72:_collision_free` | 55 | 0.0285 | 0.0323 |
| `build.py:507:_cap_edges` | 1 | 0.0010 | 0.0291 |
| `build.py:283:_outward` | 44 | 0.0015 | 0.0211 |
| `model.py:33:adjacency` | 56 | 0.0117 | 0.0193 |
| `build.py:415:_decorate` | 1 | 0.0002 | 0.0162 |
| `~:0:<built-in method numpy.array>` | 38214 | 0.0098 | 0.0098 |
| `build.py:296:_place_group` | 13 | 0.0005 | 0.0085 |
| `build.py:51:_coordinates` | 57 | 0.0016 | 0.0048 |
| `model.py:62:validate_model` | 1 | 0.0023 | 0.0047 |

### 12 × 12

| Function | Calls | Self (s) | Cumulative (s) |
|---|---:|---:|---:|
| `build.py:547:build_sheet` | 1 | 0.0005 | 0.2626 |
| `build.py:72:_collision_free` | 78 | 0.1196 | 0.1306 |
| `build.py:507:_cap_edges` | 1 | 0.0033 | 0.1030 |
| `build.py:283:_outward` | 68 | 0.0053 | 0.0731 |
| `model.py:33:adjacency` | 80 | 0.0394 | 0.0641 |
| `build.py:415:_decorate` | 1 | 0.0004 | 0.0390 |
| `~:0:<built-in method numpy.array>` | 119434 | 0.0304 | 0.0304 |
| `build.py:296:_place_group` | 13 | 0.0010 | 0.0187 |
| `build.py:51:_coordinates` | 81 | 0.0048 | 0.0147 |
| `model.py:62:validate_model` | 1 | 0.0049 | 0.0103 |

### 16 × 16

| Function | Calls | Self (s) | Cumulative (s) |
|---|---:|---:|---:|
| `build.py:547:build_sheet` | 1 | 0.0008 | 0.6486 |
| `build.py:72:_collision_free` | 104 | 0.3232 | 0.3491 |
| `build.py:507:_cap_edges` | 1 | 0.0078 | 0.2560 |
| `build.py:283:_outward` | 92 | 0.0129 | 0.1788 |
| `model.py:33:adjacency` | 104 | 0.0963 | 0.1538 |
| `build.py:415:_decorate` | 1 | 0.0006 | 0.0760 |
| `~:0:<built-in method numpy.array>` | 272794 | 0.0693 | 0.0693 |
| `build.py:296:_place_group` | 13 | 0.0018 | 0.0330 |
| `build.py:51:_coordinates` | 105 | 0.0104 | 0.0320 |
| `build.py:245:_max_matching` | 2 | 0.0004 | 0.0231 |

## Allocation profile

Largest case peak: 4.40 MiB (Python-tracked allocations). Top retained allocations:

| Location | Bytes retained | Blocks |
|---|---:|---:|
| `copy.py:202` | 490136 | 3357 |
| `model.py:38` | 398936 | 7680 |
| `copy.py:199` | 219072 | 3423 |
| `model.py:36` | 195816 | 4324 |
| `model.py:37` | 179648 | 2504 |
| `copy.py:173` | 136640 | 2440 |
| `model.py:34` | 96768 | 1070 |
| `build.py:149` | 94656 | 2958 |
| `copy.py:177` | 90912 | 2437 |
| `build.py:267` | 56128 | 877 |
| `build.py:158` | 47520 | 1980 |
| `build.py:47` | 34056 | 255 |

## Measured hotspots

The median rose 9.3× from 8 × 8 to 16 × 16 while output atoms rose 3.6×. That points to per-edge work over the full graph rather than just the linear cost of emitting atoms and bonds.

In the 16 × 16 cProfile run (`build_sheet` total 0.649 s; profiler-instrumented), `_collision_free` ran 104 times and accumulated 0.349 s (0.323 s self time). `_cap_edges` accumulated 0.256 s. Those figures overlap: edge capping invokes collision screening once per hydrogen cap, so they must not be added.
`_outward` ran 92 times (0.179 s cumulative); `adjacency` ran 104 times (0.154 s); `_coordinates` ran 105 times (0.032 s). These repeat materialization of graph-derived views while the topology/coordinates are mostly stable.

## Low-risk optimization plan

1. Batch edge-cap collision checks after creating the candidate caps, while preserving the existing minimum-image calculation, exact cutoff comparisons, and deterministic atom order. This targets the largest repeated cost. Compare pass/fail outcomes and same-seed model hashes against the current implementation before replacing it.
2. Build adjacency, coordinate arrays, and ID-to-row lookups once per stable construction phase (`_cap_edges`, then `_decorate`), and update/invalidate them only when bonds or positions change. Keep current atom/bond iteration order and arithmetic expression order to reduce floating-point or ordering drift.
3. Avoid a new spatial-index dependency initially. If profiling after batching still shows collision checks dominate at larger sizes, benchmark an internal cell-list broad phase, then apply the current exact distance test to every candidate pair near the cutoff.
4. Preserve serialized-model hashes for these recipes as compatibility checks. Exact same-seed geometry can change if traversal ordering changes random-number consumption; if output changes, treat it as a behavior change requiring explicit review rather than an automatic speed-only refactor.

This is a local microbenchmark of unrelaxed finite-sheet construction, not a performance claim for large-scale production workloads. `tracemalloc` excludes some native NumPy allocations; wall-time runs are separate from profiling runs.
