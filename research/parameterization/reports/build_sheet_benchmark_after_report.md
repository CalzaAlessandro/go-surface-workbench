# GO Builder sheet-build benchmark

Read-only benchmark of `go_builder.build.build_sheet`; source files were not changed.
Python: `python` (3.14.10) on `darwin`.
Method: 5 wall-time repetitions after warmup; isolated tracemalloc and cProfile calls; each case alarm-limited to 90 s. Recipe seed 17; fixed group counts `{"COO": 2, "COOH": 2, "OH": 4, "carbonyl": 2, "epoxide": 3}`; finite sheets, both faces.

| Cells (nx × ny) | Untrimmed sheet C | Output atoms | Output bonds | Median wall time (s) | Tracemalloc peak (MiB) |
|---:|---:|---:|---:|---:|---:|
| 8 × 8 | 256 | 303 | 403 | 0.0242 | 1.06 |
| 12 × 12 | 576 | 639 | 883 | 0.0628 | 2.92 |
| 16 × 16 | 1024 | 1103 | 1555 | 0.1352 | 4.56 |

All repeats, traced calls, and profiled calls returned the same SHA-256 model per size. Exact digests and repeat distributions are in `build_sheet_benchmark_results.json`.

## Cumulative profile

### 8 × 8

| Function | Calls | Self (s) | Cumulative (s) |
|---|---:|---:|---:|
| `build.py:573:build_sheet` | 1 | 0.0002 | 0.0525 |
| `build.py:98:_collision_free` | 55 | 0.0002 | 0.0305 |
| `build.py:72:_collision_free_vectorized` | 55 | 0.0130 | 0.0303 |
| `build.py:434:_decorate` | 1 | 0.0002 | 0.0184 |
| `build.py:526:_cap_edges` | 1 | 0.0003 | 0.0173 |
| `build.py:315:_place_group` | 13 | 0.0006 | 0.0108 |
| `model.py:33:adjacency` | 18 | 0.0035 | 0.0058 |
| `_arraysetops_impl.py:958:isin` | 371 | 0.0003 | 0.0050 |
| `model.py:62:validate_model` | 1 | 0.0023 | 0.0047 |
| `_arraysetops_impl.py:806:_isin` | 371 | 0.0026 | 0.0046 |

### 12 × 12

| Function | Calls | Self (s) | Cumulative (s) |
|---|---:|---:|---:|
| `build.py:573:build_sheet` | 1 | 0.0004 | 0.1390 |
| `build.py:98:_collision_free` | 78 | 0.0007 | 0.0866 |
| `build.py:72:_collision_free_vectorized` | 78 | 0.0392 | 0.0859 |
| `build.py:526:_cap_edges` | 1 | 0.0005 | 0.0553 |
| `build.py:434:_decorate` | 1 | 0.0004 | 0.0430 |
| `build.py:315:_place_group` | 13 | 0.0012 | 0.0230 |
| `model.py:33:adjacency` | 18 | 0.0097 | 0.0150 |
| `_arraysetops_impl.py:958:isin` | 730 | 0.0005 | 0.0103 |
| `model.py:62:validate_model` | 1 | 0.0049 | 0.0102 |
| `_arraysetops_impl.py:806:_isin` | 730 | 0.0054 | 0.0095 |

### 16 × 16

| Function | Calls | Self (s) | Cumulative (s) |
|---|---:|---:|---:|
| `build.py:573:build_sheet` | 1 | 0.0007 | 0.2930 |
| `build.py:98:_collision_free` | 104 | 0.0015 | 0.1936 |
| `build.py:72:_collision_free_vectorized` | 104 | 0.0909 | 0.1921 |
| `build.py:526:_cap_edges` | 1 | 0.0008 | 0.1324 |
| `build.py:434:_decorate` | 1 | 0.0006 | 0.0824 |
| `build.py:315:_place_group` | 13 | 0.0020 | 0.0399 |
| `model.py:33:adjacency` | 18 | 0.0175 | 0.0272 |
| `build.py:258:_max_matching` | 2 | 0.0004 | 0.0227 |
| `~:0:<method 'setdefault' of 'dict' objects>` | 313248 | 0.0213 | 0.0213 |
| `build.py:266:augment` | 10424 | 0.0059 | 0.0191 |

## Allocation profile

Largest case peak: 4.56 MiB (Python-tracked allocations). Top retained allocations:

| Location | Bytes retained | Blocks |
|---|---:|---:|
| `copy.py:202` | 490136 | 3357 |
| `model.py:38` | 398864 | 7679 |
| `copy.py:199` | 219072 | 3423 |
| `model.py:36` | 195736 | 4319 |
| `model.py:37` | 179576 | 2503 |
| `copy.py:173` | 136640 | 2440 |
| `model.py:34` | 96704 | 1069 |
| `build.py:162` | 94656 | 2958 |
| `copy.py:177` | 90912 | 2437 |
| `build.py:280` | 56000 | 875 |
| `build.py:171` | 47520 | 1980 |
| `build.py:47` | 34056 | 255 |

## Before and after

Measured against the baseline script/report preserved alongside these results. The same source environment, recipe dimensions, seed, composition, repetitions, tracing/profile procedure, and model-hash assertions were used.

| Cells | Before median (s) | After median (s) | Speedup | Before peak (MiB) | After peak (MiB) |
|---:|---:|---:|---:|---:|---:|
| 8 × 8 | 0.0492 | 0.0242 | 2.04× | 1.02 | 1.06 |
| 12 × 12 | 0.1812 | 0.0628 | 2.88× | 2.84 | 2.92 |
| 16 × 16 | 0.4667 | 0.1352 | 3.45× | 4.40 | 4.56 |

## Measured hotspots

The median rose 5.6× from 8 × 8 to 16 × 16 while output atoms rose 3.6×. In the 16 × 16 cProfile run (`build_sheet` total 0.293 s; profiler-instrumented), vectorized collision screening ran 104 times and accumulated 0.192 s (0.091 s self time). `_cap_edges` accumulated 0.132 s.
`_outward` ran 92 times (0.011 s); `adjacency` ran 18 times (0.027 s); `_coordinates` ran 20 times (0.006 s). `_cap_edges` now creates these views once for its stable sheet graph and reuses them across caps.

## Refactor and compatibility

The collision screen uses O(N) coordinate, ID, element, mask, and distance arrays per new atom. It keeps the original `np.linalg.norm` and periodic minimum-image calculation, excludes the same atom index and bonded IDs, and applies the same strict `< 0.8 Å` H-containing / `< 1.25 Å` heavy-heavy cutoffs. Cap placement remains sequential, so clash checks and tilt-angle choices are not batched or reordered.
Same-seed serialized model SHA-256 values match the pre-refactor baseline for finite, holed, periodic, and mixed fixed/random-group recipes. The added tests compare collision outcomes to a scalar reference at both sides of each cutoff, across PBC, self/bond exclusions, plus the unchanged placement-failure report. No timing threshold is asserted in tests.

## Remaining optimization options

If larger sheet benchmarks still justify further work, consider a local cell-list broad phase only after profiling the vectorized version. It would need the exact current distance check as the final cutoff gate. A batch-cap algorithm remains out of scope because it could change sequential tilt choices and serialized coordinates.

This is a local microbenchmark of unrelaxed finite-sheet construction, not a performance claim for all workloads. `tracemalloc` excludes some native NumPy allocations; wall-time runs are separate from profiling runs.
