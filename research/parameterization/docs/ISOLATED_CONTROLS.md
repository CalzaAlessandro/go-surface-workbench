# Isolated GO motif controls

The supplemental `generatedmodels_isolated_controls/` directory contains four starting models: pristine C24 and C54 honeycomb parents, plus one basal epoxide on each parent. These controls separate the isolated basal-epoxide response from the existing multi-group distance/face panel. They do not change `generatedmodels_panel/` or its 39-record catalog.

The public builder API is:

```python
from go_param.mixed_panel import export_isolated_controls

catalog = export_isolated_controls("generatedmodels_isolated_controls")
```

Each catalog record has a `model.json`, `initial.sdf`, and `initial.xyz`. A single basal epoxide is selected on a center-near basal C–C bond for which the remaining pi graph has a perfect matching. The same honeycomb scaffold, valence checks, role annotations, and constrained MMFF94 starting-geometry protocol used by the panel are applied. The pristine references use the same scaffold and geometry protocol with no functional groups.

All records are neutral, even-electron, closed-shell singlet starting ansätze. Atom charges contain only formal integer values; partial charges are unset. Geometry remains an initial constrained MMFF model, and wavefunction stability and other QM evidence are unverified. The `pairwise_distance_applicable` metadata is false, `measured_pair_motif_indices` is null, and no pair distance or self-distance is reported for a single epoxide or a pristine reference.

The C24 and C54 pristine records are direct parent references. No single basal hydroxyl is included because adding one OH to these neutral parents changes the electron/sublattice balance and would not preserve the requested closed-shell starting ansatz. Existing distant OH pairs in the main panel remain multi-group pair controls; they are not relabeled as one isolated OH.
