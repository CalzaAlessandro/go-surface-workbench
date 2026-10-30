# Atomic RESP charge context analysis

`go_param.context` summarizes fitted atomic charges in the local chemical context that produced them. It reports individual atom charges and roles, each partitioned functional group's net charge, all pairwise functional-group relationships, and background scaffold-carbon charges grouped by graph distance to the nearest functional group. It does not average charges by group or infer transferable lookup/rule/SQE parameters.

## Input

Pass a nonempty JSON list. Each item represents one model and one atomic-charge fit:

```json
[
  {
    "model_path": "models/pyrene_diol_epoxide/model.json",
    "charges_path": "calculations/.../fit_atomic_charges.json",
    "fit_evidence_path": "calculations/.../FIT_COMPLETE.json"
  }
]
```

`fit_evidence_path` is optional. Charge files need an `atoms` list in model order, with `atom_id`, `element`, and `charge_e` fields. If present, `model_id`, `formal_charge`, and per-atom `role` are checked too. The analyzer requires exact atom count, ID/order and element agreement, finite charges, and agreement between their sum and the model formal charge within six-decimal RESP output-rounding tolerance. Invalid input raises `ValueError`; it is not partially analyzed.

Run it from the project root:

```bash
PYTHONPATH=src python -m go_param.context input.json --output-dir reports/context_analysis
```

The callable API is `analyze_dataset(spec)`, where `spec` is either the JSON-list path or the parsed list. Each list entry becomes one `context_id`, even if another entry uses the same `model_id`; model counts and fit-context counts are reported separately. Three rigid rotations in one fit evidence file still produce only one context. The CLI writes `context_analysis.json`, `per_atom_charges.csv`, `functional_group_charges.csv`, `functional_group_pairs.csv`, and `background_carbon_distance_charges.csv`. JSON missing values are `null`; CSV missing values are empty fields. Empty pair tables still contain headers.

## Distance definitions

- Functional groups come from `go_param.modular.partition_model`. GRT edge caps are excluded from the functional-group list and pair combinations.
- Scaffold IDs come first from `metadata.scaffold_carbon_ids`. If absent, the analyzer uses explicit `scaffold_C` annotations in atom roles and fails if the carbon membership is ambiguous. Exocyclic `COOH_C` and `COO_C` cannot be included in the declared scaffold set; carbonyl-carbon membership follows the supplied scaffold metadata/roles because a carbonyl carbon may be part of the parent ring.
- Graph distance is the minimum count of C–C bonds on the original scaffold-carbon network between the two groups' scaffold-carbon anchors. Anchors are scaffold carbons inside the partitioned group, or directly bonded scaffold carbons if the group contains no scaffold carbon. Oxygen atoms and exocyclic carboxyl carbons are never graph nodes. The chosen atom path and anchor source are returned so the measurement is auditable. If either group has no scaffold anchor or the groups are disconnected in the scaffold graph, the distance and path are `null`.
- Geometric distance is the Euclidean distance between the centroids of the carbon atoms in the two partitioned functional-group blocks. It uses the model's coordinates in Å and does not apply periodic-boundary correction. It is reported separately from graph distance.
- Background-carbon charge rows include original scaffold carbons not assigned to a non-GRT functional group. `background_kind` labels each as `interior_scaffold` or `hydrogen_terminated_edge` (an `edge_C` role/annotation, including GRT caps). Rows retain the atom role and charge separately; the two classes are never pooled or averaged. Coverage reports counts by class. Each row also reports its nearest functional group(s) and graph distance; with no functional groups, nearest-group fields are `null`.
- Functional-group-pair rows include both group net charges and scaffold-anchor carbon charge sums alongside graph and geometric distances.

The output preserves the complete model metadata and optional fit-evidence JSON. Effective `scope` prioritizes explicit diagnostic labels from model metadata or nested fit evidence over generic protocol scope; `model_scope` and the original `protocol_scope` are also retained. `is_diagnostic` reports explicit diagnostic labels, and `reported_wavefunction_stability_validated` carries the reported boolean or `null` when absent. These fields are included in JSON and every CSV row. QM method, scope, fit status, upstream QM gate, orientation count, and source face metadata are annotations only. The analyzer never upgrades a diagnostic fit to accepted evidence. One charge file creates one analysis context; rigid rotations recorded in its fit evidence do not count as extra fits or independent samples. Face metadata is retained but not used to assign a group face.

## Scope

This is a descriptive analysis of the charges supplied; it performs no QM jobs, ESP fit, charge transfer, force-field assignment, or inference of an average atom type. Group net charges are sums of the listed atomic charges, while individual per-atom values remain available for environment comparisons. Graph distances, geometric distances, and net group charges are context descriptors, not unique observables of a particular charge or LJ parameter. Results should retain the model, conformer, method, and fit provenance when compared across the pilot dataset.
