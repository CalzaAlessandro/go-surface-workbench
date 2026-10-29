# Shared interfaces

The interchange object is a plain JSON dictionary, schema_version=1:

```
{"schema_version":1,"name":"GO","atoms":[
 {"id":1,"element":"C","name":"CX","resname":"GGG","resid":1,
  "position":[0.0,0.0,0.0],"atom_type":null,"charge":null,
  "formal_charge":0,"metadata":{}}
],"bonds":[{"a":1,"b":2,"order":1.5,"image":[0,0,0]}],
"cell":[[10,0,0],[0,10,0],[0,0,40]],"pbc":[false,false,false],
"metadata":{"seed":42,"events":[]}}
```

Coordinates and cell are Angstrom; charges in elementary charge. For a bond a->b, the physical displacement is r_b + image @ cell - r_a. Normalize a<b and flip image when swapping. Atom IDs are stable, unique integers, not positional list indices. Bond order supports aromatic 1.5, but valence must account for graphene motifs rather than naive aromatic bond sums. Periodic slab pbc is x/y only. Zero-based array indices are never atom IDs.

Geometry agent owns model.py (validation, JSON/PDB I/O), build.py and functionalize.py plus their tests. Public functions: build_sheet(recipe:dict)->dict; open_epoxides(model:dict, recipe:dict)->dict; validate_model(model:dict)->dict; load_model(path)->dict; save_model(model,path)->None; write_pdb(model,path)->None. CLI passes recipes, no argparse inside these functions. Recipe conventions should be documented by geometry agent promptly.

Topology agent owns ff.py, topology.py and tests. Public functions: prepare_amber(model:dict, output_dir, **options)->dict (write inputs, optionally run external tleap); import_amber(topology, coordinates=None, **options)->dict; export_gromacs(topology, coordinates, output_prefix, **options)->dict; gbsa_readiness(model:dict)->dict. No importing geometry module required; accepts graph above.

Charge agent owns charges.py, qm.py and tests. Public functions: audit_charge_sources(output_dir, **options)->dict; prepare_qm(model:dict, output_dir, **options)->dict; compare_esp(reference_path, candidates:dict, output_dir, **options)->dict. May expose extra focused functions. Source files under resources/reference; maintain provenance. Do not import unfinished modules of other agents.

Root owns CLI, wrappers, common execution, dependencies, docs, examples, source snapshots and integration tests. Files are disjoint; notify root before altering another owner's files. Scripts never overwrite existing results unless an explicitly implemented safe policy says otherwise. Return clear ValueError/FileExistsError; CLI formats consistent Italian/English concise status.

GO creation is local-only through `01_build_GO.py` / `go.py build`. Select exactly one of `--recipe PATH`, `--preset NAME`, or `--list-presets`. Actual creation requires `--output`; `--pdb` adds coordinate export, and explicit `--seed` overrides the selected recipe's seed. Custom recipes retain their own defaults. `--dry-run` reports a local command without creating model or scheduler files. Build rejects Slurm/resource/worker-environment flags; other stages retain their existing local/Slurm interface.

Built-in recipe lookup uses `go_builder.presets` and package resources under `go_builder/data/layerone`, independent of checkout location. Preset provenance is stored outside the scientific recipe keys in model `metadata.preset`, with the selected name, original recipe, packaged resource and content hash. Source composition and constructed chemistry must remain distinguishable; choosing a commercial-material preset is not validation of an experimental material representation.
