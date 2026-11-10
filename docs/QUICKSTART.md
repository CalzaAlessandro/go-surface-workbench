# Quickstart

Build a small GO sheet, check its connectivity and prepare Amber inputs.
Run these commands from the repository root with Python 3.12 and [uv](https://docs.astral.sh/uv/):

```bash
uv sync --frozen
uv run python examples/quickstart/run_quickstart.py --output-dir work/quickstart
```

The example uses [`examples/finite_GO.json`](../examples/finite_GO.json).
It produces an unrelaxed sheet with 179 atoms (132 C, 33 H and 14 O),
232 bonds and formal charge −1 e. The graph has one connected component
and no undercoordinated atoms.

| Output under `work/quickstart/` | Contents |
|---|---|
| `finite_GO.json`, `finite_GO.pdb` | Sheet connectivity and coordinates |
| `amber/tleap.in`, `amber/GO_gaff2.frcmod` | LEaP input and parameter modifications |
| `amber/assigned_model.json`, `amber/library_audit.json` | Atom-type and library assignments |
| `amber/manifest.json`, `amber/normalized_libraries/` | Input manifest and unit libraries |

To run the same steps individually:

```bash
uv run python 01_build_GO.py --recipe examples/finite_GO.json --output work/finite_GO.json --pdb work/finite_GO.pdb --json
uv run python check_model.py --model work/finite_GO.json --json
uv run python 03_prepare_amber.py --model work/finite_GO.json --output-dir work/amber --radii mbondi2
```

This example prepares inputs using Python dependencies only. To compile an
Amber topology, install AmberTools and pass `--run` to `03_prepare_amber.py`.
See [engine setup](ENGINE_SETUP.md) for later simulation steps.

Graph checks establish connectivity and coordination. They do not establish
electronic structure, parameter accuracy or equilibration; see
[validation status](../VALIDATION.md).
