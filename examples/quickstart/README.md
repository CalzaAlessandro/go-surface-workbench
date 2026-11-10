# Small GO sheet

From the repository root:

```bash
uv sync --frozen
uv run python examples/quickstart/run_quickstart.py --output-dir work/quickstart
```

The script uses `examples/finite_GO.json`, checks a 179-atom sheet with 232 bonds
and formal charge −1 e, and prepares Amber input templates. It runs no external
simulation engine. See the [quickstart](../../docs/QUICKSTART.md) for outputs
and individual commands.
