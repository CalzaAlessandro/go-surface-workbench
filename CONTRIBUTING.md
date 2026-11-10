# Contributing

Use Python 3.11 or later. From the repository root:

```sh
uv sync --group dev
uv run pytest -q
```

Report errors through the bug-report form. Include the command, version or
commit, Python and OS versions, a small input, the expected result and the error log.

Keep pull requests focused. When changing behavior, add a regression test and
state which tests you ran. For parameter changes, include the source and
validation evidence, and preserve complete charge, bonded, Lennard-Jones and
scaling presets. See [validation status](VALIDATION.md).

CI checks Python 3.12 on macOS/arm64, matching the platform used for the
serialized-coordinate regression baselines.
