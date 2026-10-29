# Software validation status

The private development snapshot passed 258 tests on 29 October 2026 using
Python 3.14.10 and the existing isolated GO_BUILDER environment. The suite covers
graph construction, placement failures, anti grafting, external-engine input
preparation, topology conversion and execution wrappers.

Reproduce with `uv sync --group dev` then `uv run pytest -q`. This is software
verification, not validation of new charge models, adsorption energetics or MD.

Historical external-engine logs and campaign outputs remain in the original
local GO_BUILDER validation directory; they are intentionally not versioned.
Gaussian/RESP work is a separate ongoing GO_PARAMETERIZATION campaign.
