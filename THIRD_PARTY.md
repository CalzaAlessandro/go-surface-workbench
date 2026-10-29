# Third-party and reference material

This is a private development snapshot. No project-wide license has been
selected, and no ownership claim is made over the copied reference data. Review
all source terms before any public redistribution.

## HierGO reference files

Five original reference files are retained under `resources/reference/opls/`
from the HierGO repository at commit
`60382a1a1807ef0003ffcc3903392e913b1e5e9a`: `LICENSE`, `OPLS-top.f90`,
`Oplsatom2`, `README.md`, and `atom-typing.py`. The upstream license is
GPL-3.0. `PROVENANCE.json` records exact upstream paths, URLs, and SHA-256
hashes. The charge audit reads the reference table; this snapshot does not run
the upstream programs.

## Current parameter libraries

The files under `resources/reference/current/` are user-supplied Amber/GAFF2
libraries and parameter files copied for the workbench. Their origin and
copy-time SHA-256 values are recorded in `resources/reference/SOURCES.json`.
Their individual license and redistribution status have not been established
here; do not assume they are covered by a project license.

## GOPY citation

GOPY informed the conceptual comparison of graphene-oxide construction and
functionalization workflows. No GOPY program is executed by the builder. The research audit now includes the original GO.frcmod parameter reference under research/parameterization/references/gopy; its source/terms must remain attributed to GOPY. Cite
the project when discussing that conceptual relationship: Iourarum, GOPY,
<https://github.com/Iourarum/GOPY> (repository source and tutorial links are
also listed in `docs/COMPARISON.md`).

## Runtime software

Vina, Open Babel/MGLTools, Amber/AmberTools, Gaussian, and GNINA are not bundled
here. See `docs/ENGINE_SETUP.md` for external engine expectations.
