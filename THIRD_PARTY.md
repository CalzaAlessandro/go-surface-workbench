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

## GAFF2–AIM reference preset

The published subset under `resources/reference/gaff2_aim/` is copied byte for
byte from the archived `arvpinto/GO_sheets` source at commit
`1e59f43bd6c4d8c19d0a2a208edf53ce11e00bc5`. It retains the upstream GPL-3.0
license, topology/coordinate files for C/O 1.8 and 3.0, associated force-field
files and source documentation. `PROVENANCE.json` records download URLs and
SHA-256 checksums. Cite Pinto et al., DOI `10.1021/acs.jpcb.2c06037`.
This is a published reference preset, not a validated generic assignment to
arbitrary finite GO sheets. Its charge and AIM Lennard-Jones coefficients
remain paired and unchanged.

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
