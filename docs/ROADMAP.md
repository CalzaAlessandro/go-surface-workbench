<a id="direzione-scientifica-del-progetto"></a>
# Scientific direction of the project

This is a development proposal, not a claim of absolute novelty or
a list of features already implemented. The GOPY/HierGO sources remain cited in
[COMPARISON.md](COMPARISON.md); provenance and licenses do not change when
features are added or names rewritten. The specific literature review will need to
be repeated when defining a contribution for publication.

<a id="base-già-presente"></a>
## Existing foundation

The main model is an explicit graph with stable atom IDs and periodic bond
images. Site selection is random by default, reproducible,
uses exact counts, and provides warnings for infeasible loadings. Existing features include
anti epoxide opening with permitted fragments, Amber preparation without bondbydistance,
GROMACS conversion, import of equilibrated snapshots, and docking on 3D surfaces.
Geometric support does not imply validated MD parameters for every
graft. Amber export of covalently infinite sheets remains blocked.

<a id="estensioni-che-proporrei-in-ordine"></a>
## Proposed extensions, in order

| Priority | Proposed feature | Verifiable outcome |
|---|---|---|
| 1 | Inverse recipes from XPS with uncertainty | A set of integer counts compatible with C/O and spectral components, rather than an unjustified unique composition. Ratios/fractions derived from counts, rather than duplicated count inputs. |
| 2 | Controllable spatial distribution | Dispersed or clustered groups, protected regions, independent constraints on the two faces, and correlation statistics that can be compared across samples. |
| 3 | Catalog of covalent reactions | Templates with atom mapping, atom/charge balance, and stereochemistry: epoxide opening first, then other reactions only with documented chemistry and parameters. Not a universal grafting operation without checks. |
| 4 | Parameterized building blocks and junctions | Charges and parameters with a defined chemical scope, QM protocol, local environment, and transferability tests; export blocked if a junction is missing. |
| 5 | Surface–adsorption ensembles | The same graph on a planar sheet and equilibrated snapshots; preserved IDs, accessibility of both faces, local boxes for curvature, and comparison across seeds/conformations. |
| 6 | Resumable campaigns | Recipes and seeds in a manifest, the same local/Slurm worker, resumption only from verified artifacts, and a summary of failures. |

The most useful combination for your project would be XPS recipe → spatial
ensemble → reactive functionalization → verified parameterization → docking
on equilibrated snapshots. This requires demonstrating parameter transferability and
quantifying how adsorption changes between models, in addition to producing files.

<a id="velocità-memoria-e-semplicità"></a>
## Speed, memory, and simplicity

The first measured improvement is already implemented: vectorized collision
checks and reuse of graph views during edge termination.
[PERFORMANCE.md](PERFORMANCE.md) reports speed, memory use, and limitations; outputs
remain identical in the compatibility cases tested.

The next step, if justified by larger sheets, is a local spatial
index with exact checking of candidate distances. I would avoid dense
N×N matrices and parallelism over individual atoms. Seeds or independent
models are better candidates for parallel execution, avoiding oversubscription with QM/MD engines.

The construction core already uses NumPy and imports specialist backends only
when needed. I would consider separately installable optional dependencies for
RDKit fragments, ParmEd topologies, and NetCDF trajectories. This is not yet a
package change: it must be tested with minimal installations before changing
the dependencies. A single data schema and one worker for local/Slurm execution must remain
the foundation; reducing duplication matters more than reducing line counts at any cost.
