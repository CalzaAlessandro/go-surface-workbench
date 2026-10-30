# Mixed functional-group distance panel

This panel prepares explicit finite PAH graphs and initial geometries for later
tests of charge context dependence. It contains **34 pair configurations and
five edge controls**, without QM results or fitted partial charges. Generating
this catalog does not select or submit all models for QM. A first QM stage can
select 12–18 models while keeping appropriate comparisons paired.

## Parents, chemistry and graph gate

The analytic honeycomb construction joins hexagonal cells of axial radius
`shell − 1`, merges shared vertices and retains all unique C–C edges. Shell 2
is the compact C₂₄H₁₂ coronene parent; shell 3 is C₅₄H₁₈ circumcoronene. Original
carbon IDs are one-based, deterministic, sorted by y then x, and remain stable
for every functionalization of the same parent. Hydrogen termination is the
intrinsic finite PAH boundary, not a bulk graphene environment.

Pair categories have fixed chemical counts within a parent:

| Category | Motifs | C₂₄ formula | C₅₄ formula |
|---|---|---|---|
| OH–OH | Two neutral OH | C₂₄H₁₄O₂ | C₅₄H₂₀O₂ |
| EPX–EPX | Two basal epoxides | C₂₄H₁₂O₂ | C₅₄H₁₈O₂ |
| OH–EPX | Two neutral OH and one basal epoxide | C₂₄H₁₄O₃ | C₅₄H₂₀O₃ |

Every basal functionalized carbon has three original carbon neighbors. Epoxide
motifs use an original neighboring carbon pair. No carbon is reused by another
motif. A hydroxyl pair removes one site from each graphene sublattice; an
epoxide also removes one from each. The OH–EPX category therefore includes an
explicit second, adjacent OH balancing partner. Both OH groups remain cis in
this category; the measured OH–EPX face comparison changes the epoxide face.
All three motifs and their positions are published in each model.

After all sp³ sites are removed from the π network, a bipartite perfect matching
must cover **every remaining π carbon**. Matching edges receive explicit double
bonds and all other scaffold edges single bonds. A failed matching rejects the
candidate. This gate replaces blind decrement of a preselected local double
bond. RDKit then sanitizes the full graph and adds explicit valence hydrogens;
radical atoms are rejected. Matching, balanced sublattices and even electron
count support a closed-shell starting graph. They do not establish a stable
singlet wavefunction; QM wavefunction stability remains unverified.

## Measured separation and paired faces

Graph distance is the minimum number of **original scaffold C–C edges** between
the carbon site sets of the two measured motifs. It is independent of newly
added oxygen shortcuts. Near means one edge; intermediate means two or three;
far means four or more. The actual numeric distance is always supplied. Within
each category/bin, the smallest feasible distance is selected, giving matching
numeric distances across parents wherever possible:

| Category | Near | Intermediate | Far, C₂₄ | Far, C₅₄ |
|---|---:|---:|---:|---:|
| OH–OH | 1 | 3 | 5 | 5 |
| EPX–EPX | 1 | 2 | unattainable | 4 |
| OH–EPX | 1 | 2 | 4 | 4 |

For tied graph distances, selection favors a motif centroid near the parent
center, then stable site IDs. This is a deterministic representative placement,
not exhaustive sampling of inequivalent edges or isomers. Actual motif-carbon
centroid distances and oxygen–oxygen distances in Å are stored after initial
geometry optimization. Distances can change during subsequent QM optimization.

The full requested design is two sizes × three categories × three bins × cis
and trans = 36 cells. The two C₂₄ EPX–EPX far cells cannot satisfy the distance
threshold and are explicitly listed in `omitted_cells` and `warnings` rather
than assigned a false distance. The remaining 34 cells are generated. Each
cis/trans pair shares the same sites, graph, parent, group counts and protonation;
only the relative requested face changes. Comparisons of different sizes retain
the same counts but necessarily change coverage and boundary/size context.

Cis/trans labels refer to signs of oxygen heights relative to the least-squares
plane of unfunctionalized scaffold carbons, oriented toward +z. Published heights
must agree with the requested faces and have magnitude greater than 0.35 Å.
The actual tetrahedral stereochemistry is assigned from the resulting coordinates.

## Geometry and intended use

Starting carbon positions come from the analytic honeycomb at a 1.42 Å C–C
spacing. Functional carbons and oxygen positions are displaced to the requested
faces before explicit H placement. MMFF94 then optimizes with flat-bottom
position restraints on scaffold carbons (0.35 Å) and functional oxygens (0.4 Å),
using a 100 kcal/(mol Å²) force constant. Constraints retain the intended initial
configuration; **they are not a stability calculation**. The model records
convergence or an iteration limit. Every initial bond must lie between 0.85 and
1.95 Å, coordinates must be finite, and the requested faces are rechecked.

These are constrained force-field starting geometries, not unconstrained minima,
QM geometries, equilibrated GO structures or evidence for stable epoxides. The
earlier small-pyrene epoxide pilot's C–C opening motivates checking larger finite
parents and electronic states; it does not establish instability of GO epoxides
in general. A matching-valid graph can still rearrange or reveal electronic
instability during QM. Geometry, connectivity, distances and wavefunction status
must be checked again after each actual calculation.

The default partition labels C₅₄ configurations `training_context` and C₂₄
configurations `held_out_size`. These labels are a prospective validation split,
not evidence of an existing fit. To audit context dependence, compare group
potentials/charges at matched parent, category, loading, protonation and separation
before assessing transfer to a different size. No partial-charge vector is
assigned by this module. The intended downstream modular background-carbon-zero
RESP model and shared functional-group classes require actual QM/ESP data and
separate constrained fitting.

## Edge controls and library coverage

Five separate controls are listed in `edge_controls`:

| Class | Initial parent | Explicit chemical definition |
|---|---|---|
| EPH | C₂₄H₁₂ | Epoxide with carbon attached-H counts `[0, 1]` |
| EPT | C₂₄H₁₂ | Epoxide with carbon attached-H counts `[1, 1]` |
| CAB | Existing pyrene builder | Two edge carbonyls balancing the π graph |
| COOH | Existing pyrene builder | One neutral edge carboxylic acid |
| COO | Existing pyrene builder | One edge carboxylate, formal charge −1 |

EPH/EPT are identified from actual explicit hydrogen neighbors, rather than a
residue name. EPH uses `EPH_C`, `EPH_CH`, `EPH_H` and `EPH_O`; EPT uses two
`EPT_C`, two `EPT_H` and one `EPT_O`. These roles match the modular partitioner.
They preserve the parent's edge H counts and pass the same residual π matching
gate. Basal EPX has attached-H counts `[0, 0]` and remains a distinct class.
The three non-epoxide controls reuse the existing `go_param.fragments` builder
without modifying it or its prior outputs. They are edge coverage controls,
not extra factorial distance configurations or fitted library parameters.

## Regeneration and artifacts

```sh
PYTHONPATH=src /path/to/rdkit-python -m go_param.mixed_panel generatedmodels_panel
```

Each model directory contains `model.json`, `initial.sdf` and `initial.xyz`.
`catalog.json` records models, measured separations, prospective splits, controls,
unattainable cells and warnings. Atom IDs, roles, formal charges, symmetry classes,
stereochemical tags and geometry provenance follow the existing explicit model
schema. The module invokes no QM programs, scheduler or RESP fitter.
