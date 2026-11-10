<a id="componenti-comuni-e-differenze-rispetto-a-hiergo-e-gopy"></a>
# Common components and differences from HierGO and GOPY

The comparison below divides the workflow into **seven functional components**.
This classification explains the project's choices; it does not measure
the percentage of shared code or physical accuracy.

| Component | GOPY | HierGO | GO Builder |
|---|---|---|---|
| 1. Lattice | 2D graphene generation | Replicable cells/tiles | Analytical rectangular lattice, four C atoms per cell |
| 2. Defects | Hole construction and cleanup | Vacancies, holes and topological defects | Vacancies/holes and dangling-atom cleanup; no Stone–Wales or grain boundary defects |
| 3. Oxygen-containing groups | COOH, epoxides, OH and variants | Basal/edge decoration, fresh/aged models | OH/epoxide/carbonyl/COOH/COO counts; random or fixed sites |
| 4. Molecular grafts | Dedicated PEG procedure | The compared workflow concerns GO decoration | Explicit fragments; anti opening with O retained as OH |
| 5. Periodicity | No equivalent workflow verified here | Tile construction and stitching | Bond images in the graph; Amber export of infinite lattices blocked |
| 6. MD parameters | Libraries/frcmod in the Amber workflow | Atom typing and OPLS assignment | Nine user GAFF2 libraries and GO_gaff2.frcmod; LEaP with explicit bonds |
| 7. GROMACS | Amber conversion through ParmEd | Direct OPLS workflow | Conversion of the same Amber topology and numerical round-trip comparison |

The first three components share their purpose with both tools.
Periodicity is an additional common theme with HierGO; PEG grafts
and the Amber/ParmEd workflow are additional common themes with GOPY. Algorithmic and
parameter details remain distinct.
Sources: [GOPY](https://github.com/Iourarum/GOPY),
[HierGO](https://github.com/IFM-molecular-simulation-group/HierGO).

<a id="costruzione-e-selezione-dei-siti-nel-nuovo-codice"></a>
## Construction and site selection in the new code

The engine generates atoms and bonds together. Standard groups use random
sites unless explicit IDs are supplied, with carbon-site reservation and a trial-based search.
Each trial chooses the complete arrangement and then checks the groups in
sequence; up to 128 trials are made, without enumerating every permutation. For
grafts, count and an optional id are selected; the fraction is calculated from the count.
Capacity or steric failures produce warnings and INCOMPLETE status, without partial
models. The model is an initial geometry, not an equilibrated ensemble.

GOPY randomly selects groups for removal and constructs a specific PEG
chain. Here the epoxide oxygen is retained as OH during anti opening
with a supplied fragment: this is a different transformation.
[GOPY source](https://github.com/Iourarum/GOPY/blob/master/GOPY.py).

<a id="cariche-e-parametri-origine-e-destinazione"></a>
## Charges and parameters: source and destination

In the GOPY tutorial, group parameters are derived by analogy from
AMBER99SB and graphene studies; the workflow loads libraries/frcmod, uses
bondbydistance and converts Amber to GROMACS with ParmEd. This is not our GAFF2 set.
[Original tutorial](https://github.com/Iourarum/GOPY/blob/master/Tutorial_GO_topology.md).

In the HierGO workflow, atom typing and Oplsatom2 assign tabulated types/charges;
the guide uses VMD/TopoTools and OPLS-top.f90 for the ITP, with OPLS
bonded/nonbonded parameters to be supplied. This is not a new quantum-mechanical fit for the surface.
[Official guide](https://github.com/IFM-molecular-simulation-group/HierGO/blob/master/Generating%20an%20OPLS-AA%20topology%20file%20from%20a%20pdb%20file.docx),
[OPLS generator](https://github.com/IFM-molecular-simulation-group/HierGO/blob/master/OPLS-top.f90).

In the new MD workflow, charges and types are **read from the user's nine
GAFF2 libraries**, with GO_gaff2.frcmod. LEaP receives the atoms and all bonds:
no bondbydistance. New charges are not recalculated for each sheet,
and neutralization is not imposed to achieve a desired total.

For a GO–amino acid/PEG etc. junction, intact-epoxide parameters are not
automatically reused. Types, charges and terms for the modified chemistry
are required; MD preparation remains blocked without a complete template
with documented provenance. Constructing the product geometry alone
does not resolve this limitation.

OPLS and OLD are **audit references**, not alternatives already selectable
for MD generation. OPLS charges alone are not substituted into GAFF2:
bonded potentials, Lennard–Jones terms, combination rules and 1–4 scaling factors
are also part of the model. Gaussian/MK and RESP inputs are prepared, but the common
QM/ESP benchmark has not been run. No family has been declared superior.
See [CHARGES.md](CHARGES.md).

<a id="che-cosa-è-riutilizzato-materialmente"></a>
## What is actually reused

The generator and Amber workflow do not run the GOPY/HierGO programs.
The current libraries are tracked copies of the user's files, preserved in
resources/reference/current. Five original HierGO files are archived
separately with commit/hash/license information: LICENSE, README.md, atom-typing.py,
OPLS-top.f90 and Oplsatom2. The audit reads the Oplsatom2 table and checks its hashes;
GAFF2 MD construction does not use those charges or run those programs.

Not all original tools are implemented: in particular,
Stone–Wales, mini grain boundary, multilayer stitching, dedicated bunching,
fresh/aged models and N-doping require their own extensions.
