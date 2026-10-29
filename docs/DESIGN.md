# GO Builder — agreed scope, 28 October 2026

The user authorized implementation of the proposed general GO generator in a new folder. Original models and parameter sources remain preserved. The later explicit request to replace active alignment inputs is documented separately with backups in ALIGNMENT.md. No Git is requested. XPS is an optional composition recipe, not an obligatory input. Use GOPY/HierGO as documented references and retain source/license provenance for any copied data.

## Architecture and choices

Use an explicit atom/bond graph, including bond image shifts at periodic boundaries, as the authoritative structure. Reuse the geometry principles of HierGO and modular chemical transformations of GOPY rather than repeat PDB-distance reconstruction. A minimal repair of the old sequence would preserve its reconstruction ambiguity; wholesale replacement of every upstream feature is unnecessary. Default to a finite, H-capped sheet for Amber and MM/GBSA. Make periodic graph construction explicit and gate exports/GBSA according to actual engine support.

Support reproducible rectangular honeycomb sheets, random vacancies and clustered holes, exact requested counts of OH/epoxide/carbonyl/COOH/COO, explicit face selection and anti nucleophilic epoxide opening. A user-defined nucleophile supplies an atom/bond fragment, attachment atom and a hydrogen to transfer to the oxygen; never guess an arbitrary reaction from a name. Preserve chemical transformations, counts and random seed in model metadata.

The 29 October interface simplification uses count and optional id; sites are random by default and the fraction is calculated from the achieved count. Explicit IDs must agree with count when both are supplied. Classic groups also accept count/id objects, with explicit carbons reserved before random decoration. Unmet loading raises PlacementError with INCOMPLETE status, WARNING, known capacity, best tried counts, search effort and suggestions. No partial model is saved. Random PEG placement is generalized to supplied fragments while retaining epoxide oxygen as OH; it does not silently replace OH/COOH or randomize every internal fragment torsion.

Amber is the primary parameterized representation: current GAFF2 fragments and dedicated frcmod terms, explicit LEaP bond commands (no bondbydistance), GB radii and residue masks. Export GROMACS from the same fully parameterized topology using ParmEd. Do not pretend unparameterized new junctions or periodically winding graphs are production-ready. Provide import of existing Amber structures without re-inferring bonds.

Audit current/OLD/OPLS charge sources and compare matched chemistry and finite fragment electrostatic benchmarks. The user's current charges originated from standard Amber parameterization of separate GO fragments; method, theory and RESP settings remain to be recovered, not invented. Prepare Gaussian/AmberTools reference jobs but do not claim a preferred force field without common QM reference and conformational checks. Gaussian and Amber remain system software.

## Acceptance

Tests cover seam bond images, connectedness, reproducibility, feasible/infeasible group counts, holes, anti geometry and atom/charge conservation, imported topology fidelity, charge-source provenance and same-chemistry comparison, periodic export/GBSA gates, honest parameterization failure, actual Amber/GROMACS round trips where software permits, and serial/Slurm CLI parity. Record limits separately from completed checks.
