<a id="cariche-per-il-successivo-database-dei-pollutants"></a>
# Charges for the next pollutant database

Note dated 5 November 2026, following the user's instructions. It does not modify the database,
production libraries, or parameters of calculations already started.

RESP is the reference charge model in the development of both GAFF and
GAFF2. There is no exclusive division of “RESP for GAFF, BCC for GAFF2”.
Classic AM1-BCC can be used with both. ABCG2, short for
AM1-BCC-GAFF2, is a reparameterization of BCC corrections specifically for
GAFF2, aimed primarily at experimental solvation free energies.

For the next phase of the database:

1. Evaluate **GAFF2 + ABCG2** as a candidate for small organic molecules,
   with **GAFF2 + RESP** as the reference comparison.
2. Record the charge method separately from the atom-typing/parameter method:
   `antechamber -c abcg2 -at gaff2` selects ABCG2; `-c bcc -at gaff2`
   selects classic BCC with GAFF2 types. The second combination does not
   automatically become ABCG2.
3. Check protonation, tautomers, stereochemistry, integer total charge,
   AM1 convergence, and actual coverage of the corrections. The broadest
   evidence for ABCG2 concerns neutral solutes; ions and unusual species
   require specific checks.
4. Compare molecules representative of the database families,
   including hydration and interactions with water/GO. Solvation results
   alone do not demonstrate adsorption or GBSA accuracy.
5. Preserve versions and hashes of executables, BCC tables, and the GAFF2 base,
   geometries and the RESP protocol, water model, and nonbonded conventions.
   No retroactive replacement of charges in previously studied systems.

The final choice remains open. No MD benchmarks of pollutants with ABCG2
have been run here. The GO fragment tests serve a different purpose:
they compare electrostatics against QM grids, without validating the complete
force field in solution.

Primary sources:

- [He et al., 2020, ABCG2 development](https://pmc.ncbi.nlm.nih.gov/articles/PMC7728379/).
- [He et al., 2025, extension and comparison with GAFF2/RESP](https://pmc.ncbi.nlm.nih.gov/articles/PMC11948320/).

