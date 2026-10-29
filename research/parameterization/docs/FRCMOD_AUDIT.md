# Audit of bonded and Lennard-Jones overrides

This read-only audit separates file consistency from scientific validation. `src/go_param/frcmod_audit.py` reads the archived current GAFF2, Amber, earlier GAFF2, OLD gox/rgo and GOPY frcmod files. `reports/frcmod_audit.json` preserves every parsed source term with source line, case-sensitive types, numerical values, comments and source SHA-256; `frcmod_audit_issues.csv` lists numerical defects and review flags. Original libraries and parameter files are not modified.

```sh
PYTHONPATH=src python3 -m go_param.frcmod_audit
```

The user reports that current and OLD GAFF-related frcmods were generated with parmchk and completed manually by inspecting analogous terms in Amber DAT files. This is recorded as user-reported provenance, not independently verified line by line. Unmatched parameters are therefore not presumed fabricated or QM-fitted. The new baseline controls use similarly explicit source inheritance/analogy; they do not represent methodological novelty or a QM refit.

The local source location is `references/local_snapshot/current`, not a newly generated force field. File hashes and library atom-table hashes are recorded in the JSON audit. Reports are regenerated from those files; a source change therefore produces a different recorded hash.

## What the format checks mean

Amber frcmod files contain overrides on top of a selected base force field. Empty MASS, BOND or NONBON sections do not establish missing parameters; those values may be inherited. Conversely, a clean text audit does not establish that every molecular graph path matches a base or override term. That requires the actual atom types, LEaP load order, completed topology and resolved parameters.

The parser reads MASS, BOND, ANGLE, DIHE, IMPROPER and NONBON, preserving atom-type case. It checks numerical finiteness, positive masses/bond lengths/divisors, nonnegative harmonic constants/LJ parameters, and angle/periodicity ranges. Negative proper-dihedral periodicity is a **continuation marker**, not a negative physical periodicity: components continue until a positive periodicity ends the series. Duplicate comparisons treat whole Fourier series as definitions and do not mistake legitimate multi-term torsions for conflicts. Bond, angle and proper-dihedral reversal is canonicalized. Improper order is preserved because central-atom and wildcard matching need the Amber-specific convention.

Zero torsional coefficients can be deliberate Fourier components; they are reported with context rather than declared invalid. Hydroxyl-H `ho`/`HO` LJ (0, 0) is treated as a conventional contextual zero needing a compatible base model, not an automatically invalid negative radius. Other zero LJ values are flagged for review. Generic analogy comments, large parmchk-style penalty scores, empirical estimates and attention markers identify provenance/validation priorities; penalty scores are not calibrated uncertainties or evidence that a parameter is wrong.

Units are explicit: bond constants kcal/mol/Å², bond distances Å, angle constants kcal/mol/rad², equilibrium/phase angles degrees, Fourier amplitudes kcal/mol before proper-dihedral IDIVF division, and NONBON radius **Rmin/2 in Å**, not sigma, with epsilon in kcal/mol. MASS's optional second numeric field is not a Lennard-Jones radius. The [Amber parameter-format documentation](https://archive.ambermd.org/201301/att-0490/amber4.1_manual.pdf) describes the signed continuation convention. The [ParmEd Amber parser](https://parmed.github.io/ParmEd/html/_modules/parmed/amber/parameters.html) is an independent implementation reference for modern frcmod sections.

## Actionable findings in the immutable source snapshot

1. `GO_amber.frcmod`, line 139, has the bare IMPROPER entry `CY-CX-CX-CY` with no numerical coefficients. The same incomplete record appears in historical GOPY `GO.frcmod`, line 123. This is a concrete incomplete parameter record. Its historical provenance does not make it a valid complete improper. Determine intended matching and whether it is used before preparing a corrected candidate; do not invent a default coefficient.
2. `EPH_gaff2.lib` contains both an EPH unit with lowercase `cx`, `h1`, `op` and an appended EPX unit with uppercase `CA`, `OS`. Its `!!index` lists only EPH, so the appended EPX is an orphan raw block; this static observation does not prove that `loadOff` loads or redefines EPX. Runtime LEaP testing is required to establish actual behavior. The audit retains all unit atom tables and flags the namespace mixture. Inspect the actual loaded EPX atom types and charges; do not assume the filename defines only EPH.
3. Current `GO_gaff2.frcmod` has 14 identical repeated DIHE definition groups. They include zero- and nonzero-amplitude entries; none is a conflicting coefficient in this snapshot. The duplicates warrant deterministic consolidation in a new candidate after checking effective LEaP matching. They do not by themselves establish doubled torsional energy.
4. The current GAFF2 override contains 66 differences from the earlier `GO_gaff2_old.frcmod` across MASS/BOND/ANGLE/NONBON records. These are changes between local files, not proof that either file equals installed GAFF2 or has been refitted to QM. For example, `ca-ca` is (354.25, 1.399) currently versus (461.10, 1.398) in the earlier file; `ho-oh` is (535.51, 0.973) versus (371.40, 0.973). Bond constants are kcal/mol/Å² and lengths Å. Original provenance and actual installed DAT comparisons are necessary.
5. Epoxide-specific current angles are `cx-cx-op`: 105.290 kcal/mol/rad² at 59.100°, and `cx-op-cx`: 119.010 kcal/mol/rad² at 61.760°; the earlier file has constants 92.800 and 84.300 at the same angles. Their near-60° equilibria are consistent with a three-membered ring motif, but this observation does not validate their stiffness or transferability to densely functionalized GO.
6. A separate current `cx-cx-cx` angle uses 83.510 kcal/mol/rad² and 60.000°. Any application to a non-three-membered carbon path must be checked against actual typed graph paths. An apparently reasonable ring term can be inappropriate when a custom type is reused outside its intended ring environment.
7. OLD gox/rgo each retain `!ATTENTION` on the `c2-c-oh` angle (91.182 kcal/mol/rad², 119.120°). This is an explicit source attention marker, not evidence of a missing numerical value or a quantified error. Their empty MASS/BOND/NONBON sections are overrides requiring a compatible base force field, not self-contained force fields.

## Chemistry that text cannot validate

- Basal hydroxyl and epoxy carbons are sp³ junctions within an aromatic/sp² network. Atom-type semantics and each sp²–sp³ boundary must be checked on the actual graph, especially when native GAFF labels and legacy custom aliases share names.
- Epoxide rings require ring-specific bonded terms; transferring generic ether/sp³ angles can erase the ring geometry. Coupled distortion under adjacent oxidation needs optimized fragment geometries and Hessian/normal-mode comparison.
- COOH and COO− require separate protonation and resonance treatment. Equivalent carboxylate O atoms should give consistent resolved terms and a validated geometry; a library's nominal charges do not enforce all bonded equivalences.
- Aromatic/carbonyl planarity and sp³ chirality need the intended improper centers and matching order. A wildcard or generic 1.1 kcal/mol improper cannot be endorsed merely because it parses.
- Rotation around COOH anchors, hydroxyl C–O, and mixed sp²/sp³ interfaces requires torsional scans with the correct periodicity, phase, 1–4 scaling and charge model. Electrostatic fitting does not refit these torsions or LJ interactions.

## Three candidate branches and generation gates

Prepare **new candidates**, preserving all source libraries and overrides:

1. **Native GAFF baseline:** select an explicitly versioned installed GAFF DAT, type each intended graph natively, and run parmchk2 against that base. Retain unchanged native parameters where covered and write only graph-required gaps/overrides to a new frcmod with line-level provenance. Legacy ring aliases cannot be assumed interchangeable with native types.
2. **Native GAFF2 baseline:** repeat independently with the installed/versioned GAFF2 DAT. Do not call the current local file a GAFF2 reference merely because its filename contains GAFF2. The same topology, chemical roles and fitting/validation protocol must be used for a fair comparison.
3. **AMBER99SB-compatible GOPY analogy control:** name that specific family rather than unspecified “Amber.” Use isolated, documented custom GO ring types where protein-force-field types fail to distinguish epoxide and basal sp³ chemistry. GOPY provides historical analogies; malformed entries, generic planar impropers and undefined custom-ring terms must be resolved separately. It is a provisional compatibility control, not a validated GO force field.

For scale/provenance checks, [Purtscher and Hofer (2024), corrected Table 2](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.27428) reports aromatic GAFF/GAFF2 bond constants 461.1/378.6 kcal/mol/Å² and lengths 1.3980/1.3984 Å. Its DFTB and HSEsol parameters belong to separate fitted graphitic QM/MM models. They are not generic GO replacements; harmonic normalization and the complete associated nonbonded model must be checked before transfer.

Installed DAT comparison, graph-specific atom typing, actual missing-term resolution and complete topology tests precede promotion of any candidate. QM optimization/Hessians and targeted torsion scans must decide where a refit is needed. New RESP charges alone cannot validate a new frcmod, and a text-only audit cannot select a winning force field for adsorption.

## Installed base comparison now available

The exported, hashed system sources identify GAFF 1.81 (May 2017), GAFF2 2.2.20 (March 2021), parm99.dat and ff99SB corrections. JSON `base_comparison` now compares exact/reversed case-sensitive source keys with their installed DAT entries; no exact match is **not** a claim that a parameter is missing, since wildcard/alias fallback and runtime load precedence remain to be resolved.

Installed GAFF aromatic `ca-ca` is 461.10 kcal/mol/Å² at 1.3984 Å; installed GAFF2 is 354.25 at 1.3986 Å, with `ca-ca-ca` 63.67 kcal/mol/rad² at 120.02°. The literature comparison above is not version-identical to these installed files. Native `cx` and `op` in both installed GAFF versions describe a three-membered carbon environment and three-membered-ring oxygen; GAFF2 `c6` denotes sp³ carbon in a six-membered ring. These descriptions support type-contract choices but do not establish that a single type can distinguish an intra-epoxide C–C path from a cross-epoxide path in a dense GO network.

Three new, explicitly unfitted candidate controls are documented in [NEW_FRCMOD_CANDIDATES.md](NEW_FRCMOD_CANDIDATES.md), with their exact source/term manifests and unresolved graph/context gaps under `parameters/candidates`.
