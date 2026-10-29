# Three new frcmod baseline candidates

The three new files under `parameters/candidates` are **unfitted baseline controls**, not production-validated force fields. They are generated from exact, hashed installed DAT sources rather than renamed copies of the user's current/OLD overrides. No charges or QM-fitted bonded/LJ values are invented.

```sh
PYTHONPATH=src python3 -m go_param.ff_candidates
```

| Candidate file | Source and intended namespace | Additional registration |
|---|---|---|
| `gaff_native/candidate.frcmod` | Installed GAFF 1.81, May 2017; selected native C/H/O types | Load compatible native GAFF atom types |
| `gaff2_native/candidate.frcmod` | Installed GAFF2 2.2.20, March 2021; selected native C/H/O types | Load compatible native GAFF2 atom types |
| `amber99sb_go_analogy/candidate.frcmod` | parm99.dat + ff99SB family; native Amber types plus isolated GX/GO epoxide-ring types | Source the supplied `atom_types.leap` before using GX/GO |

Each candidate includes nonempty MASS, BOND, ANGLE, DIHE, IMPROPER and NONBON selections. `manifest.json` records all source hashes, exact source line and numerical term provenance, declared role-to-type contract, parameter counts, and atom-ID-preserving contracts for the eight fragment models. `gaps.json` records unresolved required graph paths and chemistry-context collisions. `index.json` identifies all three files and their status. A role contract is a deliberate starting assumption, not proof of agreement with native antechamber typing.

The native candidates rematerialize selected source terms without changing their coefficients. Native GAFF uses `c3` for the basal OH carbon contract; native GAFF2 uses `c6` for the six-membered-ring OH carbon contract. Both use `cx`/`op` for epoxide C/O, `c` for carbonyl/carboxyl carbon, and equivalent `o` types for the two GAFF carboxylate oxygen atoms. Atom typing, aromaticity and resonance must still be reviewed per graph. These new contracts differ from the legacy use of `ca` or `c2` for certain carbonyl/carboxyl carbons and therefore are **not drop-in replacements for the old libraries**.

The Amber99SB branch is an explicitly hybrid GOPY-style analogy control, not a published all-GO AMBER99SB parameter set. It retains native `CA`, `CT`, `C`, `O`, `O2`, `OH`, `HO` and `HA` roles. Its new `GX` and `GO` types avoid overloading the existing Amber `CY` nitrile / `CZ` sp-carbon meanings or assuming legacy `CX` aliases are native. GX/GO mass and LJ values use Amber CT/OS analogies; ring bonded terms are projected exactly from GAFF1.81 `cx`/`op` entries into the isolated namespace. Every projection records its original types and source line. ff99SB's backbone corrections contain nitrogen paths and do not modify these selected C/H/O graph terms. Missing cross-environment terms remain explicit gaps; no 1.1 kcal/mol improper or zero torsion is inserted to hide a gap.

## Ring-context gate

An epoxide is a C–C–O triangle, not a three-carbon triangle. Adjacent disjoint epoxides can create a typed `cx-cx-cx` path whose outer carbons are not bonded; a native cyclopropane-like 60° angle must not be accepted for that path merely because its type signature matches. The contract gate flags this as `three_member_carbon_angle_used_outside_carbon_triangle`. It also flags a `cx-cx` (or GX–GX) bond whose endpoints do not share their ring oxygen as `cx_cx_bond_spans_distinct_epoxide_rings`. Numeric coverage and correct chemical context are separate requirements.

These context flags can require new, more specific atom types or graph-specific overrides. An unvalidated `c3-c3-c3` analogy could be a sensitivity control if separately named and benchmarked; it is not silently substituted here. Fragment QM Hessians and mixed-oxidation torsional scans must determine the final treatment.

## Literature and historical provenance

The user reports parmchk generation plus manual Amber-DAT analogy completion for the current/OLD GAFF-related frcmods. That historical workflow is retained as provenance; these new source-selected baselines do not claim a new QM fitting method. [Purtscher and Hofer (2024), corrected Table 2](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.27428) provides a useful aromatic-scale comparison, but its quoted GAFF2 constants differ from installed version 2.2.20. The HSEsol/DFTB graphite fits belong to their own QM/MM and nonbonded models and are not imported into GO. [GOPY](https://github.com/Iourarum/GOPY) supports the historical Amber-analogy context, not automatic validation of new ring terms. Published GAFF2/AIM data remain a separate basal OH/epoxy branch with limited transfer scope.

## Promotion gates

1. Validate each fragment's actual native atom types, chemical valence, aromaticity and protonation against the declared contract. Load the correct base family and registration; never load the three candidates together.
2. Run graph-specific parmchk2 and real LEaP coverage checks. The listed gaps are not padded with numbers. Confirm actual matched bond/angle/torsion/improper values, wildcard precedence, units and 1–4 scaling in a resolved topology.
3. Compare optimized geometries, Hessians/normal modes and targeted torsional scans with accepted QM states. An RHF-instability diagnostic charge fit is not a validated final electronic reference. LJ and charge choices need their own matched tests; new RESP charges do not validate bonded terms.
4. Rebuild the four definitive GO materials only after parameterization is decided. Preserve their experimental source targets and construction deviations separately. No final production model or winning force field follows from these baseline files alone.

## Portable tests and external system data

Full installed Amber DAT files are external system dependencies and are intentionally excluded from the distributable research snapshot. Five source-dependent integration tests explicitly skip when those files are absent; their skip reasons identify the external data. This is not evidence of an Amber integration pass. The source-free synthetic DAT reader, ring-context collision, static frcmod audit and orphan-index tests still run. With the pinned installed sources available, the integration tests run normally and verify source hashes/version-dependent behavior. ParmEd is required for its separate parameter-reader integration check; this technical parsing check does not replace actual LEaP/topology or scientific validation.
