<a id="audit-numerico-dei-frcmod-go"></a>
# Numerical audit of GO frcmod files

Open `numeric_report.html`. The report is separate from the integrated report and
the charge analysis. All numerical outputs can be reproduced with:

```sh
uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache \
  --python VS_POLL_HIGH_COVERAGE/.venv/bin/python \
  PARAMETERIZATION/reports/FRCMOD_COMPARISON_20261105_attempt_001/build_report.py
```

The script uses only the Python standard library. It writes exclusively to
its own folder; it reads the seven archived frcmod files, their original snapshots,
the six local DAT files, and the two GAFF2-AIM topologies. It checks the hashes of the seven files
against MANIFEST.json, their correspondence with the originals, and the hashes of the six DAT files
against the previous manifests. The reader treats `NONB` as `NONBON`, preserves
incomplete rows without filling them in, and reassembles DIHE series using the
continuation sign. It neither imports nor modifies the production parser.

- `numeric_inventory.csv`: sections, hashes, duplicates, explicit DIHE scaling.
- `numeric_terms.csv`: every valid row with its comment, units, and original line.
- `numeric_invalid_rows.csv`: every incomplete row, separate from the coefficients.
- `numeric_duplicates.csv`: repeated instances of the same signature, compared as
  complete series for DIHE. BOND/ANGLE/DIHE reversal is canonicalized; IMPROPER order is preserved.
- `numeric_improper_permutation_overlaps.csv`: additional screening with the central
  type in the third position and permutations of the three outer types; no automatic inference
  about the effect in the loader.
- `numeric_base_comparison.csv` and `numeric_base_summary.csv`: comparisons against
  all four preserved GAFF/GAFF2 bases, or the two Amber bases for the
  uppercase namespace. The comparison does not imply simultaneous loading.
- `numeric_pairwise_terms.csv` and `numeric_pairwise_summary.csv`: native signatures
  shared or present in only one override. The absence of an override is not a demonstrated
  topology gap. AMBER uppercase is not associated with GAFF lowercase.
- `numeric_current_earlier_changes.csv`: details of local differences.
- `numeric_aim_lj.csv`: σ/ε from the original topologies and conversion to Rmin/2/ε
  in Amber units, retaining all original environment classes.
- `numeric_aim_bonds.csv`: explicit bonds from the two AIM topologies with occurrence
  counts and conversion of the harmonic constant, including the factor 1/2.
- `numeric_highlights.csv`, `numeric_flags.csv`, `numeric_qm_priorities.csv`:
  key results and chemical priorities, without an accuracy ranking.
- `numeric_summary.json`: structured inventory, results, and limitations.

`same_coefficients` indicates numerical equality after normalization of DIHE
PK/IDIVF and the absolute value of the periodicity. The constant term
and the sine/cosine coefficients of the entire series are also compared. `same_to_reported_precision`
indicates compatibility within half of the coarsest last printed digit; it is not
a physical tolerance. Mass and the optional polarizability field are
distinct. The second MASS field is not a new atomic mass.

The report does not resolve wildcards, historical loading order, the GO graph/typing,
actual improper permutations, or the final LEaP topology. It therefore cannot
certify chemical coverage, accuracy, or the absence of missing terms.
Explicit SCEE/SCNB values are reported; inherited values must be checked
in the prmtop. A preset's charges, LJ, mixing rules, and 1–4 factors must
be considered together. No new QM, MD, fitting, or production update has
been performed.
