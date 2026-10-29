# LayerOne reference presets

Four local build recipes are provided in `examples/layerone/` using the corrected
integer functional-group table supplied with the LayerOne composition slides.
The recipes are deterministic at seed 17 and leave functional-group site IDs
unspecified, so site selection remains random but reproducible. The recipes are
starting graphs, not copies of the commercial powders, atom-for-atom replicas,
parameterized models, or equilibrated structures.

The thesis reports computational model O/C targets of 0.43, 0.19, 0.07, and
0.02 for GOX8-P935, REGOX8-P980, REGOX8-P995, and REGOX8-P999, respectively
([Calza, 2026, §2.3](https://doi.org/10.48676/unibo/amsdottorato/13227)). Its
90% carboxyl deprotonation is a modeling assumption. The corrected integer
recipe gives COO⁻ fractions of 84.2%, 85.7%, 80%, and 100% in those models.
Those are the actual fractions encoded below; the small integer counts cannot
all realize exactly 90%. The presets do not silently add or remove carboxylates
to force that fraction.

## Corrected source composition

The corrected slide table lists the model carbon count and integer group counts
shown below. The oxygen count follows group stoichiometry: one O per OH,
epoxide, or carbonyl; two O per COOH or COO⁻. Dividing this O count by the
listed C count gives the thesis ratios to the displayed precision.

| Preset | Source C | OH | Epoxide | Carbonyl | COOH | COO⁻ | O from groups | O/C from counts | Source vacancy label |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GOX8-P935 | 595 | 85 | 100 | 30 | 3 | 16 | 253 | 0.4252 → 0.43 | 10% |
| REGOX8-P980 | 640 | 70 | 12 | 24 | 1 | 6 | 120 | 0.1875 → 0.19 | 5% |
| REGOX8-P995 | 663 | 31 | 3 | 1 | 1 | 4 | 45 | 0.0679 → 0.07 | 2% |
| REGOX8-P999 | 671 | 5 | 3 | 3 | 0 | 1 | 13 | 0.0194 → 0.02 | 1% |

The integer table is reused as the author-provided model composition. It is not
a unique inversion of XPS peak areas: peak assignments, oxygen stoichiometry,
carboxyl deprotonation, vacancy topology, and finite sheet size require model
assumptions. The measured Abalonyx XPS spectra and the commercial product labels
are separate source records, not direct measurements of these four atomistic
recipes. `examples/layerone/source_info.json` records these distinctions and
the public source URLs.

The product suffixes must not be interpreted as oxidation percentages. The
manufacturer describes GOX8 P935 as a powder with a mesh specification; the
manufacturer data sheet reports chemical mass fractions separately. The current
REGOX8 product pages label P980, P995, and P999 by carbon content. Those
manufacturer product claims are not the thesis computational O/C ratios.

## Builder recipes and present validation

`nx` and `ny` count the current builder's rectangular honeycomb cells, four
starting sheet carbons per cell. For its default 1.42 Å bond length, the cell
edge lengths are 4.26 × `nx` Å and 2.4595 × `ny` Å. These are not necessarily
exact 4 × 4 nm sheets; the slide's 4 × 4 nm label describes the source models.

| Preset | Recipe cells | Explicit vacancies | Built C/O | Built O/C | Current status |
|---|---:|---:|---:|---:|---|
| GOX8-P935 | 10 × 16 | 64 | no output | — | Disabled: exact placement fails after 128 attempts at seed 17. The status is a current builder limitation, not evidence that the source composition is wrong. |
| REGOX8-P980 | 12 × 14 | 0 | 649 / 120 | 0.1849 | Provisional build succeeds; differs from the source's 5% vacancy label and C=640. |
| REGOX8-P995 | 12 × 14 | 13 | 632 / 45 | 0.0712 | Provisional build succeeds; builder tip pruning removes 32 additional sheet carbons beyond the 13 requested vacancies. |
| REGOX8-P999 | 10 × 17 | 7 | 638 / 13 | 0.0204 | Provisional build succeeds; builder tip pruning removes 36 additional sheet carbons beyond the 7 requested vacancies. |

Do not describe the three provisional structures as exact reproductions of the
source carbon counts or vacancy distributions. Their group counts match the
corrected table, but finite rectangular-cell dimensions, explicit vacancy
sampling, automatic dangling-tip removal, and group-placement limits change
the realized graph. The provisional structures are unrelaxed and carry no claim
of force-field suitability.

The current P935 recipe requests all corrected group counts without explicit
site IDs, with the slide's nominal 10% vacancy count for a 640-carbon starting
sheet. The builder reports adequate pre-placement site capacities but fails to
find a collision-free arrangement within 128 attempts (best observed counts:
COOH 3 and COO⁻ 2; other groups 0). Its preset is therefore listed but disabled
for generation pending placement validation. Do not reduce the source group
counts silently or claim that the literature model is impossible based on this
search result.

## Use

The local build command can list preset readiness without importing the model
builder or writing files:

```bash
python go.py build --list-presets --json
```

Create a provisional preset with its recorded seed, or override that seed:

```bash
python go.py build --preset REGOX8-P995 --output work/rego-p995.json
python go.py build --preset REGOX8-P995 --seed 23 --output work/rego-p995-seed23.json
```

The four canonical recipe JSON files and `source_info.json` are also packaged
under `src/go_builder/data/layerone/`; the examples and packaged recipes are
kept byte-identical. A built preset records the public source metadata,
source-composition values, and source recipe separately from its generated
graph. Construction geometry and counts are not presented as measured product
composition.

Custom recipe creation remains available:

```bash
python go.py build --recipe examples/layerone/REGOX8-P995.json --output work/rego-p995.json
```

The quick-creation preset interface reports readiness and refuses to write an
output for P935 while placement validation and force-field parameterization
selection are pending. The other three recipes can be generated but remain
provisional, unrelaxed graphs; definitive rebuilds await selection of the
parameterization.
