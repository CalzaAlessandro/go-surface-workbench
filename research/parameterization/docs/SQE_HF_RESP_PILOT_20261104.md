# First SQE calibration to GO HF–RESP references

This is a bounded research fit, not a new validated Amber charge preset.
Job 1468879 completed the numerical calibration; no new QM ran in that job.

## Reference, model and split

- Targets: molecular conventional HF/6-31G* Amber RESP charges on accepted
  ωB97X-D/6-311++G(d,p) geometries. RHF wavefunction stability is not claimed.
- Baseline q0: the current user block library with its chemically correct total
  charge. The new constrained modular L0 is not the fit target or baseline.
- Response: published SQE initial-charge equations, implemented and checked
  against three native ChargeFW controls. Transfers occur on actual bonds;
  total charge is conserved separately in disconnected components.
- The fitted variables are chemical-role electronegativities. Hardness, bond
  hardness and Gaussian widths are fixed in each trial. Chemical-role typing
  does not depend on arbitrary Kekulé bond orders or the choice of the charged
  oxygen in a localized COO resonance representation.
- Default: hardness 4, bond hardness 10, width scale 1, ridge 0.0001. Fixed
  elemental Gaussian widths are C 0.7 Å, O 0.6 Å and H 0.4 Å. The response
  uses distances in Å and a 1/r kernel; these are model parameters, not fitted
  physical ionization energies. The default was declared before C24 inspection.
- Objective: molecule-normalized, role-balanced squared deviations from
  molecular RESP charges, plus ridge. This differs from minimizing QM ESP.
- Split: seven C16 contexts train; seven supported C24 contexts test size and
  context transfer. EPH/EPT are absent from training and explicitly excluded
  from that test. COOH and COO have no independent held-out validation here.
- An all16 fit is a separate in-sample candidate. Eleven sensitivity settings
  are exploratory, not an independently selected optimum.

## Observed result and decision

On the same seven C24 contexts, median per-context atomic RMS charge deviation
from molecular HF–RESP decreases from 0.162055 e for the current library to
0.121557 e for the default SQE. Median HF ESP RRMS, after averaging three rigid
orientations within each context, is:

| Candidate | Median RRMS | Evaluation scope |
|---|---:|---|
| Current library | 0.350125 | Legacy baseline |
| Modular HF–RESP L0, trained C16 | 0.179167 | C24 held out |
| SQE, trained C16 | 0.233216 | C24 held out |
| SQE, fitted all16 | 0.305417 | C24 in sample |
| Molecular HF–RESP | 0.075549 | Fit reference, in sample |

SQE improves the current baseline but does not beat the recalibrated fixed
library on this potential metric. The worse all16 potential score is retained:
fitting role-balanced atomic charges need not optimize ESP, and molecular RESP
charges are not unique observables. No accuracy claim about adsorption follows
from these scores. Prefer the simpler model unless independent new contexts
establish a robust benefit. After the expansion, compare charge-target and
direct ESP/RESP-regularized objectives without tuning on the evaluation set.

## Expanded reference campaign

Array 1468993 adds 55 starting contexts: 27 C24 and 28 C54. It varies OH/EPX
loading, CAB, COOH/COO, two acids, mixed protonation and near/far combinations,
plus EPH/EPT with hydroxyls. All 55 starting graphs are connected, radical-free,
have exact requested group counts and a complete residual π matching, and pass
the modular partition check. Constrained MMFF minimizations converged; these
are starting coordinates, not evidence of a QM minimum or a stable singlet.
The bounded search failed for one requested C24 far EPX2 case; no count was
silently reduced and no proof of impossibility is asserted.

Each case first runs DFT Opt/Freq and a stability check; accepted cases proceed
to three conventional HF MK orientations and external Amber RESP. A failed
minimum/stability/reference gate remains a rejection. Charge totals are 0,
−1 or −2 e according to protonation; they are never forced to neutrality.

## APBS and the independent reference

APBS solves continuum electrostatics from assigned atomic charges and radii;
it does not derive the electronic charge distribution. Its useful next role is
comparing potential maps, salt dependence and polar solvation for candidate
models with the same geometry and controlled numerical/dielectric settings.
Such comparisons measure sensitivity unless an independent reference is
available. They cannot validate input charges merely by propagating them.

For our primary Amber charge derivation, keep HF–RESP; use whole-sheet QM as
a separate electronic/spin/size diagnostic. Do not fit Amber charges to an
already solvent-screened PB potential. APBS is complementary, not a replacement
for QM or explicit-solvent adsorption validation. Charged GO also warrants
checking dielectric/radius sensitivity and limitations of mean-field ions.
Sources: [APBS workflow](https://apbs.readthedocs.io/en/latest/using/index.html),
[solvation and limitations](https://apbs.readthedocs.io/en/latest/background.html),
[SQE method](https://doi.org/10.1186/s13321-021-00528-w).

Reproducible tables, source checksums and scheduler snapshots accompany
`reports/QM_SQE_UPDATE_20261104_attempt_001/REPORT.html`.
