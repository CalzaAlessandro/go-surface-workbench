# DFT reference level and interpretation of the HF stability flag

## Actual DFT settings

The campaign uses gas-phase omegaB97X-D/6-311++G(d,p) in Gaussian 16, with
Integral=UltraFine, SCF=(Tight,XQC), Symmetry=None and tight geometry convergence.
It performs optimization and frequencies plus a separately requested Stable
check, without automatically changing spin or following an instability.
The DFT–RESP comparison potentials use the same functional and basis.

omegaB97X-D is the original long-range-corrected hybrid with empirical damped
atom–atom dispersion, not its later D3 variants. The functional's original
study benchmarks covalent and noncovalent properties; it is a rationale for
testing a GO-fragment reference, not a validation of this specific GO model:
[Chai and Head-Gordon, 2008, DOI 10.1039/B810189B](https://web.phys.ntu.edu.tw/jdchai/Papers/p8.pdf).

The choice was intended to describe the conjugated/oxidized scaffold and
weak contacts among nearby groups with a common geometry method. The
6-311++G(d,p) basis gives split-valence triple-zeta flexibility, diffuse
functions and polarization on heavy atoms and hydrogens. Diffuse functions
are especially relevant to testing charged fragments such as COO−. No GO-
specific basis-convergence or functional-ranking benchmark has yet established
that this is an optimal level. It is computationally expensive for the C54
frequency references. The final conventional RESP target remains a separate
HF/6-31G* ESP definition.

## What the RHF flag means here

All seven HF/6-31G* probes at accepted DFT geometries completed Gaussian
normally with process exit code zero. All seven reported exactly:

`The wavefunction has an RHF -> UHF instability.`

RHF constrains alpha and beta partners to the same spatial orbitals. The
stability analysis finds an energy-lowering direction when that restriction
is relaxed. This is distinct from SCF nonconvergence or a failed geometry
optimization. It already occurs for pristine pyrene, so it cannot be assigned
solely to oxygen functionalization. A relationship to the conjugated fragment
and the limitations of restricted HF is a hypothesis, not a diagnosed unique
chemical cause or proof of the actual ground-state spin.

The triplet-like stability-vector labels/S2 values describe the perturbation
space, not a separately optimized UHF determinant or a proven triplet ground
state. No spin state or ESP density has been automatically replaced. See the
[Gaussian Stable manual](https://www.conflex.co.jp/gaussian_support/stable.php).

The research worker exits with code 4 when its additional strict stability
gate fails. That worker status must not be described as a Gaussian execution
error or automatic Amber incompatibility. Conventional restricted-HF RESP
and full electronic stability are separate questions: the
[Amber input generator](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/gcrt.c)
does not add a Stable analysis to its default Gaussian route. Retain the
restricted-HF control and assess charge/parameter sensitivity rather than
selecting UHF or DFT merely to make a strict gate pass.

The complete probe evidence is in the HF-on-DFT panel snapshot and its
HF_PANEL_SUMMARY.json. No probe computes RESP charges or selects a production
force field. See AMBER_RESP_ACCEPTANCE.md for the selection requirements.
