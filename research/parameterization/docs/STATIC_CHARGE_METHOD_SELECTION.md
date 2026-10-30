# Optional structure-specific charges: published method selection

Research finding, 2026-10-30. No new charge-correction mode is implemented.

**SQE+qp is a close published architectural match.** It combines type reference
charges with bond-domain redistribution. The 2021 study evaluated organic and
peptide sets against NPA charges; its reference parameters are fitted together
with the response parameters. They are not interchangeable with the user's
RESP-derived library values. See the [primary article](https://doi.org/10.1186/s13321-021-00528-w),
Methods and Reference charges. Native total-charge conservation does not impose
neutrality separately on each GO block. A GO/RESP-specific calibration and any
extra block constraints would need independent validation.

An existing local implementation is [ChargeFW2](https://github.com/sb-ncbr/ChargeFW2),
the core of [ACC III](https://doi.org/10.1093/nar/gkag379). Its official
`src/methods/sqeqp.cpp` was inspected: reference charges are read from the
parameter set, normalized to the requested total and combined with a
bond-incidence solution. The code does not simply start a generic QEq optimizer
from arbitrary imported partial charges. No calculations were submitted to
the public web service; source inspection only.

For this project, first benchmark a published parameter set as a separate
control. If useful, calibrate the same method family to accepted GO reference
data, defining which baseline values are held fixed and which parameters are
learned. Check identifiability, stability of the solve, size dependence,
protonation and missing atom/bond types. A new GO calibration remains a new
parameter model even when the solver is published and reused.

Classic QEq is a useful comparison rather than an automatic production default.
[Chen and Martínez's QTPIE paper](https://doi.org/10.1016/j.cplett.2007.02.065)
documents excessive long-range charge flow and dissociation behavior in older
models, and proposes transfer penalties. Static use does not make every
large-sheet extrapolation safe; test separation and flake-size effects.

For environment classes rather than a global redistribution solve,
[GAFF2–AIM](https://doi.org/10.1021/acs.jpcb.2c06037) is the relevant GO-specific
precedent. Its coupled charge/LJ model must be assessed together, with finite
edges and unsupported motifs treated explicitly. A context lookup and a
structure-specific equilibration solve are separate options.

Whichever method is selected, assign final charges once during preparation,
retain the baseline and the per-atom differences, and derive MD/MMGBSA
topologies from the same final assignment. Ordinary fixed-charge Amber can
store those values; accuracy with GAFF2 and adsorption observables still needs
testing. A classical redistribution solve should be called structure-specific
charge equilibration, rather than claiming that an embedded QM calculation ran.
