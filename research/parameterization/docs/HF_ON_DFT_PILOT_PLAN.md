# HF-on-DFT reference pilot implementation plan

> For agentic workers: execute inline with test-driven-development and verification-before-completion.

Goal: independently test closed-shell HF/6-31G* electronic stability at accepted
DFT geometries before deriving a canonical Amber RESP reference.

Spec: AMBER_RESP_ACCEPTANCE.md. User explicitly requested continuation with
Amber-consistent RESP. This adds a reference probe, not a final charge preset.

Architecture: a stdlib-only preparation/validation helper reuses the existing
DFT parent gate and Gaussian evidence parser. A separate frozen site campaign
runs Gaussian without optimizing or changing spin. Existing controllers and
running calculations remain unchanged.

Task 1: add tests for accepted DFT parent preparation and fixed-geometry HF
evidence; reject wrong method/energy, instability, geometry/identity changes,
nonfinite values, incomplete jobs and automatic Stable=Opt. Observe RED, then
implement src/go_param/hf_fixed.py; run focused tests and full research suite.

Task 2: submit one pristine pilot with a five-day limit. Preserve the source
DFT and HF evidence independently. No HF instability bypass or automatic
promotion to production. A successful HF probe still requires MK/RESP and
transfer/charge-LJ validation.

Execution record: the 11 focused tests were observed RED before implementation,
then GREEN; the full research suite passed 218 tests. Seven short electronic
probes were submitted together as array 1408585, maximum two simultaneous
tasks, instead of serially deferring the six non-pristine controls. This tests
the chemical-context scope at low additional cost and does not run new
optimizations or charge fits. Source inputs remain frozen per campaign.

Review focus: HF and DFT provenance cannot be interchanged; normal termination
does not imply stable wavefunction; a probe cannot claim HF geometry/frequency
validation; spin and chemistry stay fixed; prior numerical inputs stay intact.
