# DOC-VelocityDistortion findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/sampler/VelocityDistortion.cpp (setNMASoftModeFromHessian,
nmaKineticCorrection, previewBatScaling, applyBatScalingDrive) and the World.hpp
declarations (536-589, plus SamplerConfig NMA fields 98-108) carried full
contract documentation from SPLIT-W8. No residual edit required.

## Verified hypotheses (acceptance accounting -- Critical if wrong)
- The NMA distortion is applied at momentum-draw time, AFTER sqrt(M^-1) seeding
  is set up: the biased mixture Us = z +/- alpha*uhat shifts the white noise, then
  the SAME u = sqrt(RT) M^-1/2 Us map runs (HmcMove.cpp:314-385). Verified.
- Detailed balance is preserved by an explicit acceptance term: ke_mix = ke -
  RT ln cosh(w.mu) via nmaKineticCorrection(), subtracted in BOTH Hold_ and Hnew
  (HmcMove.cpp:405, 431, 440, 462). With alpha=0 / nullopt the correction is a
  strict 0 and H is bit-identical to plain HMC. The distortion IS accounted for
  in acceptance -- documented as the correctness seam.
- BAT-scaling drive (ScaleBendStretch family) is a POSITION map applied once
  before any MD; applyBatScalingDrive THROWS if mdSteps != 0 or distortOption !=
  ScaleBendStretch (D7 hard SHALL), and returns the D6 lnJac used in the driven
  exchange acceptance (World.hpp:544-589). Verified as a distinct family from NMA.

## Convention note
- The simtk source read per-body NMA factors from getMobodUScaleFactor but never
  populated them; the effective uScaleFactors_ default is unity, so NMA reduces
  to the plain Gaussian draw unless non-unit factors are supplied
  (World.hpp:700-704, reinitialize:323-326). Documented; not a defect, a
  documented degenerate default.

## Coverage gap
- TestNMALinearAlgebra and TestBatScalingJacobian cover the mode basis and the
  BAT Jacobian (the inputs); no direct oracle exercises the momentum-accounting
  (ke_mix) end-to-end for HMC detailed balance. Recorded as a gap.

## Assumed notes
None.
