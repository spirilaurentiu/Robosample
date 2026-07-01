# Simbody - check fixtures

The paper is mostly architectural; concrete numeric fixtures are limited but useful.

## Accuracy knob (Eq. 5)

- given accuracy request of n = 4 digits, expect alpha = 10^-4 = 0.0001 (≈ 0.01% relative error).
- given real-time use, expect alpha in the range 0.01 to 0.10 (1-10%).
- alpha plays the role of the integrator relative tolerance (rtol); there is NO separate
  absolute tolerance (atol) - scaling is folded into W and T weighting matrices.

## Performance scaling (Fig. 3)

- given a 4th-order integrator, expect number of steps (and CPU time) to scale as alpha^{-1/4}.
- fitted example curve: CPU_time ≈ 0.118 * alpha^{-1/4} seconds.
- Fixture geometry: 11 chains of 20 bodies each (revolute joints, random orientation),
  attached to a common oscillating base, with gravity and light damping; 20 s simulation;
  accuracies swept alpha = 10^-2 down to 10^-7. Absolute times are machine-dependent;
  only the alpha^{-1/4} relative scaling is the invariant to test.

## Computational complexity (asserted, testable via operation counts)

- given the recursive spatial-operator solver, expect M^{-1} v to cost O(n) (not O(n^2)/O(n^3)).
- forming Y = G M^{-1} G^T costs O(mn + m^2); RHS g_0 costs O(n + m).
- solving for lambda (QTZ / pseudoinverse) costs O(m^3); total acceleration solve O(m^3 + mn + m^2).
- QTZ (complete orthogonal factorization) is stated ~5x faster than SVD for the pseudoinverse.

## Contact model sanity fixtures

- Hertz force (Eq. 23): given circular contact sigma = 1, expect f = (4/3) R^{1/2} E* x^{3/2};
  force is nonlinear (x^{3/2}) in deformation even though materials are linear elastic.
- given identical materials (E1* = E2*), expect contact point P midway between surfaces:
  s1 = s2 = 1/2, x1 = x2 = x/2.
- given one body infinitely stiffer, expect P at the stiff (non-deforming) surface (its split fraction -> 0).
- Hunt-Crossley (Eq. 25): total normal force f_Hz + f_dissipation is clamped >= 0
  (f_dissipation = max(f_HC, -f_Hz)); no pulling/adhesive force is ever produced.
- coefficient of restitution model: e = 1 - c*v for low impact speed v (c = negated slope of
  restitution-vs-impact-speed curve near v=0); restitution is NOT a material property, c is.
- eccentricity factor sigma (Eqs. 28-29) computed via [70] approximations is accurate to
  5 decimal places; circular contact (k = a/b = 1) gives sigma = 1.

## Combining-rule consistency check (Appendix A)

- Hertz composite modulus (2/3-power rule): E* = ( E1*^{2/3} E2*^{2/3} / (E1*^{2/3} + E2*^{2/3}) )^{3/2}.
  This differs from the naive E1* E2*/(E1* + E2*).
- Elastic Foundation Model uses the LINEAR combining rule E* = E1* E2*/(E1* + E2*) instead.
- plane-strain modulus: E_i* = E_i / (1 - nu_i^2).
- deformation split fractions satisfy s1 + s2 = 1 and x1 + x2 = x.
