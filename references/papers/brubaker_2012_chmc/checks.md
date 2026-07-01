# Checks - brubaker_2012_chmc

Regression fixtures extracted from the paper. Sampler efficiency numbers (ESS,
ESS/second, timing) are hardware- and RNG-dependent and are NOT reproducible
fixtures; they are recorded for relative-ordering sanity checks only. The
deterministic, portable fixtures are the target-distribution/constraint
definitions and their parameter values.

## Special-case reductions (structural, exact)

- CHMC with $\mathcal{M}=\mathbb{R}^n$ (no constraint) == standard unconstrained HMC.
  RATTLE reduces to leapfrog. Test: with $c\equiv 0$, CHMC output distribution
  must match a reference HMC on the same $\pi$.
- CHMC with $\hat{\mathcal{H}}=\mathcal{H}$ and $L=1$ == Constrained Langevin MC.
- Unconstrained Metropolis with proposal covariance $\Sigma$ == CHMC with
  $\hat{\mathcal{H}}(p,q)=\tfrac12 p^T M^{-1}p$, $h=1$, $L=1$, $M=\Sigma^{-1}$.
- RATTLE integrator: symplectic, symmetric, order $r=2$. A symmetric integrator
  satisfies $(p,q)=\Phi_{-h}^{\mathcal{H}}(\Phi_h^{\mathcal{H}}(p,q))$ exactly
  (up to Newton-solve tolerance) - use as an integrator unit test.
- Symplecticity: $\det(F(x))^2 = 1$ where $F=\partial\Phi_h/\partial x$ - volume
  preservation check on the phase-space Jacobian.
- BvMF gauge invariance: $\pi$ with matrix $A$ and with $A+\alpha I$ must produce
  identical samples on $\mathbb{S}^{n-1}$ for any scalar $\alpha$ (since
  $q^Tq=1$). Use as a self-consistency test.

## Fixture 1 - Linearly Constrained Gaussian (Sec 4.1)

Given:
- target $\pi(q)\propto\mathcal{N}(q|\mu,\Sigma)$ subject to $c(q)=Aq-b$.
- $\mu=(0,0,0,0)^T$, $\Sigma=\mathrm{diag}(1,1,0.01,0.01)$.
- constraints $A_1=(1,1,1,1)$, $A_2=(1,1,-1,1)$, $b=(0,0)^T$.
- chains initialized at $q=(9,-9,11,-11)^T$.

Expect:
- CHMC and CLangevin converge quickly; estimated mean of first coordinate -> 0
  (the correct value).
- Because constraints are linear, this equals sampling a Gaussian in a subspace -
  the exact marginal is analytically available and can be used as ground truth.

## Fixture 2 - Bingham-von Mises-Fisher on the sphere (Sec 4.2, Table 1)

Given (Table 1 parameters):
- density $\pi(q)\propto\exp(d^Tq + q^TAq)$ on $\mathbb{S}^{n-1}$, $n=6$.
- $d=(100,0,0,0,0,0)$ (Table 1 uses $d=(100,0,0,0,0,0)$; Fig 2 caption lists a
  7-entry vector, treat the 6-entry Table-1 form as authoritative for $n=6$).
  <!-- CHECK: Fig 2 caption writes d=(100,0,0,0,0,0,0) with 7 entries but A is 6x6; likely OCR/typo, n=6. -->
- $A=\mathrm{diag}(-1000,-600,-200,200,600,1000)$.
- (Table 1 text mentions $M=2000$; steps: CHMC/CLangevin $h=1$, CMetropolis $h=0.4$.)
- Results averaged over 10 runs.

Expected efficiency table (relative ordering only; not bit-reproducible):

| Method       | $E[-\log\pi(q)]$ | ESS % | ESS/second |
|--------------|------------------|-------|------------|
| CHMC (L=4)   | -999.021         | 27.3  | 183.756    |
| CHMC (L=3)   | -998.759         | 25.4  | 217.427    |
| CHMC (L=2)   | -999.121         | 37.9  | 440.898    |
| CLangevin    | -998.757         | 33.0  | 619.339    |
| CMetropolis  | -998.82          | 3.8   | 90.1513    |
| Gibbs [Hoff] | -998.742         | 50.8  | 160.722    |

Qualitative expectations: all methods agree on $E[-\log\pi(q)]\approx-999$;
CLangevin has best ESS/second; Gibbs has highest raw ESS but loses on ESS/time;
CLangevin outperforms CHMC on this compact space; ESS decreases for $L>2$.

## Fixture 3 - Collaborative filtering (Sec 4.3, Figs 3-4)

Setup: 1M MovieLens and EachMovie, weak generalization (1 rating/user withheld);
ranks $r\in\{5,10,15\}$; mean predictions over 2000 samples; NMAE normalization
constants 1.6 (MovieLens), 1.944 (EachMovie).

MovieLens (Fig 3) RMSE / NMAE, mean $\pm$ std over ranks 5/10/15:

| Method | RMSE r5 | RMSE r10 | RMSE r15 | NMAE r5 | NMAE r10 | NMAE r15 |
|--------|---------|----------|----------|---------|----------|----------|
| HMC    | 1.577±0.39 | 2.001±0.66 | 2.306±0.25 | 0.435±0.008 | 0.465±0.016 | 0.503±0.002 |
| HMC-l  | 0.909±0.008 | 0.949±0.01 | 0.99±0.01 | 0.413±0.002 | 0.429±0.004 | 0.445±0.007 |
| CHMC   | 0.893±0.01 | 0.888±0.01 | 0.889±0.01 | 0.419±0.003 | 0.418±0.003 | 0.419±0.004 |
| CHMC-l | 0.888±0.01 | 0.881±0.01 | 0.881±0.01 | 0.418±0.004 | 0.415±0.003 | 0.416±0.002 |

EachMovie (Fig 4) RMSE / NMAE:

| Method | RMSE r5 | RMSE r10 | RMSE r15 | NMAE r5 | NMAE r10 | NMAE r15 |
|--------|---------|----------|----------|---------|----------|----------|
| HMC    | 1.153±0.002 | 1.161±0.002 | 1.204±0.018 | 0.44±0.003 | 0.44±0.003 | 0.448±0.005 |
| HMC-l  | 1.155±0.007 | 1.164±0.001 | 1.184±0.004 | 0.437±0.003 | 0.436±0.002 | 0.443±0.0015 |
| CHMC   | 1.144±0.002 | 1.121±0.001 | 1.116±0.001 | 0.444±0.003 | 0.434±0.003 | 0.432±0.002 |
| CHMC-l | 1.137±0.003 | 1.115±0.002 | 1.11±0.002 | 0.44±0.002 | 0.43±0.002 | 0.428±0.003 |

Expected: CHMC (with orthonormality constraints) beats unconstrained HMC on RMSE
in most cases and does not degrade with larger rank, whereas unconstrained HMC
overfits (RMSE grows with rank). Reference-method NMAE ranges: MovieLens
0.4342-0.3916, EachMovie 0.4520-0.4109.

## Fixture 4 - Human pose estimation (Sec 4.4)

Setup: HumanEva subject 1 walking sequence; limb-length constraints (eq:5);
density eq:6/eq:7; $L=200$ steps; mass matrix = Hessian of negative log
observations (state-dependent, gradient recomputed each step); image noise std
0-10 px; joint-to-joint error averaged over 100 frames.

Expect: CHMC/MAP clearly outperforms a projected-gradient-descent constrained
optimization baseline, which gets stuck in local optima; CHMC escapes local
minima for most frames at noise std 10.
