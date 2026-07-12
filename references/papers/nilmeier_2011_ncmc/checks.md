# Numeric fixtures - Nilmeier 2011 NCMC

All quantities in reduced units unless stated. `k_B T` is the energy unit for the bistable
dimer model.

## Bistable dimer model parameters (Eqs. 24-25)

- Double-well bonded potential (Eq. 24): given `h = 5 k_B T`, `r_0 = r_WCA`,
  `s = r_WCA / 2`, `r_WCA = 2^(1/6) * sigma`.
  - expect minima at `r = r_0` (compact) and `r = 2 r_0` (extended); barrier height `5 k_B T`.
  - expect `U_bond(r_0) = 0`, `U_bond(2 r_0) = 0`, and `U_bond(r_0 + s) = h = 5 k_B T` (top of barrier).
- WCA potential (Eq. 25): given `sigma = 3.4 Angstrom`, `epsilon = 120 k_B T`,
  `r_WCA = 2^(1/6) sigma`.
  - expect `U_WCA(r_WCA) = 0` (continuous, force zero at cutoff); `U_WCA -> +inf` as `r -> 0`;
    `U_WCA(r) = 0` for all `r >= r_WCA`.
  - expect `U_WCA` minimum value 0 attained at `r = r_WCA` (LJ shifted up by epsilon).
- Particle mass `m = 39.9 amu` (argon-like).
- WCA interaction excluded between the two bonded dimer particles.

## System / thermodynamic state

- Solvated system: 216 WCA particles total, reduced density `rho * sigma^3 = 0.96`.
- Reduced temperature `k_B T / epsilon = 0.824` for all simulations.
- Reduced time unit `tau = sqrt(sigma^2 * m / epsilon)`.
- Integrator timestep `Delta t = 0.002 tau`; collision rate `gamma = tau^-1`.
- GHMC acceptance probability at this timestep: `99.929 +/- 0.001 %`.

## Instantaneous MC move rule (Eq. 26)

- Given `r < 1.5 r_0` -> `Delta r = +r_0`.
- Given `1.5 r_0 <= r <= 3 r_0` -> `Delta r = -r_0`.
- Otherwise -> `Delta r = 0`.
- Radial Jacobian ratio `J_r = (r_new / r_old)^2`.

## Umbrella sampling (Eq. 33)

- `r_min = r_0`, `r_max = 2.05 r_0`, `K = k_B T / eta^2`, `eta = 0.3 Angstrom`.
- The `k_B T ln(r^2)` term cancels the `4 pi r^2` radial Jacobian, flattening the barrier.

## Acceptance probabilities vs switching length (dense WCA solvent)

- Instantaneous MC (`+/- r_0` proposal): acceptance `~ 10^-27` (likely underestimated).
- NCMC 1-8 steps: little to no increase over instantaneous.
- NCMC 16-1,024 steps: superlinear boost in acceptance.
- NCMC 2,048 steps: acceptance `= 12 %` (reported `gamma = 12.1 %`).
- NCMC 8,192 steps: acceptance `= 38 %`.

## Correlation times (dimer extension)

- Vacuum, MD only (500 GHMC steps/iter): `tau = 59.2` iterations.
- Vacuum, MD + instantaneous MC: `tau ~ 0.0` (uncorrelated each iteration).
- Solvent, MD only: statistical inefficiency `g ~ 600` iterations per uncorrelated sample.
- Solvent, MD + instantaneous MC: no significant reduction in `tau`.
- Solvent, MD + 2,048-step NCMC: `tau = 4.0` iterations.
- Reported `tau_MD = 299.8` iterations (solvent, MD-only correlation used in Eq. 54).

## Cost / efficiency accounting

- MD-only iteration: `T_MD = 500` force evaluations.
- MD + 2,048-step NCMC iteration: `500 + 2,048 = 2,548` force evaluations (5x cost).
- 2,048-step NCMC: 67-fold reduction in correlation time; order-of-magnitude net efficiency gain.
- Efficiency gain `E` (Eq. 23): minimum `86.9 %` of MD-alone at 128 steps (slight loss);
  plateau `~ 13x` for 2,048-4,096-step NCMC; lower gain at 8,192 steps.
- `tau_eff` diminished only when `tau_NCMC ~ tau_MD`, roughly `>= 256` switching steps.

## Eq. 54 self-consistency check (worked example in the paper)

- Given `tau_MD = 299.8` iterations and `gamma = 12.1 %` (2,048-step NCMC):
  - `tau_NCMC ~ -1 / ln(1 - 2*0.121)` (Eq. `tau_ncmc`).
  - `tau_eff = tau_MD * tau_NCMC / (tau_MD + tau_NCMC)` (Eq. 54) -> expect `tau_eff ~ 4.0`.
- Directly measured from a 10,000-iteration simulation: `tau_eff = 4.0` (agreement).

## Statistics collection (Fig. 4 protocol)

- 10,000 iterations of 2,048-step NCMC, 500 GHMC steps between NCMC trials.
- Error bars: 1,000 bootstrap resamples of the 10,000 work samples; 95% confidence intervals.
- Log mean acceptance via log-sum-exp (Eq. 31) with shift `b = max_n a_n`, `a_n = ln A(X_n)`.
