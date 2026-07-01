# Checks / fixtures - Flores 2011, RNABuilder

Benchmark numbers for regression-testing a port. Most are end-to-end folding
results (hard to reproduce exactly) but the potential-shape and sphere-parameter
fixtures are directly checkable.

## Base-pairing potential shape (eqs 5, 6, 9)

- Radial shape function g (eq:6) is continuous at r=c: short branch
  `-k·c²/(2c²) + 3k/2 = -k/2 + 3k/2 = k`; long branch `k·c/c = k`.
  Given `r = c`, expect `g = k` (the inflection-point / well-depth value).
- g at r=0: `g(0) = 3k/2`.
- Derivative g' (eq:9) at r=c is continuous: short branch `-k·c/c² = -k/c`;
  long branch `-k·c/c² = -k/c`. Given `r = c`, expect `g' = -k/c`.
- Potential (eq:5) at θ=0: `U(r,0) = g(r,k,c)·m`. At the inflection r=c, θ=0:
  `U = k·m` (per paper: "U = k" up to scaling m; Fig. 5 tick marks).
- θ dependence is quadratic everywhere; domain `-π < θ < π`.

## Sign conventions

- Given a base-pairing interaction: `κ > 0` (typically), `k < 0` (typically).

## Steric sphere parameters (Fig. 6, reduced/SelectedAtoms scheme)

- Reduced scheme default spheres: `1.75 Å` on P and C4*; `1.35 Å` on glycosidic N.
- Resulting single helix: within `1.88 Å` RMSD of an idealized helix (NAB / make_na).
- Full scheme (AllHeavyAtomSterics): `1.34 Å` spheres on all atoms except hydrogen;
  resulting helix within `1.05 Å` RMSD.
- No sterics at all: `1.96 Å` RMSD.
- Fig. 6 helix-generation params: temperature `1 K`, simulation time `10 ns`,
  forceMultiplier (m) `10`.

## P4/P6 folding application (Section 3, Fig. 9, Fig. 10)

- 160-base RNA (P4/P6 domain of *Tetrahymena* group I intron).
- 4 runs with randomized initial velocities, RMSD averaged over last ns:
  `9.3, 10.1, 11.2, 10.3 Å` (best/reported model = 10.2 Å in abstract).
- ~`6 Å` lower RMSD than best previous computational prediction.
- Prior best (NAST): `16.3 Å` RMSD, requiring `300 cpu-hours`.
- Each RNABuilder run: ~`10.5 hours` on one core (Intel Nehalem).
- Total simulation time per run: `9.6 ns`.
- P5abc folded separately, then pulled to the rest after `4.8 ns`.
- From `2.9 ns` to `6.7 ns`: base-pairing off `6 ps` / on `114 ps`, repeatedly
  (escape kinetic traps); contacts mostly converged after `6.7 ns`.
- Contacts enforced: 3 runs enforced all `166` base-pairing contacts; 1 run
  enforced `162`.
- Fig. 10 (predicted structure) params: SelectedAtoms sterics; temperature `10 K`;
  all Amber99 terms off except bond stretching scaled to `0.1` of default;
  RNABuilder base-wise forces scaled by `20`.

## Complexity & memory

- Time per step scales approximately O(n) in number of mobilities (rigid,
  flexible, and flexible+AllHeavyAtomSterics chains all ~O(n)).
- Memory: ~`226 KBytes` per residue.
- Max chain instantiated in 32-bit mode (4 GB RAM): `>= 13,000` residues.
- Rigid-body-only two-strand test (Fig. 7): rigid strands had `12` rigid-body DOFs
  (2 bodies × 6) and no internal DOFs.

## Integrator facts

- Runge-Kutta-Merson (default): variable step, `5` force evaluations per step,
  4th-order accurate trajectory, 3rd-order accurate error estimate.
- Velocity Verlet: conserves energy in fixed-step mode.
- Fig. 7 timing runs used Velocity Verlet, `1 fs` fixed steps, no thermostat, no
  MD force field.

## Software revisions used

- RNABuilder rev `284`, Molmodel rev `650`, Simbody rev `1030`.
