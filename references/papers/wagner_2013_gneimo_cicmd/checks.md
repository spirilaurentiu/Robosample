# Checks / fixtures - Wagner 2013 (GNEIMO / CICMD)

Concrete numbers an implementation can regress against. Most are simulation
outcomes (soft targets), a few are exact formula fixtures.

## Exact formula fixtures

- Thermal target (eq:1): given N generalized DOF at temperature T, expect
  `Re = 0.5 * (N - 6) * k * T`. E.g. N=100, T=300 K -> `Re = 47 * k * 300 = 14100 k`.
- Modal KE identity (eq:5): given modal velocities v, expect `Re = 0.5 * sum_k v(k)^2`
  (plain sum of squares; no mass weighting in modal coords).
- B-factor / RMSF (eq:bfactor): given RMSF, expect `B = (8*pi^2/3) * RMSF^2`
  (`8*pi^2/3 ≈ 26.319`). RMSF=1 Å -> B ≈ 26.32 Å².
- CM-nulling invariant (eq:deltaV): after applying `δ_V = -V_CM` to the base cluster,
  expect system spatial momentum `h_S = 0` (to numerical tolerance).
- Flying-ice-cube fit (eq:13): `KE_CM = c * exp(-2 ln s_t) = c / s_t^2`.

## Integrator stability (torsional dynamics, single-chain Nose-Hoover, 250 fs relaxation)

- RK4 stable up to raw timestep 16 fs; Lobatto stable up to 10 fs.
- Normalized timestep (fs per force-field evaluation): RK4 does 4 force calls/step
  -> 16 fs raw = 4 fs normalized; Lobatto does 1 force call/step -> 10 fs raw = 10 fs normalized.
- Lobatto simulations begin to fail at normalized timesteps > ~9-10 fs.
- RK4 temperature-STD roughly flat for normalized timesteps < 4 fs.
- Cartesian all-atom generally unstable for timesteps > 2 fs.
- Conclusion: Lobatto (2nd order, implicit Lobatto IIIa-b) is the most efficient
  choice (largest normalized timestep) for GNEIMO.

## Flying-ice-cube diagnostic

- Nose-Hoover thermostat: bath T = 300 K, bath relaxation constant 500 fs, RK4.
- Cluster model KE_CM at 20 fs timestep grew ~10x vs 1 fs timestep (measured at 90 ps).
- Growth rate lower with GB/SA solvation than in vacuum; far larger for cluster model
  than for all-atom Cartesian.
- Mitigation: periodically null CM velocity (eq:9/eq:deltaV).

## Equilibrium MD of crystal structures (all-torsion GNEIMO, 5 ns, 310 K, Hoover 250 fs)

- Proteins: Crambin (1CRN, 1.50 Å), Defensin (1DFN, 1.9 Å), BPTI (4PTI, 1.50 Å).
- 500 ps equilibration before 5 ns production.
- Mean backbone CRMSD < 2.5 Å for most torsional-dynamics runs; not correlated with timestep.
- Cartesian sims of crambin and BPTI drifted from crystal; Cartesian 1DFN stayed near folded.

## Ab-initio folding via dynamic clustering (GNEIMO-REXMD, CVODE Adams-Moulton)

- 12 replicas, T from 300 to 1050 K, random temperature switch every 7.5 ps.
- Secondary structure (STRIDE) frozen only when REXMD temperature > 400 K.
- Test proteins (resolved subrange): 1BDD (res 11-56), 1EON (res 7-31),
  1PRB (res 11-53), 1UBQ (res 1-35).
- Best folded structures reach molten-globule state within 4-5 Å of crystal.
- Closest backbone CRMSD to crystal: 1BDD 4.007 Å; 1EON 4.198 Å; 1PRB 3.726 Å; 1UBQ 4.325 Å.
- Population peaks (backbone CRMSD): 1BDD 5-7 Å; 1UBQ 6-10 Å; 1EON 7-8 Å; 1PRB 8-10 Å.
- Run lengths: 1PRB, 1UBQ = 3 ns/replica (36 ns total); 1BDD, 1EON = 0.3 ns/replica (3.6 ns total).
- 1BDD folding path: 12-16 Å broad interhelical contacts -> 8-11 Å incorrect 3-helix
  packings -> < 7 Å correct native topology; correct topology sampled within ~40 ns total.
