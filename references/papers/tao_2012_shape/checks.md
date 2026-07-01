# Checks / fixtures — Tao et al. 2012, SHAPE

Regression fixtures for a SHAPE implementation. Energies in kcal/mol; total
energies include a high-frequency correction based on the expectation value of
the symplectic shadow Hamiltonian. RMSD in Å.

## Algorithm invariants (unit tests)

- Linear momentum is conserved exactly by the body-fixed COM frame construction
  (no iteration needed). `given any single SHAPE correction step, expect total
  linear momentum unchanged to machine precision`.
- Single rigid body: `given typical MD Δt, expect the angular-velocity iteration
  (eq. 26) to converge to double precision in n = 3 iterations`.
- No shared atoms across rigid bodies: `expect the outer loop to converge in
  exactly 1 iteration`.
- Convergence criterion: iterate until `‖L^rig,(n)(t+Δt/2) − L^non(t+Δt/2)‖ <
  tolerance`, with `tolerance = 1e-7` used in the paper.
- Order of accuracy: `single-step coordinate error is O(Δt^3); global (fixed t)
  error is O(Δt^2)`.
- Small-angle branch: use the Taylor expansion (eq. 14) instead of Rodrigues
  (eq. 13) when ‖θ‖ is small to avoid the 1/‖θ‖ and 1/‖θ‖^2 blow-up.
- Square-root approximation: `(R^(n))^{1/2} ≈ ½(1 + R^(n))`.

## Test systems

- System a: water box, 126 TIP3P waters; 1 water treated as rigid. Simulated 1 ns.
- System b: Trp-cage (PDB 1L2Y) + 1169 waters; Trp6 side chain and Pro20 each a
  separate rigid body. 1 ns.
- System c: MMP2 protein + substrate + 12 636 waters; single rigid body of 43
  atoms from 6 residues (His288/His292/His298/Glu289 side chains, thiirane ring,
  methylene, sulfone, Zn ion). 1 ns.
- System d: nine-residue β-hairpin peptide + 290 waters; 8 peptide planes as
  individual planar rigid bodies, 7 Cα shared by adjacent planes. 100 ps, 1 fs,
  tolerance 1e-7.
- NVE ensemble; time step 1 fs (also 1.5 and 2.0 fs for time-step study).

## Table I — Total energy (kcal/mol) and standard deviation

| System / method | Δt=1.0 E | Δt=1.0 σ | Δt=1.5 E | Δt=1.5 σ | Δt=2.0 E | Δt=2.0 σ |
|---|---|---|---|---|---|---|
| a (SHAPE) | -665.49 | 0.017 | -661.78 | 0.104 | -661.44 | 0.86 |
| a (SHAKE) | -665.49 | 0.017 | -661.72 | 0.092 | -651.14 | 0.58 |
| a (no constraint) | -664.12 | 0.016 | -660.22 | 0.081 | -651.05 | 0.61 |
| b (SHAPE) | -7681.02 | 0.051 | -7631.59 | 0.26 | -7505.21 | 6.6 |
| b (no constraint) | -7783.08 | 0.040 | -7733.09 | 0.22 | -7607.02 | 3.8 |
| c (SHAPE) | -111021.16 | 2.8 | -110294.37 | 17.0 | -108361.57 | 192.8 |
| c (no constraint) | -111078.77 | 2.9 | -110341.24 | 20.2 | -108436.01 | 165.6 |

Key check: `given system a at Δt=1.0 fs, expect SHAPE total energy = -665.49 with
σ = 0.017 — identical to SHAKE (-665.49, 0.017)`. SHAPE and SHAKE are
interchangeable for this rigid water at small time step.

## Table II — RMSD (Å) of rigid bodies (average / std dev)

RMSD of rigid part vs. first frame; only rigid atoms used.

| System / method | Δt=1.0 avg | Δt=1.0 std | Δt=1.5 avg | Δt=1.5 std | Δt=2.0 avg | Δt=2.0 std |
|---|---|---|---|---|---|---|
| a (SHAPE) | 2.3e-7 | 2.2e-7 | 1.1e-6 | 9.6e-7 | 1.0e-6 | 7.5e-7 |
| a (SHAKE) | 4.2e-7 | 3.6e-7 | 7.7e-7 | 7.2e-7 | 1.3e-6 | 9.4e-7 |
| b Trp6 side chain (SHAPE) | 2.3e-7 | 9.2e-8 | 3.1e-7 | 1.2e-7 | 5.1e-7 | 2.8e-7 |
| b Pro20 (SHAPE) | 2.8e-7 | 5.5e-8 | 3.0e-7 | 4.8e-8 | 6.9e-7 | 3.4e-7 |
| c (SHAPE) | 6.7e-7 | 5.8e-8 | 6.6e-7 | 6.2e-8 | 6.9e-7 | 6.2e-8 |

Key check: `given SHAPE on any test rigid body, expect rigid-body RMSD ~1e-7 to
1e-6 Å (i.e. rigidity maintained to near machine precision)`.

## Table III — Outer-loop test, system d (β-hairpin, 100 ps, 1 fs, tol 1e-7)

| Method | Energy fluctuation (kcal/mol) | Energy drift (kcal/mol·ps) | Avg iterations | Max iterations |
|---|---|---|---|---|
| SHAPE (Cα shared) | 0.10 | 0.0034 | 18.1 | 27 |
| SHAKE (bonds, no planes) | 0.08 | 0.0028 | 13.3 | 20 |
| SHAPE non-sharing | 0.06 | 0.0020 | 1 | 1 |

Key checks:
- `given system d with shared Cα, expect SHAPE outer loop avg ~18 iterations
  (~5 more than SHAKE's 13.3), max 27`.
- `given system d with no shared atoms, expect exactly 1 outer-loop iteration`.
- `energy fluctuation column is the std dev of total energy (with high-frequency
  correction); energy drift is the least-squares slope of energy vs. time`.
