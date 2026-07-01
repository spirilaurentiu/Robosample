# Toward canonical ensemble distribution from self-guided Langevin dynamics simulation

Xiongwu Wu and Bernard R. Brooks. J. Chem. Phys. 134, 134108 (2011). doi:10.1063/1.3574397

## Abstract

This work derives a quantitative description of the conformational distribution in
self-guided Langevin dynamics (SGLD) simulations. SGLD simulations employ guiding forces
calculated from local average momentums to enhance low-frequency motion. This enhancement
in low-frequency motion dramatically accelerates conformational search efficiency, but also
induces certain perturbations in conformational distribution. Through the local averaging,
we separate properties of molecular systems into low-frequency and high-frequency portions.
The guiding force effect on the conformational distribution is quantitatively described
using these low-frequency and high-frequency properties. This quantitative relation provides
a way to convert between a canonical ensemble and a self-guided ensemble. Using example
systems, we demonstrate how to utilize the relation to obtain canonical ensemble properties
and conformational distributions from SGLD simulations. This development makes SGLD not only
an efficient approach for conformational searching, but also an accurate means for
conformational sampling.

## I. Introduction

The self-guided Langevin dynamics simulation method was developed for efficient conformational
searching so that rare events, such as protein folding and ligand binding, can be accessed with
much less computing resources. While it can accelerate slow events to an affordable time scale,
the perturbation in conformational distribution from the self-guiding force remained a major
concern. For some calculations, such as free energy simulation, conformational search efficiency
is a crucial factor to obtain convergent results, while the correct conformational distribution
is responsible for accuracy.

Because the guiding force is calculated from the so-called local averages, it has been a difficult
task to quantitatively understand the effect of the guiding force on ensemble distributions. A
common practice for SGLD simulation is to limit the guiding factor to a small range so that the
effect on conformational distribution is very small and can be neglected. To obtain correct
thermodynamic average properties, Andricioaei *et al.* proposed a Monte Carlo procedure called the
momentum-enhanced hybrid Monte Carlo method to include the benefit of the guiding force while
preserving the ensemble average properties. In dynamics simulations, the difficulty in
characterizing the guiding force effect on ensemble distributions is mainly due to the lack of a
quantitative definition of the low-frequency motion to be enhanced. This work proposes a way to
separate low-frequency and high-frequency portions of thermodynamic properties through the local
averaging procedure, and derives a quantitative relation between conformational distribution and
guiding parameters.

## II. Theory and Method

### A. The low-frequency and high-frequency properties

Thermal motion in a molecular system has a distribution of frequencies. Chemical bonds vibrate and
bend at high frequencies, while ion translation and protein folding events take a relatively long
time to happen. Low-frequency events are important for many macroscopic behaviors, such as protein
folding and binding, but often are beyond the time scale accessible by molecular simulations.

We define a low-frequency property by the so-called local average property. A local averaging
procedure, typically on force or momentum, is performed by the local averaging equation over the
most recent L points, or the most recent time period `t_L = L*δt`. Here `δt` is the time interval
between data points. L is the local averaging size and `t_L` the local average time. This average
can be approximately calculated as an *evolving average* with a constant updating of the current
value, denoted with a "~" cap. Because all local averages in this work are calculated as evolving
averages, `<p>_L` is also used to represent evolving averages. Corresponding to the low-frequency
properties, high-frequency properties are defined as the difference between instantaneous properties
and their low-frequency ones: `p - p̃`.

The local averaging suppresses high-frequency effects and emphasizes low-frequency contributions.
The local average time `t_L` determines the contribution frequency range. Rearranging the evolving
average and taking `δt -> 0` gives a first-order relaxation ODE whose solution shows that a property
at any moment provides an exponentially decaying contribution to the evolving average, with decay
rate set by `t_L`.

Using `q(t) = sin(2π ϖ t)` as an example function of frequency ϖ demonstrates the behavior: for high
frequency (`2π ϖ t_L >> 1`) the amplitude of the evolving average is inversely proportional to ϖ,
while for low frequency (`2π ϖ t_L << 1`) the evolving average approximately equals the function.
The local average time `t_L` defines the separation of high vs low frequency, corresponding to a
local averaging frequency `ϖ_L = 1/t_L`. The high-frequency portion `q - q̃` keeps high-frequency
contributions while suppressing low-frequency components.

With evolving averaging, many low-frequency properties can be obtained: low-frequency forces,
low-frequency momentums, and low-frequency potential energies (all via the same evolving-average
update), and derived quantities such as low-frequency kinetic energy and low-frequency temperature.

### B. The self-guided Langevin dynamics

Langevin dynamics (LD) is based on an equation of motion with interaction force `f_i`, a friction
term `-γ_i p_i`, and a random force `R_i` whose autocorrelation obeys the fluctuation-dissipation
relation. By adding a guiding force `g_i`, we obtain the SGLD equation of motion. The guiding force
is computed from the low-frequency momentum and includes an energy-conservation factor ξ that
cancels the energy input from the guiding force so that the net work of the guiding force is zero.

### C. Conformational distribution in SGLD

The guiding force has two types of effects. First, it enhances the low-frequency motion (measured by
the increase in low-frequency temperature) and reduces high-frequency motion via the
energy-conservation force. Second, it produces a bias in the energy surface. The SGLD partition
function is split into low-frequency and high-frequency parts weighted by a low-frequency energy
factor `λ_lf` and a high-frequency energy factor `λ_hf`, with effective temperatures `T_lf` and
`T_hf` in the two spaces. When `λ = 0`, `λ_lf = λ_hf = 1` and `T_lf = T_hf = T`, recovering the LD
canonical partition function.

`λ_lf` is computed from the projection of the total low-frequency force onto the low-frequency force
direction; `λ_hf` from the analogous projection in high-frequency space. The effective temperatures
are assumed proportional to the low- and high-frequency temperatures with proportionality constants
`C_lf` and `C_hf` estimated from an LD (or `λ=0` SGLD) run where `T_lf = T_hf = T`. Here `T̃_0` is
the reference low-frequency temperature (the low-frequency temperature at `λ=0`), which depends on
simulation conditions and `t_L`.

To avoid needing a separate `λ=0` run, `T̃_0` is estimated from the same SGLD run. The low-frequency
motion can be rewritten as a Langevin dynamics with an effective collision frequency `χ_lf γ_i`,
where `χ_lf` is the low-frequency collision factor. Because the guiding force does not affect the
random force, the fluctuation-dissipation relation gives `T̃_0 = T̃ χ_lf`. This yields the SGLD
weighting factor `w_SGLD` that reweights SGLD samples back to the canonical (LD) ensemble, and any
ensemble average is computed as a weighted average.

The reweighting scheme is based on a first-order perturbation approximation and is limited to small
differences in conformational distribution; large guiding factors make reweighting hard to converge.

### D. The self-guiding temperature

The guiding factor `λ` is an input parameter whose value is hard to choose for lack of physical
meaning. A self-guiding temperature `T_sg` is defined from the effective temperatures in the low-
and high-frequency spaces. `T_sg` provides a rough measure of conformational searching ability in
units of temperature: a SGLD run with self-guiding temperature `T_sg` has search ability comparable
to a high-temperature run at `T_sg`. For LD (`T̃ = T̃_0`), `T_sg = T`. For `λ > 0`, `T̃ > T̃_0` and
`T_sg > T`; for `λ < 0`, `T_sg < T`. When `T_sg` is too large relative to `T`, reweighting becomes
inaccurate; `λ` should balance search acceleration against reweighting accuracy.

## III. Simulation Details

A leap-frog Verlet algorithm for SGLD was implemented into CHARMM version 36 (see Appendix). SGLD
adds extra calculation only in the propagation of the equations of motion, so its cost is almost
identical to normal LD for the same number of time steps. SGLD requires additional memory to store
guiding forces and weighting-factor accumulators.

## IV. Results and Discussions

Three model systems demonstrate: (1) effect of guiding forces on conformational search, (2) effect on
conformational distribution, and (3) conversion from SGLD to LD conformational distributions.

### A. The skewed double well system

A single-particle skewed double-well system (an argon atom) with two minima of different depths along
the y-axis, forcing high-frequency motion in x-z and low-frequency motion in y. Simulated at 80 K,
`t_L = 0.2 ps`, 1 fs time step, 100 ns per run, collision frequency 10/ps. SGLD (`λ=1`,
`T_sg = 100.7 K`) shows ~10x more transitions between the wells than LD (`λ=0`, `T_sg = 80 K`).
Weighting converts the SGLD energy and y-coordinate distributions back to the `λ=0` canonical
distribution. Accuracy degrades at large `λ`.

### B. Argon fluid

500 argon atoms (Lennard-Jones, `ε = 119.8 K`, `σ = 3.405 Å`) in a cubic periodic box
(28.53 × 28.53 × 28.53 Å³), 1 fs time step, 10 ns per run, 100 K, collision frequency 1/ps.
Nonbonded interactions use rationalized-polynomial 3D isotropic periodic sum (IPS) potentials.
Weighting converges the energy distributions except for `λ > 1`. SGLD raises energies far less than a
temperature increase does; comparing 100 K and 140 K LD distributions shows little overlap, whereas
SGLD stays near the target distribution. Average potential energy vs diffusion constant plots show
SGLD increases diffusion with much smaller energy deviation than high-T LD.

### C. Alanine dipeptide

Alanine dipeptide characterized by backbone dihedrals φ (CT–N–Cα–C) and ψ (N–Cα–C–NT). CHARMM
all-atom force field, distance-dependent dielectric `4r`, 100 Å nonbonded cutoff, 2 fs time step,
SHAKE on bond lengths, 200 ns per run, frames every 2 ps, `t_L = 0.2 ps`, 300 K, collision frequency
10/ps. A transition is counted when (φ,ψ) moves from within 40° of `(-90°, -70°)` to within 40° of
`(-90°, 170°)`. SGLD guiding factors 0.2, 0.5, 1 give self-guiding temperatures 346, 458, and 1067 K
respectively. Reweighting recovers the LD φ–ψ distribution, noisier at larger `λ`.

## V. Conclusions

The conformational distribution from SGLD simulation is quantitatively described through the
low-frequency and high-frequency properties, providing a way to convert SGLD distributions to
canonical ensemble distributions. SGLD can therefore achieve both dramatically enhanced conformational
search and accurate conformational distribution.

## Appendix: SGLD simulation algorithm (leap-frog Verlet)

A leap-frog Verlet SGLD algorithm. The equation labels A1–A13 are extracted in `equations.md`.

1. **Initiate low-frequency variables:** `Ẽ_p(0) = E_p(0)`, `f̃_i(0) = 0`, `p̃_i(0) = 0`,
   `g̃_i(0) = 0`.
2. **At step t** compute interaction forces `f_i(t)`, random forces `R_i(t)` (Gaussian, zero mean,
   eq A1), and the uncorrected guiding force `g'_i(t) = λ_i γ_i p̃_i(t)`. The interaction forces must
   include any constraint force. Update the low-frequency momentum from the previous half-step
   momentum (eq A2).
3. **Compute the energy-conservation factor ξ** (eqs A3–A6) and the actual guiding force `g_i(t)`
   (eq A7).
4. **Update low-frequency variables and weighting accumulators** (low-frequency forces, potential
   energy, guiding forces, low-frequency temperature; accumulators FLF, FHF, GLF, GHF, PPLF, GPLF).
   Collision and energy factors follow (eq A8). The average potential energy is subtracted from the
   low-frequency energy in `w_SGLD` (eq A9) to avoid exponential overflow.
5. **Advance velocities** to the next half step with scaling parameter `χ_i` (eqs A10, A11), then
   **advance positions** (eq A12). If constraints are needed, apply SHAKE or semiflexible constraint
   dynamics; the constraint force is recovered from the position correction (eq A13) and must be
   included in the low-frequency force.
6. **Continue** to step 2 with `t = t + δt` until the end of the simulation.
