/**
 * @file energy_snapshot.hpp
 * @brief Energy bookkeeping and validation for HMC/Gibbs sampling, supporting
 *        both internal-coordinate (Simbody) and Cartesian (OpenMM) backends.
 *
 * -- Scientific context --------------------------------------------------------
 *
 * This struct supports Hybrid (Hamiltonian) Monte Carlo (HMC) coupled with
 * Gibbs sampling.  Two integration backends are supported:
 *
 *   • Internal coordinates (Simbody): dynamics on a curved manifold.  The
 *     Boltzmann measure picks up a metric-tensor Jacobian factor relative to
 *     the flat Cartesian measure.  Two correction terms are required:
 *
 *       1. Fixman potential  F(q) = −(½) kT ln det M(q), where M(q) is the
 *          mass-metric tensor.  Adding F restores the Cartesian invariant
 *          measure.  Reference: Fixman (1974), J. Chem. Phys. 69, 1527.
 *
 *       2. logSineSqrGamma2  = log(∏ sin²θ_i) captures the remaining
 *          geometric Jacobian from bond-angle sine factors.  Enters H with
 *          a prefactor −(½) RT.
 *
 *     The extended Hamiltonian is:
 *       H = V(q) + K(p,q) + F(q) − (½) RT · logSineSqrGamma2(q)
 *
 *   • Cartesian coordinates (OpenMM): flat Euclidean space.  No metric-tensor
 *     correction is needed.  The Jacobian is constant and cancels in all
 *     acceptance ratios:
 *       H = V(x) + K(p)        [fixman = 0, logSineSqrGamma2 = 0]
 *
 *     Use EnergySnapshot::cartesian() to construct this case; it enforces the
 *     zero values explicitly and avoids any accidental NaN in those fields.
 *
 * -- Validation philosophy -----------------------------------------------------
 *
 * The checks in this file serve TWO conceptually distinct purposes.
 * Conflating them is the root cause of the conformational-transition
 * rejection issue described below.
 *
 *   Category A — Numerical sanity (unconditional rejection)
 *     NaN / Inf values, negative kinetic energy, or bookkeeping inconsistency
 *     indicate integrator failure.  These should always trigger rejection
 *     regardless of the Monte Carlo acceptance criterion.
 *
 *   Category B — Physical plausibility heuristics (advisory only)
 *     "Suspiciously large" energy changes *might* signal a poor step size, but
 *     they can also be legitimate conformational transitions.  A proposal with
 *     large-but-finite ΔH should be handled by the Metropolis acceptance ratio
 *     min(1, exp(−β ΔH)) — not by a hard pre-filter.  Running these checks
 *     before the Metropolis step silently discards samples the chain needs for
 *     ergodicity, biasing the ensemble.
 *
 * -- The conformational-transition rejection bug -------------------------------
 *
 * The symptom — legitimate conformational transitions being rejected — is
 * caused by Category B checks (checkTermStability and hamiltonianDrift) being
 * wired into the hard-rejection path.
 *
 * A cis↔trans isomerisation, helix–coil transition, or ring flip involves
 * passing over a barrier where:
 *   • |ΔV| / (1 + |V_ref|) can vastly exceed EXPLOSION_THRESHOLD = 10
 *   • β|ΔH| can vastly exceed 2√N
 *
 * Both are physically expected at the barrier.  Correct HMC handles them via
 * Metropolis: proposals are accepted with probability exp(−β ΔH), which is
 * small but non-zero.  When a pre-filter intercepts these proposals the
 * barrier becomes impenetrable and the sampler is no longer ergodic.
 *
 * Recommendation: demote checkTermStability() and hamiltonianDrift() to
 * diagnostic warnings.  Only finite() + checkTotal() + checkKinetic() should
 * cause unconditional rejection.
 */

#pragma once

#include <SimTKcommon.h> // SimTK::Real, SimTK::NaN

#include <cmath>       // std::isfinite, std::abs, std::sqrt
#include <cstddef>     // std::size_t
#include <iostream>    // std::cerr
#include <type_traits> // std::underlying_type_t

// =============================================================================
// Strong type: DegreesOfFreedom
// =============================================================================

/**
 * @brief Strong typedef for a degree-of-freedom count.
 *
 * Using a plain `int` for `ndofs` alongside a `SimTK::Real` for
 * `beta` in the same parameter list triggers clang-tidy's
 * bugprone-easily-swappable-parameters warning because both types are
 * implicitly inter-convertible (double <-> int).  A transposed call
 * site would compile silently and produce nonsense results at runtime.
 *
 * An `enum class` wrapping `int` has the identical memory layout
 * (zero overhead) but breaks implicit conversion in both directions.  A
 * transposition now fails to compile.
 *
 * Usage at call sites:
 *   snapshot.hamiltonianDrift(ref, beta, DegreesOfFreedom{ndofs});
 *
 * Defined at file scope (not nested in EnergySnapshot) so that callers do
 * not need to qualify it through the struct name.
 */
enum class DegreesOfFreedom : int {
};

// =============================================================================
// Struct declaration
// =============================================================================

struct EnergySnapshot {
    // -- Data members ----------------------------------------------------------
    //
    // All fields default to NaN so that any use of an uninitialised snapshot
    // is caught immediately by finite().
    //
    // Consistent energy units (e.g. kJ/mol throughout) are assumed.  Mixing
    // units between fields corrupts both checkTotal() and hamiltonianDrift().

    SimTK::Real potential = SimTK::NaN;
    ///< Force-field potential energy V(q).

    SimTK::Real kinetic = SimTK::NaN;
    ///< Kinetic energy K(p, q).  Must be >= 0 for any positive-definite
    ///< mass-metric tensor.  Reduces to (1/2) p.p/m in Cartesian coordinates.

    SimTK::Real fixman = SimTK::NaN;
    ///< Fixman correction F(q) = -(1/2) kT ln det M(q).
    ///< Zero for Cartesian (OpenMM) integration -- use cartesian() factory.
    ///< Can be positive or negative depending on det M.

    SimTK::Real logSineSqrGamma2 = SimTK::NaN;
    ///< log(prod sin^2(theta_i)) -- geometric Jacobian from bond-angle sines.
    ///< Enters H with prefactor -(1/2) RT.
    ///< Zero for Cartesian (OpenMM) integration -- use cartesian() factory.
    ///< Diverges as any bond angle theta_i -> 0 or pi.

    SimTK::Real total = SimTK::NaN;
    ///< Extended Hamiltonian H = V + K + F - (1/2) RT logSineSqrGamma2.
    ///< Reduces to H = V + K for the Cartesian backend.
    ///< Stored explicitly so checkTotal() can detect bookkeeping bugs.

    // -- Class constants -------------------------------------------------------

    /**
     * Modified-relative-change threshold for checkTermStability().
     *
     * The metric |delta| / (1 + |ref|) is flagged when it exceeds this value.
     * The +1 keeps the metric well-defined when ref ~ 0.
     *
     * WARNING: 10.0 is far too tight for systems with conformational
     * transitions.  A cis-trans barrier at even 5 kJ/mol against a reference
     * near V = 0 produces a ratio >> 10 legitimately.  If used for hard
     * rejection (not just logging), raise by at least an order of magnitude
     * (>= 100) and validate against known transition rates.
     * In the Cartesian backend this threshold has no significance for the
     * fixman / logSineSqrGamma2 terms because both are always 0.
     */
    static constexpr SimTK::Real EXPLOSION_THRESHOLD = 10.0;

    /**
     * Absolute upper bound on |logSineSqrGamma2| and |fixman| used by
     * checkGeometryStability().
     *
     * 1e3 is generous for typical biomolecular systems.  In the Cartesian
     * backend both fields are 0, so checkGeometryStability() is a trivial
     * no-op there.
     */
    static constexpr SimTK::Real GEOMETRIC_LIMIT = 1e3;

    // -- Factory methods -------------------------------------------------------

    /**
     * @brief Construct a snapshot for Cartesian (OpenMM) integration.
     *
     * In flat Euclidean space the mass matrix is constant (independent of
     * coordinates), so both the Fixman correction and the bond-angle Jacobian
     * vanish identically.  Setting them to 0.0 (rather than leaving them as
     * NaN) lets the generic validation path work correctly:
     *
     *   finite()                -- 0 is finite; passes.
     *   checkTotal()            -- expected = V + K; matches total.
     *   checkGeometryStability()-- |0| < GEOMETRIC_LIMIT; trivially passes.
     *   checkTermStability()    -- fixman diff = 0; trivially passes.
     *
     * The RT parameter is not required here (the zero Jacobian terms make it
     * irrelevant), but callers must still pass the same RT to checkTotal()
     * when validating; the multiplication by 0 makes its value irrelevant.
     *
     * @param V  Potential energy from OpenMM (kJ/mol or consistent unit).
     * @param K  Kinetic energy from OpenMM (same unit).
     */
    [[nodiscard]] static auto cartesian(SimTK::Real V, SimTK::Real K) -> EnergySnapshot {
        EnergySnapshot s;
        s.potential = V;
        s.kinetic = K;
        s.fixman = 0.0;           // no metric-tensor correction in Cartesian space
        s.logSineSqrGamma2 = 0.0; // no bond-angle Jacobian in Cartesian space
        s.total = V + K;          // H = V + K; RT term vanishes because logSine = 0
        return s;
    }

    /**
     * @brief Construct a snapshot for internal-coordinate (Simbody) integration.
     *
     * All five fields are provided explicitly.  `total` is computed from the
     * canonical formula to guarantee checkTotal() passes.
     *
     * @param V          Potential energy.
     * @param K          Kinetic energy (generalised).
     * @param F          Fixman correction.
     * @param logSineSq  log(prod sin^2 theta) Jacobian term.
     * @param RT         Thermal energy k_B T (same units as V, K, F).
     */
    [[nodiscard]] static auto
    internalCoordinates(SimTK::Real V, SimTK::Real K, SimTK::Real F, SimTK::Real logSineSq, SimTK::Real RT)
        -> EnergySnapshot {
        EnergySnapshot s;
        s.potential = V;
        s.kinetic = K;
        s.fixman = F;
        s.logSineSqrGamma2 = logSineSq;
        s.total = V + K + F - (0.5 * RT * logSineSq);
        return s;
    }

    // -- Category A: numerical sanity (unconditional rejection) ----------------

    /**
     * @brief Returns true iff every energy field is finite (non-NaN, non-Inf).
     *
     * Non-finite values indicate integrator divergence or an ill-conditioned
     * geometry.  This is a hard gate: a non-finite snapshot must never reach
     * the Metropolis acceptance step.
     *
     * All five fields are checked independently so one call reports every
     * problematic field rather than stopping at the first.
     */
    [[nodiscard]] auto finite() const -> bool;

    /**
     * @brief Returns true iff `total` matches the sum of components to within
     *        floating-point rounding tolerance.
     *
     * Detects bookkeeping bugs where `total` is stale or computed differently.
     * The formula is H = V + K + F - (1/2) RT logSineSqrGamma2.
     * For Cartesian snapshots (F = logSineSqrGamma2 = 0) this reduces to
     * H = V + K, so checkTotal() works for both backends with the same path.
     *
     * BUG FIXED: original used exact `!=` comparison.  Rounding across four
     * floating-point additions means this almost always evaluates to true and
     * emits a spurious warning.  Fixed with a relative-epsilon comparison.
     *
     * @param RT  Thermal energy k_B T.  For Cartesian snapshots the RT term
     *            is multiplied by 0, so its value does not affect the result.
     */
    [[nodiscard]] auto checkTotal(SimTK::Real RT) const -> bool;

    /**
     * @brief Returns true iff the kinetic energy is non-negative.
     *
     * For any positive-definite mass tensor K = (1/2) p^T M^{-1} p >= 0.
     * A negative value indicates a bug in the momentum-draw step or an
     * ill-conditioned mass matrix.  Applies identically to both backends.
     */
    [[nodiscard]] auto checkKinetic() const -> bool;

    // -- Category B: physical-plausibility heuristics (advisory) --------------

    /**
     * @brief Returns true iff no individual energy term changed by more than
     *        EXPLOSION_THRESHOLD in the modified-relative-change sense.
     *
     * The metric |val - ref| / (1 + |ref|) is well-behaved near zero and
     * reduces to a near-absolute check when |ref| << 1.
     *
     * For Cartesian snapshots fixman and logSineSqrGamma2 are always 0, so
     * their modified-relative changes are identically 0 and trivially pass.
     * The check still provides useful diagnostics on V and K.
     *
     * RENAMED from `termExplosion`: the original returned *false* when an
     * explosion was detected -- opposite of the name's implication.  All
     * check*() methods now return true when the invariant is satisfied.
     *
     * WARNING: see class-level note on conformational-transition rejection.
     *
     * @param ref  Reference snapshot (last accepted configuration).
     */
    [[nodiscard]] auto checkTermStability(const EnergySnapshot& ref) const -> bool;

    /**
     * @brief Returns true iff both geometric correction terms are within
     *        GEOMETRIC_LIMIT in absolute value.
     *
     * Large |logSineSqrGamma2| or |fixman| signal a near-singular internal-
     * coordinate geometry (bond angle theta -> 0 or pi), where the metric
     * tensor is ill-conditioned and the integrator is unreliable regardless
     * of step size.
     *
     * For Cartesian (OpenMM) snapshots both fields are 0.0, so this check
     * always returns true -- it is a silent no-op for that backend.
     *
     * RENAMED from `geometryExplosion` for naming consistency.
     */
    [[nodiscard]] auto checkGeometryStability() const -> bool;

    /**
     * @brief Returns true iff the Hamiltonian drift beta|DeltaH| <= 2*sqrt(N) kT.
     *
     * Statistical rationale: for a harmonic system with N active degrees of
     * freedom, the variance of DeltaH across a single leapfrog step scales as
     * N*(kT)^2 by a CLT argument over independent normal modes, giving
     * sigma(DeltaH) ~ sqrt(N) kT.  A factor-of-2 buffer covers ~95% of the
     * distribution for a well-tuned integrator.
     *
     * This check applies equally to both backends.  For Cartesian integration
     * DeltaH = Delta(V + K) and the 2*sqrt(N) bound is the same CLT estimate
     * over Cartesian normal modes.
     *
     * WARNING (likely cause of conformational-transition rejection): near
     * energy barriers beta|DeltaH| routinely exceeds 2*sqrt(N).  Those
     * proposals should reach Metropolis and be accepted with probability
     * exp(-beta DeltaH), not be silently discarded here.
     *
     * -- Swappable-parameters fix ---------------------------------------------
     * The original signature placed a `SimTK::Real beta` and a `std::size_t
     * ndofs` adjacently.  clang-tidy flagged this because both types are
     * implicitly inter-convertible; a transposed call compiles silently and
     * produces nonsense.  Replaced `std::size_t` with the strong enum class
     * `DegreesOfFreedom`, which has the same layout but breaks all implicit
     * conversions, making the transposition a compile-time error.
     *
     * @param ref    Reference snapshot.
     * @param beta   Inverse temperature beta = 1/RT.
     * @param ndofs  Active degrees of freedom wrapped in DegreesOfFreedom.
     */
    [[nodiscard]] auto
    hamiltonianDrift(const EnergySnapshot& ref, SimTK::Real beta, DegreesOfFreedom ndofs) const -> bool;

    /**
     * @brief Aggregate validator: runs all checks, returns true iff all pass.
     *
     * Uses `&=` (bitwise-and-assign), not short-circuit `&&`, so that every
     * check runs even after an early failure and all warnings are emitted in
     * a single call.
     *
     * Works for both backends without any branching: for Cartesian snapshots,
     * checkGeometryStability() and the fixman/logSine sub-checks inside
     * checkTermStability() are trivially satisfied because those fields are 0.
     *
     * WARNING: because Category B checks are included, a `false` return does
     * not necessarily mean the proposal is unphysical.  Log which checks
     * failed and consult the per-method notes.
     *
     * @param ref    Reference snapshot (last accepted configuration).
     * @param RT     Thermal energy k_B T.
     * @param ndofs  Active degrees of freedom.
     */
    [[nodiscard]] auto validate(const EnergySnapshot& ref, SimTK::Real RT, DegreesOfFreedom ndofs) const
        -> bool;
};
