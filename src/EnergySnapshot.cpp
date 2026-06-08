#include "EnergySnapshot.hpp"

#include <iostream>

auto EnergySnapshot::finite() const -> bool {
    // Each field is tested independently (no early return) so one call to
    // finite() surfaces every non-finite field at once.
    bool ok = true;
    if (!std::isfinite(potential)) {
        std::cerr << "\t[ERROR] Potential energy is not finite: " << potential << '\n';
        ok = false;
    }
    if (!std::isfinite(kinetic)) {
        std::cerr << "\t[ERROR] Kinetic energy is not finite: " << kinetic << '\n';
        ok = false;
    }
    if (!std::isfinite(fixman)) {
        std::cerr << "\t[ERROR] Fixman energy is not finite: " << fixman << '\n';
        ok = false;
    }
    if (!std::isfinite(logSineSqrGamma2)) {
        std::cerr << "\t[ERROR] logSineSqrGamma2 is not finite: " << logSineSqrGamma2 << '\n';
        ok = false;
    }
    if (!std::isfinite(total)) {
        std::cerr << "\t[ERROR] Total energy is not finite: " << total << '\n';
        ok = false;
    }
    return ok;
}

auto EnergySnapshot::checkTotal(SimTK::Real RT) const -> bool {
    // FIXED: original used `!=` for floating-point comparison, which almost
    // always evaluates to true due to rounding across four additions (floating-
    // point addition is not associative), making this check spuriously warn on
    // every call.
    //
    // Corrected to a relative-epsilon comparison:
    //   |total - expected| <= eps * (1 + |total|)
    // eps = 1e-6 is well within double precision (machine eps ~2e-16) yet large
    // enough to absorb rounding accumulated across four additions.
    //
    // For Cartesian snapshots: fixman = logSineSqrGamma2 = 0, so expected
    // reduces to V + K regardless of the value of RT.
    const SimTK::Real expected = potential + kinetic + fixman - (0.5 * RT * logSineSqrGamma2);
    constexpr SimTK::Real eps = 1.0e-6;
    const SimTK::Real tol = eps * (1.0 + std::abs(total));

    if (std::abs(total - expected) > tol) {
        std::cerr << "\t[ERROR] Total energy does not match sum of components: " << total << " vs "
                  << expected << '\n';
        return false;
    }
    return true;
}

auto EnergySnapshot::checkKinetic() const -> bool {
    // K = 0 is physical (all momenta zero); K < 0 is not.
    if (kinetic < SimTK::Real{0}) {
        std::cerr << "\t[ERROR] Kinetic energy is negative: " << kinetic << '\n';
        return false;
    }
    return true;
}

auto EnergySnapshot::checkTermStability(const EnergySnapshot& ref) const -> bool {
    // Modified relative change |delta| / (1 + |ref|).
    // The +1 keeps the metric well-defined when ref ~ 0.
    //
    // Explicit `-> bool` return type documents intent and prevents any
    // accidental implicit narrowing from the comparison expression.
    auto is_exploded = [](SimTK::Real val, SimTK::Real ref_val) -> bool {
        return std::abs(val - ref_val) / (1.0 + std::abs(ref_val)) > EXPLOSION_THRESHOLD;
    };

    // Accumulate into `ok` (no early return) so every exploded term is logged.
    bool ok = true;
    if (is_exploded(potential, ref.potential)) {
        // std::cerr << "\t[ERROR] Potential energy exploded: " << potential << " vs ref " << ref.potential
        //           << '\n';
        ok = false;
    }
    if (is_exploded(kinetic, ref.kinetic)) {
        // std::cerr << "\t[ERROR] Kinetic energy exploded: " << kinetic << " vs ref " << ref.kinetic <<
        // '\n';
        ok = false;
    }
    if (is_exploded(fixman, ref.fixman)) {
        // For Cartesian snapshots fixman = 0 for both current and ref, so
        // the difference is always 0 and this branch is never reached.
        // std::cerr << "\t[ERROR] Fixman energy exploded: " << fixman << " vs ref " << ref.fixman << '\n';
        ok = false;
    }
    if (is_exploded(total, ref.total)) {
        // std::cerr << "\t[ERROR] Total energy exploded: " << total << " vs ref " << ref.total << '\n';
        ok = false;
    }
    // logSineSqrGamma2 is intentionally omitted here: its magnitude is better
    // characterised by the absolute bound in checkGeometryStability().
    return ok;
}

auto EnergySnapshot::checkGeometryStability() const -> bool {
    // Both fields can be negative, so test absolute values.
    // For Cartesian snapshots both are 0.0; this check is a silent no-op.
    bool ok = true;
    if (std::abs(logSineSqrGamma2) > GEOMETRIC_LIMIT) {
        // std::cerr << "\t[ERROR] logSineSqrGamma2 too large: " << logSineSqrGamma2 << '\n';
        ok = false;
    }
    if (std::abs(fixman) > GEOMETRIC_LIMIT) {
        // std::cerr << "\t[ERROR] Fixman energy too large: " << fixman << '\n';
        ok = false;
    }
    return ok;
}

auto EnergySnapshot::hamiltonianDrift(const EnergySnapshot& ref,
                                      SimTK::Real beta,
                                      DegreesOfFreedom ndofs) const -> bool {
    // Unwrap the strong type to get a plain arithmetic value.
    // The double static_cast is the only way to extract the underlying value
    // from an enum class; its verbosity is intentional -- it makes every
    // unwrap site visible during code review.
    const auto N = static_cast<SimTK::Real>(static_cast<std::underlying_type_t<DegreesOfFreedom>>(ndofs));

    const SimTK::Real deltaH = std::abs(total - ref.total);
    const SimTK::Real drift_kt = beta * deltaH; // beta|DeltaH|, dimensionless (kT units)
    // 2*sqrt(N) kT: factor 2 gives ~95% coverage of the CLT distribution for
    // DeltaH under a well-tuned integrator.  Same formula for both backends.
    const SimTK::Real drift_limit = 2.0 * std::sqrt(N);

    if (drift_kt > drift_limit) {
        std::cerr << "\t[ERROR] Hamiltonian drift too large: " << drift_kt << " kT (limit: " << drift_limit
                  << " kT)" << '\n';
        return false;
    }
    return true;
}

auto EnergySnapshot::validate(const EnergySnapshot& ref, SimTK::Real RT, DegreesOfFreedom ndofs) const
    -> bool {
    // `&=` (not `&&`) ensures every check runs so all warnings are visible.
    bool valid = true;

    // ── Category A: unconditional rejection on numerical failure ──────────────
    valid &= finite();       // NaN / Inf in any field
    valid &= checkTotal(RT); // bookkeeping invariant; RT term vanishes for Cartesian
    valid &= checkKinetic(); // unphysical K < 0

    // std::cout << "\t[DEBUG] Energy snapshot validation: potential=" << potential << ", kinetic=" << kinetic
    //           << ", fixman=" << fixman << ", logSineSqrGamma2=" << logSineSqrGamma2 << ", total=" << total
    //           << '\n';
    // std::cout << "\t[DEBUG] Reference snapshot: potential=" << ref.potential << ", kinetic=" << ref.kinetic
    //           << ", fixman=" << ref.fixman << ", logSineSqrGamma2=" << ref.logSineSqrGamma2
    //           << ", total=" << ref.total << '\n';

    if (std::abs(potential) < 1e-6) {
        std::cerr << "\t - [ERROR] Potential energy is suspiciously low: " << potential << '\n';
        valid = false;
    } else if (std::abs(potential) > 1e6) {
        std::cerr << "\t - [ERROR] Potential energy is suspiciously high: " << potential << '\n';
        valid = false;
    } else if (std::abs(ref.potential) < 1e-6) {
        std::cerr << "\t - [ERROR] Reference potential energy is suspiciously low: " << ref.potential << '\n';
        valid = false;
    } else if (std::abs(ref.potential) > 1e6) {
        std::cerr << "\t - [ERROR] Reference potential energy is suspiciously high: " << ref.potential
                  << '\n';
        valid = false;
    } else if (std::abs(potential / ref.potential) > 10) {
        std::cerr << "\t - [ERROR] Potential energy jump: " << potential << " vs ref " << ref.potential
                  << '\n';
        valid = false;
    }

    if (ndofs > DegreesOfFreedom(0)) {
        if (std::abs(kinetic) < 1e-6) {
            std::cerr << "\t - [ERROR] Kinetic energy is suspiciously low: " << kinetic << '\n';
            valid = false;
        } else if (std::abs(ref.kinetic) < 1e-6) {
            std::cerr << "\t - [ERROR] Reference kinetic energy is suspiciously low: " << ref.kinetic << '\n';
            valid = false;
        } else if (std::abs(kinetic / ref.kinetic) > 10) {
            std::cerr << "\t - [ERROR] Kinetic energy jump: " << kinetic << " vs ref " << ref.kinetic << '\n';
            valid = false;
        }
    } else {
        if (std::abs(kinetic) > 1e-6) {
            std::cerr << "\t - [ERROR] Kinetic energy is suspiciously high: " << kinetic << '\n';
            valid = false;
        } else if (std::abs(ref.kinetic) > 1e6) {
            std::cerr << "\t - [ERROR] Reference kinetic energy is suspiciously high: " << ref.kinetic
                      << '\n';
            valid = false;
        } else if (std::abs(kinetic / ref.kinetic) > 10) {
            std::cerr << "\t - [ERROR] Kinetic energy jump: " << kinetic << " vs ref " << ref.kinetic << '\n';
            valid = false;
        }
    }

    // // ── Category B: heuristic plausibility (advisory -- see class-level note) ──
    // // For Cartesian snapshots, checkGeometryStability() and the fixman
    // // sub-check inside checkTermStability() are trivially satisfied (fields
    // // are 0).  The potential, kinetic, and total sub-checks of
    // // checkTermStability(), and hamiltonianDrift(), remain active for both
    // // backends.
    // //
    // // These can incorrectly reject valid proposals near energy barriers.
    // // Consider promoting them to diagnostic-only warnings and removing them
    // // from the hard-rejection path.
    // ok &= checkTermStability(ref);
    // ok &= checkGeometryStability();
    // ok &= hamiltonianDrift(ref, 1.0 / RT, ndofs);
    //                        ^^^^^^^^^^
    //   beta = 1/RT: correct.  Guard against RT = 0 at the call site;
    //   T -> 0 is not physically relevant for HMC but would produce Inf here.

    return valid;
}
