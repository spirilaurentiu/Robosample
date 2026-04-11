#include "EnergySnapshot.hpp"

bool EnergySnapshot::finite() const {
    bool ok = true;
    if (!std::isfinite(potential)) {
        std::cout << "\t[WARNING] Potential energy is not finite: " << potential << std::endl;
        ok = false;
    }
    if (!std::isfinite(kinetic)) {
        std::cout << "\t[WARNING] Kinetic energy is not finite: " << kinetic << std::endl;
        ok = false;
    }
    if (!std::isfinite(fixman)) {
        std::cout << "\t[WARNING] Fixman energy is not finite: " << fixman << std::endl;
        ok = false;
    }
    if (!std::isfinite(logSineSqrGamma2)) {
        std::cout << "\t[WARNING] logSineSqrGamma2 is not finite: " << logSineSqrGamma2 << std::endl;
        ok = false;
    }
    if (!std::isfinite(total)) {
        std::cout << "\t[WARNING] Total energy is not finite: " << total << std::endl;
        ok = false;
    }
    return ok;
}

bool EnergySnapshot::checkTotal(SimTK::Real RT) const {
    if (total != potential + kinetic + fixman - (0.5 * RT * logSineSqrGamma2)) {
        std::cout << "\t[WARNING] Total energy does not match sum of components: " << total << " vs "
                  << (potential + kinetic + fixman - (0.5 * RT * logSineSqrGamma2)) << std::endl;
        return false;
    }
    return true;
}

bool EnergySnapshot::checkKinetic() const {
    if (kinetic < 0) {
        std::cout << "\t[WARNING] Kinetic energy is negative: " << kinetic << std::endl;
        return false;
    }
    return true;
}

bool EnergySnapshot::termExplosion(const EnergySnapshot& ref) const {
    auto is_exploded = [](SimTK::Real val, SimTK::Real ref_val) {
        return std::abs(val - ref_val) / (1.0 + std::abs(ref_val)) > EXPLOSION_THRESHOLD;
    };

    if (is_exploded(potential, ref.potential)) {
        std::cout << "\t[WARNING] Potential energy exploded: " << potential << " vs ref " << ref.potential
                  << std::endl;
        return false;
    }
    if (is_exploded(kinetic, ref.kinetic)) {
        std::cout << "\t[WARNING] Kinetic energy exploded: " << kinetic << " vs ref " << ref.kinetic
                  << std::endl;
        return false;
    }
    if (is_exploded(fixman, ref.fixman)) {
        std::cout << "\t[WARNING] Fixman energy exploded: " << fixman << " vs ref " << ref.fixman
                  << std::endl;
        return false;
    }
    if (is_exploded(total, ref.total)) {
        std::cout << "\t[WARNING] Total energy exploded: " << total << " vs ref " << ref.total << std::endl;
        return false;
    }
    return true;
}

bool EnergySnapshot::geometryExplosion() const {
    // if (std::abs(logSineSqrGamma2) > GEOMETRIC_LIMIT) {
    //     std::cout << "\t[WARNING] logSineSqrGamma2 too large: " << logSineSqrGamma2 << std::endl;
    //     return false;
    // }
    // if (std::abs(fixman) > GEOMETRIC_LIMIT) {
    //     std::cout << "\t[WARNING] Fixman energy too large: " << fixman << std::endl;
    //     return false;
    // }
    return true;
}

bool EnergySnapshot::hamiltonianDrift(const EnergySnapshot& ref, SimTK::Real beta, std::size_t ndofs) const {
    const SimTK::Real deltaH = std::abs(total - ref.total);
    const SimTK::Real drift_kt = beta * deltaH;
    const SimTK::Real drift_limit = 2.0 * std::sqrt(static_cast<SimTK::Real>(ndofs));

    if (drift_kt > drift_limit) {
        std::cout << "\t[WARNING] Hamiltonian drift too large: " << drift_kt << " kT (limit: " << drift_limit
                  << " kT)" << std::endl;
        return false;
    }
    return true;
}

bool EnergySnapshot::validate(const EnergySnapshot& ref, SimTK::Real RT, std::size_t ndofs) const {
    bool ok = true;
    ok &= finite();
    ok &= checkTotal(RT);
    ok &= checkKinetic();
    // ok &= termExplosion(ref);
    // ok &= geometryExplosion();
    // ok &= hamiltonianDrift(ref, 1 / RT, ndofs);
    return ok;
}
