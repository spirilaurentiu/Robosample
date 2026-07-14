// ============================================================================
//  VelocityDistortion.cpp - World's momentum-draw / driven-map distortion
//  concerns.
//
//  Relocated verbatim from World.cpp (SPLIT-W8, pure code motion): two
//  distortion concerns that sit around the momentum draw and the
//  driven-exchange position map:
//   1. NMA Route-B velocity distortion -- setNMASoftModeFromHessian (computes
//      the softest internal-coordinate mode and stores the per-DOF scale
//      factors) and nmaKineticCorrection (the RT ln cosh(w*mu) term that
//      restores detailed balance for the biased-mixture momentum draw).
//   2. BAT-scaling drive (docs/specs/replica-exchange-nonequilibrium-work.md)
//      -- previewBatScaling (side-effect-free preview + D6 Jacobian) and
//      applyBatScalingDrive (commit to this world's geometry, D7 hard SHALL
//      guards).
//  Both write only the distortion members (uScaleFactors_, nmaBias_,
//  lastDistortJacobianDetLog_, lastNScaled_).
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "NMA.hpp"
#include "robo/world/detail/nma_debug.hpp"

using robo::Real;

// ----------------------------------------------------------------------------
//  BAT-scaling drive (docs/specs/replica-exchange-nonequilibrium-work.md
//  B4/D5/D6/D7). Stage 2a: drive + Jacobian only.
// ----------------------------------------------------------------------------
World::BatScalingResult World::previewBatScaling(const std::vector<robo::Vec3>& atomPosIn,
                                                  double s,
                                                  const std::unordered_map<int, double>& anchorR,
                                                  const std::unordered_map<int, double>& anchorTheta) const {
    BatScalingResult out;
    out.atomPos = robo::applyBatScaling(model_, atomPosIn, s, anchorR, anchorTheta, out.nScaled, out.lnJac);
    return out;
}

void World::applyBatScalingDrive(double s,
                                 const std::unordered_map<int, double>& anchorR,
                                 const std::unordered_map<int, double>& anchorTheta) {
    // D7 hard SHALL: a driven world runs NO post-scale MD (the deterministic-
    // map limit the exact two-endpoint work acceptance relies on, INV-8).
    if (sampler_.mdSteps != 0) {
        throw std::logic_error(
            "World::applyBatScalingDrive: mdSteps must be 0 for a driven (ScaleBendStretch) "
            "world (D7/INV-8) -- got mdSteps=" + std::to_string(sampler_.mdSteps));
    }
    if (!sampler_.distortOption.has_value() || *sampler_.distortOption != DistortOption::ScaleBendStretch) {
        throw std::logic_error(
            "World::applyBatScalingDrive: called on a world without "
            "distortOption=DistortOption.ScaleBendStretch (caller error, not a silent no-op)");
    }

    std::vector<robo::Vec3> current(model_.numAtoms);
    std::copy(state_.atomPosG(), state_.atomPosG() + model_.numAtoms, current.begin());

    const BatScalingResult result = previewBatScaling(current, s, anchorR, anchorTheta);
    lastNScaled_ = result.nScaled;
    lastDistortJacobianDetLog_ = result.lnJac;

    // Re-fit this world's internal frames/q to the scaled geometry (same path
    // every Gibbs-block handoff uses, World.hpp:930) so the driven endpoint is
    // immediately readable via getAtomsLocationsInGround() (D7: x^tau = x').
    setAtomsLocationsInGround(result.atomPos);
}

double World::setNMASoftModeFromHessian(const std::vector<double>& atomPosGFlat, double h, double zeroTol) {
    if (cartesian_) {
        return 0.0; // Route B is an internal-coordinate move; Cartesian world has nu==1 body
    }
    // Position the world at the minimized geometry. This rebuilds the rigid-body
    // frames (recomputeGeometry) and sets q<-0, u<-0, then realizePosition. The
    // minimum is therefore at q0==0 in these frames -- exactly where we want H.
    const int nA = model_.numAtoms;
    std::vector<robo::Vec3> pos(static_cast<std::size_t>(nA));
    for (int a = 0; a < nA; ++a) {
        pos[a] = robo::Vec3(atomPosGFlat[3 * a], atomPosGFlat[3 * a + 1], atomPosGFlat[3 * a + 2]);
    }
    setAtomsLocationsInGround(pos);

    robo::RouteBNMA nma = robo::computeRouteBNMA(model_, state_, h, zeroTol);

    // This is the hand-off: a length-nu vector here makes reinitialize() use it
    // (its all-ones fallback only fires when the size != nu).
    uScaleFactors_ = nma.uScaleFactors;

    return (nma.softMode >= 0) ? nma.eigval[nma.softMode] : 0.0;
}

double World::nmaKineticCorrection() {
    // RT*ln cosh(w.mu), w = M^(1/2) u / sqrt(RT), mu = nmaBias_. ke_mix = ke - this
    // equals -RT*ln g(u|q) for the symmetric biased-mixture momentum draw (up to the
    // same q-independent constant the standard kinetic term drops), so subtracting it
    // from both Hold_ and Hnew makes the Route B move detailed-balanced. Returns 0
    // (strict no-op) for ordinary HMC, so it never perturbs the standard path.
    if (sampler_.distortOption != DistortOption::NMA) {
        return 0.0;
    }
    const int nu = model_.nu;
    if (static_cast<int>(nmaBias_.size()) != nu) {
        return 0.0; // bias not yet built (no reinitialize this block) => no correction
    }
    // multiplyBySqrtM needs the articulated-body inertias at the current config; it
    // uses local scratch, so it does NOT disturb V_GB (the velocities calcKineticEnergy
    // just read). Position is already realized at every call site.
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    std::vector<Real> sqrtMu(nu);
    RobotEngine::multiplyBySqrtM(model_, state_, state_.u(), sqrtMu.data()); // M^(1/2) u
    const Real inv = Real(1) / std::sqrt(RT_);
    Real wDotMu = 0.0;
    for (int i = 0; i < nu; ++i) {
        wDotMu += (sqrtMu[i] * inv) * nmaBias_[i];
    }
    // Numerically stable ln cosh(x) = |x| + log1p(exp(-2|x|)) - ln 2.
    const Real ax = std::abs(wDotMu);
    const Real lnCosh = ax + std::log1p(std::exp(-2.0 * ax)) - std::log(2.0);
    if (nmaDebugEnabled()) {
        // w.mu via multiplyBySqrtM; at the START this must match the mu.Us printed
        // by reinitialize (a live check that sqrt(M) inverts sqrt(M^-1) on the draw).
        std::cout << "[nma]   w.mu=" << wDotMu << " (via M^1/2), RT*ln cosh(w.mu)=" << (RT_ * lnCosh) << "\n"
                  << std::flush;
    }
    return RT_ * lnCosh;
}
