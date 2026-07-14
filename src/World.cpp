// ============================================================================
//  World.cpp - SimTK-free World.
//
//  CHANGES in this revision:
//   * boostMDSteps removed (was dead).
//   * Fixman potential + logSineSqrGamma2 wired into the acceptance Hamiltonian
//     of torsional worlds, summed over ALL bodies / ALL root bodies (the old
//     code applied logSineSqrGamma2 to molecule 0 only -- a multi-molecule bug).
//   * Rigid-kick docking move (RigidKick): a symmetric Cartesian proposal,
//     Metropolised on dU, optionally followed by MD-HMC relaxation.
//   * Root mobility is supplied per-world via buildModel's rootMobilities arg.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

#include "engine_helpers.hpp" // EngineHelpers::rotationToQuaternion/quatToRotation (SPLIT-W5 hoist)

using robo::Real;
using robo::Rotation;
using robo::Transform;
using robo::Vec3;

namespace {
constexpr double kBoltzmann_kJ = 0.0083144626; // kJ/mol/K
} // namespace

World::World(int index, bool cartesian, std::uint32_t seed)
    : index_(index)
    , cartesian_(cartesian)
    , bridge_(model_)
    , rng_(static_cast<std::uint64_t>(seed) ^ (0x9e3779b97f4a7c15ULL * (index + 1))) {
}

World& World::add_sampler(double timeStep,
                          int mdSteps,
                          AcceptRejectMode mode,
                          bool useNuts,
                          double sphereFactor,
                          std::optional<bool> useFixman,
                          bool alwaysKick,
                          double clashThreshold,
                          int maxInitialKickTries,
                          std::optional<DistortOption> distortOption,
                          double nmaBiasScale) {
    sampler_.timeStep = timeStep;
    sampler_.mdSteps = mdSteps;
    sampler_.acceptRejectMode = mode;
    sampler_.useNuts = useNuts;
    sampler_.sphereFactor = (sphereFactor > 0.0) ? sphereFactor : 1.0;
    sampler_.alwaysKick = alwaysKick;
    sampler_.clashThreshold = (clashThreshold > 0.0) ? clashThreshold : 1.0e4;
    sampler_.maxInitialKickTries = (maxInitialKickTries > 0) ? maxInitialKickTries : 0;
    sampler_.distortOption = distortOption; // nullopt => no velocity distortion
    sampler_.nmaBiasScale = nmaBiasScale;   // NMA Route B bias magnitude (thermal sigmas)

    // D7 hard SHALL (docs/specs/replica-exchange-nonequilibrium-work.md,
    // INV-8): a BAT-scaling-driven world runs ZERO post-scale MD -- the
    // two-endpoint work acceptance is exact only in the pure deterministic-
    // scaling-map limit. Fail loud at construction, not silently truncate
    // like the timeStep==0/mdSteps==0 "pure proposal mode" branch below
    // (that branch is a convenience default for OTHER move types; here a
    // nonzero mdSteps is a caller error, not a mode to auto-correct).
    if (distortOption.has_value() && *distortOption == DistortOption::ScaleBendStretch && mdSteps != 0) {
        throw std::invalid_argument(
            "World::add_sampler: distort_option=ScaleBendStretch requires mdSteps=0 (D7/INV-8 hard "
            "SHALL) -- got mdSteps=" + std::to_string(mdSteps));
    }

    if (sampler_.timeStep == 0.0 || sampler_.mdSteps == 0) {
        // Pure proposal world: no internal dynamics. For docking this means a
        // rigid teleport only -- placement energy is the sole criterion.
        // timeStep=0 on a Cartesian world is the "OpenMM HMC reference" mode:
        // the docking world places the ligand, then the Cartesian world runs
        // its own integrator each round from that starting pose.
        sampler_.timeStep = 0.0;
        sampler_.mdSteps = 0;
        std::fprintf(stderr,
                     "[world %d] timeStep=0 / mdSteps=0: pure proposal mode "
                     "(no internal dynamics; %s)\n",
                     index_,
                     docking_ ? "rigid teleport only" : "OpenMM HMC reference mode");
    } else if (sampler_.timeStep > 0.005) {
        std::fprintf(stderr,
                     "[world %d] WARNING: timeStep=%.4g ps is large for all-atom MD; "
                     "the integrator may go NON-FINITE every step. "
                     "Typical stable values are 0.001-0.002 ps.\n",
                     index_,
                     sampler_.timeStep);
    }
    // AUTO default: Fixman ON for non-Cartesian (torsional + docking) worlds, OFF
    // for Cartesian (flat space => constant metric => no correction). The
    // orientational Jacobian (useOrientationJacobian) stays at its struct default
    // (OFF) -- it is the Euler-angle factor and is wrong for quaternion roots; see
    // SamplerConfig and World::calcLogSineSqrGamma2.
    sampler_.useFixman = useFixman.has_value() ? (*useFixman && !cartesian_) : !cartesian_;
    sampler_.moveType = docking_ ? MoveType::RigidKick : MoveType::MdHmc;
    return *this;
}

void World::setTemperature(double T) {
    temperature_ = T;
    RT_ = kBoltzmann_kJ * T;
    beta_ = (RT_ > 0) ? 1.0 / RT_ : 0.0;
}

void World::setTimeStep(double dt) {
    sampler_.timeStep = dt;
}

void World::setMdSteps(int n) {
    sampler_.mdSteps = n;
}

void World::setAcceptRejectMode(AcceptRejectMode mode) {
    sampler_.acceptRejectMode = mode;
}

// ----------------------------------------------------------------------------
//  Kinetic-metric preconditioning (fictitious mass; sampling only).
//  Scales the spatial inertia used ONLY in the proposal (momentum draw, KE,
//  Fixman ln det M) -> raises the stable dt ~sqrt(scale) for the scaled body
//  with ZERO configurational bias (RobotModel::bodyMassScale documents why).
//  Call AFTER buildModel(). Re-sizes defensively if the model array is missing
//  (e.g. a Cartesian world that took the compact build path).
// ----------------------------------------------------------------------------
void World::setBodyMassScale(int body, double scale) {
    if ((int)model_.bodyMassScale.size() != model_.numBodies) {
        model_.bodyMassScale.assign(model_.numBodies, robo::Real(1));
    }
    if (body >= 1 && body < model_.numBodies) {
        model_.bodyMassScale[body] = robo::Real(scale);
    }
}

void World::setMassScaleByJoint(JointType jt, double scale) {
    if ((int)model_.bodyMassScale.size() != model_.numBodies) {
        model_.bodyMassScale.assign(model_.numBodies, robo::Real(1));
    }
    for (int b = 1; b < model_.numBodies; ++b) {
        if (model_.bodyJoint[b] == jt) {
            model_.bodyMassScale[b] = robo::Real(scale);
        }
    }
}

void World::setReversibilityCheck(int interval) {
    sampler_.reversibilityCheckInterval = (interval > 0) ? interval : 0;
    if (sampler_.reversibilityCheckInterval > 0) {
        // Typically set from add_*_world() BEFORE add_sampler(), so mdSteps and
        // timeStep are not known here yet; the per-round "[rev]" line logs the
        // actual step count and dt at run time.
        std::fprintf(stderr,
                     "[world %d] reversibility probe ENABLED: every %d round(s) "
                     "(non-destructive; residual logged per round; see THEORY 5.7).\n",
                     index_,
                     sampler_.reversibilityCheckInterval);
    }
}
