#include "OpenMM.hpp"

#include <chrono>
#include <vector>

#include "HarmonicAngleForce.h"
#include "HarmonicBondForce.h"
#include "Integrator.h"
#include "TopologyElements.hpp"

#if USE_CPU
#    include "../openmm/platforms/cpu/include/CpuPlatform.h"
#elif USE_REFERENCE
#    include "../openmm/platforms/reference/include/ReferencePlatform.h"
#elif USE_OPENCL
#    include "../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif USE_CUDA
#    include "../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

class MTSIntegrator : public OpenMM::CustomIntegrator {
    public:
    MTSIntegrator(SimTK::Real stepSize, std::vector<std::pair<int, int>> groups)
        : OpenMM::CustomIntegrator(stepSize) {
        if (groups.empty()) {
            throw std::invalid_argument("No force groups specified");
        }

        // Match Python: sort ascending by substep count
        std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

        addPerDofVariable("x1", 0);
        addUpdateContextState();
        createSubsteps(1, groups);
        addConstrainVelocities();
    }

    private:
    void createSubsteps(int parentSubsteps, const std::vector<std::pair<int, int>>& groups) {
        auto [group, substeps] = groups[0];

        if (substeps % parentSubsteps != 0 || substeps / parentSubsteps < 1) {
            throw std::invalid_argument("Substeps for each group must be a multiple of the previous group");
        }
        if (group < 0 || group > 31) {
            throw std::invalid_argument("Force group must be between 0 and 31");
        }

        int stepsPerParentStep = substeps / parentSubsteps;
        std::string n = std::to_string(substeps);
        std::string g = std::to_string(group);

        for (int i = 0; i < stepsPerParentStep; ++i) {
            // Half kick with this group's forces
            addComputePerDof("v", "v+0.5*(dt/" + n + ")*f" + g + "/m");

            if (groups.size() == 1) {
                // Innermost group: do the position update
                addComputePerDof("x", "x+(dt/" + n + ")*v");
                addComputePerDof("x1", "x");
                addConstrainPositions();
                addComputePerDof("v", "v+(x-x1)/(dt/" + n + ")");
                addConstrainVelocities();
            } else {
                // Recurse into faster groups
                createSubsteps(substeps, {groups.begin() + 1, groups.end()});
            }

            // Second half kick
            addComputePerDof("v", "v+0.5*(dt/" + n + ")*f" + g + "/m");
        }
    }
};

auto OPENMM::initialize(const std::vector<std::vector<int>>& worlds,
                        const SystemTopology& systemTopology,
                        const ForceFieldParams& ffParams,
                        const SimulationSettings& simSettings) -> bool {
    // Instantiate
    OPENMM& omm = get();

    // Copy data
    omm.atomMbxByWorld = worlds;
    omm.systemTopology = systemTopology;
    omm.ffParams = ffParams;
    omm.simSettings = simSettings;

    omm.numWorlds = static_cast<int>(worlds.size());
    omm.numAtomsInRigidBodiesByWorld.resize(omm.numWorlds);
    for (int wIx = 0; wIx < omm.numWorlds; ++wIx) {
        for (auto mbx : worlds[wIx]) {
            omm.numAtomsInRigidBodiesByWorld[wIx][mbx]++;
        }
    }

    // Set up variables
    omm.numAtoms = systemTopology.atoms.size();
    omm.ommAtomsPositionsCache = std::vector<OpenMM::Vec3>(omm.numAtoms);
    // omm.ommAtomsPositionsCacheOld = std::vector<OpenMM::Vec3>(omm.numAtoms);
    // omm.simbodyAtomsPositionsCache = std::vector<SimTK::Vec3>(omm.numAtoms);

    // Allocate OpenMM system and add particles to it
    omm.system = std::make_unique<OpenMM::System>();

    for (const auto& atom : systemTopology.atoms) {
        const auto moleculeIndex = atom.identity.moleculeIndex;
        const bool isRoot = atom.connectivity.root;
        const bool isWeld = systemTopology.rootMobilities[moleculeIndex] == SimTK::RootMobility::Weld;

        if (isRoot && isWeld) {
            omm.system->addParticle(0.0); // massless root
        } else {
            omm.system->addParticle(atom.physics.massInDaltons);
        }
    }

    const int FAST_GROUP = 0;
    const int SLOW_GROUP = 1;

    // Nonbonded force
    auto* nonbondedForce = omm.createNonbondedForce(systemTopology.atoms,
                                                    systemTopology.scaling14s,
                                                    systemTopology.exclusions,
                                                    ffParams,
                                                    ffParams.hasNBfix);
    // nonbondedForce->setForceGroup(SLOW_GROUP);
    omm.nonbondedForceIndex = omm.system->addForce(nonbondedForce);

    // GBSA OBC Force
    // implicitSolventKappa
    if (ffParams.useGBSAOBC2) {
        auto* gbsaOBCForce = omm.createGBSAOBCForce(systemTopology.atoms, ffParams);
        // gbsaOBCForce->setForceGroup(SLOW_GROUP);
        omm.system->addForce(gbsaOBCForce);
    }

    // 1-4 VdW Correction (Bond-based so we can target specific pairs)
    if (ffParams.hasNBfix) {
        auto* customNonbondedForce = omm.createCustomNonbondedForce(systemTopology.atoms,
                                                                    systemTopology.scaling14s,
                                                                    systemTopology.exclusions,
                                                                    ffParams,
                                                                    nonbondedForce->getUseSwitchingFunction(),
                                                                    nonbondedForce->getSwitchingDistance());
        // customNonbondedForce->setForceGroup(SLOW_GROUP);
        omm.system->addForce(customNonbondedForce);
    }

    // Add bonds
    auto* harmonicBondForce = omm.createHarmonicBondForce();
    // harmonicBondForce->setForceGroup(FAST_GROUP);
    omm.harmonicBondForceIndex = omm.system->addForce(harmonicBondForce);

    // Add angles
    auto* harmonicAngleForce = omm.createHarmonicAngleForce();
    // harmonicAngleForce->setForceGroup(SLOW_GROUP);
    omm.harmonicAngleForceIndex = omm.system->addForce(harmonicAngleForce);

    // For AMBER style torsions, OpenMM does not care if they are proper or improper
    // Note that OpenMM only supports only periodic torsions by default (AMBER style), so we need to handle
    // improper harmonic torsions differently
    auto* periodicTorsionForce = omm.createPeriodicTorsionForce();
    // periodicTorsionForce->setForceGroup(SLOW_GROUP);
    omm.periodicTorsionForceIndex = omm.system->addForce(periodicTorsionForce);

    // Harmonic improper torsions (CHARMM) are handled using a CustomTorsionForce since OpenMM does not
    // support them natively
    if (!systemTopology.harmonicImproperTorsions.empty()) {
        auto* improperTorsionForce = omm.createImproperHarmonicTorsionForce();
        // improperTorsionForce->setForceGroup(SLOW_GROUP);
        omm.improperHarmonicTorsionForceIndex = omm.system->addForce(improperTorsionForce);
    }

    // Add correction map torsions (CMAPs) force
    if (!systemTopology.cmapGrids.empty()) {
        auto* cmapTorsionForce = omm.createCMAPTorsionForce();
        // cmapTorsionForce->setForceGroup(SLOW_GROUP);
        omm.cmapTorsionForceIndex = omm.system->addForce(cmapTorsionForce);
    }

    // Add Urey-Bradley Potential
    if (!systemTopology.ureyBradleys.empty()) {
        auto* ubForce = omm.createUreyBradleyForce();
        // ubForce->setForceGroup(SLOW_GROUP);
        omm.ureyBradleyForceIndex = omm.system->addForce(ubForce);
    }

    // Create the integrator
    std::vector<std::pair<int, int>> groups = {{SLOW_GROUP, 1}, {FAST_GROUP, 4}};
    omm.integrator = std::make_unique<OpenMM::VerletIntegrator>(0.001);

#if USE_CPU
    OpenMM::Platform* platform = new OpenMM::CpuPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CPU";

#elif USE_REFERENCE
    OpenMM::Platform* platform = new OpenMM::ReferencePlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "Reference";

#elif USE_OPENCL
    OpenMM::Platform* platform = new OpenMM::OpenCLPlatform();
    OpenMM::Platform::registerPlatform(platform);
    platform->setPropertyDefaultValue("Precision", "mixed");
    constexpr auto PLATFORM_NAME = "OpenCL";

#elif USE_CUDA
    OpenMM::Platform* platform = new OpenMM::CudaPlatform();
    OpenMM::Platform::registerPlatform(platform);
    platform->setPropertyDefaultValue("Precision", "mixed");
    constexpr auto PLATFORM_NAME = "CUDA";
#endif

    try {
        omm.context = std::make_unique<OpenMM::Context>(*omm.system, *omm.integrator, *platform);

        const SimTK::Real speed = omm.context->getPlatform().getSpeed();
        std::cout << "Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed "
                  << speed << "\n";

    } catch (const std::exception& e) {
        std::cerr << "[ERROR]: Failed to create OpenMM Context.\n";
        std::cerr << "[ERROR]: " << e.what() << "\n";
        return false;
    }

    // All OpenMM components initialized successfully
    omm.initialized = true;

    // // Set initial positions
    // for (const auto& atom : systemTopology.atoms) {
    //     omm.ommAtomsPositionsCache[atom.identity.globalIndex] =
    //         OpenMM::Vec3(atom.position[0], atom.position[1], atom.position[2]);
    // }
    // omm.context->setPositions(omm.ommAtomsPositionsCache);

    // Log initialization
    std::cout << "[INFO] Initialized OpenMM. Using version " << platform->getOpenMMVersion() << ".\n";

    return true;
}

void OPENMM::setActiveForceGroup(int forceGroupIndex) {
    // ensureInitialized();

    // // Update harmonic bond parameters for the active force group
    // auto& harmonicBondForceRef = system->getForce(harmonicBondForceIndex);
    // auto* harmonicBondForce = dynamic_cast<OpenMM::HarmonicBondForce*>(&harmonicBondForceRef);

    // for (const auto& bond : harmonicBondParamsByWorld[forceGroupIndex]) {
    //     harmonicBondForce->setBondParameters(bond.index, bond.particle1, bond.particle2, bond.length,
    //     bond.k);
    // }
    // harmonicBondForce->updateParametersInContext(*context);

    // // Update harmonic angle parameters for the active force group
    // auto& harmonicAngleForceRef = system->getForce(harmonicAngleForceIndex);
    // auto* harmonicAngleForce = dynamic_cast<OpenMM::HarmonicAngleForce*>(&harmonicAngleForceRef);

    // for (const auto& angle : harmonicAngleParamsByWorld[forceGroupIndex]) {
    //     harmonicAngleForce->setAngleParameters(angle.index,
    //                                            angle.particle1,
    //                                            angle.particle2,
    //                                            angle.particle3,
    //                                            angle.angle,
    //                                            angle.k);
    // }
    // harmonicAngleForce->updateParametersInContext(*context);

    // // Update periodic torsion parameters for the active force group
    // auto& periodicTorsionForceRef = system->getForce(periodicTorsionForceIndex);
    // auto* periodicTorsionForce = dynamic_cast<OpenMM::PeriodicTorsionForce*>(&periodicTorsionForceRef);

    // for (const auto& torsion : periodicTorsionParamsByWorld[forceGroupIndex]) {
    //     periodicTorsionForce->setTorsionParameters(torsion.index,
    //                                                torsion.particle1,
    //                                                torsion.particle2,
    //                                                torsion.particle3,
    //                                                torsion.particle4,
    //                                                torsion.periodicity,
    //                                                torsion.phase,
    //                                                torsion.k);
    // }
    // periodicTorsionForce->updateParametersInContext(*context);

    // // Update nonbonded parameters for the active force group
    // auto& nonbondedForceRef = system->getForce(nonbondedForceIndex);
    // auto* nonbondedForce = dynamic_cast<OpenMM::NonbondedForce*>(&nonbondedForceRef);

    // for (const auto& nb : nonbondedParamsByWorld[forceGroupIndex]) {
    //     nonbondedForce->setParticleParameters(nb.index, nb.charge, nb.sigma, nb.epsilon);
    // }
    // nonbondedForce->updateParametersInContext(*context);
}

void OPENMM::evaluateEnergiesFromPositionCache(SimTK::Real& newPotentialEnergy,
                                               SimTK::Real& newKineticEnergy) {
    ensureInitialized();

    // Intentional value capture: relies on C++17 guaranteed copy elision.
    const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
    newPotentialEnergy = state.getPotentialEnergy();
    newKineticEnergy = state.getKineticEnergy();
}

auto OPENMM::integrateTrajectory(int steps, SimTK::Real timeStepInPicoseconds) -> bool {
    ensureInitialized();

    integrator->setStepSize(timeStepInPicoseconds);

    // Try to integrate
    bool success = true;
    try {
        integrator->step(steps);
    } catch (const std::exception& e) {
        // // Restore old positions in case of integration failure
        // ommAtomsPositionsCache = ommAtomsPositionsCacheOld;
        // context->setPositions(ommAtomsPositionsCache);
        success = false;
    }

    // Intentional value capture: relies on C++17 guaranteed copy elision.
    const auto state =
        context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Velocities,
                          enforcePeriodicBox);
    const auto& positions = state.getPositions();
    const auto& velocities = state.getVelocities();

    potentialEnergy = state.getPotentialEnergy();
    kineticEnergy = state.getKineticEnergy();

    // // Copy positions and velocities
    // simbodyAtomsPositionsCache.resize(positions.size());
    // simbodyAtomsVelocitiesCache.resize(velocities.size());
    // ommAtomsPositionsCache = positions;
    // ommAtomsVelocitiesCache = velocities;

    // for (size_t i = 0; i < positions.size(); ++i) {
    //     const auto& pos = positions[i];
    //     simbodyAtomsPositionsCache[i] = {pos[0], pos[1], pos[2]};
    // }

    // for (size_t i = 0; i < velocities.size(); ++i) {
    //     const auto& vel = velocities[i];
    //     simbodyAtomsVelocitiesCache[i] = {vel[0], vel[1], vel[2]};
    // }

    return success;
}

auto OPENMM::integrateTrajectory(std::vector<OpenMM::Vec3>& positions,
                                 int direction,
                                 std::vector<OpenMM::Vec3>& velocities,
                                 int steps,
                                 SimTK::Real timeStepInPicoseconds,
                                 SimTK::Real& potentialEnergy,
                                 SimTK::Real& kineticEnergy) -> bool {
    ensureInitialized();

    const auto initialState = context->getState(OpenMM::State::Energy);
    const auto initialPE = initialState.getPotentialEnergy();
    const auto initialKE = initialState.getKineticEnergy();
    const auto initialTotalEnergy = initialPE + initialKE;

    // // Set positions
    // context->setPositions(positions);

    // // Negate velocities for backward integration
    // if (direction == -1) {
    //     for (auto& vel : velocities) {
    //         vel = -vel;
    //     }
    // }
    // context->setVelocities(velocities);

    // Integrate
    bool success = true;
    try {
        integrator->setStepSize(timeStepInPicoseconds);
        integrator->step(steps);
    } catch (const std::exception& e) {
        success = false;
    }

    // Restore directions of velocities
    if (direction == -1) {
        // Get new velocities
        const auto state = context->getState(OpenMM::State::Velocities);
        auto newVelocities = state.getVelocities();
        for (auto& vel : newVelocities) {
            vel = -vel;
        }
        context->setVelocities(newVelocities);
    }

    // Get new positions, energies, and velocities
    const auto state =
        context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Velocities,
                          enforcePeriodicBox);
    positions = state.getPositions();
    velocities = state.getVelocities();
    potentialEnergy = state.getPotentialEnergy();
    kineticEnergy = state.getKineticEnergy();

    const auto totalEnergy = potentialEnergy + kineticEnergy;

    return success;
}

void OPENMM::updatePositionsCache(const NonBondedMappings& nonBondedMappings,
                                  const SimTK::Vector_<SimTK::Vec3>& inclAtomPos_G) {
    // SimTK::Vec3 is SimTK::Real[3] and OpenMM::Vec3 is {SimTK::Real x,y,z} — same layout.
    // Verified at compile time; lets us treat both sides as flat SimTK::Real* and
    // skip the OpenMM::Vec3 constructor, exposing plain scalar stores to the vectorizer.
    static_assert(sizeof(SimTK::Vec3) == 3 * sizeof(SimTK::Real), "SimTK::Vec3 layout changed");
    static_assert(sizeof(OpenMM::Vec3) == 3 * sizeof(SimTK::Real), "OpenMM::Vec3 layout changed");

    const int N = static_cast<int>(nonBondedMappings.dummAtomIndex.size());
    const int* __restrict__ dAIxArr = nonBondedMappings.dummAtomIndex.data();
    const int* __restrict__ iaxArr = nonBondedMappings.includedAtomIndex.data();

    // Flat views — no Vec3 abstraction in the hot loop
    const auto* __restrict__ src = reinterpret_cast<const SimTK::Real*>(&inclAtomPos_G[0]);
    auto* __restrict__ dst = reinterpret_cast<SimTK::Real*>(ommAtomsPositionsCache.data());

// Each iteration is fully independent (gather-scatter with disjoint writes).
// schedule(static) avoids dynamic overhead; if() guard skips thread-launch
// cost for tiny systems where serial is faster.
#pragma omp parallel for schedule(static) if (N > 512)
    for (int i = 0; i < N; ++i) {
        const int srcIx = iaxArr[i] * 3;
        const int dstIx = dAIxArr[i] * 3;
        dst[dstIx] = src[srcIx];
        dst[dstIx + 1] = src[srcIx + 1];
        dst[dstIx + 2] = src[srcIx + 2];
    }

    context->setPositions(ommAtomsPositionsCache);
}

void OPENMM::evaluateForcesFromPositionsCache(const NonBondedMappings& nonBondedMappings,
                                              const SimTK::Vector_<SimTK::Vec3>& inclAtomStation_G,
                                              SimTK::Vector_<SimTK::SpatialVec>& inclBodyForces_G) const {
    ensureInitialized();

    const auto state = context->getState(OpenMM::State::Forces, enforcePeriodicBox);
    const auto& forces = state.getForces();

    static_assert(sizeof(OpenMM::Vec3) == 3 * sizeof(SimTK::Real));
    static_assert(sizeof(SimTK::Vec3) == 3 * sizeof(SimTK::Real));
    static_assert(sizeof(SimTK::SpatialVec) == 6 * sizeof(SimTK::Real));

    const int N = static_cast<int>(nonBondedMappings.dummAtomIndex.size());

    // SoA pointers — no struct stride in the hot loop
    const int* __restrict__ dAIxArr = nonBondedMappings.dummAtomIndex.data();
    const int* __restrict__ iaxArr = nonBondedMappings.includedAtomIndex.data();
    const int* __restrict__ ibxArr = nonBondedMappings.bodyIndex.data();

    // Flat SimTK::Real views — kills Vec3 constructor overhead, lets the compiler
    // see plain loads and cross-product as scalar FMAs
    const auto* __restrict__ frc = reinterpret_cast<const SimTK::Real*>(forces.data());
    const auto* __restrict__ sta = reinterpret_cast<const SimTK::Real*>(&inclAtomStation_G[0]);

    // inclBodyForces_G is small (numBodies << N) and stays hot in L1/L2.
    // Scatter writes to it are cache hits — no atomics or body-sort needed.
    // SpatialVec layout: [0]=torque (3d), [1]=force (3d)
    auto* __restrict__ bfrc = reinterpret_cast<SimTK::Real*>(&inclBodyForces_G[0]);

    // Prefetch distance: tune to hide load latency (cache line = 64 bytes = ~2.6 Vec3s)
    constexpr int PF = 8;

    for (int i = 0; i < N; ++i) {
        // Software prefetch — hides the gather latency for random reads into
        // forces and stations since dAIx/iax can jump anywhere in a large array
        if (i + PF < N) {
            __builtin_prefetch(frc + (dAIxArr[i + PF] * 3), 0, 1);
            __builtin_prefetch(sta + (iaxArr[i + PF] * 3), 0, 1);
        }

        const int fi = dAIxArr[i] * 3;
        const int si = iaxArr[i] * 3;
        const int bi = ibxArr[i] * 6; // SpatialVec = 6 doubles

        const SimTK::Real f0 = frc[fi];
        const SimTK::Real f1 = frc[fi + 1];
        const SimTK::Real f2 = frc[fi + 2];
        const SimTK::Real s0 = sta[si];
        const SimTK::Real s1 = sta[si + 1];
        const SimTK::Real s2 = sta[si + 2];

        // Cross product s % f, written as FMAs explicitly
        // Without this, compiler may emit separate mul+sub+add without -ffast-math
        bfrc[bi] += std::fma(s1, f2, -(s2 * f1));     // torque x
        bfrc[bi + 1] += std::fma(s2, f0, -(s0 * f2)); // torque y
        bfrc[bi + 2] += std::fma(s0, f1, -(s1 * f0)); // torque z
        bfrc[bi + 3] += f0;                           // force x
        bfrc[bi + 4] += f1;                           // force y
        bfrc[bi + 5] += f2;                           // force z
    }
}

auto OPENMM::evaluatePotentialEnergyFromPositionsCache() const -> SimTK::Real {
    ensureInitialized();
    const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
    return state.getPotentialEnergy();
}

auto OPENMM::computePeriodicBoxVectors_Context(SimTK::Real a_length,
                                               SimTK::Real b_length,
                                               SimTK::Real c_length,
                                               SimTK::Real alpha,
                                               SimTK::Real beta,
                                               SimTK::Real gamma)
    -> std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> {
    ensureInitialized();

    const SimTK::Real TOL = 1e-6;

    // // Convert angles from degrees to radians
    // alpha = SimTK::Deg2Rad * alpha;
    // beta  = SimTK::Deg2Rad * beta;
    // gamma = SimTK::Deg2Rad * gamma;

    // Compute the box vectors
    OpenMM::Vec3 a(a_length, 0.0, 0.0);

    OpenMM::Vec3 b(b_length * std::cos(gamma), b_length * std::sin(gamma), 0.0);

    SimTK::Real cx = c_length * std::cos(beta);
    SimTK::Real cy = c_length * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    SimTK::Real cz = std::sqrt(c_length * c_length - cx * cx - cy * cy);

    OpenMM::Vec3 c(cx, cy, cz);

    // Zero out small components
    for (int i = 0; i < 3; i++) {
        if (std::abs(a[i]) < TOL) {
            a[i] = 0.0;
        }
        if (std::abs(b[i]) < TOL) {
            b[i] = 0.0;
        }
        if (std::abs(c[i]) < TOL) {
            c[i] = 0.0;
        }
    }

    // Reduced form (OpenMM requirement)
    if (b[1] != 0.0) {
        c -= b * std::round(c[1] / b[1]);
    }
    if (a[0] != 0.0) {
        c -= a * std::round(c[0] / a[0]);
    }
    if (a[0] != 0.0) {
        b -= a * std::round(b[0] / a[0]);
    }

    return std::make_tuple(a, b, c);
}

auto OPENMM::computeIntraRigidPairs(int wIx, int minBondedDistance) -> std::vector<std::pair<int, int>> {
    std::vector<std::vector<int>> adj(numAtoms);
    for (const auto& bond : systemTopology.bonds) {
        int src = bond.globalIndices[0];
        int dst = bond.globalIndices[1];
        adj[src].push_back(dst);
        adj[dst].push_back(src);
    }

    // Get rigid bodies for this world
    const auto& rigidBodies = atomMbxByWorld[wIx];

    std::unordered_map<int, std::vector<int>> atomsByRigidBody;
    for (int idx = 0; idx < static_cast<int>(rigidBodies.size()); ++idx) {
        atomsByRigidBody[rigidBodies[idx]].push_back(idx);
    }

    // BFS must stop BEFORE enqueuing nodes at minBondedDistance
    // e.g. minBondedDistance=4: bonded = {src, 1-2, 1-3, 1-4}, returns 1-5+
    const int stopDepth = minBondedDistance - 1;

    std::vector<std::pair<int, int>> result;

    for (const auto& [bodyId, atomList] : atomsByRigidBody) {
        for (int src : atomList) {
            std::unordered_set<int> bonded;
            std::queue<std::pair<int, int>> bfsQueue;
            bonded.insert(src);
            bfsQueue.emplace(src, 0);

            while (!bfsQueue.empty()) {
                auto [cur, depth] = bfsQueue.front();
                bfsQueue.pop();

                if (depth == stopDepth) {
                    continue; // don't enqueue next level
                }

                for (int nbr : adj[cur]) {
                    if (!isSameRigidBody(wIx, src, nbr)) {
                        continue;
                    }
                    if (bonded.insert(nbr).second) {
                        bfsQueue.emplace(nbr, depth + 1);
                    }
                }
            }

            for (int dst : atomList) {
                if (dst <= src) {
                    continue;
                }
                if (bonded.find(dst) == bonded.end()) {
                    result.emplace_back(src, dst);
                }
            }
        }
    }

    return result;
}

[[nodiscard]] auto OPENMM::createNonbondedForce(const std::vector<RoboAtom>& atoms,
                                                const std::vector<Scaling14>& scaling14s,
                                                const std::vector<Exclusion>& exclusions,
                                                const ForceFieldParams& ffParams,
                                                bool hasNBfix) -> OpenMM::NonbondedForce* {
    auto* nonbondedForce = new OpenMM::NonbondedForce();
    nonbondedForce->setCutoffDistance(ffParams.nonbondedCutoffInNm);

    switch (ffParams.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
            nonbondedForce->setUseDispersionCorrection(false);
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
            if (ffParams.useGBSAOBC2) {
                nonbondedForce->setReactionFieldDielectric(1.0);
                nonbondedForce->setUseDispersionCorrection(false);
            } else {
                nonbondedForce->setUseDispersionCorrection(true);
            }
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method");
    }

    nonbondedParamsByWorld.resize(numWorlds);

    // Add all particles
    for (const auto& atom : atoms) {
        int index = 0;
        if (hasNBfix) {
            index = nonbondedForce->addParticle(atom.physics.chargeInE, 1.0, 0.0);
        } else {
            index = nonbondedForce->addParticle(atom.physics.chargeInE,
                                                atom.physics.sigmaInNm,
                                                atom.physics.vdwWellDepthInKJ);
        }

        for (int wIx = 0; wIx < numWorlds; ++wIx) {
            nonbondedParamsByWorld[wIx].push_back({atom.identity.globalIndex,
                                                   atom.physics.chargeInE,
                                                   atom.physics.sigmaInNm,
                                                   atom.physics.vdwWellDepthInKJ});
        }
    }

    // Add 1-4 scalings
    for (const auto& scaling14 : scaling14s) {
        nonbondedForce->addException(scaling14.atom1GlobalIndex,
                                     scaling14.atom4GlobalIndex,
                                     scaling14.chargeProduct,
                                     scaling14.sigma,
                                     scaling14.epsilon);
    }

    // Add 1-2 and 1-3 exclusions
    for (const auto& exclusion : exclusions) {
        // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely
        // omitted from force and energy calculations
        nonbondedForce->addException(exclusion.atom1GlobalIndex, exclusion.atom2GlobalIndex, 0.0, 0.1, 0.0);
    }

    for (int wIx = 0; wIx < numWorlds; ++wIx) {
        const auto intraRigidPairs = computeIntraRigidPairs(wIx, 4);
        for (const auto& [particle1, particle2] : intraRigidPairs) {
            nonbondedParamsByWorld[wIx][particle1].charge = 0.0;
            nonbondedParamsByWorld[wIx][particle1].sigma = 1.0;
            nonbondedParamsByWorld[wIx][particle1].epsilon = 0.0;

            nonbondedParamsByWorld[wIx][particle2].charge = 0.0;
            nonbondedParamsByWorld[wIx][particle2].sigma = 1.0;
            nonbondedParamsByWorld[wIx][particle2].epsilon = 0.0;
        }
    }

    return nonbondedForce;
}

[[nodiscard]] auto OPENMM::createGBSAOBCForce(const std::vector<RoboAtom>& atoms,
                                              const ForceFieldParams& ffParams) -> OpenMM::GBSAOBCForce* {
    OpenMM::GBSAOBCForce::NonbondedMethod gbsaForceMethod;
    switch (ffParams.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            gbsaForceMethod = OpenMM::GBSAOBCForce::NoCutoff;
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            gbsaForceMethod = OpenMM::GBSAOBCForce::CutoffNonPeriodic;
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method for GBSAOBCForce");
    }

    auto* force = new OpenMM::GBSAOBCForce();

    force->setSolventDielectric(ffParams.gbsaSolventDielectric);
    force->setSoluteDielectric(ffParams.gbsaSoluteDielectric);
    force->setNonbondedMethod(gbsaForceMethod);
    force->setCutoffDistance(ffParams.nonbondedCutoffInNm);

    for (const auto& atom : atoms) {
        force->addParticle(atom.physics.chargeInE, atom.physics.solventRadiusInNm, atom.physics.screen);
    }

    return force;
}

[[nodiscard]] auto OPENMM::createCustomNonbondedForce(const std::vector<RoboAtom>& atoms,
                                                      const std::vector<Scaling14>& scaling14s,
                                                      const std::vector<Exclusion>& exclusions,
                                                      const ForceFieldParams& ffParams,
                                                      bool useSwitchingFunction,
                                                      SimTK::Real switchingDistance)
    -> OpenMM::CustomNonbondedForce* {
    auto* force = new OpenMM::CustomNonbondedForce(
        "(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2);");
    force->addTabulatedFunction(
        "acoef",
        new OpenMM::Discrete2DFunction(ffParams.numTypes, ffParams.numTypes, ffParams.aCoef));
    force->addTabulatedFunction(
        "bcoef",
        new OpenMM::Discrete2DFunction(ffParams.numTypes, ffParams.numTypes, ffParams.bCoef));
    force->addPerParticleParameter("type");

    force->setCutoffDistance(ffParams.nonbondedCutoffInNm);
    force->setUseSwitchingFunction(useSwitchingFunction);
    force->setSwitchingDistance(switchingDistance);

    switch (ffParams.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            force->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            force->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method");
    }

    for (const auto& atom : atoms) {
        force->addParticle({static_cast<SimTK::Real>(atom.identity.nonbondedIndex)});
    }

    for (const auto& scaling14 : scaling14s) {
        force->addExclusion(scaling14.atom1GlobalIndex, scaling14.atom4GlobalIndex);
    }

    // Exclude 1-2 and 1-3 interactions
    for (const auto& exclusion : exclusions) {
        // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely
        // omitted from force and energy calculations
        force->addExclusion(exclusion.atom1GlobalIndex, exclusion.atom2GlobalIndex);
    }

    return force;
}

[[nodiscard]] auto OPENMM::createHarmonicBondForce() -> OpenMM::HarmonicBondForce* {
    auto* force = new OpenMM::HarmonicBondForce();
    harmonicBondParamsByWorld.resize(numWorlds);

    for (const auto& bond : systemTopology.bonds) {
        const auto particle1 = bond.globalIndices[0];
        const auto particle2 = bond.globalIndices[1];
        const auto length = bond.nominalLengthInNm;

        // OpenMM defines the harmonic bond potential as 0.5 * k * (r-r0)^2, so we need to multiply by 2 to
        // get the correct stiffness
        auto stiffness = bond.stiffnessInKJPerNmSq * 2;

        // Add bond and save the interaction index
        const auto index = force->addBond(particle1, particle2, length, stiffness);

        for (int wIx = 0; wIx < numWorlds; ++wIx) {
            if (isSameRigidBody(wIx, particle1, particle2)) {
                stiffness = 0.0;
            }

            harmonicBondParamsByWorld[wIx].push_back({index, particle1, particle2, length, stiffness});
        }
    }

    return force;
}

[[nodiscard]] auto OPENMM::createHarmonicAngleForce() -> OpenMM::HarmonicAngleForce* {
    auto* force = new OpenMM::HarmonicAngleForce();
    harmonicAngleParamsByWorld.resize(numWorlds);

    for (const auto& angle : systemTopology.angles) {
        const auto particle1 = angle.globalIndices[0];
        const auto particle2 = angle.globalIndices[1];
        const auto particle3 = angle.globalIndices[2];
        const auto angleInRad = angle.nominalAngleInDeg * SimTK::DuMM::Deg2Rad;

        // OpenMM defines the harmonic angle potential as 0.5 * k * (theta-theta0)^2, so we need to multiply
        // by 2 to get the correct stiffness
        auto stiffness = angle.stiffnessInKJPerRadSq * 2;

        const auto index = force->addAngle(particle1, particle2, particle3, angleInRad, stiffness);

        for (int wIx = 0; wIx < numWorlds; ++wIx) {
            if (isSameRigidBody(wIx, particle1, particle2) && isSameRigidBody(wIx, particle2, particle3)) {
                stiffness = 0.0;
            }

            harmonicAngleParamsByWorld[wIx].push_back(
                {index, particle1, particle2, particle3, angleInRad, stiffness});
        }
    }
    return force;
}

[[nodiscard]] auto OPENMM::createPeriodicTorsionForce() -> OpenMM::PeriodicTorsionForce* {
    auto* force = new OpenMM::PeriodicTorsionForce();
    periodicTorsionParamsByWorld.resize(numWorlds);

    for (const auto& torsion : systemTopology.periodicTorsions) {
        const auto particle1 = torsion.globalIndices[0];
        const auto particle2 = torsion.globalIndices[1];
        const auto particle3 = torsion.globalIndices[2];
        const auto particle4 = torsion.globalIndices[3];

        for (const auto& term : torsion.terms) {
            const auto periodicity = term.periodicity;
            if (periodicity <= 0) {
                continue;
            }
            const auto phase = term.phaseDeg * SimTK::DuMM::Deg2Rad;
            auto stiffness = term.amplitudeKJ;

            const auto index =
                force->addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, stiffness);

            for (int wIx = 0; wIx < numWorlds; ++wIx) {
                if (!torsion.improper) {
                    // Skip propers where the middle two atoms are in the same rigid body
                    if (isSameRigidBody(wIx, particle2, particle3)) {
                        stiffness = 0.0;
                    }
                } else {
                    // Skip improper torsions where all four atoms are in the same rigid body
                    if (isSameRigidBody(wIx, particle1, particle2)
                        && isSameRigidBody(wIx, particle2, particle3)
                        && isSameRigidBody(wIx, particle3, particle4)) {
                        stiffness = 0.0;
                    }
                }

                periodicTorsionParamsByWorld[wIx].push_back(
                    {index, particle1, particle2, particle3, particle4, periodicity, phase, stiffness});
            }
        }
    }

    return force;
}

[[nodiscard]] auto OPENMM::createImproperHarmonicTorsionForce() -> OpenMM::CustomTorsionForce* {
    // Create force expression for harmonic improper torsions (CHARMM style)
    std::stringstream strStream;
    strStream << std::setprecision(17) << "k*min(dtheta, 2*" << SimTK::Pi
              << "-dtheta)^2; dtheta = abs(theta-theta0)";

    auto* force = new OpenMM::CustomTorsionForce(strStream.str());
    force->addPerTorsionParameter("k");
    force->addPerTorsionParameter("theta0");

    for (const auto& torsion : systemTopology.harmonicImproperTorsions) {
        const auto atom1GlobalIndex = torsion.globalIndices[0];
        const auto atom2GlobalIndex = torsion.globalIndices[1];
        const auto atom3GlobalIndex = torsion.globalIndices[2];
        const auto atom4GlobalIndex = torsion.globalIndices[3];

        const std::vector<SimTK::Real> params = {torsion.stiffnessInKJPerRadSq, torsion.nominalAngleInRad};
        force->addTorsion(atom1GlobalIndex, atom2GlobalIndex, atom3GlobalIndex, atom4GlobalIndex, params);
    }

    return force;
}

[[nodiscard]] auto OPENMM::createCMAPTorsionForce() -> OpenMM::CMAPTorsionForce* {
    auto* force = new OpenMM::CMAPTorsionForce();
    force->setUsesPeriodicBoundaryConditions(enforcePeriodicBox);

    for (const auto& grid : systemTopology.cmapGrids) {
        force->addMap(grid.size, grid.energy);
    }

    for (const auto& cmapTorsion : systemTopology.cmapTorsions) {
        force->addTorsion(cmapTorsion.mapIndex,
                          cmapTorsion.torsionAAtom1GlobalIndex,
                          cmapTorsion.torsionAAtom2GlobalIndex,
                          cmapTorsion.torsionAAtom3GlobalIndex,
                          cmapTorsion.torsionAAtom4GlobalIndex,
                          cmapTorsion.torsionBAtom1GlobalIndex,
                          cmapTorsion.torsionBAtom2GlobalIndex,
                          cmapTorsion.torsionBAtom3GlobalIndex,
                          cmapTorsion.torsionBAtom4GlobalIndex);
    }

    return force;
}

[[nodiscard]] auto OPENMM::createUreyBradleyForce() -> OpenMM::HarmonicBondForce* {
    auto* force = new OpenMM::HarmonicBondForce();

    for (const auto& term : systemTopology.ureyBradleys) {
        force->addBond(term.atom1GlobalIndex,
                       term.atom3GlobalIndex,
                       term.nominalLengthInNm,
                       term.stiffnessInKJPerNmSq * 2);
    }

    return force;
}
