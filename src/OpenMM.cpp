#include "OpenMM.hpp"

#include <chrono>
#include <vector>

#if USE_CPU
#    include "../Molmodel/src/gbsa/cpuObcInterface.h"
#    include "../openmm/platforms/cpu/include/CpuPlatform.h"
#elif USE_REFERENCE
#    include "../openmm/platforms/reference/include/ReferencePlatform.h"
#elif USE_OPENCL
#    include "../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif USE_CUDA
#    include "../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

void ForceGroup::initialize(int forceGroupIndex,
                            const std::vector<int>& rigidBodies,
                            const SystemTopology& systemTopology,
                            const ForceFieldParams& ffParams,
                            const SimulationSettings& simSettings) {
    std::cout << "Initializing OpenMM force group " << forceGroupIndex << ":\n";

    OpenMM::NonbondedForce::NonbondedMethod nonbondedForceMethod;
    OpenMM::CustomNonbondedForce::NonbondedMethod customNonbondedForceMethod;
    OpenMM::GBSAOBCForce::NonbondedMethod gbsaForceMethod;

    switch (ffParams.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            nonbondedForceMethod = OpenMM::NonbondedForce::NoCutoff;
            customNonbondedForceMethod = OpenMM::CustomNonbondedForce::NoCutoff;
            gbsaForceMethod = OpenMM::GBSAOBCForce::NoCutoff;
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            nonbondedForceMethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
            customNonbondedForceMethod = OpenMM::CustomNonbondedForce::CutoffNonPeriodic;
            gbsaForceMethod = OpenMM::GBSAOBCForce::CutoffNonPeriodic;
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method");
    }

    // Save force group index
    fg = forceGroupIndex;

    // Instantiate the thermostat with adjusted temperature
    auto* thermostat =
        new OpenMM::AndersenThermostat(simSettings.thermostatTemperatureInK, simSettings.collisionFrequency);
    thermostat->setRandomNumberSeed(simSettings.seed);
    registerForce(thermostat);

    // Nonbonded forces
    auto* nonbondedForce = new OpenMM::NonbondedForce();
    nonbondedForce->setNonbondedMethod(nonbondedForceMethod);
    nonbondedForce->setCutoffDistance(ffParams.nonbondedCutoffInNm);

    if (nonbondedForceMethod == OpenMM::NonbondedForce::CutoffNonPeriodic) {
        nonbondedForce->setUseDispersionCorrection(true);
    } else {
        nonbondedForce->setUseDispersionCorrection(false);
    }

    // 1-4 VdW Correction (Bond-based so we can target specific pairs)
    OpenMM::CustomNonbondedForce* customNonbonded = nullptr;
    if (ffParams.hasNBfix) {
        customNonbonded = new OpenMM::CustomNonbondedForce(
            "(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2);");
        customNonbonded->addTabulatedFunction(
            "acoef",
            new OpenMM::Discrete2DFunction(ffParams.numTypes, ffParams.numTypes, ffParams.aCoef));
        customNonbonded->addTabulatedFunction(
            "bcoef",
            new OpenMM::Discrete2DFunction(ffParams.numTypes, ffParams.numTypes, ffParams.bCoef));
        customNonbonded->addPerParticleParameter("type");

        customNonbonded->setNonbondedMethod(customNonbondedForceMethod);
        customNonbonded->setCutoffDistance(ffParams.nonbondedCutoffInNm);

        customNonbonded->setUseSwitchingFunction(nonbondedForce->getUseSwitchingFunction());
        customNonbonded->setSwitchingDistance(nonbondedForce->getSwitchingDistance());
    }

    // GBSA OBC Force
    // implicitSolventKappa
    OpenMM::GBSAOBCForce* GBSAOBCForce = nullptr;
    if (ffParams.useGBSAOBC2) {
        GBSAOBCForce = new OpenMM::GBSAOBCForce();
        GBSAOBCForce->setSolventDielectric(ffParams.gbsaSolventDielectric);
        GBSAOBCForce->setSoluteDielectric(ffParams.gbsaSoluteDielectric);
        GBSAOBCForce->setNonbondedMethod(gbsaForceMethod);
        GBSAOBCForce->setCutoffDistance(ffParams.nonbondedCutoffInNm);

        /*
         * The following comment is copy-pasted from the OpenMM documentation for GBSAOBCForce (note that we
         * only support CutoffNonPeriodic): When using GBSAOBCForce, the System should also include a
         * NonbondedForce, and both objects must specify identical charges for all particles. Otherwise, the
         * results will not be correct. Furthermore, if the nonbonded method is set to CutoffNonPeriodic or
         * CutoffPeriodic, you should call setReactionFieldDielectric(1.0) on the NonbondedForce to turn off
         * the reaction field approximation, which does not produce correct results when combined with GBSA.
         */
        if (nonbondedForceMethod == OpenMM::NonbondedForce::CutoffNonPeriodic) {
            nonbondedForce->setReactionFieldDielectric(1.0);
            nonbondedForce->setUseDispersionCorrection(false);
        }
    }

    // Compute rigid body sizes - only multi-atom bodies cause terms to be skipped
    std::unordered_map<int, int> rigidBodySizes;
    for (int rb : rigidBodies) {
        rigidBodySizes[rb]++;
    }

    auto inSameMultiAtomRigidBody = [&](int a, int b) {
        return rigidBodies[a] == rigidBodies[b] && rigidBodySizes.at(rigidBodies[a]) > 1;
    };

    // Add atoms
    for (const auto& atom : systemTopology.atoms) {
        const SimTK::Real charge = atom.physics.chargeInE;
        const SimTK::Real sigma = atom.physics.sigmaInNm;
        const SimTK::Real epsilon = atom.physics.vdwWellDepthInKJ;

        if (!std::isfinite(charge)) {
            throw std::runtime_error("Non-finite charge for atom " + atom.identity.uniqueAtomName + ": "
                                     + std::to_string(charge));
        }

        if (!std::isfinite(sigma)) {
            throw std::runtime_error("Non-finite sigma for atom " + atom.identity.uniqueAtomName + ": "
                                     + std::to_string(sigma));
        }

        if (!std::isfinite(epsilon)) {
            throw std::runtime_error("Non-finite epsilon for atom " + atom.identity.uniqueAtomName + ": "
                                     + std::to_string(epsilon));
        }

        if (ffParams.hasNBfix) {
            nonbondedForce->addParticle(charge, 1.0, 0.0);
            customNonbonded->addParticle({static_cast<SimTK::Real>(atom.identity.nonbondedIndex)});
        } else {
            nonbondedForce->addParticle(charge, sigma, epsilon);
        }

        if (ffParams.useGBSAOBC2) {
            const SimTK::Real solventRadiusInNm = atom.physics.solventRadiusInNm;
            const bool validSolventRadius = std::isfinite(solventRadiusInNm) && solventRadiusInNm > 0.0;
            if (!validSolventRadius) {
                throw std::runtime_error("Invalid solvent radius for atom " + atom.identity.uniqueAtomName
                                         + ": " + std::to_string(solventRadiusInNm));
            }

            const SimTK::Real screen = atom.physics.screen;
            const bool validScreen = std::isfinite(screen) && screen >= 0.0 && screen <= 1.0;
            if (!validScreen) {
                throw std::runtime_error("Invalid screen value for atom " + atom.identity.uniqueAtomName
                                         + ": " + std::to_string(screen));
            }

            GBSAOBCForce->addParticle(charge, solventRadiusInNm, screen);
        }
    }

    if (GBSAOBCForce != nullptr) {
        registerForce(GBSAOBCForce);
    }

    int numAddedExceptions = 0;

    for (const auto& scaling14 : systemTopology.scaling14s) {
        if (!std::isfinite(scaling14.chargeProduct)) {
            throw std::runtime_error("Scaling 1-4 between atoms " + std::to_string(scaling14.atom1GlobalIndex)
                                     + " and " + std::to_string(scaling14.atom4GlobalIndex)
                                     + " has non-finite charge product: "
                                     + std::to_string(scaling14.chargeProduct));
        }

        if (!std::isfinite(scaling14.sigma)) {
            throw std::runtime_error("Scaling 1-4 between atoms " + std::to_string(scaling14.atom1GlobalIndex)
                                     + " and " + std::to_string(scaling14.atom4GlobalIndex)
                                     + " has non-finite sigma: " + std::to_string(scaling14.sigma));
        }

        if (!std::isfinite(scaling14.epsilon)) {
            throw std::runtime_error("Scaling 1-4 between atoms " + std::to_string(scaling14.atom1GlobalIndex)
                                     + " and " + std::to_string(scaling14.atom4GlobalIndex)
                                     + " has non-finite epsilon: " + std::to_string(scaling14.epsilon));
        }

        nonbondedForce->addException(scaling14.atom1GlobalIndex,
                                     scaling14.atom4GlobalIndex,
                                     scaling14.chargeProduct,
                                     scaling14.sigma,
                                     scaling14.epsilon);

        if (ffParams.hasNBfix) {
            customNonbonded->addExclusion(scaling14.atom1GlobalIndex, scaling14.atom4GlobalIndex);
        }
    }

    // Exclude 1-2 and 1-3 interactions
    for (const auto& exclusion : systemTopology.exclusions) {
        // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely
        // omitted from force and energy calculations
        nonbondedForce->addException(exclusion.atom1GlobalIndex, exclusion.atom2GlobalIndex, 0.0, 0.1, 0.0);
        ++numAddedExceptions;

        if (ffParams.hasNBfix) {
            customNonbonded->addExclusion(exclusion.atom1GlobalIndex, exclusion.atom2GlobalIndex);
        }
    }

    std::vector<std::pair<int, int>> intraRigidPairs;

    // Build adjacency list
    std::vector<std::vector<int>> adj(systemTopology.atoms.size());
    for (const auto& b : systemTopology.bonds) {
        int a = b.globalIndices[0];
        int c = b.globalIndices[1];
        adj[a].push_back(c);
        adj[c].push_back(a);
    }

    // Group atoms by rigid body
    std::unordered_map<int, std::vector<int>> rbAtoms;
    for (int i = 0; i < rigidBodies.size(); ++i) {
        rbAtoms[rigidBodies[i]].push_back(i);
    }

    for (const auto& [rb, atomList] : rbAtoms) {
        for (int a : atomList) {
            std::queue<std::pair<int, int>> q;
            std::unordered_set<int> visited;

            q.push({a, 0});
            visited.insert(a);

            while (!q.empty()) {
                auto [v, depth] = q.front();
                q.pop();

                if (depth == 4) {
                    continue;
                }

                for (int nb : adj[v]) {
                    if (rigidBodies[nb] != rb) {
                        continue;
                    }
                    if (visited.insert(nb).second) {
                        q.push({nb, depth + 1});
                    }
                }
            }

            for (int b : atomList) {
                if (b <= a) {
                    continue;
                }
                if (!visited.count(b)) {
                    intraRigidPairs.emplace_back(a, b);
                }
            }
        }
    }

    for (const auto& pair : intraRigidPairs) {
        // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be
        // completelyomitted from force and energy calculations
        nonbondedForce->addException(pair.first, pair.second, 0.0, 0.1, 0.0);
        ++numAddedExceptions;

        if (ffParams.hasNBfix) {
            customNonbonded->addExclusion(pair.first, pair.second);
        }
    }

    std::cout << "\tAdded " << numAddedExceptions << " nonbonded exceptions.\n";

    // Register non-bonded forces
    registerForce(nonbondedForce);
    if (ffParams.hasNBfix) {
        registerForce(customNonbonded);
    }

    // Add bonds
    int numAddedBonds = 0;
    auto* harmonicBondStretch = new OpenMM::HarmonicBondForce();

    for (const auto& bond : systemTopology.bonds) {
        const auto& g = bond.globalIndices;

        // Don't add bonds between atoms in the same rigid body
        if (inSameMultiAtomRigidBody(g[0], g[1])) {
            continue;
        }

        harmonicBondStretch->addBond(g[0], g[1], bond.nominalLengthInNm, bond.stiffnessInKJPerNmSq * 2);
        ++numAddedBonds;
    }

    registerForce(harmonicBondStretch);
    std::cout << "\tAdded " << numAddedBonds << " bonds.\n";

    // Add angles
    int numAddedAngles = 0;
    auto* harmonicAngleForce = new OpenMM::HarmonicAngleForce();

    for (const auto& angle : systemTopology.angles) {
        const auto& g = angle.globalIndices;

        // Skip angles where all three atoms are in the same rigid body
        if (inSameMultiAtomRigidBody(g[0], g[1]) && inSameMultiAtomRigidBody(g[1], g[2])) {
            continue;
        }

        harmonicAngleForce->addAngle(g[0],
                                     g[1],
                                     g[2],
                                     angle.nominalAngleInDeg * SimTK::DuMM::Deg2Rad,
                                     angle.stiffnessInKJPerRadSq * 2);
        ++numAddedAngles;
    }

    registerForce(harmonicAngleForce);
    std::cout << "\tAdded " << numAddedAngles << " angles.\n";

    // Define torsions
    // Note that OpenMM only supports only periodic torsions by default (AMBER style), so we need to handle
    // improper harmonic torsions differently For AMBER style torsions, OpenMM does not care if they are
    // proper or improper
    int numAddedPeriodicTorsions = 0;
    auto* periodicTorsionForce = new OpenMM::PeriodicTorsionForce();

    for (const auto& t : systemTopology.periodicTorsions) {
        const auto& g = t.globalIndices;
        for (int i = 0; i < t.numTerms; ++i) {
            const auto& term = t.terms[i];

            // Skip propers where the middle two atoms are in the same rigid body
            if (!t.improper && inSameMultiAtomRigidBody(g[1], g[2])) {
                continue;
            }

            // Skip improper torsions where all four atoms are in the same rigid body
            if (t.improper && inSameMultiAtomRigidBody(g[0], g[1]) && inSameMultiAtomRigidBody(g[1], g[2])
                && inSameMultiAtomRigidBody(g[2], g[3])) {
                continue;
            }

            periodicTorsionForce->addTorsion(g[0],
                                             g[1],
                                             g[2],
                                             g[3],
                                             term.periodicity,
                                             term.phaseDeg * SimTK::DuMM::Deg2Rad,
                                             term.amplitudeKJ);
            ++numAddedPeriodicTorsions;
        }
    }
    std::cout << "\tAdded " << numAddedPeriodicTorsions << " periodic torsions.\n";
    registerForce(periodicTorsionForce);

    // Harmonic improper torsions (CHARMM) are handled using a CustomTorsionForce since OpenMM does not
    // support them natively
    if (!systemTopology.harmonicImproperTorsions.empty()) {
        std::stringstream ss;
        ss << std::setprecision(17) << "k*min(dtheta, 2*" << SimTK::Pi
           << "-dtheta)^2; dtheta = abs(theta-theta0)";

        auto* improperTorsionForce = new OpenMM::CustomTorsionForce(ss.str());
        improperTorsionForce->addPerTorsionParameter("k");
        improperTorsionForce->addPerTorsionParameter("theta0");

        int numAddedImproperTorsions = 0;

        for (const auto& t : systemTopology.harmonicImproperTorsions) {
            const int a1 = t.globalIndices[0];
            const int a2 = t.globalIndices[1];
            const int a3 = t.globalIndices[2];
            const int a4 = t.globalIndices[3];

            // Improper torsions where all four atoms are in the same rigid body won't change, so we can skip
            // them
            if (inSameMultiAtomRigidBody(a1, a2) && inSameMultiAtomRigidBody(a2, a3)
                && inSameMultiAtomRigidBody(a3, a4)) {
                continue;
            }

            std::vector<double> params = {t.stiffnessInKJPerRadSq, t.nominalAngleInRad};
            improperTorsionForce->addTorsion(a1, a2, a3, a4, params);
            ++numAddedImproperTorsions;
        }
        registerForce(improperTorsionForce);
        std::cout << "\tAdded " << numAddedImproperTorsions << " harmonic improper torsions.\n";
    }

    // Add correction map torsions (CMAPs) force
    if (!systemTopology.cmapGrids.empty()) {
        auto* cmapTorsionForce = new OpenMM::CMAPTorsionForce();
        cmapTorsionForce->setUsesPeriodicBoundaryConditions(false);

        for (const auto& grid : systemTopology.cmapGrids) {
            cmapTorsionForce->addMap(grid.size, grid.energy);
        }

        int numAddedCmapTorsions = 0;

        for (const auto& cmapTorsion : systemTopology.cmapTorsions) {
            const bool phi_fixed = inSameMultiAtomRigidBody(cmapTorsion.torsionAAtom1GlobalIndex,
                                                            cmapTorsion.torsionAAtom2GlobalIndex)
                                   && inSameMultiAtomRigidBody(cmapTorsion.torsionAAtom2GlobalIndex,
                                                               cmapTorsion.torsionAAtom3GlobalIndex)
                                   && inSameMultiAtomRigidBody(cmapTorsion.torsionAAtom3GlobalIndex,
                                                               cmapTorsion.torsionAAtom4GlobalIndex);
            const bool psi_fixed = inSameMultiAtomRigidBody(cmapTorsion.torsionBAtom1GlobalIndex,
                                                            cmapTorsion.torsionBAtom2GlobalIndex)
                                   && inSameMultiAtomRigidBody(cmapTorsion.torsionBAtom2GlobalIndex,
                                                               cmapTorsion.torsionBAtom3GlobalIndex)
                                   && inSameMultiAtomRigidBody(cmapTorsion.torsionBAtom3GlobalIndex,
                                                               cmapTorsion.torsionBAtom4GlobalIndex);
            if (phi_fixed && psi_fixed) {
                continue;
            }

            cmapTorsionForce->addTorsion(cmapTorsion.mapIndex,
                                         cmapTorsion.torsionAAtom1GlobalIndex,
                                         cmapTorsion.torsionAAtom2GlobalIndex,
                                         cmapTorsion.torsionAAtom3GlobalIndex,
                                         cmapTorsion.torsionAAtom4GlobalIndex,
                                         cmapTorsion.torsionBAtom1GlobalIndex,
                                         cmapTorsion.torsionBAtom2GlobalIndex,
                                         cmapTorsion.torsionBAtom3GlobalIndex,
                                         cmapTorsion.torsionBAtom4GlobalIndex);
            ++numAddedCmapTorsions;
        }

        std::cout << "\tAdded " << numAddedCmapTorsions << " CMAP torsions.\n";
        registerForce(cmapTorsionForce);
    }

    // Add Urey-Bradley Potential
    if (!systemTopology.ureyBradleys.empty()) {
        auto* ubForce = new OpenMM::HarmonicBondForce();
        int numAddedUreyBradleys = 0;

        for (const auto& term : systemTopology.ureyBradleys) {
            if (inSameMultiAtomRigidBody(term.atom1GlobalIndex, term.atom3GlobalIndex)) {
                continue;
            }

            ubForce->addBond(term.atom1GlobalIndex,
                             term.atom3GlobalIndex,
                             term.nominalLengthInNm,
                             term.stiffnessInKJPerNmSq * 2);
            ++numAddedUreyBradleys;
        }

        std::cout << "\tAdded " << numAddedUreyBradleys << " Urey-Bradley terms.\n";
        registerForce(ubForce);
    }
}

bool OPENMM::initialize(const std::vector<std::vector<int>>& worlds,
                        const SystemTopology& systemTopology,
                        const ForceFieldParams& ffParams,
                        const SimulationSettings& simSettings) {
    // Instantiate
    OPENMM& omm = get();

    // Set up variables
    omm.numAtoms = systemTopology.atoms.size();
    omm.ommAtomsPositionsCache = std::vector<OpenMM::Vec3>(omm.numAtoms);
    omm.ommAtomsPositionsCacheOld = std::vector<OpenMM::Vec3>(omm.numAtoms);
    omm.simbodyAtomsPositionsCache = std::vector<SimTK::Vec3>(omm.numAtoms);

    // Allocate OpenMM system and add particles to it
    omm.system = std::make_unique<OpenMM::System>();

    // Add atoms
    for (const auto& atom : systemTopology.atoms) {
        const SimTK::Real charge = atom.physics.chargeInE;
        const SimTK::Real sigma = atom.physics.sigmaInNm;
        const SimTK::Real epsilon = atom.physics.vdwWellDepthInKJ;

        if (!std::isfinite(charge)) {
            throw std::runtime_error("Non-finite charge for atom " + atom.identity.uniqueAtomName + ": "
                                     + std::to_string(charge));
        }

        if (!std::isfinite(sigma)) {
            throw std::runtime_error("Non-finite sigma for atom " + atom.identity.uniqueAtomName + ": "
                                     + std::to_string(sigma));
        }

        if (!std::isfinite(epsilon)) {
            throw std::runtime_error("Non-finite epsilon for atom " + atom.identity.uniqueAtomName + ": "
                                     + std::to_string(epsilon));
        }

        if (atom.connectivity.root) {
            omm.system->addParticle(0.0); // massless root
        } else {
            omm.system->addParticle(atom.physics.massInDaltons);
        }
    }

    // Add a force group for each world
    for (int forceGroupIndex = 0; forceGroupIndex < worlds.size(); ++forceGroupIndex) {
        omm.forceGroups.emplace_back();
        omm.forceGroups.back().initialize(forceGroupIndex,
                                          worlds[forceGroupIndex],
                                          systemTopology,
                                          ffParams,
                                          simSettings);

        // Gather all forces from all force groups into the system
        // Force groups indices have already been set
        for (const auto& force : omm.forceGroups.back().getForces()) {
            omm.system->addForce(force);
        }
    }

    // Create the integrator
    omm.integrator = std::make_unique<OpenMM::VerletIntegrator>(0.0007);

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

        const double speed = omm.context->getPlatform().getSpeed();
        std::cout << "Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed "
                  << speed << "\n";

    } catch (const std::exception& e) {
        std::cerr << "[ERROR]: Failed to create OpenMM Context.\n";
        std::cerr << "[ERROR]: " << e.what() << "\n";
        return false;
    }

    // All OpenMM components initialized successfully
    omm.initialized = true;

    // Set initial positions
    for (const auto& atom : systemTopology.atoms) {
        omm.ommAtomsPositionsCache[atom.identity.globalIndex] =
            OpenMM::Vec3(atom.position[0], atom.position[1], atom.position[2]);
    }
    omm.context->setPositions(omm.ommAtomsPositionsCache);

    // Log initialization
    std::cout << "[INFO] Initialized OpenMM. Using version " << platform->getOpenMMVersion() << ".\n";

    return true;
}

void OPENMM::evaluateEnergiesFromPositionCache(SimTK::Real& newPotentialEnergy,
                                               SimTK::Real& newKineticEnergy) {
    ensureInitialized();

    // Intentional value capture: relies on C++17 guaranteed copy elision.
    //
    const auto state =
        context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1 << activeForceGroupIndex);
    newPotentialEnergy = state.getPotentialEnergy();
    newKineticEnergy = state.getKineticEnergy();
}

auto OPENMM::integrateTrajectory(int steps, SimTK::Real timeStepInPicoseconds) -> bool {
    ensureInitialized();

    integrator->setStepSize(timeStepInPicoseconds);

    // Try to integrate
    bool success = true;
    try {
        integrator->setIntegrationForceGroups(1 << activeForceGroupIndex);
        integrator->step(steps);
    } catch (const std::exception& e) {
        // Restore old positions in case of integration failure
        ommAtomsPositionsCache = ommAtomsPositionsCacheOld;
        context->setPositions(ommAtomsPositionsCache);
        success = false;
    }

    // Intentional value capture: relies on C++17 guaranteed copy elision.
    const auto state =
        context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Velocities,
                          enforcePeriodicBox,
                          1 << activeForceGroupIndex);
    const auto& positions = state.getPositions();
    const auto& velocities = state.getVelocities();

    // DO NOT UNCOMMENT THIS - THIS IS A BUG AND BREAK SOMETHING SOMEWHERE Store energies
    potentialEnergy = state.getPotentialEnergy();
    kineticEnergy = state.getKineticEnergy();

    // Copy positions and velocities
    simbodyAtomsPositionsCache.resize(positions.size());
    simbodyAtomsVelocitiesCache.resize(velocities.size());

    for (size_t i = 0; i < positions.size(); ++i) {
        const auto& pos = positions[i];
        simbodyAtomsPositionsCache[i] = {pos[0], pos[1], pos[2]};
    }

    for (size_t i = 0; i < velocities.size(); ++i) {
        const auto& vel = velocities[i];
        simbodyAtomsVelocitiesCache[i] = {vel[0], vel[1], vel[2]};
    }

    return success;
}

void OPENMM::updatePositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings,
                                  const SimTK::Vector_<SimTK::Vec3>& inclAtomPos_G) {
    // TODO This loop transforms from SimTK::Vec3 to OpenMM::Vec3
    // When you have the time and consider that this is worth the effort, you may force OpenMM to accept
    // SimTK::Vec3 directly
    // 1. Context::setPositions which calls
    // 2. ContextImpl::setPositions which calls (for CUDA):
    // 3. CudaUpdateStateDataKernel::setPositions
    for (const auto& mapping : nonBondedMappings) {
        const std::size_t dAIx = mapping.dummAtomIndex;
        const std::size_t iax = mapping.includedAtomIndex;
        const SimTK::Vec3& pos_G = inclAtomPos_G[iax];

        ommAtomsPositionsCache[dAIx] = OpenMM::Vec3(pos_G[0], pos_G[1], pos_G[2]);
    }

    context->setPositions(ommAtomsPositionsCache);
}

void OPENMM::evaluateForcesFromPositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings,
                                              const SimTK::Vector_<SimTK::Vec3>& inclAtomStation_G,
                                              SimTK::Vector_<SimTK::SpatialVec>& inclBodyForces_G) const {
    ensureInitialized();

    // Intentional value capture: relies on C++17 guaranteed copy elision
    const auto state =
        context->getState(OpenMM::State::Forces, enforcePeriodicBox, 1 << activeForceGroupIndex);
    const auto& forces = state.getForces();

    // Map forces from atoms to bodies
    for (const auto& mapping : nonBondedMappings) {
        const std::size_t dAIx = mapping.dummAtomIndex;
        const std::size_t iax = mapping.includedAtomIndex;
        const std::size_t ibx = mapping.bodyIndex;

        const SimTK::Vec3 simForce(forces[dAIx][0], forces[dAIx][1], forces[dAIx][2]);
        inclBodyForces_G[ibx] += SimTK::SpatialVec(inclAtomStation_G[iax] % simForce, simForce);
    }
}

auto OPENMM::evaluatePotentialEnergyFromPositionsCache() const -> SimTK::Real {
    ensureInitialized();
    const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
    return state.getPotentialEnergy();
}

auto OPENMM::computePeriodicBoxVectors_Context(double a_length,
                                               double b_length,
                                               double c_length,
                                               double alpha,
                                               double beta,
                                               double gamma)
    -> std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> {
    ensureInitialized();

    const double TOL = 1e-6;

    // // Convert angles from degrees to radians
    // alpha = SimTK::Deg2Rad * alpha;
    // beta  = SimTK::Deg2Rad * beta;
    // gamma = SimTK::Deg2Rad * gamma;

    // Compute the box vectors
    OpenMM::Vec3 a(a_length, 0.0, 0.0);

    OpenMM::Vec3 b(b_length * std::cos(gamma), b_length * std::sin(gamma), 0.0);

    double cx = c_length * std::cos(beta);
    double cy = c_length * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    double cz = std::sqrt(c_length * c_length - cx * cx - cy * cy);

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
