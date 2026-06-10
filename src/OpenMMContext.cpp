#include "OpenMMContext.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>


#if USE_CPU
#    include "../openmm/platforms/cpu/include/CpuPlatform.h"
#elif USE_REFERENCE
#    include "../openmm/platforms/reference/include/ReferencePlatform.h"
#elif USE_OPENCL
#    include "../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif USE_CUDA
#    include "../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

class MTSIntegrator : public OpenMMContext::CustomIntegrator {
    public:
    MTSIntegrator(double stepSize, std::vector<std::pair<int, int>> groups)
        : OpenMMContext::CustomIntegrator(stepSize) {
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

auto OpenMMContext::initialize(const SystemTopology& systemTopology) -> bool {
    // Instantiate
    OpenMMContext& omm = get();

    // Set up variables
    numAtoms = systemTopology.numAtoms;

    // Allocate OpenMM system and add particles to it
    system = std::make_unique<OpenMMContext::System>();

    for (int i = 0; i < systemTopology.numAtoms; ++i) {
        // const bool isRoot = atom.connectivity.root;
        // const auto moleculeIndex = atom.identity.moleculeIndex;
        // const bool isWeld = systemTopology.rootMobilities[moleculeIndex] == SimTK::RootMobility::Weld;

        system->addParticle(systemTopology.atomsMass[i]);

        // if (isRoot && isWeld) {
        //     system->addParticle(0.0); // massless root
        // } else {
        //     system->addParticle(atom.physics.massInDaltons);
        // }
    }

    // Nonbonded force
    auto* nonbondedForce = createNonbondedForce(systemTopology);
    system->addForce(nonbondedForce);

    // GBSA OBC Force
    // implicitSolventKappa
    if (systemTopology.useGBSAOBC2) {
        auto* gbsaOBCForce = createGBSAOBCForce(systemTopology);
        system->addForce(gbsaOBCForce);
    }

    // 1-4 VdW Correction (Bond-based so we can target specific pairs)
    if (systemTopology.hasNBfix) {
        auto* customNonbondedForce = createCustomNonbondedForce(systemTopology);
        system->addForce(customNonbondedForce);
    }

    // Add bonds
    auto* harmonicBondForce = createHarmonicBondForce(systemTopology);
    system->addForce(harmonicBondForce);

    // Add angles
    auto* harmonicAngleForce = createHarmonicAngleForce(systemTopology);
    system->addForce(harmonicAngleForce);

    // For AMBER style torsions, OpenMM does not care if they are proper or improper
    // Note that OpenMM only supports only periodic torsions by default (AMBER style), so we need to handle
    // improper harmonic torsions differently
    auto* periodicTorsionForce = createPeriodicTorsionForce(systemTopology);
    system->addForce(periodicTorsionForce);

    // Harmonic improper torsions (CHARMM) are handled using a CustomTorsionForce since OpenMM does not
    // support them natively
    if (systemTopology.numHarmonicTorsions > 0) {
        auto* improperTorsionForce = createImproperHarmonicTorsionForce(systemTopology);
        system->addForce(improperTorsionForce);
    }

    // // Add correction map torsions (CMAPs) force
    // if (systemTopology.cmapGridSize > 0) {
    //     auto* cmapTorsionForce = createCMAPTorsionForce();
    //     system->addForce(cmapTorsionForce);
    // }

    // // Add Urey-Bradley Potential
    // if (systemTopology.numUreyBradley > 0) {
    //     auto* ubForce = createUreyBradleyForce();
    //     system->addForce(ubForce);
    // }

    // Create the integrator
    // std::vector<std::pair<int, int>> groups = {{SLOW_GROUP, 1}, {FAST_GROUP, 4}};
    integrator = std::make_unique<OpenMMContext::VerletIntegrator>(0.001);

#if USE_CPU
    OpenMMContext::Platform* platform = new OpenMMContext::CpuPlatform();
    OpenMMContext::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CPU";

#elif USE_REFERENCE
    OpenMMContext::Platform* platform = new OpenMMContext::ReferencePlatform();
    OpenMMContext::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "Reference";

#elif USE_OPENCL
    OpenMMContext::Platform* platform = new OpenMMContext::OpenCLPlatform();
    OpenMMContext::Platform::registerPlatform(platform);
    platform->setPropertyDefaultValue("Precision", "mixed");
    constexpr auto PLATFORM_NAME = "OpenCL";

#elif USE_CUDA
    OpenMMContext::Platform* platform = new OpenMMContext::CudaPlatform();
    OpenMMContext::Platform::registerPlatform(platform);
    platform->setPropertyDefaultValue("Precision", "mixed");
    constexpr auto PLATFORM_NAME = "CUDA";
#endif

    try {
        context = std::make_unique<OpenMMContext::Context>(*system, *integrator, *platform);

        const double speed = context->getPlatform().getSpeed();
        std::cout << "Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed "
                  << speed << "\n";

    } catch (const std::exception& e) {
        std::cerr << "[ERROR]: Failed to create OpenMM Context.\n";
        std::cerr << "[ERROR]: " << e.what() << "\n";
        return false;
    }

    // All OpenMM components initialized successfully
    initialized = true;

    // Log initialization
    std::cout << "[INFO] Initialized OpenMM. Using version " << platform->getOpenMMVersion() << ".\n";

    return true;
}

auto OpenMMContext::integrateTrajectory(int steps, double timeStepInPicoseconds) -> bool {
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
    const auto state = context->getState(OpenMMContext::State::Positions | OpenMMContext::State::Energy
                                             | OpenMMContext::State::Velocities,
                                         enforcePeriodicBox);
    const auto& positions = state.getPositions();
    const auto& velocities = state.getVelocities();

    potentialEnergy = state.getPotentialEnergy();
    kineticEnergy = state.getKineticEnergy();

    return success;
}

void OpenMMContext::evaluateForcesFromPositionsCache(const std::vector<OpenMMContext::Vec3>& positions,
                                                     std::vector<OpenMMContext::Vec3>& outForces) const {
    ensureInitialized();

    context->setPositions(positions);
    const auto state = context->getState(OpenMMContext::State::Forces, enforcePeriodicBox);
    outForces = state.getForces();
}

auto OpenMMContext::computePeriodicBoxVectors_Context(double a_length,
                                                      double b_length,
                                                      double c_length,
                                                      double alpha,
                                                      double beta,
                                                      double gamma)
    -> std::tuple<OpenMMContext::Vec3, OpenMMContext::Vec3, OpenMMContext::Vec3> {
    ensureInitialized();

    const double TOL = 1e-6;

    // // Convert angles from degrees to radians
    // alpha = SimTK::Deg2Rad * alpha;
    // beta  = SimTK::Deg2Rad * beta;
    // gamma = SimTK::Deg2Rad * gamma;

    // Compute the box vectors
    OpenMMContext::Vec3 a(a_length, 0.0, 0.0);

    OpenMMContext::Vec3 b(b_length * std::cos(gamma), b_length * std::sin(gamma), 0.0);

    double cx = c_length * std::cos(beta);
    double cy = c_length * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    double cz = std::sqrt(c_length * c_length - cx * cx - cy * cy);

    OpenMMContext::Vec3 c(cx, cy, cz);

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

auto OpenMMContext::createNonbondedForce(const SystemTopology& systemTopology)
    -> OpenMMContext::NonbondedForce* {
    auto* nonbondedForce = new OpenMMContext::NonbondedForce();
    nonbondedForce->setCutoffDistance(systemTopology.nonbondedCutoff);

    switch (systemTopology.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            nonbondedForce->setNonbondedMethod(OpenMMContext::NonbondedForce::NoCutoff);
            nonbondedForce->setUseDispersionCorrection(false);
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            nonbondedForce->setNonbondedMethod(OpenMMContext::NonbondedForce::CutoffNonPeriodic);
            if (systemTopology.useGBSAOBC2) {
                nonbondedForce->setReactionFieldDielectric(1.0);
                nonbondedForce->setUseDispersionCorrection(false);
            } else {
                nonbondedForce->setUseDispersionCorrection(true);
            }
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method");
    }

    // Add all particles
    for (int index = 0; index < systemTopology.numAtoms; ++index) {
        if (systemTopology.hasNBfix) {
            index = nonbondedForce->addParticle(systemTopology.atomsCharge[index], 1.0, 0.0);
        } else {
            index = nonbondedForce->addParticle(systemTopology.atomsCharge[index],
                                                systemTopology.atomsSigma[index],
                                                systemTopology.atomsEpsilon[index]);
        }
    }

    // Add 1-4 scalings
    for (int index = 0; index < systemTopology.numScaling14; ++index) {
        nonbondedForce->addException(systemTopology.scaling14I[index],
                                     systemTopology.scaling14L[index],
                                     systemTopology.scaling14ChargeProduct[index],
                                     systemTopology.scaling14Sigma[index],
                                     systemTopology.scaling14Epsilon[index]);
    }

    // Add 1-2 and 1-3 exclusions
    for (int index = 0; index < systemTopology.numExclusions; ++index) {
        // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely
        // omitted from force and energy calculations
        nonbondedForce->addException(systemTopology.exclusionI[index],
                                     systemTopology.exclusionJ[index],
                                     0.0,
                                     0.1,
                                     0.0);
    }

    return nonbondedForce;
}

auto OpenMMContext::createGBSAOBCForce(const SystemTopology& systemTopology) -> OpenMMContext::GBSAOBCForce* {
    OpenMMContext::GBSAOBCForce::NonbondedMethod gbsaForceMethod;
    switch (systemTopology.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            gbsaForceMethod = OpenMMContext::GBSAOBCForce::NoCutoff;
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            gbsaForceMethod = OpenMMContext::GBSAOBCForce::CutoffNonPeriodic;
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method for GBSAOBCForce");
    }

    auto* force = new OpenMMContext::GBSAOBCForce();

    force->setSolventDielectric(systemTopology.gbsaSolventDielectric);
    force->setSoluteDielectric(systemTopology.gbsaSoluteDielectric);
    force->setNonbondedMethod(gbsaForceMethod);
    force->setCutoffDistance(systemTopology.nonbondedCutoff);

    for (int index = 0; index < systemTopology.numAtoms; ++index) {
        force->addParticle(systemTopology.atomsCharge[index],
                           systemTopology.atomsRadius[index],
                           systemTopology.atomsScreen[index]);
    }

    return force;
}

auto OpenMMContext::createCustomNonbondedForce(const SystemTopology& systemTopology)
    -> OpenMMContext::CustomNonbondedForce* {
    throw std::runtime_error("not implemented");
    // -> OpenMMContext::CustomNonbondedForce* {
    // auto* force = new OpenMMContext::CustomNonbondedForce(
    //     "(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2);");
    // force->addTabulatedFunction(
    //     "acoef",
    //     new OpenMMContext::Discrete2DFunction(systemTopology.numNBTypes, systemTopology.numNBTypes,
    //     systemTopology.aCoef));
    // force->addTabulatedFunction(
    //     "bcoef",
    //     new OpenMMContext::Discrete2DFunction(systemTopology.numNBTypes, systemTopology.numNBTypes,
    //     systemTopology.bCoef));
    // force->addPerParticleParameter("type");

    // force->setCutoffDistance(systemTopology.nonbondedCutoff);
    // force->setUseSwitchingFunction(false);
    // force->setSwitchingDistance(0);

    // switch (systemTopology.nonbondedMethod) {
    //     case NonbondedMethod::NoCutoff:
    //         force->setNonbondedMethod(OpenMMContext::CustomNonbondedForce::NoCutoff);
    //         break;
    //     case NonbondedMethod::CutoffNonPeriodic:
    //         force->setNonbondedMethod(OpenMMContext::CustomNonbondedForce::CutoffNonPeriodic);
    //         break;
    //     default:
    //         throw std::invalid_argument("Unsupported nonbonded method");
    // }

    // for (const auto& atom : atoms) {
    //     force->addParticle({static_cast<double>(atom.identity.nonbondedIndex)});
    // }

    // for (const auto& scaling14 : scaling14s) {
    //     force->addExclusion(scaling14.atom1GlobalIndex, scaling14.atom4GlobalIndex);
    // }

    // // Exclude 1-2 and 1-3 interactions
    // for (const auto& exclusion : exclusions) {
    //     // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely
    //     // omitted from force and energy calculations
    //     force->addExclusion(exclusion.atom1GlobalIndex, exclusion.atom2GlobalIndex);
    // }

    // return force;
}

auto OpenMMContext::createHarmonicBondForce(const SystemTopology& systemTopology)
    -> OpenMMContext::HarmonicBondForce* {
    auto* force = new OpenMMContext::HarmonicBondForce();

    for (int index = 0; index < systemTopology.numBonds; ++index) {
        const auto particle1 = systemTopology.bondsI[index];
        const auto particle2 = systemTopology.bondsJ[index];
        const auto length = systemTopology.bondsEquilibrium[index];

        // OpenMM defines the harmonic bond potential as 0.5 * k * (r-r0)^2, so we need to multiply by 2 to
        // get the correct stiffness
        auto stiffness = systemTopology.bondsStiffness[index] * 2;

        // Add bond and save the interaction index
        force->addBond(particle1, particle2, length, stiffness);
    }

    return force;
}

auto OpenMMContext::createHarmonicAngleForce(const SystemTopology& systemTopology)
    -> OpenMMContext::HarmonicAngleForce* {
    auto* force = new OpenMMContext::HarmonicAngleForce();

    for (int index = 0; index < systemTopology.numAngles; ++index) {
        const auto particle1 = systemTopology.anglesI[index];
        const auto particle2 = systemTopology.anglesJ[index];
        const auto particle3 = systemTopology.anglesK[index];
        const auto angleInRad = systemTopology.anglesEquilibrium[index];

        // OpenMM defines the harmonic angle potential as 0.5 * k * (theta-theta0)^2, so we need to
        // multiply by 2 to get the correct stiffness
        auto stiffness = systemTopology.anglesStiffness[index] * 2;
        force->addAngle(particle1, particle2, particle3, angleInRad, stiffness);
    }
    return force;
}

auto OpenMMContext::createPeriodicTorsionForce(const SystemTopology& systemTopology)
    -> OpenMMContext::PeriodicTorsionForce* {
    auto* force = new OpenMMContext::PeriodicTorsionForce();

    for (int index = 0; index < systemTopology.numPeriodicTorsions; ++index) {
        const auto particle1 = systemTopology.periodicTorsionsI[index];
        const auto particle2 = systemTopology.periodicTorsionsJ[index];
        const auto particle3 = systemTopology.periodicTorsionsK[index];
        const auto particle4 = systemTopology.periodicTorsionsL[index];
        const auto periodicity = systemTopology.periodicTorsionsN[index];
        const auto phase = systemTopology.periodicTorsionsPhase[index];
        auto stiffness = systemTopology.periodicTorsionsStiffness[index];

        force->addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, stiffness);
    }

    return force;
}

auto OpenMMContext::createImproperHarmonicTorsionForce(const SystemTopology& systemTopology)
    -> OpenMMContext::CustomTorsionForce* {
    throw std::runtime_error("not implemented");

    // // Create force expression for harmonic improper torsions (CHARMM style)
    // std::stringstream strStream;
    // strStream << std::setprecision(17) << "k*min(dtheta, 2*" << SimTK::Pi
    //           << "-dtheta)^2; dtheta = abs(theta-theta0)";

    // auto* force = new OpenMMContext::CustomTorsionForce(strStream.str());
    // force->addPerTorsionParameter("k");
    // force->addPerTorsionParameter("theta0");

    // for (const auto& torsion : systemTopology.harmonicImproperTorsions) {
    //     const auto atom1GlobalIndex = torsion.globalIndices[0];
    //     const auto atom2GlobalIndex = torsion.globalIndices[1];
    //     const auto atom3GlobalIndex = torsion.globalIndices[2];
    //     const auto atom4GlobalIndex = torsion.globalIndices[3];

    //     const std::vector<double> params = {torsion.stiffnessInKJPerRadSq, torsion.nominalAngleInRad};
    //     force->addTorsion(atom1GlobalIndex, atom2GlobalIndex, atom3GlobalIndex, atom4GlobalIndex, params);
    // }

    // return force;
}

auto OpenMMContext::createCMAPTorsionForce(const SystemTopology& systemTopology)
    -> OpenMMContext::CMAPTorsionForce* {
    return nullptr;

    // auto* force = new OpenMMContext::CMAPTorsionForce();
    // force->setUsesPeriodicBoundaryConditions(enforcePeriodicBox);

    // for (int index = 0; index < systemTopology.cmapGridSize; ++index) {
    //     force->addMap(grid.size, grid.energy);
    // }

    // for (const auto& cmapTorsion : systemTopology.cmapTorsions) {
    //     force->addTorsion(cmapTorsion.mapIndex,
    //                       cmapTorsion.torsionAAtom1GlobalIndex,
    //                       cmapTorsion.torsionAAtom2GlobalIndex,
    //                       cmapTorsion.torsionAAtom3GlobalIndex,
    //                       cmapTorsion.torsionAAtom4GlobalIndex,
    //                       cmapTorsion.torsionBAtom1GlobalIndex,
    //                       cmapTorsion.torsionBAtom2GlobalIndex,
    //                       cmapTorsion.torsionBAtom3GlobalIndex,
    //                       cmapTorsion.torsionBAtom4GlobalIndex);
    // }

    // return force;
}

auto OpenMMContext::createUreyBradleyForce(const SystemTopology& systemTopology)
    -> OpenMMContext::HarmonicBondForce* {
    // auto* force = new OpenMMContext::HarmonicBondForce();

    // for (const auto& term : systemTopology.ureyBradleys) {
    //     force->addBond(term.atom1GlobalIndex,
    //                    term.atom3GlobalIndex,
    //                    term.nominalLengthInNm,
    //                    term.stiffnessInKJPerNmSq * 2);
    // }

    // return force;

    return nullptr;
}
