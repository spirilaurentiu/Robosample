#include "bridge/OpenMMSystemBuilder.hpp"

#include "OpenMMContext.hpp"       // OpenMMContext::isPeriodic
#include "bridge/ForceFactory.hpp" // robo::forcefactory::create*Force
#include "bridge/MTSIntegrator.hpp"

#include "../../openmm/platforms/reference/include/ReferencePlatform.h"

#if USE_CPU
#    include "../../openmm/platforms/cpu/include/CpuPlatform.h"
#elif USE_OPENCL
#    include "../../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif USE_CUDA
#    include "../../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Force-group policy constants (see OpenMMSystemBuilder::build's addForce
// lambda): the only consumer of these, so they live here rather than as
// OpenMMContext members.
constexpr int kMtsSlowGroup = 0; // Nonbonded, GBSA, NBFIX
constexpr int kMtsFastGroup = 1; // bonds, angles, torsions, CMAP, UB

} // namespace

auto OpenMMSystemBuilder::build(const SystemTopology& systemTopology,
                                bool useMTS,
                                int mtsInnerSubsteps,
                                bool separateForceGroups,
                                AlchemyForceFactory& alchemyFactory) -> OpenMMSystemBuildResult {
    OpenMMSystemBuildResult result;
    result.numAtoms = static_cast<std::size_t>(systemTopology.numAtoms);

    result.system = std::make_unique<OpenMM::System>();
    for (int i = 0; i < systemTopology.numAtoms; ++i) {
        result.system->addParticle(systemTopology.atomsMass[i]);
    }

    // ------------------------------------------------------------------------
    //  Periodic box (explicit solvent). For any periodic method the box vectors
    //  MUST be set on the System BEFORE the Context is constructed -- PME builds
    //  its reciprocal-space grid from them at context-creation time. The vectors
    //  arrive already reduced (lower-triangular) from ParmEd, so we copy them in
    //  verbatim. OpenMM also requires the cutoff to be at most half the smallest
    //  box width; check it here so the failure is legible instead of a deep
    //  OpenMM assertion.
    // ------------------------------------------------------------------------
    if (OpenMMContext::isPeriodic(systemTopology.nonbondedMethod)) {
        if (systemTopology.boxVectors.size() != 9) {
            throw std::runtime_error("A periodic nonbonded method (CutoffPeriodic/Ewald/PME) requires "
                                     "box_vectors of length 9 (three reduced lattice vectors). Got "
                                     + std::to_string(systemTopology.boxVectors.size()) + ".");
        }
        const auto& bv = systemTopology.boxVectors;
        const OpenMM::Vec3 a(bv[0], bv[1], bv[2]);
        const OpenMM::Vec3 b(bv[3], bv[4], bv[5]);
        const OpenMM::Vec3 c(bv[6], bv[7], bv[8]);
        // Diagonal entries are the box widths for reduced (lower-triangular) vectors.
        const double minWidth = std::min({bv[0], bv[4], bv[8]});
        if (systemTopology.nonbondedCutoff > 0.5 * minWidth) {
            throw std::runtime_error("nonbonded_cutoff (" + std::to_string(systemTopology.nonbondedCutoff)
                                     + " nm) exceeds half the smallest box width (" + std::to_string(minWidth)
                                     + " nm). Lower the cutoff or enlarge the box.");
        }
        result.system->setDefaultPeriodicBoxVectors(a, b, c);
        std::cout << "[INFO] Periodic box set: a=" << bv[0] << " b=" << bv[4] << " c=" << bv[8]
                  << " nm; method="
                  << (systemTopology.nonbondedMethod == NonbondedMethod::PME     ? "PME"
                      : systemTopology.nonbondedMethod == NonbondedMethod::Ewald ? "Ewald"
                                                                                 : "CutoffPeriodic")
                  << ", cutoff=" << systemTopology.nonbondedCutoff << " nm.\n";
    }

    // ------------------------------------------------------------------------
    //  Virtual sites (extra points). Declared BEFORE the Context is created so
    //  OpenMM treats each EP as a dependent particle: the integrator does not
    //  move it, computeVirtualSites() places it from its parents, and its force
    //  is redistributed onto the parents. Without this an EP (mass 0) is frozen
    //  in space and detaches from its molecule as the molecule moves.
    // ------------------------------------------------------------------------
    result.hasVirtualSites = systemTopology.numVirtualSites > 0;
    for (int i = 0; i < systemTopology.numVirtualSites; ++i) {
        result.system->setVirtualSite(systemTopology.vsSite[i],
                                      new OpenMM::ThreeParticleAverageSite(systemTopology.vsAtom1[i],
                                                                           systemTopology.vsAtom2[i],
                                                                           systemTopology.vsAtom3[i],
                                                                           systemTopology.vsWeight1[i],
                                                                           systemTopology.vsWeight2[i],
                                                                           systemTopology.vsWeight3[i]));
    }
    if (result.hasVirtualSites) {
        std::cout << "[INFO] Declared " << systemTopology.numVirtualSites
                  << " virtual site(s) (3-particle average).\n";
    }

    // Safety net: every massless particle MUST be a declared virtual site. A
    // massless real particle is silently frozen by the integrator (this is
    // exactly how an undeclared water EP detaches and corrupts the electrostatics
    // without any crash). Fail loudly instead.
    {
        std::vector<bool> isVS(static_cast<std::size_t>(systemTopology.numAtoms), false);
        for (int i = 0; i < systemTopology.numVirtualSites; ++i) {
            isVS[static_cast<std::size_t>(systemTopology.vsSite[i])] = true;
        }
        for (int a = 0; a < systemTopology.numAtoms; ++a) {
            if (systemTopology.atomsMass[a] == 0.0 && !isVS[static_cast<std::size_t>(a)]) {
                throw std::runtime_error(
                    "Particle " + std::to_string(a)
                    + " is massless but is not a declared virtual site. Massless real particles are "
                      "frozen by the integrator (an extra point would detach from its molecule). "
                      "Declare it as a virtual site (populate system_topology.vs*), or give it mass.");
            }
        }
    }

    // GBSA needs implicit-solvent radii. If GBSA is requested but every radius is
    // zero (not populated by the prmtop reader), GBSA would be silently wrong, so
    // skip it with a warning rather than emit a bogus solvation energy.
    bool gbsaUsable = systemTopology.useGBSAOBC2;
    // Implicit (GBSA) and explicit (periodic) solvent are mutually exclusive --
    // a periodic box means real waters carry the solvation, so GBSA must be off.
    if (gbsaUsable && OpenMMContext::isPeriodic(systemTopology.nonbondedMethod)) {
        std::cerr << "[WARN] use_gbsa_obc2 is set together with a periodic nonbonded method; "
                     "GBSA (implicit) is incompatible with explicit solvent and will be skipped.\n";
        gbsaUsable = false;
    }
    if (gbsaUsable) {
        bool anyRadius = false;
        for (int i = 0; i < systemTopology.numAtoms; ++i) {
            if (systemTopology.atomsRadius[i] > 0.0) {
                anyRadius = true;
                break;
            }
        }
        if (!anyRadius) {
            std::cerr << "[WARN] use_gbsa_obc2 is set but all GBSA radii are 0; skipping GBSAOBCForce. "
                         "Populate atoms_radius / atoms_screen in the prmtop reader.\n";
            gbsaUsable = false;
        }
    }

    // Force-group assignment:
    //  * MTS on  -> slow forces to group 0, fast forces to group 1 (2 tiers).
    //  * MTS off, separateForceGroups -> one sequential group per force.
    //  * otherwise -> all in group 0.
    int nextForceGroup = 0;
    auto addForce = [&](OpenMM::Force* force, const std::string& name, bool slow) {
        if (force == nullptr) {
            return;
        }
        if (useMTS) {
            const int grp = slow ? kMtsSlowGroup : kMtsFastGroup;
            force->setForceGroup(grp);
            result.forceGroupLabels.emplace_back(grp, name);
        } else if (separateForceGroups) {
            if (nextForceGroup > 31) {
                throw std::runtime_error("Cannot give every force its own group: OpenMM supports at most "
                                         "32 force groups (0-31).");
            }
            force->setForceGroup(nextForceGroup);
            result.forceGroupLabels.emplace_back(nextForceGroup, name);
            ++nextForceGroup;
        }
        result.system->addForce(force);
    };

    // Slow (long-range) forces -> outer tier. Capture the main NonbondedForce: the
    // PME path mutates it in place (charge offsets) for alchemy.
    auto* mainNonbonded = robo::forcefactory::createNonbondedForce(systemTopology);
    addForce(mainNonbonded, "NonbondedForce", /*slow*/ true);
    if (gbsaUsable) {
        addForce(robo::forcefactory::createGBSAOBCForce(systemTopology), "GBSAOBCForce", /*slow*/ true);
    }
    if (systemTopology.hasNBfix) {
        addForce(robo::forcefactory::createCustomNonbondedForce(systemTopology), "CustomNonbondedForce",
                /*slow*/ true);
    }
    // NCMC: per-molecule intermolecular decoupling.
    if (alchemyFactory.enabled()) {
        if (systemTopology.hasNBfix) {
            throw std::runtime_error("NCMC alchemy + NBFIX is not supported in v1.");
        }
        if (OpenMMContext::isPeriodic(systemTopology.nonbondedMethod)) {
            // Explicit solvent (PME/Ewald/CutoffPeriodic): see AlchemyForceFactory::createAlchemyDecouplingForces.
            auto [soft, hard] = alchemyFactory.createAlchemyDecouplingForces(systemTopology, mainNonbonded);
            addForce(soft, "AlchemySoftcoreAxR", /*slow*/ true);
            addForce(hard, "AlchemyIntraAxA", /*slow*/ true);
            alchemyFactory.setAlchemyForce(nullptr); // lambda is driven through mainNonbonded's offset
        } else {
            // Vacuum/implicit: exact linear A x rest scaling (unchanged).
            auto* alchemyForce = alchemyFactory.createAlchemyCorrectionForce(systemTopology);
            alchemyFactory.setAlchemyForce(alchemyForce);
            addForce(alchemyForce, "AlchemyCorrection", /*slow*/ true);
        }
    }

    // Fast (bonded) forces -> inner tier.
    addForce(robo::forcefactory::createHarmonicBondForce(systemTopology), "HarmonicBondForce", /*slow*/ false);
    addForce(robo::forcefactory::createHarmonicAngleForce(systemTopology), "HarmonicAngleForce",
            /*slow*/ false);
    addForce(robo::forcefactory::createPeriodicTorsionForce(systemTopology), "PeriodicTorsionForce",
            /*slow*/ false);
    if (systemTopology.numHarmonicTorsions > 0) {
        addForce(robo::forcefactory::createImproperHarmonicTorsionForce(systemTopology), "CustomTorsionForce",
                /*slow*/ false);
    }
    if (systemTopology.cmapGridSize > 0) {
        addForce(robo::forcefactory::createCMAPTorsionForce(systemTopology), "CMAPTorsionForce",
                /*slow*/ false);
    }
    if (systemTopology.numUreyBradley > 0) {
        addForce(robo::forcefactory::createUreyBradleyForce(systemTopology), "HarmonicBondForce",
                /*slow*/ false);
    }

    if (useMTS) {
        // group 0 (slow) once per outer step, group 1 (fast) mtsInnerSubsteps times.
        std::vector<std::pair<int, int>> tiers{{kMtsSlowGroup, 1}, {kMtsFastGroup, mtsInnerSubsteps}};
        result.integrator = std::make_unique<MTSIntegrator>(0.001, tiers);
        std::cout << "[INFO] OpenMM MTS (r-RESPA) integrator: slow group every step, fast group x"
                  << mtsInnerSubsteps << ".\n";
    } else {
        result.integrator = std::make_unique<OpenMM::VerletIntegrator>(0.001);
    }

#if USE_CPU
    auto* platform = new OpenMM::CpuPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CPU";
#elif USE_REFERENCE
    auto* platform = new OpenMM::ReferencePlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "Reference";
#elif USE_OPENCL
    auto* platform = new OpenMM::OpenCLPlatform();
    OpenMM::Platform::registerPlatform(platform);
    platform->setPropertyDefaultValue("Precision", "mixed");
    constexpr auto PLATFORM_NAME = "OpenCL";
#elif USE_CUDA
    auto* platform = new OpenMM::CudaPlatform();
    OpenMM::Platform::registerPlatform(platform);
    platform->setPropertyDefaultValue("Precision", "mixed");
    constexpr auto PLATFORM_NAME = "CUDA";
#endif

    try {
        result.context = std::make_unique<OpenMM::Context>(*result.system, *result.integrator, *platform);
        const double speed = result.context->getPlatform().getSpeed();
        std::cout << "Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed "
                  << speed << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[ERROR]: Failed to create OpenMM Context.\n";
        std::cerr << "[ERROR]: " << e.what() << "\n";
        result.success = false;
        return result;
    }

    result.openMMVersion = platform->getOpenMMVersion();
    result.success = true;
    return result;
}
