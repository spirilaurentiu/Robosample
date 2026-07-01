#include "OpenMMContext.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
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

namespace {
constexpr double kPi = 3.14159265358979323846;
} // namespace

// Multiple-timestep (MTS) r-RESPA Verlet integrator. groups = (force group,
// substeps); slower groups have fewer substeps. Built by initialize() when MTS
// is enabled. (Reversible + symplectic -> safe inside HMC acceptance.)
class MTSIntegrator : public OpenMM::CustomIntegrator {
    public:
    MTSIntegrator(double stepSize, std::vector<std::pair<int, int>> groups)
        : OpenMM::CustomIntegrator(stepSize) {
        if (groups.empty()) {
            throw std::invalid_argument("No force groups specified");
        }
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
        const int stepsPerParentStep = substeps / parentSubsteps;
        const std::string n = std::to_string(substeps);
        const std::string g = std::to_string(group);
        for (int i = 0; i < stepsPerParentStep; ++i) {
            addComputePerDof("v", "v+0.5*(dt/" + n + ")*f" + g + "/m");
            if (groups.size() == 1) {
                addComputePerDof("x", "x+(dt/" + n + ")*v");
                addComputePerDof("x1", "x");
                addConstrainPositions();
                addComputePerDof("v", "v+(x-x1)/(dt/" + n + ")");
                addConstrainVelocities();
            } else {
                createSubsteps(substeps, {groups.begin() + 1, groups.end()});
            }
            addComputePerDof("v", "v+0.5*(dt/" + n + ")*f" + g + "/m");
        }
    }
};

auto OpenMMContext::initialize(const SystemTopology& systemTopology) -> bool {
    numAtoms = static_cast<std::size_t>(systemTopology.numAtoms);
    forceGroupLabels.clear();

    system = std::make_unique<OpenMM::System>();
    for (int i = 0; i < systemTopology.numAtoms; ++i) {
        system->addParticle(systemTopology.atomsMass[i]);
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
    if (isPeriodic(systemTopology.nonbondedMethod)) {
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
        system->setDefaultPeriodicBoxVectors(a, b, c);
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
    hasVirtualSites = systemTopology.numVirtualSites > 0;
    for (int i = 0; i < systemTopology.numVirtualSites; ++i) {
        system->setVirtualSite(systemTopology.vsSite[i],
                               new OpenMM::ThreeParticleAverageSite(systemTopology.vsAtom1[i],
                                                                    systemTopology.vsAtom2[i],
                                                                    systemTopology.vsAtom3[i],
                                                                    systemTopology.vsWeight1[i],
                                                                    systemTopology.vsWeight2[i],
                                                                    systemTopology.vsWeight3[i]));
    }
    if (hasVirtualSites) {
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
    if (gbsaUsable && isPeriodic(systemTopology.nonbondedMethod)) {
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
            forceGroupLabels.emplace_back(grp, name);
        } else if (separateForceGroups) {
            if (nextForceGroup > 31) {
                throw std::runtime_error("Cannot give every force its own group: OpenMM supports at most "
                                         "32 force groups (0-31).");
            }
            force->setForceGroup(nextForceGroup);
            forceGroupLabels.emplace_back(nextForceGroup, name);
            ++nextForceGroup;
        }
        system->addForce(force);
    };

    // Slow (long-range) forces -> outer tier. Capture the main NonbondedForce: the
    // PME path mutates it in place (charge offsets) for alchemy.
    auto* mainNonbonded = createNonbondedForce(systemTopology);
    addForce(mainNonbonded, "NonbondedForce", /*slow*/ true);
    if (gbsaUsable) {
        addForce(createGBSAOBCForce(systemTopology), "GBSAOBCForce", /*slow*/ true);
    }
    if (systemTopology.hasNBfix) {
        addForce(createCustomNonbondedForce(systemTopology), "CustomNonbondedForce", /*slow*/ true);
    }
    // NCMC: per-molecule intermolecular decoupling.
    if (alchemyEnabled) {
        if (systemTopology.hasNBfix) {
            throw std::runtime_error("NCMC alchemy + NBFIX is not supported in v1.");
        }
        if (isPeriodic(systemTopology.nonbondedMethod)) {
            // Explicit solvent (PME/Ewald/CutoffPeriodic): see createAlchemyDecouplingForces.
            auto [soft, hard] = createAlchemyDecouplingForces(systemTopology, mainNonbonded);
            addForce(soft, "AlchemySoftcoreAxR", /*slow*/ true);
            addForce(hard, "AlchemyIntraAxA", /*slow*/ true);
            alchemyForce = nullptr; // lambda is driven through mainNonbonded's offset
        } else {
            // Vacuum/implicit: exact linear A x rest scaling (unchanged).
            alchemyForce = createAlchemyCorrectionForce(systemTopology);
            addForce(alchemyForce, "AlchemyCorrection", /*slow*/ true);
        }
    }

    // Fast (bonded) forces -> inner tier.
    addForce(createHarmonicBondForce(systemTopology), "HarmonicBondForce", /*slow*/ false);
    addForce(createHarmonicAngleForce(systemTopology), "HarmonicAngleForce", /*slow*/ false);
    addForce(createPeriodicTorsionForce(systemTopology), "PeriodicTorsionForce", /*slow*/ false);
    if (systemTopology.numHarmonicTorsions > 0) {
        addForce(createImproperHarmonicTorsionForce(systemTopology), "CustomTorsionForce", /*slow*/ false);
    }
    if (systemTopology.cmapGridSize > 0) {
        addForce(createCMAPTorsionForce(systemTopology), "CMAPTorsionForce", /*slow*/ false);
    }
    if (systemTopology.numUreyBradley > 0) {
        addForce(createUreyBradleyForce(systemTopology), "HarmonicBondForce", /*slow*/ false);
    }

    if (useMTS) {
        // group 0 (slow) once per outer step, group 1 (fast) mtsInnerSubsteps times.
        std::vector<std::pair<int, int>> tiers{{kMtsSlowGroup, 1}, {kMtsFastGroup, mtsInnerSubsteps}};
        integrator = std::make_unique<MTSIntegrator>(0.001, tiers);
        std::cout << "[INFO] OpenMM MTS (r-RESPA) integrator: slow group every step, fast group x"
                  << mtsInnerSubsteps << ".\n";
    } else {
        integrator = std::make_unique<OpenMM::VerletIntegrator>(0.001);
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
        context = std::make_unique<OpenMM::Context>(*system, *integrator, *platform);
        const double speed = context->getPlatform().getSpeed();
        std::cout << "Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed "
                  << speed << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[ERROR]: Failed to create OpenMM Context.\n";
        std::cerr << "[ERROR]: " << e.what() << "\n";
        return false;
    }

    initialized = true;
    std::cout << "[INFO] Initialized OpenMM. Using version " << platform->getOpenMMVersion() << ".\n";
    return true;
}

auto OpenMMContext::computePotentialEnergy(const std::vector<OpenMM::Vec3>& positions) -> double {
    ensureInitialized();
    context->setPositions(positions);
    if (hasVirtualSites) {
        context->computeVirtualSites();
    }
    const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
    potentialEnergy = state.getPotentialEnergy();
    return potentialEnergy;
}

auto OpenMMContext::computePotentialEnergyByGroup(const std::vector<OpenMM::Vec3>& positions)
    -> std::pair<double, std::vector<OpenMMContext::ForceGroupEnergy>> {
    ensureInitialized();
    context->setPositions(positions);
    if (hasVirtualSites) {
        context->computeVirtualSites();
    }
    const auto totalState = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
    potentialEnergy = totalState.getPotentialEnergy();

    std::vector<ForceGroupEnergy> breakdown;
    if ((separateForceGroups || useMTS) && !forceGroupLabels.empty()) {
        breakdown.reserve(forceGroupLabels.size());
        for (const auto& [group, name] : forceGroupLabels) {
            const int groupMask = 1 << group;
            const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, groupMask);
            breakdown.push_back(ForceGroupEnergy{group, name, state.getPotentialEnergy()});
        }
    } else {
        breakdown.push_back(ForceGroupEnergy{0, "All", potentialEnergy});
    }
    return {potentialEnergy, breakdown};
}

auto OpenMMContext::integrateTrajectory(int steps, double timeStepInPicoseconds) -> bool {
    ensureInitialized();
    integrator->setStepSize(timeStepInPicoseconds);
    bool success = true;
    try {
        integrator->step(steps);
    } catch (const std::exception&) {
        success = false;
    }
    const auto state =
        context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Velocities,
                          enforcePeriodicBox);
    potentialEnergy = state.getPotentialEnergy();
    kineticEnergy = state.getKineticEnergy();
    return success;
}

void OpenMMContext::evaluateForcesFromPositionsCache(const std::vector<OpenMM::Vec3>& positions,
                                                     std::vector<OpenMM::Vec3>& outForces) const {
    ensureInitialized();
    context->setPositions(positions);
    if (hasVirtualSites) {
        context->computeVirtualSites();
    }
    const auto state = context->getState(OpenMM::State::Forces, enforcePeriodicBox);
    outForces = state.getForces();
}

auto OpenMMContext::computePeriodicBoxVectors_Context(double a_length,
                                                      double b_length,
                                                      double c_length,
                                                      double alpha,
                                                      double beta,
                                                      double gamma)
    -> std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> {
    ensureInitialized();
    constexpr double TOL = 1e-6;
    OpenMM::Vec3 a(a_length, 0.0, 0.0);
    OpenMM::Vec3 b(b_length * std::cos(gamma), b_length * std::sin(gamma), 0.0);
    const double cx = c_length * std::cos(beta);
    const double cy = c_length * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    const double cz = std::sqrt(c_length * c_length - cx * cx - cy * cy);
    OpenMM::Vec3 c(cx, cy, cz);
    for (int i = 0; i < 3; ++i) {
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

auto OpenMMContext::createNonbondedForce(const SystemTopology& systemTopology) -> OpenMM::NonbondedForce* {
    auto* nonbondedForce = new OpenMM::NonbondedForce();
    nonbondedForce->setCutoffDistance(systemTopology.nonbondedCutoff);
    switch (systemTopology.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
            nonbondedForce->setUseDispersionCorrection(false);
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
            if (systemTopology.useGBSAOBC2) {
                nonbondedForce->setReactionFieldDielectric(1.0);
                nonbondedForce->setUseDispersionCorrection(false);
            } else {
                nonbondedForce->setUseDispersionCorrection(true);
            }
            break;
        case NonbondedMethod::CutoffPeriodic:
            // Periodic reaction-field cutoff. Explicit solvent => isotropic
            // long-range dispersion correction on.
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
            nonbondedForce->setUseDispersionCorrection(true);
            break;
        case NonbondedMethod::Ewald:
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::Ewald);
            nonbondedForce->setEwaldErrorTolerance(systemTopology.ewaldErrorTolerance);
            nonbondedForce->setUseDispersionCorrection(true);
            break;
        case NonbondedMethod::PME:
            // Particle-Mesh Ewald: the standard explicit-solvent electrostatics.
            // The box was already set on the System above. Exceptions/1-4 pairs
            // added below are PME-aware (OpenMM applies the reciprocal-space
            // correction for excluded pairs automatically).
            nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::PME);
            nonbondedForce->setEwaldErrorTolerance(systemTopology.ewaldErrorTolerance);
            nonbondedForce->setUseDispersionCorrection(true);
            break;
        default:
            throw std::invalid_argument("Unsupported nonbonded method");
    }
    for (int index = 0; index < systemTopology.numAtoms; ++index) {
        if (systemTopology.hasNBfix) {
            nonbondedForce->addParticle(systemTopology.atomsCharge[index], 1.0, 0.0);
        } else {
            nonbondedForce->addParticle(systemTopology.atomsCharge[index],
                                        systemTopology.atomsSigma[index],
                                        systemTopology.atomsEpsilon[index]);
        }
    }
    for (int index = 0; index < systemTopology.numScaling14; ++index) {
        nonbondedForce->addException(systemTopology.scaling14I[index],
                                     systemTopology.scaling14L[index],
                                     systemTopology.scaling14ChargeProduct[index],
                                     systemTopology.scaling14Sigma[index],
                                     systemTopology.scaling14Epsilon[index]);
    }
    for (int index = 0; index < systemTopology.numExclusions; ++index) {
        nonbondedForce->addException(systemTopology.exclusionI[index],
                                     systemTopology.exclusionJ[index],
                                     0.0,
                                     0.1,
                                     0.0);
    }
    return nonbondedForce;
}

auto OpenMMContext::createGBSAOBCForce(const SystemTopology& systemTopology) -> OpenMM::GBSAOBCForce* {
    OpenMM::GBSAOBCForce::NonbondedMethod gbsaForceMethod = OpenMM::GBSAOBCForce::NoCutoff;
    switch (systemTopology.nonbondedMethod) {
        case NonbondedMethod::NoCutoff:
            gbsaForceMethod = OpenMM::GBSAOBCForce::NoCutoff;
            break;
        case NonbondedMethod::CutoffNonPeriodic:
            gbsaForceMethod = OpenMM::GBSAOBCForce::CutoffNonPeriodic;
            break;
        default:
            // CutoffPeriodic/Ewald/PME mean explicit solvent -- GBSA is implicit
            // and must never be combined with it. initialize() already skips GBSA
            // for periodic methods, so reaching here is a logic error.
            throw std::invalid_argument("GBSA (implicit solvent) is incompatible with a periodic "
                                        "nonbonded method (explicit solvent).");
    }
    auto* force = new OpenMM::GBSAOBCForce();
    force->setSolventDielectric(systemTopology.gbsaSolventDielectric); // default 78.5
    force->setSoluteDielectric(systemTopology.gbsaSoluteDielectric);   // default 1.0
    // Match OpenMM's app-layer implicit solvent. `createSystem(implicitSolvent=
    // OBC2)` builds a CustomGBForce with the POLAR GB term only -- it has NO
    // nonpolar surface-area term. The built-in GBSAOBCForce, by contrast, adds
    // an ACE SA term by default (surfaceAreaEnergy = 2.25936 kJ/mol/nm^2). That
    // term is exactly the ~15 kJ/mol that made our implicit-solvent energy
    // diverge from the OpenMM reference, so zero it to reproduce app.OBC2.
    force->setSurfaceAreaEnergy(0.0);
    force->setNonbondedMethod(gbsaForceMethod);
    force->setCutoffDistance(systemTopology.nonbondedCutoff);
    for (int index = 0; index < systemTopology.numAtoms; ++index) {
        force->addParticle(systemTopology.atomsCharge[index],
                           systemTopology.atomsRadius[index],
                           systemTopology.atomsScreen[index]);
    }
    return force;
}

auto OpenMMContext::createCustomNonbondedForce(const SystemTopology& /*systemTopology*/)
    -> OpenMM::CustomNonbondedForce* {
    throw std::runtime_error("createCustomNonbondedForce (NBFIX) not implemented yet");
}

auto OpenMMContext::createAlchemyCorrectionForce(const SystemTopology& sys) -> OpenMM::CustomNonbondedForce* {
    // Total [begin,end) x rest pair energy = standard + (lambda_inter-1)*standard
    //                                      = lambda_inter * standard.
    // Lorentz-Berthelot combining, matching OpenMM NonbondedForce defaults. All
    // 1-4/exclusion pairs are intramolecular, so the A x rest interaction group
    // carries no exceptions and needs no exclusion list (scales to assemblies).
    const double ONE_4PI_EPS0 = 138.935456; // kJ*nm/(mol*e^2)
    auto* f =
        new OpenMM::CustomNonbondedForce("(lambda_inter - 1)*(4*eps*((sig/r)^12 - (sig/r)^6) + k*q1*q2/r);"
                                         "eps=sqrt(eps1*eps2); sig=0.5*(sig1+sig2)");
    f->addGlobalParameter("lambda_inter", 1.0);
    f->addGlobalParameter("k", ONE_4PI_EPS0);
    f->addPerParticleParameter("q");
    f->addPerParticleParameter("sig");
    f->addPerParticleParameter("eps");
    for (int i = 0; i < sys.numAtoms; ++i) {
        std::vector<double> p{sys.atomsCharge[i], sys.atomsSigma[i], sys.atomsEpsilon[i]};
        f->addParticle(p);
    }
    if (sys.nonbondedMethod == NonbondedMethod::CutoffNonPeriodic) {
        f->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        f->setCutoffDistance(sys.nonbondedCutoff);
    } else {
        f->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
    }
    std::set<int> aSet, restSet;
    for (int i = 0; i < sys.numAtoms; ++i) {
        if (i >= alchemyBegin && i < alchemyEnd) {
            aSet.insert(i);
        } else {
            restSet.insert(i);
        }
    }
    f->addInteractionGroup(aSet, restSet); // A x rest ONLY
    return f;
}

void OpenMMContext::enableAlchemy(int atomBegin, int atomEnd) {
    alchemyEnabled = true;
    alchemyBegin = atomBegin;
    alchemyEnd = atomEnd;
}

void OpenMMContext::setAlchemicalLambda(double lambdaInter) {
    if (!alchemyEnabled) {
        return;
    }
    ensureInitialized();
    context->setParameter("lambda_inter", lambdaInter);
}

auto OpenMMContext::createHarmonicBondForce(const SystemTopology& systemTopology)
    -> OpenMM::HarmonicBondForce* {
    auto* force = new OpenMM::HarmonicBondForce();
    for (int index = 0; index < systemTopology.numBonds; ++index) {
        const auto particle1 = systemTopology.bondsI[index];
        const auto particle2 = systemTopology.bondsJ[index];
        const auto length = systemTopology.bondsEquilibrium[index];
        const auto stiffness = systemTopology.bondsStiffness[index] * 2.0;
        force->addBond(particle1, particle2, length, stiffness);
    }
    return force;
}

auto OpenMMContext::createHarmonicAngleForce(const SystemTopology& systemTopology)
    -> OpenMM::HarmonicAngleForce* {
    auto* force = new OpenMM::HarmonicAngleForce();
    for (int index = 0; index < systemTopology.numAngles; ++index) {
        const auto particle1 = systemTopology.anglesI[index];
        const auto particle2 = systemTopology.anglesJ[index];
        const auto particle3 = systemTopology.anglesK[index];
        const auto angleInRad = systemTopology.anglesEquilibrium[index];
        const auto stiffness = systemTopology.anglesStiffness[index] * 2.0;
        force->addAngle(particle1, particle2, particle3, angleInRad, stiffness);
    }
    return force;
}

auto OpenMMContext::createPeriodicTorsionForce(const SystemTopology& systemTopology)
    -> OpenMM::PeriodicTorsionForce* {
    auto* force = new OpenMM::PeriodicTorsionForce();
    for (int index = 0; index < systemTopology.numPeriodicTorsions; ++index) {
        force->addTorsion(systemTopology.periodicTorsionsI[index],
                          systemTopology.periodicTorsionsJ[index],
                          systemTopology.periodicTorsionsK[index],
                          systemTopology.periodicTorsionsL[index],
                          systemTopology.periodicTorsionsN[index],
                          systemTopology.periodicTorsionsPhase[index],
                          systemTopology.periodicTorsionsStiffness[index]);
    }
    return force;
}

auto OpenMMContext::createImproperHarmonicTorsionForce(const SystemTopology& systemTopology)
    -> OpenMM::CustomTorsionForce* {
    std::ostringstream expr;
    expr << std::setprecision(17) << "k*min(dtheta, 2*" << kPi << "-dtheta)^2; dtheta=abs(theta-theta0)";
    auto* force = new OpenMM::CustomTorsionForce(expr.str());
    force->addPerTorsionParameter("k");
    force->addPerTorsionParameter("theta0");
    for (int index = 0; index < systemTopology.numHarmonicTorsions; ++index) {
        const std::vector<double> params = {systemTopology.harmonicTorsionsStiffness[index],
                                            systemTopology.harmonicTorsionsPhase[index]};
        force->addTorsion(systemTopology.harmonicTorsionsI[index],
                          systemTopology.harmonicTorsionsJ[index],
                          systemTopology.harmonicTorsionsK[index],
                          systemTopology.harmonicTorsionsL[index],
                          params);
    }
    return force;
}

auto OpenMMContext::createCMAPTorsionForce(const SystemTopology& systemTopology)
    -> OpenMM::CMAPTorsionForce* {
    auto* force = new OpenMM::CMAPTorsionForce();
    const int res = systemTopology.cmapGridSize;
    const int gridPoints = res * res;
    const int numGrids = (res > 0) ? static_cast<int>(systemTopology.cmapGridEnergy.size()) / gridPoints : 0;
    for (int g = 0; g < numGrids; ++g) {
        const int offset = g * gridPoints;
        const std::vector<double> slice(systemTopology.cmapGridEnergy.begin() + offset,
                                        systemTopology.cmapGridEnergy.begin() + offset + gridPoints);
        force->addMap(res, slice);
    }
    const int numTorsions = static_cast<int>(systemTopology.cmapTorsionMapIndex.size());
    for (int index = 0; index < numTorsions; ++index) {
        force->addTorsion(systemTopology.cmapTorsionMapIndex[index],
                          systemTopology.cmapTorsionA1[index],
                          systemTopology.cmapTorsionA2[index],
                          systemTopology.cmapTorsionA3[index],
                          systemTopology.cmapTorsionA4[index],
                          systemTopology.cmapTorsionB1[index],
                          systemTopology.cmapTorsionB2[index],
                          systemTopology.cmapTorsionB3[index],
                          systemTopology.cmapTorsionB4[index]);
    }
    return force;
}

auto OpenMMContext::createUreyBradleyForce(const SystemTopology& systemTopology)
    -> OpenMM::HarmonicBondForce* {
    auto* force = new OpenMM::HarmonicBondForce();
    for (int index = 0; index < systemTopology.numUreyBradley; ++index) {
        force->addBond(systemTopology.ureyBradleyI[index],
                       systemTopology.ureyBradleyK[index],
                       systemTopology.ureyBradleyEquilibrium[index],
                       systemTopology.ureyBradleyStiffness[index] * 2.0);
    }
    return force;
}

auto OpenMMContext::createAlchemyDecouplingForces(const SystemTopology& sys, OpenMM::NonbondedForce* main)
    -> std::pair<OpenMM::CustomNonbondedForce*, OpenMM::CustomNonbondedForce*> {
    constexpr double kSoftcoreAlpha = 0.5; // Beutler soft-core; standard value

    // (1) MAIN (PME) force. Electrostatics: charge(lambda)=lambda_inter*q via a
    // parameter offset (base set to 0, scale = q) -- this is the ONLY route that
    // scales the reciprocal-space sum correctly. It scales A's charge against
    // everything, so intra-A electrostatics are ANNIHILATED for lambda<1 (a
    // documented departure from pure decoupling; exact at lambda=1, hence
    // unbiased -- the protocol is guidance only, acceptance is full H at lambda=1).
    // Sterics: zero A's epsilon so MAIN computes no LJ involving A; rebuilt below.
    // 1-4/exclusion exceptions are intramolecular and left untouched.
    main->addGlobalParameter("lambda_inter", 1.0);
    for (int i = alchemyBegin; i < alchemyEnd; ++i) {
        double q = 0.0, sig = 0.0, eps = 0.0;
        main->getParticleParameters(i, q, sig, eps);
        main->setParticleParameters(i, 0.0, sig, 0.0);
        main->addParticleParameterOffset("lambda_inter", i, q, 0.0, 0.0);
    }

    std::set<int> aSet, restSet;
    for (int i = 0; i < sys.numAtoms; ++i) {
        (i >= alchemyBegin && i < alchemyEnd ? aSet : restSet).insert(i);
    }

    // (2) Soft-core A x rest LJ, scaled by lambda_inter. lambda=1 -> exact LJ;
    // lambda=0 -> 0; finite for all r at lambda<1 (no overlap singularity).
    auto* soft = new OpenMM::CustomNonbondedForce("lambda_inter*4*eps*(1/(d*d) - 1/d);"
                                                  "d = alpha*(1 - lambda_inter) + (r/sig)^6;"
                                                  "eps = sqrt(eps1*eps2); sig = 0.5*(sig1 + sig2)");
    soft->addGlobalParameter("lambda_inter", 1.0);
    soft->addGlobalParameter("alpha", kSoftcoreAlpha);
    soft->addPerParticleParameter("sig");
    soft->addPerParticleParameter("eps");
    for (int i = 0; i < sys.numAtoms; ++i) {
        soft->addParticle({sys.atomsSigma[i], sys.atomsEpsilon[i]});
    }
    soft->addInteractionGroup(aSet, restSet); // A x rest only; no intermolecular exceptions exist

    // (3) Hard intra-A LJ (lambda-independent). Restores the intra-solute LJ that
    // (1) removed from MAIN, so lambda=1 reproduces the unmodified field. Excludes
    // every intra-A pair MAIN carries as an exception (1-2/1-3 and 1-4) so they are
    // not double-counted.
    auto* hard = new OpenMM::CustomNonbondedForce(
        "4*eps*((sig/r)^12 - (sig/r)^6); eps = sqrt(eps1*eps2); sig = 0.5*(sig1 + sig2)");
    hard->addPerParticleParameter("sig");
    hard->addPerParticleParameter("eps");
    for (int i = 0; i < sys.numAtoms; ++i) {
        hard->addParticle({sys.atomsSigma[i], sys.atomsEpsilon[i]});
    }
    hard->addInteractionGroup(aSet, aSet); // A x A only
    const auto bothAlch = [&](int a, int b) {
        return a >= alchemyBegin && a < alchemyEnd && b >= alchemyBegin && b < alchemyEnd;
    };
    for (int k = 0; k < sys.numScaling14; ++k) {
        if (bothAlch(sys.scaling14I[k], sys.scaling14L[k])) {
            hard->addExclusion(sys.scaling14I[k], sys.scaling14L[k]);
        }
    }
    for (int k = 0; k < sys.numExclusions; ++k) {
        if (bothAlch(sys.exclusionI[k], sys.exclusionJ[k])) {
            hard->addExclusion(sys.exclusionI[k], sys.exclusionJ[k]);
        }
    }

    // Match MAIN's LJ treatment (cutoff under any periodic method).
    for (auto* f : {soft, hard}) {
        if (isPeriodic(sys.nonbondedMethod)) {
            f->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        } else if (sys.nonbondedMethod == NonbondedMethod::CutoffNonPeriodic) {
            f->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        } else {
            f->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
        }
        if (f->getNonbondedMethod() != OpenMM::CustomNonbondedForce::NoCutoff) {
            f->setCutoffDistance(sys.nonbondedCutoff);
        }
    }
    return {soft, hard};
}
