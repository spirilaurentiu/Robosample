#include "bridge/ForceFactory.hpp"

#include "OpenMMContext.hpp" // OpenMMContext::isPeriodic (the shared NonbondedMethod predicate)

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace robo::forcefactory {

auto createNonbondedForce(const SystemTopology& systemTopology) -> OpenMM::NonbondedForce* {
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

auto createGBSAOBCForce(const SystemTopology& systemTopology) -> OpenMM::GBSAOBCForce* {
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
    // Replicate native OpenMM's default implicit solvent. `AmberPrmtopFile
    // .createSystem(implicitSolvent=OBC2)` with no salt builds this same built-in
    // GBSAOBCForce and, because its default sasaMethod is 'ACE', leaves the ACE
    // nonpolar surface-area term on (surfaceAreaEnergy = 2.25936 kJ/mol/nm^2). Set
    // it explicitly (rather than relying on OpenMM::GBSAOBCForce's own constructor
    // default) so this is not silently dependent on that default staying
    // 2.25936 across vendored-OpenMM versions -- zeroing/losing it would drop the
    // ~15 kJ/mol ACE term and diverge from the native OpenMM reference.
    force->setSurfaceAreaEnergy(2.25936);
    force->setNonbondedMethod(gbsaForceMethod);
    force->setCutoffDistance(systemTopology.nonbondedCutoff);
    for (int index = 0; index < systemTopology.numAtoms; ++index) {
        force->addParticle(systemTopology.atomsCharge[index],
                           systemTopology.atomsRadius[index],
                           systemTopology.atomsScreen[index]);
    }
    return force;
}

auto createCustomNonbondedForce(const SystemTopology& sys) -> OpenMM::CustomNonbondedForce* {
    // P1: the NBFIX A/B-coefficient table must be well-formed before it is fed to
    // OpenMM::Discrete2DFunction (which does no bounds checking of its own).
    if (sys.numNBTypes <= 0) {
        throw std::runtime_error("createCustomNonbondedForce (NBFIX): numNBTypes must be > 0");
    }
    const std::size_t expectedSize =
        static_cast<std::size_t>(sys.numNBTypes) * static_cast<std::size_t>(sys.numNBTypes);
    if (sys.aCoef.size() != expectedSize || sys.bCoef.size() != expectedSize) {
        throw std::runtime_error("createCustomNonbondedForce (NBFIX): aCoef/bCoef size must equal "
                                 "numNBTypes^2");
    }
    if (static_cast<int>(sys.atomsNonbondedIndex.size()) != sys.numAtoms) {
        throw std::runtime_error("createCustomNonbondedForce (NBFIX): atomsNonbondedIndex.size() must "
                                 "equal numAtoms");
    }
    for (int index = 0; index < sys.numAtoms; ++index) {
        const int type = sys.atomsNonbondedIndex[index];
        if (type < 0 || type >= sys.numNBTypes) {
            throw std::runtime_error("createCustomNonbondedForce (NBFIX): atomsNonbondedIndex["
                                     + std::to_string(index) + "] out of range [0, numNBTypes)");
        }
    }
    // P2: 12-6-4 (LENNARD_JONES_CCOEF) is not represented in SystemTopology; there
    // is no field to check, so reaching here with such data would silently drop the
    // C/r^4 term. Nothing to guard against today (gfcd has none) -- see spec OQ1.

    auto* force = new OpenMM::CustomNonbondedForce(
        "(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2);");
    force->addTabulatedFunction("acoef",
                                new OpenMM::Discrete2DFunction(sys.numNBTypes, sys.numNBTypes, sys.aCoef));
    force->addTabulatedFunction("bcoef",
                                new OpenMM::Discrete2DFunction(sys.numNBTypes, sys.numNBTypes, sys.bCoef));
    force->addPerParticleParameter("type");
    for (int index = 0; index < sys.numAtoms; ++index) {
        force->addParticle({double(sys.atomsNonbondedIndex[index])});
    }

    // Method/cutoff branch mirroring createNonbondedForce (CustomNonbondedForce has
    // no Ewald/PME reciprocal-space option of its own, so every periodic method
    // collapses to CutoffPeriodic here -- same pattern as createAlchemyDecouplingForces).
    if (OpenMMContext::isPeriodic(sys.nonbondedMethod)) {
        force->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        force->setCutoffDistance(sys.nonbondedCutoff);
        force->setUseLongRangeCorrection(true);
    } else if (sys.nonbondedMethod == NonbondedMethod::CutoffNonPeriodic) {
        force->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        force->setCutoffDistance(sys.nonbondedCutoff);
    } else {
        force->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
    }

    // I2: exclude every pair the main NonbondedForce already accounts for (1-4
    // exceptions with CHAMBER-specific sigma/epsilon, plus zeroed 1-2/1-3
    // exclusions) so LJ is never double-counted between the two forces.
    addStandardExclusions(force, sys);
    return force;
}

void addStandardExclusions(OpenMM::CustomNonbondedForce* force, const SystemTopology& sys) {
    // scaling14 (1-4) and exclusion (1-2/1-3) pairs are disjoint -- the main
    // NonbondedForce adds both as exceptions without duplicate-key errors, so the
    // same two loops here never double-add a pair.
    for (int k = 0; k < sys.numScaling14; ++k) {
        force->addExclusion(sys.scaling14I[k], sys.scaling14L[k]);
    }
    for (int k = 0; k < sys.numExclusions; ++k) {
        force->addExclusion(sys.exclusionI[k], sys.exclusionJ[k]);
    }
}

auto createHarmonicBondForce(const SystemTopology& systemTopology) -> OpenMM::HarmonicBondForce* {
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

auto createHarmonicAngleForce(const SystemTopology& systemTopology) -> OpenMM::HarmonicAngleForce* {
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

auto createPeriodicTorsionForce(const SystemTopology& systemTopology) -> OpenMM::PeriodicTorsionForce* {
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

auto createImproperHarmonicTorsionForce(const SystemTopology& systemTopology) -> OpenMM::CustomTorsionForce* {
    std::ostringstream expr;
    expr << std::setprecision(17) << "k*min(dtheta, 2*" << M_PI << "-dtheta)^2; dtheta=abs(theta-theta0)";
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

auto createCMAPTorsionForce(const SystemTopology& systemTopology) -> OpenMM::CMAPTorsionForce* {
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

auto createUreyBradleyForce(const SystemTopology& systemTopology) -> OpenMM::HarmonicBondForce* {
    auto* force = new OpenMM::HarmonicBondForce();
    for (int index = 0; index < systemTopology.numUreyBradley; ++index) {
        force->addBond(systemTopology.ureyBradleyI[index],
                       systemTopology.ureyBradleyK[index],
                       systemTopology.ureyBradleyEquilibrium[index],
                       systemTopology.ureyBradleyStiffness[index] * 2.0);
    }
    return force;
}

} // namespace robo::forcefactory
