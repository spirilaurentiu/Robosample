#include "bridge/AlchemyForceFactory.hpp"

#include "OpenMMContext.hpp"       // OpenMMContext::isPeriodic
#include "bridge/ForceFactory.hpp" // robo::forcefactory::addStandardExclusions

// ONE_4PI_EPS0 (Coulomb constant, kJ*nm/(mol*e^2)): reachable only through the
// reference-platform utility header, matching the same unconditional include
// OpenMMContext.cpp itself carries for the identical reason.
#include "../../openmm/platforms/reference/include/ReferencePlatform.h"

#include <algorithm>
#include <vector>

void AlchemyForceFactory::enableAlchemy(const std::vector<int>& atomIndices) {
    alchemyEnabled_ = true;
    alchemyAtoms_ = atomIndices;
    std::sort(alchemyAtoms_.begin(), alchemyAtoms_.end());
    alchemyAtoms_.erase(std::unique(alchemyAtoms_.begin(), alchemyAtoms_.end()), alchemyAtoms_.end());
    alchemyAtomSet_ = std::set<int>(alchemyAtoms_.begin(), alchemyAtoms_.end());
}

void AlchemyForceFactory::enableAlchemy(int atomBegin, int atomEnd) {
    std::vector<int> atomIndices;
    for (int i = atomBegin; i < atomEnd; ++i) {
        atomIndices.push_back(i);
    }
    enableAlchemy(atomIndices);
}

auto AlchemyForceFactory::createAlchemyCorrectionForce(const SystemTopology& sys)
    -> OpenMM::CustomNonbondedForce* {
    // Total [begin,end) x rest pair energy = standard + (lambda_inter-1)*standard
    //                                      = lambda_inter * standard.
    // Lorentz-Berthelot combining, matching OpenMM NonbondedForce defaults. All
    // 1-4/exclusion pairs are intramolecular, so the A x rest interaction group
    // carries no exceptions and needs no exclusion list (scales to assemblies).
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
        (alchemyAtomSet_.count(i) ? aSet : restSet).insert(i);
    }
    f->addInteractionGroup(aSet, restSet); // A x rest ONLY
    // Share the main NonbondedForce's exclusion list so the CPU platform accepts
    // the Context (see robo::forcefactory::addStandardExclusions). Energy-neutral:
    // no excluded pair is an A x rest pair.
    robo::forcefactory::addStandardExclusions(f, sys);
    return f;
}

auto AlchemyForceFactory::createAlchemyDecouplingForces(const SystemTopology& sys, OpenMM::NonbondedForce* main)
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
    for (int i : alchemyAtoms_) {
        double q = 0.0, sig = 0.0, eps = 0.0;
        main->getParticleParameters(i, q, sig, eps);
        main->setParticleParameters(i, 0.0, sig, 0.0);
        main->addParticleParameterOffset("lambda_inter", i, q, 0.0, 0.0);
    }

    std::set<int> aSet, restSet;
    for (int i = 0; i < sys.numAtoms; ++i) {
        (alchemyAtomSet_.count(i) ? aSet : restSet).insert(i);
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
    // Both custom forces share the main NonbondedForce's full exclusion list so
    // the CPU platform accepts the Context (see robo::forcefactory::addStandardExclusions).
    // This is energy-neutral on both platforms: for `hard` (A x A) the only excluded
    // pairs that fall inside the group are the intra-A 1-2/1-3/1-4 exceptions --
    // exactly the pairs MAIN carries as exceptions and that must not be
    // double-counted -- while the rest-involving exclusions are never A x A; for
    // `soft` (A x rest) no excluded (intramolecular) pair is ever an A x rest pair.
    robo::forcefactory::addStandardExclusions(soft, sys);
    robo::forcefactory::addStandardExclusions(hard, sys);

    // Match MAIN's LJ treatment (cutoff under any periodic method).
    for (auto* f : {soft, hard}) {
        if (OpenMMContext::isPeriodic(sys.nonbondedMethod)) {
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
