#include "OpenMMContext.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bridge/OpenMMSystemBuilder.hpp"

#if USE_CUDA
#    include <cstdlib>
#endif

auto OpenMMContext::initialize(const SystemTopology& systemTopology) -> bool {
    OpenMMSystemBuildResult result =
        OpenMMSystemBuilder::build(systemTopology, useMTS, mtsInnerSubsteps, separateForceGroups, alchemyFactory_);

    // These are populated regardless of success (mirrors the pre-split
    // initialize(), which had already written the corresponding members before
    // the OpenMM::Context construction try/catch below could fail).
    numAtoms = result.numAtoms;
    forceGroupLabels = std::move(result.forceGroupLabels);
    hasVirtualSites = result.hasVirtualSites;
    system = std::move(result.system);
    integrator = std::move(result.integrator);
    if (!result.success) {
        return false;
    }
    context = std::move(result.context);
    initialized = true;

    // Opt-in fused CUDA robot-kinematics pipeline. Env is the default control; a
    // Context/World setter (setCudaKinematics) may override it afterward. No-op unless
    // built with USE_CUDA (cudaKinematicsAvailable() is false otherwise).
#if USE_CUDA
    if (const char* env = std::getenv("ROBO_CUDA_KINEMATICS")) {
        cudaKinematicsEnabled_ = (std::string(env) == "1");
        std::cout << "[INFO] ROBO_CUDA_KINEMATICS=" << env << " -> fused CUDA robot kinematics "
                  << (cudaKinematicsEnabled_ ? "ENABLED" : "disabled") << ".\n";
    }
#endif

    std::cout << "[INFO] Initialized OpenMM. Using version " << result.openMMVersion << ".\n";
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

void OpenMMContext::enableAlchemy(const std::vector<int>& atomIndices) {
    alchemyFactory_.enableAlchemy(atomIndices);
}

void OpenMMContext::enableAlchemy(int atomBegin, int atomEnd) {
    alchemyFactory_.enableAlchemy(atomBegin, atomEnd);
}

void OpenMMContext::setAlchemicalLambda(double lambdaInter) {
    if (!alchemyFactory_.enabled()) {
        return;
    }
    ensureInitialized();
    context->setParameter("lambda_inter", lambdaInter);
}
