#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "OpenMM.h"

/**
 * @brief Multiple-timestep (r-RESPA) Verlet integrator built on OpenMM's
 *        `CustomIntegrator`.
 *
 * The outer step is split into nested substeps by force group: fast groups
 * (larger substep count) are evaluated more often than slow groups, following
 * the reversible reference-system propagator (r-RESPA) recursion. The scheme is
 * time-reversible and symplectic, so it preserves the sampled distribution and
 * is safe as the inner propagator of an HMC acceptance test.
 *
 * @note Built by `OpenMMSystemBuilder::build` only when MTS is enabled
 *       (`OpenMMContext::setMTS`, reachable from Python via
 *       `Context::setMTS`); it is disabled by default and has no dedicated test.
 */
class MTSIntegrator : public OpenMM::CustomIntegrator {
    public:
    /**
     * @brief Constructs the nested r-RESPA substep program.
     *
     * @param[in] stepSize  Outer step size [ps].
     * @param[in] groups    `(force group, substeps)` pairs; sorted internally by
     *                      substep count so slower groups (fewer substeps) form
     *                      the outer tiers and faster groups the inner. Each
     *                      group's substep count must be an integer multiple of
     *                      the enclosing group's.
     * @throws std::invalid_argument  if @p groups is empty, any substep count is
     *         not a multiple of its parent's, or a force group is outside `[0, 31]`.
     */
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
