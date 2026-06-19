#include "Constraints.hpp"

namespace robo {

auto ConstraintSet::mapAtomForcesToGeneralizedForces(const RobotModel& model,
                                                     const RobotState& state,
                                                     const Vec3* atomForce,
                                                     Real* tauOut) -> void {
    const Transform* xGroundBody = state.X_GB();
    const Vec3* posGround = state.atomPosG();
    const SpatialVec* matrixH = state.H();

    std::vector<SpatialVec> bias(static_cast<std::size_t>(model.numBodies), SpatialVec(Vec3(0), Vec3(0)));

    for (int atom = 0; atom < model.numAtoms; ++atom) {
        const int body = model.atomBody[atom];
        const Vec3& force = atomForce[atom];
        const Vec3 lever = posGround[atom] - xGroundBody[body].p();
        bias[body].linear += force;
        bias[body].angular += (lever % force);
    }

    // Inward sweep (children have higher ids, so high->low visits them first).
    for (int body = model.numBodies - 1; body >= 1; --body) {
        const int uOff = model.bodyUIndex[body];
        const int dof = model.bodyNU[body];
        for (int col = 0; col < dof; ++col) {
            tauOut[uOff + col] = dot(matrixH[uOff + col], bias[body]);
        }
        const int parent = model.bodyParent[body];
        if (parent >= 1) {
            const Vec3 lever = xGroundBody[body].p() - xGroundBody[parent].p();
            bias[parent].angular += bias[body].angular + (lever % bias[body].linear);
            bias[parent].linear += bias[body].linear;
        }
    }
}

// RATTLE: project u so that G u = 0. Requires realizeVelocity() current.
auto ConstraintSet::enforceVelocityConstraints(const RobotModel& model, RobotState& state) const -> void {
    const int numC = numConstraints();
    if (numC == 0) {
        return;
    }
    const int numU = model.nu;

    const Transform* xGroundBody = state.X_GB();
    const Vec3* posGround = state.atomPosG();
    const SpatialVec* velGroundBody = state.V_GB();
    Real* velU = state.u();

    std::vector<std::vector<Real>> jacobianT(static_cast<std::size_t>(numC), std::vector<Real>(numU, 0));
    std::vector<std::vector<Real>> minvJacobianT(static_cast<std::size_t>(numC), std::vector<Real>(numU, 0));
    std::vector<Real> rhs(static_cast<std::size_t>(numC), 0);
    std::vector<Vec3> atomForce(static_cast<std::size_t>(model.numAtoms), Vec3(0));

    for (int con = 0; con < numC; ++con) {
        const DistanceConstraint& bond = distance[con];
        const Vec3 vecAB = posGround[bond.atomA] - posGround[bond.atomB];

        const int bodyA = model.atomBody[bond.atomA];
        const int bodyB = model.atomBody[bond.atomB];
        const Vec3 stationA = posGround[bond.atomA] - xGroundBody[bodyA].p();
        const Vec3 stationB = posGround[bond.atomB] - xGroundBody[bodyB].p();
        const Vec3 velA = velGroundBody[bodyA].linear + (velGroundBody[bodyA].angular % stationA);
        const Vec3 velB = velGroundBody[bodyB].linear + (velGroundBody[bodyB].angular % stationB);
        rhs[con] = dot(vecAB, velA - velB);

        for (Vec3& force : atomForce) {
            force = Vec3(0);
        }
        atomForce[bond.atomA] = vecAB;
        atomForce[bond.atomB] = Vec3(0) - vecAB;
        mapAtomForcesToGeneralizedForces(model, state, atomForce.data(), jacobianT[con].data());
        RobotEngine::multiplyByMInv(model, state, jacobianT[con].data(), minvJacobianT[con].data());
    }

    std::vector<Real> multiplier = solveCoupling(jacobianT, minvJacobianT, rhs, numC, numU);
    for (int con = 0; con < numC; ++con) {
        for (int idx = 0; idx < numU; ++idx) {
            velU[idx] -= multiplier[con] * minvJacobianT[con][idx];
        }
    }
}

auto ConstraintSet::applyVelSpaceIncrementToQ(const RobotModel& model,
                                              RobotState& state,
                                              const Real* deltaVel) -> void {
    Real* coordQ = state.q();
    for (int body = 1; body < model.numBodies; ++body) {
        const int qOff = model.bodyQIndex[body];
        const int uOff = model.bodyUIndex[body];
        const int dof = model.bodyNU[body];
        if (model.bodyJoint[body] == JointType::Free) {
            const Quat quat(coordQ[qOff], coordQ[qOff + 1], coordQ[qOff + 2], coordQ[qOff + 3]);
            const Vec3 deltaOmega(deltaVel[uOff], deltaVel[uOff + 1], deltaVel[uOff + 2]);
            const Quat deltaQuat = Quat::angVelToQdot(quat, deltaOmega);
            for (int comp = 0; comp < 4; ++comp) {
                coordQ[qOff + comp] += deltaQuat.elems[static_cast<std::size_t>(comp)];
            }
            coordQ[qOff + 4] += deltaVel[uOff + 3];
            coordQ[qOff + 5] += deltaVel[uOff + 4];
            coordQ[qOff + 6] += deltaVel[uOff + 5];
        } else {
            for (int col = 0; col < dof; ++col) {
                coordQ[qOff + col] += deltaVel[uOff + col];
            }
        }
    }
    RobotEngine::normalizeQuaternions(model, state);
}

// ----------------------------------------------------------------------------
//  calcConstraintLogDet : ln det( G M^-1 G^T ) for the loop-closure Fixman term.
//
//  Spiridon & Minh 2017 (JCTC 13:4649), Eqs. 2-3:
//    rho(phi_f) ~ |M_{N_f}|^{1/2} e^{-beta U}   and   U_F = (1/2) RT ln(|M_{N_f}|/|M_3N|).
//  |M_{N_f}| is the mass-metric determinant of the constrained system's FLEXIBLE
//  coordinates. With m loop-closure (RATTLE) constraints, projecting the seeded
//  momenta onto the velocity-constraint surface G M^-1 p = 0 puts a factor
//  det(G M^-1 G^T)^{-1/2} into that marginal, i.e.
//        |M_{N_f}| = |M_tree| / det(G M^-1 G^T),
//  so the Fixman potential acquires a  -(1/2) RT ln det(G M^-1 G^T)  term on top
//  of the tree term (1/2) RT ln|M_tree| that calcLogDetM already supplies.
//
//  This assembles G M^-1 G^T with exactly the same Jacobian convention as
//  enforceVelocityConstraints (so the two are guaranteed consistent), then reads
//  off ln|det| from the same factorization the SHAKE/RATTLE solver already uses.
//  The zero rhs means we do no useful solve -- we only want the determinant.
//
//  Returns 0 when there are no loop closures (acyclic molecule, or both ring
//  atoms welded into one rigid body): then m == 0 and the term vanishes exactly,
//  so this change is a guaranteed no-op on every acyclic system -- including the
//  C4 chain, butane, and alanine-dipeptide cases the original paper validated.
//
//  Requires realizeArticulatedBodyInertias() to be current (multiplyByMInv reads
//  the articulated quantities). calcFixman() guarantees this by calling it first.
// ----------------------------------------------------------------------------
auto ConstraintSet::calcConstraintLogDet(const RobotModel& model, RobotState& state) const -> Real {
    const int numC = numConstraints();
    if (numC == 0) {
        return Real(0); // no loop closures -> no correction (acyclic case)
    }
    const int numU = model.nu;
    const Vec3* posGround = state.atomPosG();

    std::vector<std::vector<Real>> jacobianT(static_cast<std::size_t>(numC), std::vector<Real>(numU, 0));
    std::vector<std::vector<Real>> minvJacobianT(static_cast<std::size_t>(numC), std::vector<Real>(numU, 0));
    std::vector<Vec3> atomForce(static_cast<std::size_t>(model.numAtoms), Vec3(0));

    // Build one row of G (= (1/2) dC/dq, the same convention RATTLE uses) and
    // M^-1 G^T per constraint -- identical assembly to enforceVelocityConstraints.
    for (int con = 0; con < numC; ++con) {
        const DistanceConstraint& bond = distance[con];
        const Vec3 vecAB = posGround[bond.atomA] - posGround[bond.atomB];
        for (Vec3& force : atomForce) {
            force = Vec3(0);
        }
        atomForce[bond.atomA] = vecAB;
        atomForce[bond.atomB] = Vec3(0) - vecAB;
        mapAtomForcesToGeneralizedForces(model, state, atomForce.data(), jacobianT[con].data());
        RobotEngine::multiplyByMInv(model, state, jacobianT[con].data(), minvJacobianT[con].data());
    }

    // Harvest ln|det(G M^-1 G^T)| from the factorization; the solve itself is
    // discarded (zero rhs).
    Real logAbsDet = Real(0);
    std::vector<Real> zeroRhs(static_cast<std::size_t>(numC), Real(0));
    (void)solveCoupling(jacobianT, minvJacobianT, zeroRhs, numC, numU, &logAbsDet);
    return logAbsDet;
}

auto ConstraintSet::solveCoupling(const std::vector<std::vector<Real>>& jacobianT,
                                  const std::vector<std::vector<Real>>& minvJacobianT,
                                  const std::vector<Real>& rhs,
                                  int numC,
                                  int numU,
                                  Real* logAbsDetOut) -> std::vector<Real> {
    std::vector<std::vector<Real>> coupling(static_cast<std::size_t>(numC), std::vector<Real>(numC, 0));
    for (int row = 0; row < numC; ++row) {
        for (int col = 0; col < numC; ++col) {
            Real acc = 0;
            for (int idx = 0; idx < numU; ++idx) {
                acc += jacobianT[row][idx] * minvJacobianT[col][idx];
            }
            coupling[row][col] = acc;
        }
    }
    return solveSmallSpd(coupling, rhs, logAbsDetOut);
}

auto ConstraintSet::solveSmallSpd(std::vector<std::vector<Real>> matA,
                                  std::vector<Real> vecB,
                                  Real* logAbsDetOut) -> std::vector<Real> {
    const int dim = static_cast<int>(vecB.size());
    // ln|det matA| = sum over pivot columns of ln|pivot|. The pivot at step `col`
    // (matA[col][col] after the row swap, before it eliminates the other rows) is
    // the same value the LU factorization would use, so the product of pivots is
    // +/- det. Partial-pivot row swaps only flip the sign; matA = G M^-1 G^T is
    // symmetric positive-definite, so |det| = product of |pivots| and the sign is
    // irrelevant. A degenerate (near-zero) pivot is floored so the log stays finite.
    if (logAbsDetOut != nullptr) {
        *logAbsDetOut = Real(0);
    }
    for (int col = 0; col < dim; ++col) {
        int pivot = col;
        for (int row = col + 1; row < dim; ++row) {
            if (std::abs(matA[row][col]) > std::abs(matA[pivot][col])) {
                pivot = row;
            }
        }
        std::swap(matA[pivot], matA[col]);
        std::swap(vecB[pivot], vecB[col]);
        const Real diag = matA[col][col];
        if (logAbsDetOut != nullptr) {
            *logAbsDetOut += std::log(std::max(std::abs(diag), Real(1e-300)));
        }
        if (std::abs(diag) < Real(1e-300)) {
            continue; // degenerate constraint row
        }
        for (int row = 0; row < dim; ++row) {
            if (row == col) {
                continue;
            }
            const Real factor = matA[row][col] / diag;
            for (int idx = col; idx < dim; ++idx) {
                matA[row][idx] -= factor * matA[col][idx];
            }
            vecB[row] -= factor * vecB[col];
        }
    }
    std::vector<Real> result(static_cast<std::size_t>(dim), 0);
    for (int idx = 0; idx < dim; ++idx) {
        result[idx] = (std::abs(matA[idx][idx]) < Real(1e-300)) ? Real(0) : vecB[idx] / matA[idx][idx];
    }
    return result;
}

} // namespace robo