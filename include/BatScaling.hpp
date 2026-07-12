#pragma once

// ============================================================================
//  BatScaling -- the deterministic BAT (bond/angle) scaling drive and its
//  Cartesian log-Jacobian (docs/specs/replica-exchange-nonequilibrium-work.md
//  B4, D5, D6, D7). STAGE 2a: the drive + Jacobian ONLY -- no acceptance, no
//  WORK accumulation, no run-type dispatch (Stage 2b).
//
//  D5 selection: a scaled DOF exists only for Slider/Cylinder (bond length r)
//  and BendStretch/SphericalCoords (bond length r AND bond angle theta);
//  Torsion/Rigid/root joints (Cartesian/Ball/FreeLine/Free) contribute
//  nothing. r = dist(zI,zJ), theta = angle(zI,zJ,zK), read from the
//  RobotModel z-matrix (World::buildModel populates zI/zJ/zK/zL/bodyZRow).
//
//  D6 Jacobian: lnJac = (J(x') - J(x^0)) + N_scaled*ln(s), where
//  J(x) = sum over scaled bodies [2 ln r + ln sin(theta)] (only the
//  coordinate(s) that body's joint type actually exposes -- Slider/Cylinder
//  contribute "2 ln r" only, no angle term). This is the CORRECTED
//  composition (REFUTES the original's `J_ini + J_scale - J_fin`, F7): see
//  D6 for the full derivation and the numeric counter-example.
//
//  D7 / the deterministic map: each scaled body's subtree (itself + every
//  descendant body) is moved as ONE RIGID ASSEMBLY so its own (r,theta)
//  relative to its PARENT become (s*r - (s-1)*mu_r, s*theta - (s-1)*mu_theta)
//  -- first a "bend" (rigid rotation about the axis perpendicular to the
//  (zI,zJ,zK) plane, pivoting at zJ, by the angle DELTA -- undefined for
//  Slider/Cylinder, which have no angle DOF), then a "stretch" (rigid
//  translation along the resulting bond direction by the DELTA in r). Because
//  internal (z-matrix-relative) BAT coordinates are invariant under a rigid
//  motion of their own ancestor chain, processing bodies in TOPOLOGICAL order
//  (parent before child, the RobotModel body-index convention) makes each
//  body's own scaling read the correct (already-cascaded) input geometry --
//  this is the property the >=3-body-chain V5 case exists to catch.
// ============================================================================

#include <unordered_map>
#include <vector>

#include "RobotModel.hpp"
#include "robot_math.hpp"

namespace robo {

// Which BAT coordinate(s) a JointType's mobilizer exposes AND the D5 rule
// scales. Torsion/azimuth and root DOFs are never scaled (spec D5 table).
struct ScaledDofFlags {
    bool r = false;
    bool theta = false;
};

[[nodiscard]] constexpr ScaledDofFlags batScaledDofs(JointType jt) {
    switch (jt) {
        case JointType::Slider:
            return {true, false}; // r only
        case JointType::Cylinder:
            return {true, false}; // r only (rotation DOF is TORSION, excluded)
        case JointType::BendStretch:
            return {true, true}; // theta, r
        case JointType::SphericalCoords:
            return {true, true}; // zenith theta, radius r (azimuth excluded)
        default:
            return {false, false};
    }
}

// N_scaled (D5/RQ-2): sum of scaled bond+angle DOF contributions over every
// body `model` exposes (Slider/Cylinder=1, BendStretch/SphericalCoords=2,
// everything else 0). This is the STATIC (topology-only) count; the count of
// DOFs actually displaced by applyBatScaling can be lower in a degenerate
// geometry (see applyBatScaling's outNScaled) -- lnJac SHALL use the latter
// (the count of what the map actually did), which is what
// applyBatScaling/getDistortJacobianDetLog use, not this static count.
[[nodiscard]] int countScaledDofsStatic(const RobotModel& model);

// J(x) = sum over scaled bodies [2 ln r + ln sin(theta)] (D6), evaluated on a
// GIVEN Cartesian snapshot `atomPos` (not necessarily the model's own live
// RobotState -- callers pass pre-/post-scale copies). Reads r/theta from the
// z-matrix (RobotModel::zI/zJ/zK, World::buildModel). A body whose z-row is
// missing (-1, e.g. Rigid) or whose zJ/zK ancestor is absent (root-adjacent,
// flagged limitation -- see applyBatScaling) contributes 0.
[[nodiscard]] double calcBatVolumeLogJac(const RobotModel& model, const std::vector<Vec3>& atomPos);

// The deterministic BAT-scaling drive (B4, D7) + its D6 Jacobian, computed
// together so outLnJac is guaranteed consistent with what the map actually
// applied (outNScaled). Pure function: `atomPosIn` is not mutated; the
// scaled copy is returned. `anchorR`/`anchorTheta` key the (D2/INV-9) shared
// scaling anchor by the scaled body's zI (placed) atom index -- a missing
// key defaults to 0.0 (an anchor-independent identity for the Jacobian; see
// the header's INV-9 note below).
//
// INV-9 (D2 revision 2, corrected per reviewer S1): the affine map
// q' = s*q - (s-1)*mu has LOCAL Jacobian dq'/dq = s REGARDLESS of mu (mu is
// a pure translation of the pivot; det dS = s^N is anchor-independent). But
// `outLnJac` returned here is NOT that local slope alone -- it is
// (J(x')-J(x0)) + N_scaled*ln(s), evaluated at the ACTUAL produced geometry
// x' = f(x0, s, mu). Since r1 = s*r0 - (s-1)*mu_r (and similarly for theta1)
// depends on mu, x' itself depends on mu, so J(x') -- and hence the numeric
// value of `outLnJac` this call returns -- DOES depend on mu. What is
// anchor-independent is only the map's local determinant s^N, not this
// function's output for a GIVEN input. The involution property
// M_{1/s,mu} o M_{s,mu} = id (and therefore correctionTerm == 1, D2) holds
// ONLY when both legs share the SAME mu -- callers implementing INV-9
// (Context's shared, state-independent running-mean anchor) SHALL pass the
// SAME frozen anchorR/anchorTheta maps to both partners' drives in one round.
//
// FLAGGED LIMITATION (RobotModel/BAT-API ambiguity, coder checkpoint): a
// scaled body whose PARENT is itself Ground-adjacent (bodyParent[bodyParent
// [b]] == Ground) has no zK -- the angle DOF is silently skipped for that
// one body (its `theta` contribution is omitted from N_scaled/lnJac, kept
// self-consistent: the map does not move it either). This is a real, valid
// molecular configuration (a BendStretch/SphericalCoords body whose parent is
// itself a Free/Torsion-rooted body) that the current engine's z-matrix has
// no fallback reference atom for; production configs SHOULD avoid stacking a
// scaled joint directly on a root-adjacent body until a fallback (e.g. a
// sibling atom) is added.
std::vector<Vec3> applyBatScaling(const RobotModel& model,
                                  const std::vector<Vec3>& atomPosIn,
                                  double s,
                                  const std::unordered_map<int, double>& anchorR,
                                  const std::unordered_map<int, double>& anchorTheta,
                                  int& outNScaled,
                                  double& outLnJac);

// ----------------------------------------------------------------------------
//  Shared/global scaling anchor (INV-9, D2 revision 2)
// ----------------------------------------------------------------------------
// The affine scaling anchor `mu` SHALL be a single, STATE-INDEPENDENT value
// shared by both partners of a swap pair and frozen across their forward and
// reverse drive (INV-9): only then is the paired map an exact involution
// (M_{1/s,mu} o M_{s,mu} = id). BatAnchorStats is the Context-owned (NOT
// per-ThermodynamicState) accumulator this requires: a single running mean
// per scaled BAT coordinate (keyed by the scaled body's zI atom index, same
// key applyBatScaling reads), updated from equilibrium samples across ALL
// thermodynamic states, and frozen into an immutable Snapshot once per round
// -- BOTH partners of a swap pair SHALL read the SAME Snapshot object.
class BatAnchorStats {
    public:
    // Welford's online mean update: one sample per scaled body's r (and
    // theta, if that body's joint scales an angle) read from `atomPos`.
    void accumulate(const RobotModel& model, const std::vector<Vec3>& atomPos);

    struct Snapshot {
        std::unordered_map<int, double> meanR;
        std::unordered_map<int, double> meanTheta;
    };
    // Immutable copy of the CURRENT running means. Take exactly ONE snapshot
    // per round and pass it to every drive that round (INV-9).
    [[nodiscard]] Snapshot snapshot() const;

    void reset();

    private:
    struct RunningMean {
        long long n = 0;
        double mean = 0.0;
    };
    std::unordered_map<int, RunningMean> meanR_;
    std::unordered_map<int, RunningMean> meanTheta_;
};

} // namespace robo
