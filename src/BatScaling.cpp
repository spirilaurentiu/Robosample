#include "BatScaling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace robo {

namespace {

// Rodrigues' rotation formula: rotate vector v about the UNIT axis by angle
// (radians), right-hand rule.
Vec3 rotateAboutAxis(const Vec3& v, const Vec3& axisUnit, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return (v * c) + ((axisUnit % v) * s) + (axisUnit * (dot(axisUnit, v) * (1.0 - c)));
}

// Collect every atom belonging to `rootBody` and all of its descendant
// bodies (the RIGID subtree a scaled body's own bend/stretch moves).
void collectSubtreeAtoms(const RobotModel& model, int rootBody, std::vector<int>& out) {
    std::vector<int> stack{rootBody};
    while (!stack.empty()) {
        const int b = stack.back();
        stack.pop_back();
        for (int a = model.bodyAtomsBeg[b]; a < model.bodyAtomsEnd[b]; ++a) {
            out.push_back(model.bodyAtoms[a]);
        }
        for (int c = model.bodyChildrenBeg[b]; c < model.bodyChildrenEnd[b]; ++c) {
            stack.push_back(model.bodyChildren[c]);
        }
    }
}

double lookupOr(const std::unordered_map<int, double>& m, int key, double fallback) {
    const auto it = m.find(key);
    return (it != m.end()) ? it->second : fallback;
}

// One body's D5-scaled BAT coordinate value(s), read geometrically from a
// GIVEN atom-position snapshot via the z-matrix (RobotModel::zI/zJ/zK).
// Shared by calcBatVolumeLogJac and BatAnchorStats::accumulate so both read
// EXACTLY the same (r, theta) definition applyBatScaling's forward map uses.
struct BatCoordSample {
    int zI = -1;
    bool hasR = false;
    double r = 0.0;
    bool hasTheta = false;
    double theta = 0.0;
};

BatCoordSample readBatCoord(const RobotModel& model, const std::vector<Vec3>& atomPos, int b) {
    BatCoordSample out;
    const auto sel = batScaledDofs(model.bodyJoint[b]);
    if (!sel.r && !sel.theta) {
        return out;
    }
    const int zrow = model.bodyZRow[b];
    if (zrow < 0) {
        return out;
    }
    const int zI = model.zI[zrow];
    const int zJ = model.zJ[zrow];
    const int zK = model.zK[zrow];
    if (zJ < 0) {
        return out; // no parent atom -- cannot happen for a legal scaled joint
    }
    out.zI = zI;
    const Vec3 v1 = atomPos[zI] - atomPos[zJ];
    const double r = v1.norm();
    if (!(r > 0.0) || std::isnan(r)) {
        return out; // degenerate zero-length bond
    }
    if (sel.r) {
        out.hasR = true;
        out.r = r;
    }
    if (sel.theta && zK >= 0) {
        const Vec3 v2 = atomPos[zK] - atomPos[zJ];
        const double r2 = v2.norm();
        if (r2 > 0.0 && !std::isnan(r2)) {
            double cosTheta = dot(v1, v2) / (r * r2);
            cosTheta = std::min(1.0, std::max(-1.0, cosTheta));
            const double theta = std::acos(cosTheta);
            if (!std::isnan(theta)) {
                out.hasTheta = true;
                out.theta = theta;
            }
        }
    }
    return out;
}

// Result of one body's "bend" step: the (possibly rotated) bond direction
// and whether the angle DOF was actually displaced (false for the
// degenerate zK<0 / collinear-axis skips, self-consistent with
// outNScaled/outLnJac not counting them).
struct BendResult {
    Vec3 bondDir;
    bool displaced = false;
};

// The "bend" step of applyBatScaling's D7 map: rigid rotation of the scaled
// body's subtree about the axis perpendicular to the (zI,zJ,zK) plane,
// pivoting at zJ, by the angle DELTA between the scaled theta1 and theta0.
// Callers gate this on sel.theta (undefined for Slider/Cylinder, which have
// no angle DOF). Mutates atomPos for every atom in `subtree`.
//
// FP note (INV-7, dedup scope): this reads theta INLINE, not via
// readBatCoord -- readBatCoord's cosTheta = dot(v1,v2)/(r*r2) and this
// cosTheta0 = dot(bondDir, ez) on pre-normalized vectors agree
// mathematically but round differently, so routing this through
// readBatCoord would perturb dtheta and the rotated positions (a B3
// bitwise violation). See the ticket's FP note for the full argument.
BendResult applyBend(std::vector<Vec3>& atomPos,
                     const Vec3& posJ,
                     int zK,
                     int zI,
                     double s,
                     const std::unordered_map<int, double>& anchorTheta,
                     const std::vector<int>& subtree,
                     Vec3 bondDir) {
    BendResult result{bondDir, false};
    if (zK < 0) {
        return result;
    }
    const Vec3 v2 = atomPos[zK] - posJ;
    const double n2 = v2.norm();
    if (!(n2 > 1e-9) || std::isnan(n2)) {
        return result;
    }
    const Vec3 ez = v2 / n2;
    double cosTheta0 = dot(bondDir, ez);
    cosTheta0 = std::min(1.0, std::max(-1.0, cosTheta0));
    const double theta0 = std::acos(cosTheta0);
    Vec3 axis = ez % bondDir; // perpendicular to the (zI,zJ,zK) plane
    const double axisNorm = axis.norm();
    if (!(axisNorm > 1e-9) || std::isnan(theta0)) {
        // (zI-zJ) parallel/antiparallel to (zK-zJ) -- the bend
        // axis is undefined at this exact geometry; skip the angle
        // DOF for this body only (self-consistent: not counted in
        // outNScaled either, since the map did not move it).
        return result;
    }
    axis = axis / axisNorm;
    const double muTheta = lookupOr(anchorTheta, zI, 0.0);
    const double theta1 = (s * theta0) - ((s - 1.0) * muTheta);
    // A bond angle's physical domain is (0, pi); std::sin(theta1)
    // <= 0 outside it would make calcBatVolumeLogJac's ln(sin
    // theta1) term undefined. FAIL LOUD here (Rule 11) instead
    // of the silent sinTheta<=0 skip calcBatVolumeLogJac used to
    // apply on its own -- that guard, by itself, made the map
    // MOVE the atom while the Jacobian silently dropped the
    // term, corrupting lnJac (caught by this file's V5 FD gate:
    // an aggressive (s, anchor) pair pushed theta1 past pi).
    // Callers (Stage 2b) SHALL treat this as an invalid/
    // automatic-reject drive, not retry with a different s.
    if (!(theta1 > 0.0 && theta1 < Pi) || std::isnan(theta1)) {
        throw std::domain_error(
            "applyBatScaling: scaled angle theta1=" + std::to_string(theta1)
            + " is outside the physical domain (0, pi) for body zI=" + std::to_string(zI)
            + " (theta0=" + std::to_string(theta0) + ", s=" + std::to_string(s)
            + ", mu_theta=" + std::to_string(muTheta) + ")");
    }
    const double dtheta = theta1 - theta0;
    for (int a : subtree) {
        const Vec3 rel = atomPos[a] - posJ;
        atomPos[a] = posJ + rotateAboutAxis(rel, axis, dtheta);
    }
    result.bondDir = rotateAboutAxis(bondDir, axis, dtheta);
    result.displaced = true;
    return result;
}

// The "stretch" step of applyBatScaling's D7 map: rigid translation of the
// scaled body's subtree along the (possibly bend-rotated) bond direction by
// the DELTA in r. Every scaled joint type has this DOF (D5). Mutates
// atomPos for every atom in `subtree`. Returns whether the DOF was
// displaced (always true on success; a domain violation throws instead of
// returning false, matching applyBatScaling's fail-loud contract).
//
// `r0` is the STRETCH radial read, bitwise-identical to
// readBatCoord(...).r (same atomPos[zI]-atomPos[zJ] subtraction, same
// .norm()) -- callers SHALL source it from readBatCoord so the map and the
// Cartesian log-Jacobian (calcBatVolumeLogJac -> readBatCoord) read the
// same r (INV-7).
bool applyStretch(std::vector<Vec3>& atomPos,
                  double r0,
                  int zI,
                  double s,
                  const std::unordered_map<int, double>& anchorR,
                  const std::vector<int>& subtree,
                  const Vec3& bondDir) {
    const double muR = lookupOr(anchorR, zI, 0.0);
    const double r1 = (s * r0) - ((s - 1.0) * muR);
    // Bond length's physical domain is r1 > 0 (calcBatVolumeLogJac's
    // ln(r1) term is undefined otherwise) -- fail loud (Rule 11), same
    // rationale as the theta1 domain guard above.
    if (!(r1 > 0.0) || std::isnan(r1)) {
        throw std::domain_error(
            "applyBatScaling: scaled bond length r1=" + std::to_string(r1)
            + " is outside the physical domain (r1 > 0) for body zI=" + std::to_string(zI)
            + " (r0=" + std::to_string(r0) + ", s=" + std::to_string(s)
            + ", mu_r=" + std::to_string(muR) + ")");
    }
    const double dr = r1 - r0;
    const Vec3 shift = bondDir * dr;
    for (int a : subtree) {
        atomPos[a] += shift;
    }
    return true;
}

} // namespace

int countScaledDofsStatic(const RobotModel& model) {
    int n = 0;
    for (int b = 1; b < model.numBodies; ++b) {
        const auto sel = batScaledDofs(model.bodyJoint[b]);
        n += (sel.r ? 1 : 0) + (sel.theta ? 1 : 0);
    }
    return n;
}

double calcBatVolumeLogJac(const RobotModel& model, const std::vector<Vec3>& atomPos) {
    double J = 0.0;
    for (int b = 1; b < model.numBodies; ++b) {
        const BatCoordSample sample = readBatCoord(model, atomPos, b);
        if (sample.hasR) {
            J += 2.0 * std::log(sample.r);
        }
        if (sample.hasTheta) {
            const double sinTheta = std::sin(sample.theta);
            if (sinTheta > 0.0 && !std::isnan(sinTheta)) {
                J += std::log(sinTheta);
            }
        }
    }
    return J;
}

std::vector<Vec3> applyBatScaling(const RobotModel& model,
                                  const std::vector<Vec3>& atomPosIn,
                                  double s,
                                  const std::unordered_map<int, double>& anchorR,
                                  const std::unordered_map<int, double>& anchorTheta,
                                  int& outNScaled,
                                  double& outLnJac) {
    std::vector<Vec3> atomPos = atomPosIn;
    const double J0 = calcBatVolumeLogJac(model, atomPos);
    int nScaled = 0;

    for (int b = 1; b < model.numBodies; ++b) {
        const auto sel = batScaledDofs(model.bodyJoint[b]);
        if (!sel.r && !sel.theta) {
            continue;
        }
        const int zrow = model.bodyZRow[b];
        if (zrow < 0) {
            continue;
        }
        const int zI = model.zI[zrow];
        const int zJ = model.zJ[zrow];
        const int zK = model.zK[zrow];
        if (zJ < 0) {
            continue; // defensive: a legal scaled joint always has a parent atom
        }

        const Vec3 posJ = atomPos[zJ];
        Vec3 bondVec = atomPos[zI] - posJ;
        // Bitwise-identical to readBatCoord(...).r (same atomPos[zI]-
        // atomPos[zJ] subtraction, same .norm()) -- sourcing it from
        // readBatCoord makes the map and the Cartesian log-Jacobian
        // (calcBatVolumeLogJac -> readBatCoord) read the same r (INV-7).
        // This does NOT replace applyBatScaling's own domain guard below
        // with readBatCoord's r > 0.0 guard -- that guard stays here.
        const double r0 = readBatCoord(model, atomPos, b).r;
        if (!(r0 > 1e-9) || std::isnan(r0)) {
            continue; // degenerate zero-length bond (flagged edge case)
        }
        Vec3 bondDir = bondVec / r0;

        std::vector<int> subtree;
        collectSubtreeAtoms(model, b, subtree);

        // ---- bend (theta), if this joint scales an angle -------------------
        if (sel.theta) {
            const BendResult bend = applyBend(atomPos, posJ, zK, zI, s, anchorTheta, subtree, bondDir);
            bondDir = bend.bondDir;
            if (bend.displaced) {
                ++nScaled;
            }
        }

        // ---- stretch (r): every scaled joint type has this ------------------
        if (sel.r) {
            if (applyStretch(atomPos, r0, zI, s, anchorR, subtree, bondDir)) {
                ++nScaled;
            }
        }
    }

    outNScaled = nScaled;
    const double J1 = calcBatVolumeLogJac(model, atomPos);
    outLnJac = (J1 - J0) + (static_cast<double>(nScaled) * std::log(s));
    return atomPos;
}

// ----------------------------------------------------------------------------
//  BatAnchorStats (INV-9 shared/global anchor)
// ----------------------------------------------------------------------------
void BatAnchorStats::accumulate(const RobotModel& model, const std::vector<Vec3>& atomPos) {
    for (int b = 1; b < model.numBodies; ++b) {
        const BatCoordSample sample = readBatCoord(model, atomPos, b);
        if (sample.hasR) {
            RunningMean& rm = meanR_[sample.zI];
            rm.n += 1;
            rm.mean += (sample.r - rm.mean) / static_cast<double>(rm.n);
        }
        if (sample.hasTheta) {
            RunningMean& tm = meanTheta_[sample.zI];
            tm.n += 1;
            tm.mean += (sample.theta - tm.mean) / static_cast<double>(tm.n);
        }
    }
}

BatAnchorStats::Snapshot BatAnchorStats::snapshot() const {
    Snapshot out;
    for (const auto& kv : meanR_) {
        out.meanR[kv.first] = kv.second.mean;
    }
    for (const auto& kv : meanTheta_) {
        out.meanTheta[kv.first] = kv.second.mean;
    }
    return out;
}

void BatAnchorStats::reset() {
    meanR_.clear();
    meanTheta_.clear();
}

} // namespace robo
