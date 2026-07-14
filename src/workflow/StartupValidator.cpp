#include "Context.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "OpenMMContext.hpp"
#include "PeriodicBox.hpp"

// ----------------------------------------------------------------------------
//  checkStartupGeometry -- refuse to start (or warn) on an unminimized system.
//  Set the env var ROBO_ALLOW_BAD_START=1 to downgrade the hard error to a
//  warning (e.g. if you intend to run a Cartesian relaxation world first).
// ----------------------------------------------------------------------------
void Context::checkStartupGeometry() {
    const int n = systemTopology.numAtoms;
    if (n <= 0) {
        return;
    }

    // Excluded pairs: not just 1-2 bonds, but the FULL intramolecular
    // non-interacting set -- the OpenMM exclusion list (1-2 and 1-3) and the 1-4
    // (scaling14) pairs. These are bonded geometry, never clashes. Using only
    // bonds is wrong for 4-point water: the extra point sits ~0.078 nm from each
    // hydrogen (a 1-3 pair that OpenMM excludes), so a bonds-only filter reports
    // two phantom "clashes" per water (2*Nwater of them). The clash scan is meant
    // to catch INTERMOLECULAR overlaps from a bad/unminimized placement, not the
    // internal geometry of a virtual-site water. Keyed as min*N+max in a hash set.
    std::unordered_set<long long> excluded;
    auto pairKey = [n](int i, int j) -> long long {
        if (i > j) {
            std::swap(i, j);
        }
        return static_cast<long long>(i) * n + j;
    };
    for (int k = 0; k < systemTopology.numBonds; ++k) {
        excluded.insert(pairKey(systemTopology.bondsI[k], systemTopology.bondsJ[k]));
    }
    for (int k = 0; k < systemTopology.numExclusions; ++k) {
        excluded.insert(pairKey(systemTopology.exclusionI[k], systemTopology.exclusionJ[k]));
    }
    for (int k = 0; k < systemTopology.numScaling14; ++k) {
        excluded.insert(pairKey(systemTopology.scaling14I[k], systemTopology.scaling14L[k]));
    }

    // Virtual sites / massless particles (e.g. a 4-point water's extra point) have
    // no Lennard-Jones term -- no steric presence -- so they cannot clash and are
    // skipped entirely. (Belt-and-suspenders with the exclusion set above.)
    std::vector<bool> isVirtual(static_cast<std::size_t>(n), false);
    for (int a = 0; a < n; ++a) {
        if (systemTopology.atomsMass[a] == 0.0) {
            isVirtual[static_cast<std::size_t>(a)] = true;
        }
    }

    // Hard-clash distance. Two non-bonded heavy/H atoms closer than this are in
    // r^-12 overlap. 0.08 nm (0.8 A) is well inside any real contact (vdW
    // contacts are >= ~0.2 nm) yet above bonded H distances we already excluded.
    constexpr double kClashNm = 0.08;

    // O(N^2) is fine: this runs once. (For very large systems a grid would help.)
    int nClash = 0;
    double minNonbondedNm = 1e30;
    std::pair<int, int> worst{-1, -1};
    const auto& X = systemTopology.atomsX;
    const auto& Y = systemTopology.atomsY;
    const auto& Z = systemTopology.atomsZ;

    // Under explicit solvent the box is periodic, so a "distance" must be the
    // MINIMUM-IMAGE distance: two atoms on opposite faces are actually neighbours.
    // box_vectors are reduced (lower-triangular) a=(ax,0,0) b=(bx,by,0) c=(cx,cy,cz);
    // wrap the displacement by subtracting whole lattice vectors in c,b,a order
    // (the same order OpenMM reduces them). Without this the scan reports phantom
    // clashes (or misses real cross-boundary ones) on a solvated box.
    const bool periodic =
        OpenMMContext::isPeriodic(systemTopology.nonbondedMethod) && systemTopology.boxVectors.size() == 9;
    const auto& bv = systemTopology.boxVectors;
    auto minImage = [&](double& dx, double& dy, double& dz) {
        if (!periodic) {
            return;
        }
        // single source of truth: robo::pbc::minimumImage (PeriodicBox.hpp).
        robo::pbc::minimumImage(dx, dy, dz, bv.data());
    };

    bool anyNaN = false;
    for (int a = 0; a < n; ++a) {
        if (!std::isfinite(X[a]) || !std::isfinite(Y[a]) || !std::isfinite(Z[a])) {
            anyNaN = true;
        }
    }
    for (int a = 0; a < n && !anyNaN; ++a) {
        if (isVirtual[static_cast<std::size_t>(a)]) {
            continue;
        }
        for (int b = a + 1; b < n; ++b) {
            if (isVirtual[static_cast<std::size_t>(b)]) {
                continue;
            }
            if (excluded.count(pairKey(a, b)) != 0) {
                continue;
            }
            double dx = X[a] - X[b], dy = Y[a] - Y[b], dz = Z[a] - Z[b];
            minImage(dx, dy, dz);
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (d < minNonbondedNm) {
                minNonbondedNm = d;
                worst = {a, b};
            }
            if (d < kClashNm) {
                ++nClash;
            }
        }
    }

    const double pe = calcOpenMMPotentialEnergy();
    const bool peBad = !std::isfinite(pe) || pe > 1.0e4; // +1e4 kJ/mol is already pathological

    std::fprintf(stderr,
                 "[init] startup geometry check: nAtoms=%d  initial PE=%.4g kJ/mol  "
                 "min non-bonded distance=%.4f nm  hard clashes(<%.2f nm)=%d\n",
                 n,
                 pe,
                 (minNonbondedNm > 1e29 ? 0.0 : minNonbondedNm),
                 kClashNm,
                 nClash);

    if (!anyNaN && nClash == 0 && !peBad) {
        return; // healthy start
    }

    std::string msg = "Context::initialize: the input structure is NOT usable as-is.\n";
    if (anyNaN) {
        msg += "  * coordinates contain NaN/Inf.\n";
    }
    if (nClash > 0 && worst.first >= 0) {
        const auto nm = [&](int i) {
            return (i < (int)systemTopology.atomsUniqueName.size()) ? systemTopology.atomsUniqueName[i]
                                                                    : std::to_string(i);
        };
        char buf[256];
        std::snprintf(buf,
                      sizeof(buf),
                      "  * %d steric clash(es): non-bonded atoms closer than %.2f nm "
                      "(worst: %s -- %s at %.4f nm).\n",
                      nClash,
                      kClashNm,
                      nm(worst.first).c_str(),
                      nm(worst.second).c_str(),
                      minNonbondedNm);
        msg += buf;
    }
    if (peBad) {
        char buf[160];
        std::snprintf(buf,
                      sizeof(buf),
                      "  * initial potential energy is %.4g kJ/mol (clashing/unminimized).\n",
                      pe);
        msg += buf;
    }
    msg += "  The docking world welds the receptor RIGID, so it cannot relax these\n"
           "  clashes -- the run would loop forever rejecting kicks. Energy-MINIMIZE\n"
           "  the structure first (e.g. tleap/sander/OpenMM LocalEnergyMinimizer, or a\n"
           "  Cartesian relaxation world before the docking world). To proceed anyway,\n"
           "  set ROBO_ALLOW_BAD_START=1.\n";

    const char* allow = std::getenv("ROBO_ALLOW_BAD_START");
    if (allow != nullptr && allow[0] != '0' && allow[0] != '\0') {
        std::fprintf(stderr, "[init] WARNING (continuing, ROBO_ALLOW_BAD_START set):\n%s", msg.c_str());
        return;
    }
    throw std::runtime_error(msg);
}
