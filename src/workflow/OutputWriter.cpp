#include "Context.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "DCDWriter.hpp"

namespace {

// Convert the three reduced lattice vectors (nm, lower-triangular, row-major
// a,b,c) into a CHARMM/DCD Box: side lengths in Angstrom and angles in degrees.
dcd::Box boxFromReducedVectors(const std::vector<double>& bv) {
    dcd::Box box;
    if (bv.size() != 9) {
        return box; // default 1 A cube
    }
    const double ax = bv[0];
    const double bx = bv[3], by = bv[4];
    const double cx = bv[6], cy = bv[7], cz = bv[8];
    const double a = ax;
    const double b = std::sqrt(bx * bx + by * by);
    const double c = std::sqrt(cx * cx + cy * cy + cz * cz);
    constexpr double kRad2Deg = 57.29577951308232;
    const double cosGamma = (a > 0 && b > 0) ? (ax * bx) / (a * b) : 0.0;
    const double cosBeta = (a > 0 && c > 0) ? (ax * cx) / (a * c) : 0.0;
    const double cosAlpha = (b > 0 && c > 0) ? (bx * cx + by * cy) / (b * c) : 0.0;
    constexpr double kNm2Ang = 10.0;
    box.sideA = a * kNm2Ang;
    box.sideB = b * kNm2Ang;
    box.sideC = c * kNm2Ang;
    box.angleAlpha = std::acos(std::max(-1.0, std::min(1.0, cosAlpha))) * kRad2Deg;
    box.angleBeta = std::acos(std::max(-1.0, std::min(1.0, cosBeta))) * kRad2Deg;
    box.angleGamma = std::acos(std::max(-1.0, std::min(1.0, cosGamma))) * kRad2Deg;
    return box;
}

} // namespace

void Context::writeOutputs(int replica, int round, bool verbose) {
    writeOutputsCore(replica, round, verbose, replicaCoords_[replica], temperatures_[replica]);
}

// Pure extraction from the former writeOutputs body -- `idx` used to be
// `replica`/`replicaCoords_[replica]`/`temperatures_[replica]` verbatim.
// RunREX (the label-swap driver) reuses this on Replica-owned coordinates,
// indexed by THERMODYNAMIC STATE rather than replica-object identity (see
// the RunREX doc comment in Context.hpp).
void Context::writeOutputsCore(int idx, int round, bool verbose, const std::vector<robo::Vec3>& coords, double T) {
    const double pe = openmmPotential(coords);
    const std::string path = baseName + "." + std::to_string(idx) + ".csv";
    std::ofstream out(path, std::ios::app);
    if (out) {
        out << round << "," << idx << "," << T << "," << pe << "\n";
    }
    if (verbose) {
        std::printf("[rex] round=%d replica=%d T=%.1f PE=%.4f kJ/mol\n", round, idx, T, pe);
    }

    if (idx >= 0 && idx < static_cast<int>(dcdWriters_.size())) {
        const int n = systemTopology.numAtoms;
        const auto& perm = systemTopology.atomsPrmtopIndex;
        const bool havePerm = (static_cast<int>(perm.size()) == n);
        dcdScratch_.resize(static_cast<std::size_t>(3 * n));

        // Whole-molecule periodic imaging, applied ONLY to the output copy.
        // replicaCoords_ stays unwrapped (contiguous per molecule) so the next
        // setAtomsLocationsInGround/frame rebuild is unaffected. Each molecule is
        // shifted by integer lattice vectors so its center of mass lands in the
        // primary cell, then translated rigidly -- bonds never straddle a face.
        const auto& bv = systemTopology.boxVectors;
        const int numMol = systemTopology.numMolecules;
        const bool haveRanges = OpenMMContext::isPeriodic(systemTopology.nonbondedMethod) && bv.size() == 9
                                && static_cast<int>(systemTopology.atomsBegin.size()) == numMol
                                && static_cast<int>(systemTopology.atomsEnd.size()) == numMol;

        auto scatter = [&](int a, double sx, double sy, double sz) {
            const int p = havePerm ? perm[a] : a;
            dcdScratch_[3 * p + 0] = (coords[a][0] + sx) * 10.0;
            dcdScratch_[3 * p + 1] = (coords[a][1] + sy) * 10.0;
            dcdScratch_[3 * p + 2] = (coords[a][2] + sz) * 10.0;
        };

        if (!haveRanges) {
            // Non-periodic (or missing ranges): write coordinates verbatim.
            for (int a = 0; a < n; ++a) {
                scatter(a, 0.0, 0.0, 0.0);
            }
        } else {
            // Reduced lower-triangular box: a=(bv0,0,0) b=(bv3,bv4,0) c=(bv6,bv7,bv8).
            // Wrap in c -> b -> a order (same convention as the clash-scan minImage).
            for (int m = 0; m < numMol; ++m) {
                const int beg = systemTopology.atomsBegin[m];
                const int end = systemTopology.atomsEnd[m];
                if (beg >= end) {
                    continue;
                }

                // Mass-weighted COM in engine order (molecule is intact here).
                double cx = 0.0, cy = 0.0, cz = 0.0, mtot = 0.0;
                const int nMass = static_cast<int>(systemTopology.atomsMass.size());
                for (int a = beg; a < end; ++a) {
                    const double mass = (a < nMass) ? systemTopology.atomsMass[a] : 1.0;
                    cx += mass * coords[a][0];
                    cy += mass * coords[a][1];
                    cz += mass * coords[a][2];
                    mtot += mass;
                }
                if (mtot > 0.0) {
                    cx /= mtot;
                    cy /= mtot;
                    cz /= mtot;
                } else { // all-massless (e.g. pure virtual sites): use first atom
                    cx = coords[beg][0];
                    cy = coords[beg][1];
                    cz = coords[beg][2];
                }

                // Accumulate the rigid shift that brings the COM into [0, L).
                double sx = 0.0, sy = 0.0, sz = 0.0;
                const double nc = (bv[8] != 0.0) ? std::floor(cz / bv[8]) : 0.0;
                cx -= nc * bv[6];
                cy -= nc * bv[7];
                sx -= nc * bv[6];
                sy -= nc * bv[7];
                sz -= nc * bv[8];
                const double nb = (bv[4] != 0.0) ? std::floor(cy / bv[4]) : 0.0;
                cx -= nb * bv[3];
                sx -= nb * bv[3];
                sy -= nb * bv[4];
                const double na = (bv[0] != 0.0) ? std::floor(cx / bv[0]) : 0.0;
                sx -= na * bv[0];

                for (int a = beg; a < end; ++a) {
                    scatter(a, sx, sy, sz);
                }
            }
        }

        dcdWriters_[idx].append(dcdScratch_, boxFromReducedVectors(systemTopology.boxVectors));
    }
}

// docs/specs/reaction-force-monitoring.md Sec.4: one CSV per replica (the
// filename carries the replica, matching the per-replica .dcd), comma-
// delimited, appended as frames are produced. No-op on an empty row set (a
// non-write round never reaches here, and a reporter world with zero
// interesting bodies would otherwise write an empty file with only a header).
// Always the same 10 columns regardless of which term(s) the reporter world
// summed into `force`/`torque` (docs/specs/reaction-force-monitoring.md
// Sec.1.2) -- that choice is made once, at World::enableReactionReporter,
// and is invisible to the CSV format.
void Context::writeReactionRows(int replica, int frame, const std::vector<ReactionSample>& rows) {
    if (rows.empty()) {
        return;
    }
    const std::string path = baseName + "." + std::to_string(replica) + ".reactions.csv";
    std::ofstream out(path, std::ios::app);
    if (!out) {
        return;
    }
    if (out.tellp() == 0) {
        out << "# frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz\n";
    }
    for (const auto& s : rows) {
        out << frame << ',' << replica << ',' << s.bodyIdx << ',' << s.atomIdx << ',' << s.force[0] << ','
            << s.force[1] << ',' << s.force[2] << ',' << s.torque[0] << ',' << s.torque[1] << ',' << s.torque[2]
            << '\n';
    }
}
