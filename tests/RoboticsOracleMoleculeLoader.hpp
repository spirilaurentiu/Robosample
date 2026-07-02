#pragma once

// ============================================================================
//  RoboticsOracleMoleculeLoader.hpp -- Scope B (docs/specs/robotics-oracle-
//  differential.md §4.2/§4.4) loader for the PORT's own SystemTopology,
//  dumped by tests/fixtures/robotics_oracle_molecules/_generate_port_topology.py
//  from the SAME context.load_amber path loader_differential uses.
//
//  Deliberately minimal: it fills ONLY the SystemTopology fields
//  World::buildModel (src/World.cpp) and Context::buildFlexibilities
//  (src/Context.cpp) actually read (atoms, bonds, molecule ranges, root
//  mobilities) -- Scope B needs the DECOMPOSITION, never OpenMM energetics,
//  so the force-field arrays are never populated (left at their
//  SystemTopology default, which World::buildModel never touches).
// ============================================================================

#include "RoboticsOracleJson.hpp"
#include "TopologyElements.hpp"
#include "cnpy.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace robotics_oracle_loader {

inline SystemTopology loadPortSystemTopology(const std::string& fixturesDir, const std::string& caseName) {
    const cnpy::npz_t npz = cnpy::npz_load(fixturesDir + "/" + caseName + ".systopo.npz");
    auto need = [&](const std::string& name) -> const cnpy::NpyArray& {
        const auto it = npz.find(name);
        if (it == npz.end()) {
            throw std::runtime_error("robotics_oracle_molecules: npz '" + caseName + "' missing array '" + name
                                     + "'");
        }
        return it->second;
    };

    SystemTopology sys;
    sys.numMolecules = static_cast<int>(need("num_molecules").data<std::int32_t>()[0]);
    sys.numAtoms = static_cast<int>(need("num_atoms").data<std::int32_t>()[0]);
    sys.numBonds = static_cast<int>(need("num_bonds").data<std::int32_t>()[0]);

    const auto copyI32 = [&](const std::string& name, int n) {
        const cnpy::NpyArray& a = need(name);
        if (static_cast<int>(a.num_vals()) != n) {
            throw std::runtime_error("robotics_oracle_molecules: '" + name + "' has wrong extent");
        }
        const std::int32_t* d = a.data<std::int32_t>();
        return std::vector<int>(d, d + n);
    };
    const auto copyF64 = [&](const std::string& name, int n) {
        const cnpy::NpyArray& a = need(name);
        if (static_cast<int>(a.num_vals()) != n) {
            throw std::runtime_error("robotics_oracle_molecules: '" + name + "' has wrong extent");
        }
        const double* d = a.data<double>();
        return std::vector<double>(d, d + n);
    };

    sys.atomsBegin = copyI32("atoms_begin", sys.numMolecules);
    sys.atomsEnd = copyI32("atoms_end", sys.numMolecules);
    {
        const std::vector<int> rm = copyI32("root_mobilities", sys.numMolecules);
        sys.rootMobilities.resize(rm.size());
        for (std::size_t i = 0; i < rm.size(); ++i) {
            sys.rootMobilities[i] = static_cast<JointType>(rm[i]);
        }
    }

    sys.atomsX = copyF64("atoms_x", sys.numAtoms);
    sys.atomsY = copyF64("atoms_y", sys.numAtoms);
    sys.atomsZ = copyF64("atoms_z", sys.numAtoms);
    sys.atomsMass = copyF64("atoms_mass", sys.numAtoms);
    sys.atomsNumBondsInvolved = copyI32("atoms_num_bonds_involved", sys.numAtoms);
    sys.atomsPrmtopIndex = copyI32("atoms_prmtop_index", sys.numAtoms);

    sys.bondsI = copyI32("bonds_i", sys.numBonds);
    sys.bondsJ = copyI32("bonds_j", sys.numBonds);
    {
        const cnpy::NpyArray& a = need("bonds_ring_closing");
        if (static_cast<int>(a.num_vals()) != sys.numBonds) {
            throw std::runtime_error("robotics_oracle_molecules: 'bonds_ring_closing' has wrong extent");
        }
        const bool* d = a.data<bool>();
        sys.bondsRingClosing.assign(d, d + sys.numBonds);
    }

    return sys;
}

// ============================================================================
//  Clone-side numeric dump loader (docs/specs/robotics-oracle-differential.md
//  §4.2 Scope-B build path, resolved follow-up 4 / §6 staged comparison):
//  reads the <case>.moldyn.npz + <case>.moldyn.manifest.json pair written by
//  tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py
//  (Robosample/src/RoboticsOracleMoleculeDump.cpp's dump, via the additive
//  World::dump_robotics_oracle_molecule pybind binding). One Simbody body's
//  worth of §7 quantities per MoleculeOracleBody; one state (rest/random) per
//  MoleculeOracleState; the whole case (all states, constant body layout) is
//  MoleculeOracleCase.
// ============================================================================
struct MoleculeOracleBody {
    std::vector<int> atomPrmtopIndices; // sorted; the §5 correspondence key
    // Frame-invariant per-atom Ground position anchor (docs/specs/robotics-
    // oracle-differential.md refined §6): flattened xyz, PARALLEL to
    // atomPrmtopIndices -- atomPosG[3*i + {0,1,2}] is atom
    // atomPrmtopIndices[i]'s Ground location.
    std::vector<double> atomPosG;
    int nq = 0;
    int nu = 0;
    std::vector<double> q;
    std::vector<double> u;
    std::vector<double> udot;
    double X_GB_R[9]{}; // Scope-A-only anchor (§6 NOTE); NOT an admissible
                         // cross-engine anchor in Scope B -- kept for diagnostics.
    double X_GB_p[3]{};
    double V_GB_ang[3]{};
    double V_GB_lin[3]{};
    double A_GB_ang[3]{};
    double A_GB_lin[3]{};
    // Zero-force mobilizer reaction transmitted to this body, Ground-expressed,
    // [angular=torque; linear=force] (docs/specs/robotics-oracle-reactions.md
    // §4). reactionBo (at body origin) is the convention-free cross-engine
    // anchor; reactionMo (at outboard frame origin Mo) has a convention-gated
    // angular part. This is the quantity World::calcSpatialForces stores.
    double reactionBo_ang[3]{};
    double reactionBo_lin[3]{};
    double reactionMo_ang[3]{};
    double reactionMo_lin[3]{};
    // Per-body min-eig(D_b) (docs/specs/singular-dof-fixman.md); +inf for a
    // 0-dof body (RoboticsOracleMoleculeDump.hpp::BodyDump::minEigD). Lets
    // the numeric differential guard V_GB.angular/A_GB.angular PER-BODY
    // (phantom/near-singular hinge -> gauge-unobservable axis direction)
    // instead of tree-wide.
    double minEigD = std::numeric_limits<double>::infinity();
};

struct MoleculeOracleState {
    std::string label; // "rest" | "random"
    std::vector<MoleculeOracleBody> bodies; // parallel to MoleculeOracleCase::bodyLayout
    double logDetM = 0;
    double kineticEnergy = 0;
    double normUdot = 0;
    double minEigD = 0;
};

struct MoleculeOracleCase {
    std::string caseName;
    std::string moleculeClass; // "rigid" | "regular" | "cyclic"
    std::vector<MoleculeOracleState> states;
    // The CLONE's own ACTUAL ring-closing (cotree) bond set, keyed by raw
    // prmtopIndex bond endpoints (§5 cross-engine atom identity), read back
    // from `context.system_topology.bonds` by `_generate_molecule_oracle.py`
    // AFTER Molmodel built the molecule -- not merely what was requested via
    // the shared-tree override. Empty for classes with no ring closures
    // (rigid/regular). Shared-tree fix (docs/specs/robotics-oracle-
    // differential.md Scope B §6): lets the port-side test assert SET
    // equality (not just count) against its own SystemTopology.bondsRingClosing.
    std::vector<std::pair<int, int>> ringClosingBondPrmtopPairs;
};

inline MoleculeOracleCase loadMoleculeOracleCase(const std::string& fixturesDir, const std::string& caseName) {
    const JsonValue manifest = loadManifestJsonFile(fixturesDir + "/" + caseName + ".moldyn.manifest.json");
    const cnpy::npz_t npz = cnpy::npz_load(fixturesDir + "/" + caseName + ".moldyn.npz");

    auto need = [&](const std::string& name) -> const cnpy::NpyArray& {
        const auto it = npz.find(name);
        if (it == npz.end()) {
            throw std::runtime_error("robotics_oracle_molecules: '" + caseName + ".moldyn.npz' missing array '"
                                     + name + "'");
        }
        return it->second;
    };
    auto copyF64 = [&](const std::string& name, std::size_t n) {
        const cnpy::NpyArray& a = need(name);
        if (a.num_vals() != n) {
            throw std::runtime_error("robotics_oracle_molecules: '" + name + "' has wrong extent (expected "
                                     + std::to_string(n) + ", got " + std::to_string(a.num_vals()) + ")");
        }
        const double* d = a.data<double>();
        return std::vector<double>(d, d + n);
    };
    auto scalarF64 = [&](const std::string& name) { return copyF64(name, 1)[0]; };

    MoleculeOracleCase out;
    out.caseName = caseName;
    out.moleculeClass = manifest.at("molecule_class").strVal;

    if (manifest.has("ring_closing_bond_prmtop_pairs")) {
        for (const JsonValue& pairVal : manifest.at("ring_closing_bond_prmtop_pairs").arrVal) {
            if (pairVal.arrVal.size() != 2) {
                throw std::runtime_error("robotics_oracle_molecules: '" + caseName
                                         + "' ring_closing_bond_prmtop_pairs entry has wrong arity");
            }
            out.ringClosingBondPrmtopPairs.emplace_back(static_cast<int>(pairVal.arrVal[0].numVal),
                                                         static_cast<int>(pairVal.arrVal[1].numVal));
        }
    }

    const int numBodies = static_cast<int>(manifest.at("num_bodies").numVal);
    const JsonValue& bodiesSpec = manifest.at("bodies");
    if (static_cast<int>(bodiesSpec.arrVal.size()) != numBodies) {
        throw std::runtime_error("robotics_oracle_molecules: '" + caseName + "' manifest num_bodies mismatch");
    }

    for (const JsonValue& stateNameVal : manifest.at("states").arrVal) {
        const std::string label = stateNameVal.strVal;
        const int s = static_cast<int>(out.states.size());
        MoleculeOracleState st;
        st.label = label;
        st.logDetM = scalarF64("state" + std::to_string(s) + "_logDetM");
        st.kineticEnergy = scalarF64("state" + std::to_string(s) + "_KE");
        st.normUdot = scalarF64("state" + std::to_string(s) + "_normUdot");
        st.minEigD = scalarF64("state" + std::to_string(s) + "_minEigD");

        st.bodies.resize(static_cast<std::size_t>(numBodies));
        for (int b = 0; b < numBodies; ++b) {
            const JsonValue& bSpec = bodiesSpec.arrVal[static_cast<std::size_t>(b)];
            MoleculeOracleBody& bd = st.bodies[static_cast<std::size_t>(b)];
            bd.nq = static_cast<int>(bSpec.at("nq").numVal);
            bd.nu = static_cast<int>(bSpec.at("nu").numVal);
            bd.atomPrmtopIndices.reserve(bSpec.at("atom_prmtop_indices").arrVal.size());
            for (const JsonValue& a : bSpec.at("atom_prmtop_indices").arrVal) {
                bd.atomPrmtopIndices.push_back(static_cast<int>(a.numVal));
            }

            const std::string prefix = "state" + std::to_string(s) + "_body" + std::to_string(b) + "_";
            bd.q = copyF64(prefix + "q", static_cast<std::size_t>(bd.nq));
            bd.u = copyF64(prefix + "u", static_cast<std::size_t>(bd.nu));
            bd.udot = copyF64(prefix + "udot", static_cast<std::size_t>(bd.nu));
            bd.atomPosG = copyF64(prefix + "atomPosG", bd.atomPrmtopIndices.size() * 3);
            const std::vector<double> R = copyF64(prefix + "X_GB_R", 9);
            std::copy(R.begin(), R.end(), bd.X_GB_R);
            const std::vector<double> p = copyF64(prefix + "X_GB_p", 3);
            std::copy(p.begin(), p.end(), bd.X_GB_p);
            const std::vector<double> vAng = copyF64(prefix + "V_GB_ang", 3);
            std::copy(vAng.begin(), vAng.end(), bd.V_GB_ang);
            const std::vector<double> vLin = copyF64(prefix + "V_GB_lin", 3);
            std::copy(vLin.begin(), vLin.end(), bd.V_GB_lin);
            const std::vector<double> aAng = copyF64(prefix + "A_GB_ang", 3);
            std::copy(aAng.begin(), aAng.end(), bd.A_GB_ang);
            const std::vector<double> aLin = copyF64(prefix + "A_GB_lin", 3);
            std::copy(aLin.begin(), aLin.end(), bd.A_GB_lin);
            const std::vector<double> rBoAng = copyF64(prefix + "reactionBo_ang", 3);
            std::copy(rBoAng.begin(), rBoAng.end(), bd.reactionBo_ang);
            const std::vector<double> rBoLin = copyF64(prefix + "reactionBo_lin", 3);
            std::copy(rBoLin.begin(), rBoLin.end(), bd.reactionBo_lin);
            const std::vector<double> rMoAng = copyF64(prefix + "reactionMo_ang", 3);
            std::copy(rMoAng.begin(), rMoAng.end(), bd.reactionMo_ang);
            const std::vector<double> rMoLin = copyF64(prefix + "reactionMo_lin", 3);
            std::copy(rMoLin.begin(), rMoLin.end(), bd.reactionMo_lin);
            bd.minEigD = scalarF64(prefix + "minEigD");
        }
        out.states.push_back(std::move(st));
    }

    return out;
}

} // namespace robotics_oracle_loader
