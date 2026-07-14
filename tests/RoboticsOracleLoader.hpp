#pragma once

// ============================================================================
//  RoboticsOracleLoader.hpp -- runtime loader for the robotics-oracle
//  fixtures (docs/specs/robotics-oracle-differential.md §9,
//  docs/specs/robotics-oracle-data-provenance.md §4/§6): reads a case's
//  <name>.npz (cnpy) + <name>.manifest.json (hand-rolled JSON, see below)
//  from ROBOTICS_ORACLE_FIXTURE_DIR and reassembles the SAME
//  OracleCase/OracleMultiCase/OracleAggregateCase PODs
//  (fixtures/robotics_oracle/RoboticsOracleTypes.hpp) that used to be
//  `#include`d as generated C++ constants. TestRoboticsOracle.cpp's staged
//  comparison (runOracleCase/runOracleMultiCase/runOracleAggregateCase) is
//  therefore UNCHANGED -- this file only changes where the struct's data
//  comes from, never what is compared (the migration's hard constraint).
//
//  Included by tests/support/RoboticsOracleRunners.hpp (TEST-005), which is
//  in turn included by every TestRoboticsOracle{SingleState,MultiSystem,
//  Aggregate,Fuzz}.cpp split binary and TestRoboticsOracleMolecule{Symbolic,
//  Numeric}.cpp -- one binary per TU, so everything here is `inline`/
//  header-only without ODR risk (the `inline` functions already tolerate
//  multiple TUs).
//
//  JSON: the manifest is small, fixed-schema, machine-written-by-us-only
//  JSON (§9: "small and human-reviewable"). A full JSON library dependency
//  buys nothing a ~150-line recursive-descent parser doesn't already give
//  (Rule 2) -- this parser handles the JSON subset the generator emits
//  (objects/arrays/strings/numbers/bool/null), nothing more.
// ============================================================================

#include "RobotModel.hpp"
#include "RoboticsOracleJson.hpp"
#include "cnpy.hpp"
#include "fixtures/robotics_oracle/RoboticsOracleTypes.hpp"

#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace robotics_oracle_loader {

// JsonValue/JsonParser/readWholeTextFile/loadManifestJson/internString now
// live in RoboticsOracleJson.hpp (shared with Scope B's
// RoboticsOracleMoleculeLoader.hpp, same namespace) -- pure extraction, see
// that header's banner; nothing below this line changed.

// ---------------------------------------------------------------------------
//  §9 no-silent-gap guard: the manifest's array-name list must cover every
//  quantity this schema_version's OracleState/OracleBodyOutput/
//  OracleAggregateState is expected to carry. A newly added struct field
//  whose writer forgot to emit it (or whose fixtures are stale) fails HERE,
//  loudly, before any numeric comparison -- not by silently reading zeros.
// ---------------------------------------------------------------------------
inline const std::vector<std::string>& expectedSingleStateFields(int schemaVersion) {
    static const std::vector<std::string> v1 = {
        "q", "u", "bodyForceTorque", "bodyForceForce", "mobilityForce", "X_GB_R", "X_GB_p", "X_FM_R", "X_FM_p",
        "V_GB_ang", "V_GB_lin", "qdot", "P_J", "P_F", "P_M", "PPlus_J", "PPlus_F", "PPlus_M", "DI", "G_ang",
        "G_lin", "minEigD", "Z_ang", "Z_lin", "ZPlus_ang", "ZPlus_lin", "eps", "udot", "A_GB_ang", "A_GB_lin",
        "Mdense", "logDetM", "reactionBoAng", "reactionBoLin", "reactionMoAng", "reactionMoLin"};
    if (schemaVersion == 1) {
        return v1;
    }
    throw std::runtime_error("robotics_oracle: no expected single-state field list for schema_version " +
                              std::to_string(schemaVersion));
}
inline const std::vector<std::string>& expectedBodyOutputFields(int schemaVersion) {
    static const std::vector<std::string> v1 = {
        "X_GB_R", "X_GB_p", "X_FM_R", "X_FM_p", "V_GB_ang", "V_GB_lin", "qdot", "P_J", "P_F", "P_M", "PPlus_J",
        "PPlus_F", "PPlus_M", "DI", "G_ang", "G_lin", "minEigD", "Z_ang", "Z_lin", "ZPlus_ang", "ZPlus_lin",
        "eps", "udot", "A_GB_ang", "A_GB_lin", "reactionBoAng", "reactionBoLin", "reactionMoAng",
        "reactionMoLin"};
    if (schemaVersion == 1) {
        return v1;
    }
    throw std::runtime_error("robotics_oracle: no expected body-output field list for schema_version " +
                              std::to_string(schemaVersion));
}
inline const std::vector<std::string>& expectedMultiSystemFields(int schemaVersion) {
    static const std::vector<std::string> v1 = {"q",  "u", "bodyForceTorque",   "bodyForceForce",
                                                 "groundForceTorque", "groundForceForce", "mobilityForce",
                                                 "Mdense", "logDetM"};
    if (schemaVersion == 1) {
        return v1;
    }
    throw std::runtime_error("robotics_oracle: no expected multi-system field list for schema_version " +
                              std::to_string(schemaVersion));
}
inline const std::vector<std::string>& expectedAggregateStateFields(int schemaVersion) {
    static const std::vector<std::string> v1 = {
        "q", "u", "bodyForceTorque", "bodyForceForce", "logDetM", "totalKE", "normUdot", "minEigD",
        "reportX_GB_R", "reportX_GB_p", "reportV_GB_ang", "reportV_GB_lin", "reportA_GB_ang", "reportA_GB_lin"};
    if (schemaVersion == 1) {
        return v1;
    }
    throw std::runtime_error("robotics_oracle: no expected aggregate-state field list for schema_version " +
                              std::to_string(schemaVersion));
}
// §8.3 fuzz batch (added under the SAME schema_version=1: a new case_kind is
// schema-additive -- it introduces no new field on any EXISTING struct, so
// no existing fixture's expected-field list changes; only "fuzz" fixtures
// are checked against this list, §9's no-silent-gap guard applied to the new
// kind exactly like every other kind).
inline const std::vector<std::string>& expectedFuzzStateFields(int schemaVersion) {
    static const std::vector<std::string> v1 = {"q",       "u",        "bodyForceTorque", "bodyForceForce",
                                                 "logDetM", "totalKE",  "normUdot",        "minEigD",
                                                 "udot",    "A_GB_ang", "A_GB_lin"};
    if (schemaVersion == 1) {
        return v1;
    }
    throw std::runtime_error("robotics_oracle: no expected fuzz-state field list for schema_version " +
                              std::to_string(schemaVersion));
}

inline void checkNoSilentGap(const JsonValue& manifest, const std::string& caseName) {
    const int schemaVersion = static_cast<int>(manifest.at("schema_version").numVal);
    const std::string kind = manifest.at("case_kind").strVal;

    std::vector<std::string> present;
    present.reserve(manifest.at("arrays").arrVal.size());
    for (const JsonValue& a : manifest.at("arrays").arrVal) {
        present.push_back(a.at("name").strVal);
    }
    std::sort(present.begin(), present.end());
    auto has = [&present](const std::string& name) {
        return std::binary_search(present.begin(), present.end(), name);
    };
    auto requireField = [&](const std::string& name) {
        if (!has(name)) {
            throw std::runtime_error("robotics_oracle: fixture '" + caseName + "' (schema_version " +
                                      std::to_string(schemaVersion) + ") is missing expected array '" + name +
                                      "' -- regenerate the fixtures (gen_robotics_oracle)");
        }
    };

    const auto& states = manifest.at("states").arrVal;
    if (kind == "single") {
        for (std::size_t i = 0; i < states.size(); ++i) {
            const std::string prefix = "state" + std::to_string(i) + "_";
            for (const std::string& f : expectedSingleStateFields(schemaVersion)) {
                requireField(prefix + f);
            }
        }
    } else if (kind == "multi") {
        const int numBodies = static_cast<int>(manifest.at("model").at("numBodies").numVal);
        for (std::size_t i = 0; i < states.size(); ++i) {
            const std::string sp = "state" + std::to_string(i) + "_";
            for (const std::string& f : expectedMultiSystemFields(schemaVersion)) {
                requireField(sp + f);
            }
            for (int b = 0; b < numBodies; ++b) {
                const std::string bp = sp + "body" + std::to_string(b) + "_";
                for (const std::string& f : expectedBodyOutputFields(schemaVersion)) {
                    requireField(bp + f);
                }
            }
        }
    } else if (kind == "aggregate") {
        for (std::size_t i = 0; i < states.size(); ++i) {
            const std::string sp = "state" + std::to_string(i) + "_";
            for (const std::string& f : expectedAggregateStateFields(schemaVersion)) {
                requireField(sp + f);
            }
        }
    } else if (kind == "fuzz") {
        for (std::size_t i = 0; i < states.size(); ++i) {
            const std::string sp = "state" + std::to_string(i) + "_";
            for (const std::string& f : expectedFuzzStateFields(schemaVersion)) {
                requireField(sp + f);
            }
        }
    } else {
        throw std::runtime_error("robotics_oracle: fixture '" + caseName + "' has unknown case_kind '" + kind +
                                  "'");
    }
}

// ---------------------------------------------------------------------------
//  npz array -> struct field helpers.
// ---------------------------------------------------------------------------
inline void copyNpzArr(const cnpy::npz_t& npz, const std::string& name, double* out, std::size_t n) {
    const auto it = npz.find(name);
    if (it == npz.end()) {
        throw std::runtime_error("robotics_oracle: npz missing array '" + name + "'");
    }
    const cnpy::NpyArray& a = it->second;
    if (a.fortran_order && a.shape.size() > 1) {
        // The generator (Robosample/tools/gen_robotics_oracle.cpp) always
        // writes 'fortran_order': False, so this can't currently trigger --
        // but a Fortran-ordered multi-dim array's flat byte layout does NOT
        // match the row-major copy below, and would silently transpose data
        // instead of erroring. Fail loud rather than assume.
        throw std::runtime_error("robotics_oracle: npz array '" + name +
                                  "' is Fortran-ordered with rank > 1; row-major copy would silently "
                                  "transpose it");
    }
    if (a.num_vals() != n) {
        throw std::runtime_error("robotics_oracle: npz array '" + name + "' has " +
                                  std::to_string(a.num_vals()) + " values, expected " + std::to_string(n));
    }
    if (n > 0) {
        const double* d = a.data<double>();
        std::copy(d, d + n, out);
    }
}
inline double scalarFromNpz(const cnpy::npz_t& npz, const std::string& name) {
    double v = 0;
    copyNpzArr(npz, name, &v, 1);
    return v;
}
// Per-body nu from the npz itself (the "eps" array's true extent) rather
// than re-deriving it from JointType -- keeps the loader purely data-driven
// (matches whatever the generator actually wrote, not a second hand-kept
// joint table that could drift from it).
inline int nuFromNpz(const cnpy::npz_t& npz, const std::string& prefix) {
    const auto it = npz.find(prefix + "eps");
    if (it == npz.end()) {
        throw std::runtime_error("robotics_oracle: npz missing '" + prefix + "eps' (needed to size nu)");
    }
    return it->second.shape.empty() ? 0 : static_cast<int>(it->second.shape[0]);
}
inline int nqFromNpz(const cnpy::npz_t& npz, const std::string& prefix) {
    const auto it = npz.find(prefix + "qdot");
    if (it == npz.end()) {
        throw std::runtime_error("robotics_oracle: npz missing '" + prefix + "qdot' (needed to size nq)");
    }
    return it->second.shape.empty() ? 0 : static_cast<int>(it->second.shape[0]);
}

inline void copyJsonNumArr(const JsonValue& arr, double* out, int n) {
    if (static_cast<int>(arr.arrVal.size()) != n) {
        throw std::runtime_error("robotics_oracle: manifest JSON array has wrong length");
    }
    for (int i = 0; i < n; ++i) {
        out[static_cast<std::size_t>(i)] = arr.arrVal[static_cast<std::size_t>(i)].numVal;
    }
}

inline robotics_oracle::BodyModelSpec bodySpecFromJson(const JsonValue& b) {
    robotics_oracle::BodyModelSpec out{};
    out.parent = static_cast<int>(b.at("parent").numVal);
    out.joint = static_cast<int>(b.at("joint").numVal);
    copyJsonNumArr(b.at("X_PF_R"), out.X_PF_R, 9);
    copyJsonNumArr(b.at("X_PF_p"), out.X_PF_p, 3);
    copyJsonNumArr(b.at("X_BM_R"), out.X_BM_R, 9);
    copyJsonNumArr(b.at("X_BM_p"), out.X_BM_p, 3);
    out.mass = b.at("mass").numVal;
    copyJsonNumArr(b.at("com_B"), out.com_B, 3);
    copyJsonNumArr(b.at("unitInertia_B"), out.unitInertia_B, 6);
    return out;
}

inline void loadStateArrays(const cnpy::npz_t& npz, const std::string& p, robotics_oracle::OracleState& out) {
    const std::size_t nq = static_cast<std::size_t>(out.nq);
    const std::size_t nu = static_cast<std::size_t>(out.nu);
    copyNpzArr(npz, p + "q", out.q, nq);
    copyNpzArr(npz, p + "u", out.u, nu);
    copyNpzArr(npz, p + "bodyForceTorque", out.bodyForceTorque, 3);
    copyNpzArr(npz, p + "bodyForceForce", out.bodyForceForce, 3);
    copyNpzArr(npz, p + "mobilityForce", out.mobilityForce, nu);

    copyNpzArr(npz, p + "X_GB_R", out.X_GB_R, 9);
    copyNpzArr(npz, p + "X_GB_p", out.X_GB_p, 3);
    copyNpzArr(npz, p + "X_FM_R", out.X_FM_R, 9);
    copyNpzArr(npz, p + "X_FM_p", out.X_FM_p, 3);

    copyNpzArr(npz, p + "V_GB_ang", out.V_GB_ang, 3);
    copyNpzArr(npz, p + "V_GB_lin", out.V_GB_lin, 3);
    copyNpzArr(npz, p + "qdot", out.qdot, nq);

    copyNpzArr(npz, p + "P_J", out.P_J, 6);
    copyNpzArr(npz, p + "P_F", out.P_F, 9);
    copyNpzArr(npz, p + "P_M", out.P_M, 6);
    copyNpzArr(npz, p + "PPlus_J", out.PPlus_J, 6);
    copyNpzArr(npz, p + "PPlus_F", out.PPlus_F, 9);
    copyNpzArr(npz, p + "PPlus_M", out.PPlus_M, 6);
    copyNpzArr(npz, p + "DI", out.DI, nu * nu);
    copyNpzArr(npz, p + "G_ang", &out.G_ang[0][0], nu * 3);
    copyNpzArr(npz, p + "G_lin", &out.G_lin[0][0], nu * 3);
    out.minEigD = scalarFromNpz(npz, p + "minEigD");

    copyNpzArr(npz, p + "Z_ang", out.Z_ang, 3);
    copyNpzArr(npz, p + "Z_lin", out.Z_lin, 3);
    copyNpzArr(npz, p + "ZPlus_ang", out.ZPlus_ang, 3);
    copyNpzArr(npz, p + "ZPlus_lin", out.ZPlus_lin, 3);
    copyNpzArr(npz, p + "eps", out.eps, nu);
    copyNpzArr(npz, p + "udot", out.udot, nu);
    copyNpzArr(npz, p + "A_GB_ang", out.A_GB_ang, 3);
    copyNpzArr(npz, p + "A_GB_lin", out.A_GB_lin, 3);

    copyNpzArr(npz, p + "Mdense", out.Mdense, nu * nu);
    out.logDetM = scalarFromNpz(npz, p + "logDetM");
    copyNpzArr(npz, p + "reactionBoAng", out.reactionBoAng, 3);
    copyNpzArr(npz, p + "reactionBoLin", out.reactionBoLin, 3);
    copyNpzArr(npz, p + "reactionMoAng", out.reactionMoAng, 3);
    copyNpzArr(npz, p + "reactionMoLin", out.reactionMoLin, 3);
}

inline void loadBodyOutputArrays(const cnpy::npz_t& npz, const std::string& p,
                                 robotics_oracle::OracleBodyOutput& out, int nuI, int nqI) {
    const std::size_t nu = static_cast<std::size_t>(nuI);
    const std::size_t nq = static_cast<std::size_t>(nqI);
    copyNpzArr(npz, p + "X_GB_R", out.X_GB_R, 9);
    copyNpzArr(npz, p + "X_GB_p", out.X_GB_p, 3);
    copyNpzArr(npz, p + "X_FM_R", out.X_FM_R, 9);
    copyNpzArr(npz, p + "X_FM_p", out.X_FM_p, 3);
    copyNpzArr(npz, p + "V_GB_ang", out.V_GB_ang, 3);
    copyNpzArr(npz, p + "V_GB_lin", out.V_GB_lin, 3);
    copyNpzArr(npz, p + "qdot", out.qdot, nq);
    copyNpzArr(npz, p + "P_J", out.P_J, 6);
    copyNpzArr(npz, p + "P_F", out.P_F, 9);
    copyNpzArr(npz, p + "P_M", out.P_M, 6);
    copyNpzArr(npz, p + "PPlus_J", out.PPlus_J, 6);
    copyNpzArr(npz, p + "PPlus_F", out.PPlus_F, 9);
    copyNpzArr(npz, p + "PPlus_M", out.PPlus_M, 6);
    copyNpzArr(npz, p + "DI", out.DI, nu * nu);
    copyNpzArr(npz, p + "G_ang", &out.G_ang[0][0], nu * 3);
    copyNpzArr(npz, p + "G_lin", &out.G_lin[0][0], nu * 3);
    out.minEigD = scalarFromNpz(npz, p + "minEigD");
    copyNpzArr(npz, p + "Z_ang", out.Z_ang, 3);
    copyNpzArr(npz, p + "Z_lin", out.Z_lin, 3);
    copyNpzArr(npz, p + "ZPlus_ang", out.ZPlus_ang, 3);
    copyNpzArr(npz, p + "ZPlus_lin", out.ZPlus_lin, 3);
    copyNpzArr(npz, p + "eps", out.eps, nu);
    copyNpzArr(npz, p + "udot", out.udot, nu);
    copyNpzArr(npz, p + "A_GB_ang", out.A_GB_ang, 3);
    copyNpzArr(npz, p + "A_GB_lin", out.A_GB_lin, 3);
    copyNpzArr(npz, p + "reactionBoAng", out.reactionBoAng, 3);
    copyNpzArr(npz, p + "reactionBoLin", out.reactionBoLin, 3);
    copyNpzArr(npz, p + "reactionMoAng", out.reactionMoAng, 3);
    copyNpzArr(npz, p + "reactionMoLin", out.reactionMoLin, 3);
}

inline void loadAggregateStateArrays(const cnpy::npz_t& npz, const std::string& p,
                                     robotics_oracle::OracleAggregateState& out, int numBodies) {
    const std::size_t nq = static_cast<std::size_t>(out.nq);
    const std::size_t nu = static_cast<std::size_t>(out.nu);
    copyNpzArr(npz, p + "q", out.q, nq);
    copyNpzArr(npz, p + "u", out.u, nu);
    copyNpzArr(npz, p + "bodyForceTorque", &out.bodyForceTorque[0][0], static_cast<std::size_t>(numBodies) * 3);
    copyNpzArr(npz, p + "bodyForceForce", &out.bodyForceForce[0][0], static_cast<std::size_t>(numBodies) * 3);
    out.logDetM = scalarFromNpz(npz, p + "logDetM");
    out.totalKE = scalarFromNpz(npz, p + "totalKE");
    out.normUdot = scalarFromNpz(npz, p + "normUdot");
    out.minEigD = scalarFromNpz(npz, p + "minEigD");
    copyNpzArr(npz, p + "reportX_GB_R", out.reportX_GB_R, 9);
    copyNpzArr(npz, p + "reportX_GB_p", out.reportX_GB_p, 3);
    copyNpzArr(npz, p + "reportV_GB_ang", out.reportV_GB_ang, 3);
    copyNpzArr(npz, p + "reportV_GB_lin", out.reportV_GB_lin, 3);
    copyNpzArr(npz, p + "reportA_GB_ang", out.reportA_GB_ang, 3);
    copyNpzArr(npz, p + "reportA_GB_lin", out.reportA_GB_lin, 3);
}

// ---------------------------------------------------------------------------
//  Public entry points -- one per case shape (§4.1/§9), mirroring the
//  generator's writeCase/writeMultiCase/writeAggregateCase exactly.
// ---------------------------------------------------------------------------
inline robotics_oracle::OracleCase loadCase(const std::string& fixturesDir, const std::string& caseName) {
    const JsonValue manifest = loadManifestJson(fixturesDir, caseName);
    checkNoSilentGap(manifest, caseName);
    const cnpy::npz_t npz = cnpy::npz_load(fixturesDir + "/" + caseName + ".npz");

    robotics_oracle::OracleCase c{};
    c.name = internString(caseName);

    const JsonValue& model = manifest.at("model");
    const robotics_oracle::BodyModelSpec bm = bodySpecFromJson(model.at("bodies").arrVal.at(0));
    std::memcpy(c.X_PF_R, bm.X_PF_R, sizeof(c.X_PF_R));
    std::memcpy(c.X_PF_p, bm.X_PF_p, sizeof(c.X_PF_p));
    std::memcpy(c.X_BM_R, bm.X_BM_R, sizeof(c.X_BM_R));
    std::memcpy(c.X_BM_p, bm.X_BM_p, sizeof(c.X_BM_p));
    c.mass = bm.mass;
    std::memcpy(c.com_B, bm.com_B, sizeof(c.com_B));
    std::memcpy(c.unitInertia_B, bm.unitInertia_B, sizeof(c.unitInertia_B));
    c.nq = static_cast<int>(model.at("nq").numVal);
    c.nu = static_cast<int>(model.at("nu").numVal);

    const auto& states = manifest.at("states").arrVal;
    c.numStates = static_cast<int>(states.size());
    for (int i = 0; i < c.numStates; ++i) {
        robotics_oracle::OracleState& s = c.states[i];
        s.label = internString(states[static_cast<std::size_t>(i)].strVal);
        s.nq = c.nq;
        s.nu = c.nu;
        loadStateArrays(npz, "state" + std::to_string(i) + "_", s);
    }
    return c;
}

inline robotics_oracle::OracleMultiCase loadMultiCase(const std::string& fixturesDir,
                                                       const std::string& caseName) {
    const JsonValue manifest = loadManifestJson(fixturesDir, caseName);
    checkNoSilentGap(manifest, caseName);
    const cnpy::npz_t npz = cnpy::npz_load(fixturesDir + "/" + caseName + ".npz");

    robotics_oracle::OracleMultiCase c{};
    c.name = internString(caseName);

    const JsonValue& model = manifest.at("model");
    c.numBodies = static_cast<int>(model.at("numBodies").numVal);
    for (int b = 0; b < c.numBodies; ++b) {
        c.bodies[b] = bodySpecFromJson(model.at("bodies").arrVal.at(static_cast<std::size_t>(b)));
    }

    const auto& states = manifest.at("states").arrVal;
    c.numStates = static_cast<int>(states.size());
    for (int i = 0; i < c.numStates; ++i) {
        robotics_oracle::OracleMultiState& s = c.states[i];
        const std::string sp = "state" + std::to_string(i) + "_";
        s.label = internString(states[static_cast<std::size_t>(i)].strVal);
        s.nq = static_cast<int>(model.at("nq").numVal);
        s.nu = static_cast<int>(model.at("nu").numVal);

        copyNpzArr(npz, sp + "q", s.q, static_cast<std::size_t>(s.nq));
        copyNpzArr(npz, sp + "u", s.u, static_cast<std::size_t>(s.nu));
        copyNpzArr(npz, sp + "bodyForceTorque", &s.bodyForceTorque[0][0], static_cast<std::size_t>(c.numBodies) * 3);
        copyNpzArr(npz, sp + "bodyForceForce", &s.bodyForceForce[0][0], static_cast<std::size_t>(c.numBodies) * 3);
        copyNpzArr(npz, sp + "groundForceTorque", s.groundForceTorque, 3);
        copyNpzArr(npz, sp + "groundForceForce", s.groundForceForce, 3);
        copyNpzArr(npz, sp + "mobilityForce", s.mobilityForce, static_cast<std::size_t>(s.nu));
        for (int b = 0; b < c.numBodies; ++b) {
            const std::string bp = sp + "body" + std::to_string(b) + "_";
            const int bodyNu = nuFromNpz(npz, bp);
            const int bodyNq = nqFromNpz(npz, bp);
            loadBodyOutputArrays(npz, bp, s.body[b], bodyNu, bodyNq);
        }
        copyNpzArr(npz, sp + "Mdense", s.Mdense, static_cast<std::size_t>(s.nu) * static_cast<std::size_t>(s.nu));
        s.logDetM = scalarFromNpz(npz, sp + "logDetM");
    }
    return c;
}

inline robotics_oracle::OracleAggregateCase loadAggregateCase(const std::string& fixturesDir,
                                                               const std::string& caseName) {
    const JsonValue manifest = loadManifestJson(fixturesDir, caseName);
    checkNoSilentGap(manifest, caseName);
    const cnpy::npz_t npz = cnpy::npz_load(fixturesDir + "/" + caseName + ".npz");

    robotics_oracle::OracleAggregateCase c{};
    c.name = internString(caseName);

    const JsonValue& model = manifest.at("model");
    c.numBodies = static_cast<int>(model.at("numBodies").numVal);
    for (int b = 0; b < c.numBodies; ++b) {
        c.bodies[b] = bodySpecFromJson(model.at("bodies").arrVal.at(static_cast<std::size_t>(b)));
    }
    const int reportBody = model.has("reportBody") ? static_cast<int>(model.at("reportBody").numVal) : 0;

    const auto& states = manifest.at("states").arrVal;
    c.numStates = static_cast<int>(states.size());
    for (int i = 0; i < c.numStates; ++i) {
        robotics_oracle::OracleAggregateState& s = c.states[i];
        s.label = internString(states[static_cast<std::size_t>(i)].strVal);
        s.reportBody = reportBody;
        s.nq = static_cast<int>(model.at("nq").numVal);
        s.nu = static_cast<int>(model.at("nu").numVal);
        loadAggregateStateArrays(npz, "state" + std::to_string(i) + "_", s, c.numBodies);
    }
    return c;
}

// §8.3 fuzz batch: same shape as loadAggregateCase, plus the per-body
// udot/A_GB end-products and the rng_seed/resample_count provenance fields
// (§9: "the RNG seed for fuzz cases").
inline robotics_oracle::OracleFuzzCase loadFuzzCase(const std::string& fixturesDir, const std::string& caseName) {
    const JsonValue manifest = loadManifestJson(fixturesDir, caseName);
    checkNoSilentGap(manifest, caseName);
    const cnpy::npz_t npz = cnpy::npz_load(fixturesDir + "/" + caseName + ".npz");

    robotics_oracle::OracleFuzzCase c{};
    c.name = internString(caseName);

    const JsonValue& model = manifest.at("model");
    c.numBodies = static_cast<int>(model.at("numBodies").numVal);
    for (int b = 0; b < c.numBodies; ++b) {
        c.bodies[b] = bodySpecFromJson(model.at("bodies").arrVal.at(static_cast<std::size_t>(b)));
    }

    if (manifest.at("rng_seed").type == JsonValue::Type::Number) {
        c.rngSeed = static_cast<long long>(manifest.at("rng_seed").numVal);
    }
    if (manifest.has("resample_count")) {
        c.resampleCount = static_cast<int>(manifest.at("resample_count").numVal);
    }

    const auto& states = manifest.at("states").arrVal;
    c.numStates = static_cast<int>(states.size());
    if (c.numStates > robotics_oracle::kMaxFuzzStates) {
        throw std::runtime_error("robotics_oracle: fuzz fixture '" + caseName + "' has " +
                                  std::to_string(c.numStates) + " states, exceeding kMaxFuzzStates (" +
                                  std::to_string(robotics_oracle::kMaxFuzzStates) + ")");
    }
    for (int i = 0; i < c.numStates; ++i) {
        robotics_oracle::OracleFuzzState& s = c.states[i];
        const std::string sp = "state" + std::to_string(i) + "_";
        s.label = internString(states[static_cast<std::size_t>(i)].strVal);
        s.nq = static_cast<int>(model.at("nq").numVal);
        s.nu = static_cast<int>(model.at("nu").numVal);

        copyNpzArr(npz, sp + "q", s.q, static_cast<std::size_t>(s.nq));
        copyNpzArr(npz, sp + "u", s.u, static_cast<std::size_t>(s.nu));
        copyNpzArr(npz, sp + "bodyForceTorque", &s.bodyForceTorque[0][0], static_cast<std::size_t>(c.numBodies) * 3);
        copyNpzArr(npz, sp + "bodyForceForce", &s.bodyForceForce[0][0], static_cast<std::size_t>(c.numBodies) * 3);
        s.logDetM = scalarFromNpz(npz, sp + "logDetM");
        s.totalKE = scalarFromNpz(npz, sp + "totalKE");
        s.normUdot = scalarFromNpz(npz, sp + "normUdot");
        s.minEigD = scalarFromNpz(npz, sp + "minEigD");
        copyNpzArr(npz, sp + "udot", s.udot, static_cast<std::size_t>(s.nu));
        copyNpzArr(npz, sp + "A_GB_ang", &s.A_GB_ang[0][0], static_cast<std::size_t>(c.numBodies) * 3);
        copyNpzArr(npz, sp + "A_GB_lin", &s.A_GB_lin[0][0], static_cast<std::size_t>(c.numBodies) * 3);
    }
    return c;
}

} // namespace robotics_oracle_loader
