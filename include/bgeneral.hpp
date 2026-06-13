#pragma once

// Standard Input/Output
#include <cstdio>
#include <iomanip>
#include <iostream>

// Containers
#include <array>
#include <bitset>
#include <deque>
#include <forward_list>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Algorithms & Numerics
#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <random>
#include <valarray>

// Strings & Text
#include <cctype>
#include <cstring>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>

// Utilities & Modern C++ Features
#include <any>
#include <chrono>
#include <execution>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

// System, Filesystem & Resources
#include <filesystem>
#include <fstream>
#include <sys/resource.h>
#include <system_error>

// Multithreading
#include <atomic>
#include <condition_variable>
#include <future>
#include <mutex>
#include <thread>

// Debugging & Limits
#include <cassert>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <stdexcept>

// #ifndef __DRILLING__
// #define __DRILLING__
// #endif

#include "Molmodel.h"
#include "SimTKcommon.h"
#include "Simbody.h"
#include "pcg_random.hpp"

enum class ReplicaMixingScheme : std::uint8_t {
    All = 0,
    Neighboring = 1
};

enum class RUN_TYPE : std::uint8_t {
    Default = 0,
    REMC,      // Replica Exchange Monte Carlo
    RENEMC,    // Replica Exchange Non-Equilibrium Monte Carlo
    RENE,      // Replica Exchange Non-Equilibrium
    REBASONTOP // Replica
};

enum class TopologyRangeType : std::uint8_t {
    Atom = 0,
    Bond,
    Angle,
    PeriodicTorsion,
    ImproperHarmonicTorsion,
    NofTopologyRangeTypes
};

/*
 * Sampling
 */
enum class AcceptRejectMode : std::uint8_t {
    AlwaysAccept = 0,   // MD
    MetropolisHastings, // MCMC
};

enum struct ThermostatName : std::uint8_t { // Thermostats
    None = 0,
    Andersen,
    Berendsen,
    Langevin,
    NoseHoover
};

/*
 * Simulation
 */
enum struct IntegratorType : std::uint8_t { // Integrators
    Empty = 0,
    Verlet,
    Euler,
    Euler2,
    CPodes,
    RungeKutta,
    RungeKutta2,
    RungeKutta3,
    RungeKuttaFeldberg,
    BendStretch,
    OpenMMVelocityVerlet,
    BoundWalk,
    BoundHMC,
    StationsTask,
    NofIntegrators
};

enum struct PositionsPerturbMethod : std::uint8_t {
    Empty = 0,
    BendStretch1,
    BendStretch2,
    BendStretch3,
    BendStretch4,
    BendStretch5,
    BendStretch6,
    NofPositionsPerturbMethod
};

enum struct VelocitiesPerturbMethod : std::uint8_t {
    ToTemperature = 0,
    ToZero,
    NofVelocitiesPerturbMethod
};

enum struct ForcesPerturbMethod : std::uint8_t {
    Empty = 0,
    NotImplemented,
    NofForcesPerturbMethod
};


enum struct NMAOptions : std::uint8_t {
    Empty = 0,
    NMAAltSign,
    NMABernoulli,
    NMAGauss,
    NMAGaussScale,
    NMALenpert,
    NMAFinal,
    NofNMAOptions
};

/*
 * The type of distribution to draw a random number from.
 */
enum struct GmolRandDistributionType : std::uint8_t {
    Uniform = 0,
    Normal
};

// Samplers
enum struct SamplerName : std::uint8_t {
    Empty = 0,
    MC,
    HMC,
    LAHMC
};

enum struct JointType : std::uint8_t {
    Linear = 0,
    Angular180,
    Angular360,
    QuaternionA,
    QuaternionB,
    QuaternionC,
    QuaternionD
};

struct BondFlexibility {
    BondFlexibility() = default;

    int globalIndex1 = -1;
    int globalIndex2 = -1;
    std::string uniqueAtomName1;
    std::string uniqueAtomName2;
    SimTK::BondMobility::Mobility mobility = SimTK::BondMobility::Default;
};

/**
 * @brief Computes the Log-Sum-Exp (LSE) of two values.
 * This is a numerically stable way to calculate log(exp(left) + exp(right)).
 * Fixes: Implicit bool conversion and variable naming requirements.
 */
[[nodiscard]] auto calculateLogSumExp2(SimTK::Real leftValue, SimTK::Real rightValue) -> SimTK::Real;

/**
 * @brief Constructs a new frame F_out with:
 *        1) Origin translated to G_v1
 *        2) X-axis aligned with the vector from F1 origin to G_v1 (expressed in G)
 *        3) Orientation constrained so that F1’s local direction toward v1 is aligned with +X
 *        4) Secondary rotation about X to align Y-axis using a dihedral constraint
 *
 * @param gTransform_F1 Transform from frame F1 to global frame G (G_X_F1)
 * @param gPoint_v1     Point v1 expressed in frame G
 * @return              Transform from F1 to the new frame F_out
 *
 * @details
 * The function constructs a local frame F_out attached at G_v1 such that:
 *  - Its X-axis points from F1 origin toward v1 (in G space)
 *  - Its orientation is fully defined by resolving the remaining rotational
 *    degree of freedom using a dihedral constraint involving the Y-axis
 *
 * This is a rigid-body alignment operation combining translation + rotation.
 */
[[nodiscard]] auto alignFlipAndTranslateFrameAlongXAxis(const SimTK::Transform& gTransform_F1, const SimTK::Vec3& gPoint_v1) -> SimTK::Transform;

void PrintMat33(const SimTK::Mat33& matrix, int decimalPlaces, const std::string& header = "unknown");

void PrintTransform(const SimTK::Transform& transform, int decimalPlaces, const std::string& header = "unknown");

/**
 * @brief Calculates the angle in radians between vectors (pos1-pos0) and (pos2-pos0).
 * Uses std::clamp to prevent NaN results from floating-point drift.
 */
[[nodiscard]] auto calculateAngleInRad(const SimTK::Vec3& pos0, const SimTK::Vec3& pos1, const SimTK::Vec3& pos2) -> SimTK::Real;

/**
 * @brief Calculates the dihedral angle in radians between four positions.
 * * Uses the Praxeolitic formula for high numerical stability.
 * Formula: atan2( |b1| * b0 · (b1 × b2), (b0 × b1) · (b1 × b2) )
 */
[[nodiscard]] auto calculateDihedralInRad(const SimTK::Vec3& pos0, const SimTK::Vec3& pos1, const SimTK::Vec3& pos2, const SimTK::Vec3& pos3) -> SimTK::Real;

/** * @brief Squared magnitude of a vector.
 * Uses std::transform_reduce for better compiler optimization/vectorization.
 */
[[nodiscard]] auto calculateMagSq(const std::vector<SimTK::Real>& vec) -> SimTK::Real;

/**
 * @brief Normalizes a vector in-place.
 * @return A reference to the modified vector.
 */
void normalizeInPlace(std::vector<SimTK::Real>& inputVector);

/**
 * @brief Multiplies a source vector by a scalar and stores the result in the destination.
 * Clang-Tidy: Avoids swapping by clearly distinguishing 'source' and 'destination'.
 */
void multiplyByScalar(const std::vector<SimTK::Real>& sourceVector, SimTK::Real scalarValue, std::vector<SimTK::Real>& destinationVector);

/**
 * @brief Numerically stable computation of log(sin^2(pitch))
 * using a Taylor expansion near pitch = 0 for smoothness.
 *
 * This avoids log(0) and ensures continuous derivatives,
 * useful for energy/gradient computations.
 */
[[nodiscard]] auto safeLogSineSqr(SimTK::Real pitch) -> SimTK::Real;
