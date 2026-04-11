#include "bgeneral.hpp"

namespace {
// Define canonical X axis for F1 frame as a constant to avoid repeated construction
const SimTK::Vec3 f1Axis_X{1.0, 0.0, 0.0};
const SimTK::UnitVec3 f1UnitAxis_X(f1Axis_X);

// Define axes in F2 space
const SimTK::Vec3 f2Axis_Y{0.0, 1.0, 0.0};
const SimTK::Vec3 f2Origin{0.0, 0.0, 0.0};

const SimTK::UnitVec3 kUnitX{1.0, 0.0, 0.0};
} // namespace

/**
 * @brief Computes the Log-Sum-Exp (LSE) of two values.
 * This is a numerically stable way to calculate log(exp(left) + exp(right)).
 * Fixes: Implicit bool conversion and variable naming requirements.
 */
[[nodiscard]] auto calculateLogSumExp2(SimTK::Real leftValue, SimTK::Real rightValue) -> SimTK::Real {
    const SimTK::Real maxValue = std::max(leftValue, rightValue);

    // Handle infinities explicitly
    const int infinityCount = std::isinf(maxValue);
    if (infinityCount == -1 || infinityCount == 1) {
        return maxValue;
    }

    const SimTK::Real exponentSum = std::exp(leftValue - maxValue) + std::exp(rightValue - maxValue);

    return maxValue + std::log(exponentSum);
}

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
[[nodiscard]] auto alignFlipAndTranslateFrameAlongXAxis(const SimTK::Transform& gTransform_F1,
                                                        const SimTK::Vec3& gPoint_v1) -> SimTK::Transform {
    // Vector from F1 origin to v1, expressed in G
    const SimTK::Vec3 gVec_F1ToV1 = (gPoint_v1 - gTransform_F1.p());

    // Re-express that vector in F1 coordinates
    const SimTK::Vec3 f1Vec_F1ToV1 = ~(gTransform_F1.R()) * gVec_F1ToV1;

    // Normalize safely via UnitVec3
    const SimTK::UnitVec3 f1Dir_ToV1(f1Vec_F1ToV1);

    // Compute rotation angle between vectors
    const SimTK::Real cosTheta = SimTK::dot(f1Dir_ToV1, f1UnitAxis_X);

    // Numerical safety
    const SimTK::Real clampedCosTheta = std::clamp(cosTheta, SimTK::Real(-1.0), SimTK::Real(1.0));
    const SimTK::Angle theta = std::acos(clampedCosTheta);

    // Rotation axis (handle degeneracy implicitly via UnitVec3)
    const SimTK::UnitVec3 rotAxis_F1(SimTK::cross(f1Dir_ToV1, f1UnitAxis_X));

    // First transform: align F1 direction to X axis and translate to v1
    const SimTK::Transform f1_X_f2(SimTK::Rotation((-theta) + SimTK::Pi, rotAxis_F1), f1Vec_F1ToV1);

    // Invert to express quantities in F2
    const SimTK::Transform f2_X_f1 = ~f1_X_f2;

    // Vector from F2 origin to F1 origin expressed in F2
    const SimTK::Vec3& f2Vec_F1 = f2_X_f1.p();

    // F1 X axis expressed in F2
    const SimTK::Vec3 f2Axis_F1X = f2_X_f1.R() * kUnitX;

    // Resolve remaining rotational DOF using dihedral angle
    const SimTK::Angle dihedralAngle = calculateDihedralInRad(f2Axis_Y, f2Origin, f2Vec_F1, f2Axis_F1X);

    // Final correction rotation about X axis
    const SimTK::Transform f2_X_f3(SimTK::Rotation(dihedralAngle, kUnitX));

    // Compose final transform
    const SimTK::Transform f1_X_f3 = f1_X_f2 * f2_X_f3;

    return f1_X_f3;
}

void PrintMat33(const SimTK::Mat33& matrix, int decimalPlaces, const std::string& header) {
    std::cout << header << "\n";
    std::cout << std::setw(6 + decimalPlaces) << std::fixed << std::setprecision(decimalPlaces);
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            std::cout << matrix(i, k) << " ";
        }
        std::cout << "\n";
    }
}

/**
 * @brief Calculates the angle in radians between vectors (pos1-pos0) and (pos2-pos0).
 * Uses std::clamp to prevent NaN results from floating-point drift.
 */
[[nodiscard]] auto calculateAngleInRad(const SimTK::Vec3& pos0,
                                       const SimTK::Vec3& pos1,
                                       const SimTK::Vec3& pos2) -> SimTK::Real {
    const SimTK::Vec3 vec10 = pos1 - pos0;
    const SimTK::Vec3 vec20 = pos2 - pos0;

    const SimTK::Real magSquared = vec10.normSqr() * vec20.normSqr();

    // Handle edge case where points are coincident to avoid division by zero
    if (magSquared <= SimTK::TinyReal) {
        return 0.0;
    }

    // Dot product divided by the product of magnitudes
    const SimTK::Real cosTheta = SimTK::dot(vec10, vec20) / std::sqrt(magSquared);

    // C++17 std::clamp prevents acos(1.0000000000001) -> NaN
    return std::acos(std::clamp<SimTK::Real>(cosTheta, -1.0, 1.0));
}

/**
 * @brief Calculates the dihedral angle in radians between four positions.
 * * Uses the Praxeolitic formula for high numerical stability.
 * Formula: atan2( |b1| * b0 * (b1 x b2), (b0 x b1) * (b1 x b2) )
 */
[[nodiscard]] auto calculateDihedralInRad(const SimTK::Vec3& pos0,
                                          const SimTK::Vec3& pos1,
                                          const SimTK::Vec3& pos2,
                                          const SimTK::Vec3& pos3) -> SimTK::Real {
    // Vector segments between the four points
    const SimTK::Vec3 vec10 = pos1 - pos0;
    const SimTK::Vec3 vec21 = pos2 - pos1;
    const SimTK::Vec3 vec32 = pos3 - pos2;

    // Normals to the two planes defined by (pos0, pos1, pos2) and (pos1, pos2, pos3)
    // SimTK '%' operator is the cross product
    const SimTK::Vec3 normalLeft = vec10 % vec21;
    const SimTK::Vec3 normalRight = vec21 % vec32;

    // psin = |vec21| * (normalLeft * vec32)
    const SimTK::Real psin = vec21.norm() * SimTK::dot(normalLeft, vec32);

    // pcos = normalLeft * normalRight
    const SimTK::Real pcos = SimTK::dot(normalLeft, normalRight);
    return std::atan2(psin, pcos);
}


/** * @brief Squared magnitude of a vector.
 * Uses std::transform_reduce for better compiler optimization/vectorization.
 */
[[nodiscard]] auto calculateMagSq(const std::vector<SimTK::Real>& vec) -> SimTK::Real {
    // std::transform_reduce is the modern, more optimizable C++17 version of inner_product
    return std::transform_reduce(std::execution::unseq, // Allow vectorization
                                 vec.begin(),
                                 vec.end(),
                                 SimTK::Real{0.0},
                                 std::plus<>(),
                                 [](const SimTK::Real val) {
                                     return val * val;
                                 });
}

/**
 * @brief Normalizes a vector in-place.
 * @return A reference to the modified vector.
 */
void normalizeInPlace(std::vector<SimTK::Real>& inputVector) {
    // C++17 transform_reduce is faster and more precise than a manual loop
    const SimTK::Real normSquared = std::transform_reduce(std::execution::unseq,
                                                          inputVector.begin(),
                                                          inputVector.end(),
                                                          0.0,
                                                          std::plus<>(),
                                                          [](const SimTK::Real value) -> SimTK::Real {
                                                              return value * value;
                                                          });

    const SimTK::Real magnitude = std::sqrt(normSquared);

    if (magnitude > 1e-15) { // Check against a small epsilon
        const SimTK::Real inverseMagnitude = 1.0 / magnitude;
        std::for_each(std::execution::unseq,
                      inputVector.begin(),
                      inputVector.end(),
                      [inverseMagnitude](SimTK::Real& value) {
                          value *= inverseMagnitude;
                      });
    }
}

/**
 * @brief Multiplies a source vector by a scalar and stores the result in the destination.
 * Clang-Tidy: Avoids swapping by clearly distinguishing 'source' and 'destination'.
 */
void multiplyByScalar(const std::vector<SimTK::Real>& sourceVector,
                      SimTK::Real scalarValue,
                      std::vector<SimTK::Real>& destinationVector) {
    // Ensure the destination is the correct size before operating
    if (destinationVector.size() != sourceVector.size()) {
        destinationVector.resize(sourceVector.size());
    }

    std::transform(std::execution::unseq,
                   sourceVector.begin(),
                   sourceVector.end(),
                   destinationVector.begin(),
                   [scalarValue](const SimTK::Real value) -> SimTK::Real {
                       return value * scalarValue;
                   });
}

/**
 * @brief Numerically stable computation of log(sin^2(pitch))
 * using a Taylor expansion near pitch = 0 for smoothness.
 *
 * This avoids log(0) and ensures continuous derivatives,
 * useful for energy/gradient computations.
 */
[[nodiscard]] auto safeLogSineSqr(SimTK::Real pitch) -> SimTK::Real {
    // Threshold for small angles to prevent log(0)
    constexpr SimTK::Real deltaThreshold = 1e-6;
    const SimTK::Real absolutePitch = std::abs(pitch);

    if (absolutePitch < deltaThreshold) {
        // Use series expansion: log(sin^2(x)) ≈ 2log|x| - (x^2 / 3)
        // Explicit parentheses used to ensure operation order and clarity.
        const SimTK::Real deltaSquared = (deltaThreshold * deltaThreshold);
        const SimTK::Real pitchSquared = (pitch * pitch);

        // Taylor-based smoothing near zero
        return (2.0 * std::log(deltaThreshold)) - ((pitchSquared - deltaSquared) / (3.0 * deltaSquared));
    }

    // No 'else' after return (Clang-Tidy: readability-else-after-return)
    const SimTK::Real sineValue = std::sin(pitch);

    // log(sin^2(x)) is equivalent to 2 * log(|sin(x)|)
    return 2.0 * std::log(std::abs(sineValue));
}
