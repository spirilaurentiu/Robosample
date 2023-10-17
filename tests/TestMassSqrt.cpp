#include "Context.hpp"
#include <random>
#include <vector>

constexpr auto DECIMAL_PLACES = 5;

bool areEqual(SimTK::Real a, SimTK::Real b) {
    SimTK::Real epsilon = std::pow(10, -DECIMAL_PLACES);
    return std::abs(a - b) < epsilon;
}

int main() {
    // this has cartesian coordinates (one fully flexible world)
    Context c;
    if (!c.initializeFromFile("inp.2but.mono")) {
        std::cout << "Failed to initialize from inp.2but.mono" << std::endl;
        return -1;
    }

    const auto state = c.getWorld(0)->integ->getState();
    const auto nu = state.getNU();

    c.getWorld(0)->compoundSystem->realize(state);

    // get the original matrices
    SimTK::Matrix M(nu, nu), MInv(nu, nu), MInvSqrt(nu, nu);
    c.getWorld(0)->matter->calcM(state, M);
    c.getWorld(0)->matter->calcMInv(state, MInv);
    c.getWorld(0)->matter->calcMInvSqrt(state, MInvSqrt);

    // check that the inverse multiplied twice yields the identity matrix
    SimTK::Matrix testMInvSq0 = M * MInv;
    SimTK::Matrix testMInvSq1 = MInv * M;
    for (int i = 0; i < testMInvSq0.nrow(); i++) {
        for (int j = 0; i < testMInvSq0.ncol(); i++) {
            if (i == j) {
                if (!areEqual(testMInvSq0(i, j), 1) || !areEqual(testMInvSq1(i, j), 1) || !areEqual(testMInvSq0(i, j), testMInvSq1(i, j))) {
                    std::cout << "Inverse mass matrix does not yield the identity matrix when squared." << std::endl;
                    std::cout << "Expected 1, got testMInvSq0 = " << testMInvSq0(i, j) << ", testMInvSq1 = " << testMInvSq1(i, j) << std::endl;
                    return 1;
                }
            }
            else {
                if (!areEqual(testMInvSq0(i, j), 0) || !areEqual(testMInvSq1(i, j), 0) || !areEqual(testMInvSq0(i, j), testMInvSq1(i, j))) {
                    std::cout << "Inverse mass matrix does not yield the identity matrix when squared." << std::endl;
                    std::cout << "Expected 0, got testMInvSq0 = " << testMInvSq0(i, j) << ", testMInvSq1 = " << testMInvSq1(i, j) << std::endl;
                    return 1;
                }
            }
        }
    }

    // check the inverse sqrt of mass matrix
    SimTK::Matrix testMInvSqrtSq = MInvSqrt * MInvSqrt.transpose();
    for (int i = 0; i < testMInvSqrtSq.nrow(); i++) {
        for (int j = 0; i < testMInvSqrtSq.ncol(); i++) {
            if (!areEqual(testMInvSqrtSq(i, j), MInv(i, j))) {
                std::cout << "Inverse square root of matrix does not yield the inverse mass matrix when squared." << std::endl;
                std::cout << "Expected 1, got " << testMInvSqrtSq(i, j) << std::endl;
                return 1;
            }
        }
    }

	return 0;
}