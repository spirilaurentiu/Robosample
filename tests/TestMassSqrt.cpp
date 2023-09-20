#include "Context.hpp"
#include <random>
#include <vector>

constexpr int DOF = 10;
constexpr auto DECIMAL_PLACES = 5;

bool areEqual(SimTK::Real a, SimTK::Real b) {
    SimTK::Real epsilon = std::pow(10, -DECIMAL_PLACES);
    return std::abs(a - b) < epsilon;
}

int main() {
    // // Fill all state 19k bits of Mersenne Twister
	// std::vector<uint32_t> RandomData(624);
	// std::random_device Source;
	// std::generate(RandomData.begin(), RandomData.end(), std::ref(Source));
	// std::seed_seq SeedSeq(RandomData.begin(), RandomData.end());

    // // randm source
	// std::mt19937 engine(SeedSeq);
	// std::uniform_real_distribution<> dist(0.0, 1.0);

    // // create an empty mass matrix (diagonal matrix)
    // std::vector<SimTK::Real> A(DOF * DOF, 0);
    // for(int i = 0; i < DOF; i++)
    //     A[i * DOF + i] = dist(engine);

    // int M = 0, info = 0;
    // SimTK::Real W[DOF], Z[DOF * DOF], WORK[26 * DOF];
    // int ISUPPZ[2 * DOF], IWORK[10 * DOF];

    // dsyevr_('V', 'A', 'L', DOF, &A[0], DOF, 0, 0, 0, 0, dlamch_('S'), M, W, Z, DOF, ISUPPZ, WORK, 26 * DOF, IWORK, 10 * DOF, info);

    // // std::reverse(W, W + DOF);
    // // std::reverse(Z, Z + DOF * DOF);

    // SimTK::Mat<DOF, DOF> D, U;
    // for(int i = 0; i < DOF; i++)
    // {
    //     for(int j = 0; j < DOF; j++)
    //     {
    //         if(i == j) D[i][j] = sqrt(W[i]);
    //         else D[i][j] = 0;

    //         U[i][j] = Z[i * DOF + j];
    //     }
    // }
    
    // auto sqrtDI = U.transpose() * D * U;
    // auto test = sqrtDI * sqrtDI;

    // // for(int i = 0; i < DOF; i++)
    // // {
    // //     for(int j = 0; j < DOF; j++)
    // //         std::cout << A[i * DOF + j] << " ";
    // //     std::cout << std::endl;
    // // }
    // // std::cout << sqrtDI << std::endl;
    // // std::cout << sqrtDI * sqrtDI << std::endl;

    // for (std::size_t i = 0; i < DOF; i++) {
    //     for (std::size_t j = 0; j < DOF; j++) {
    //         if (!areEqual(A[i * DOF + j], test[i][j])) {
    //             std::cout << "Expected " << A[i] << ", got " << test[i][j];
    //             return 1;
    //         }
    //     }
    // }

    // return 0;


    Context c;
    if (!c.initializeFromFile("inp.2but.mono")) {
        std::cout << "Failed to initialize from inp.2but.mono" << std::endl;
        return -1;
    }

    const auto state = c.getWorld(0)->integ->getState();
    const auto nu = state.getNU();

    // c.getWorld(0)->realizeTopology();
    // c.getWorld(0)->realize();


    SimTK::Matrix M(nu, nu);
    c.getWorld(0)->matter->calcM(state, M);
    
    std::cout << "got here" << std::endl;
    SimTK::Matrix MInv(nu, nu);
	c.getWorld(0)->matter->calcMInv(state, MInv);
    std::cout << "got here" << std::endl;
    SimTK::Matrix MInvSq = MInv * MInv;

    for (std::size_t i = 0; i < DOF; i++) {
        for (std::size_t j = 0; j < DOF; j++) {
            if (!areEqual(M[i][j], MInvSq[i][j])) {
                std::cout << "Expected " << M[i][j] << ", got " << MInvSq[i][j];
                return 1;
            }
        }
    }

	// // Compute expected sqrt
	// auto m = matrixSqrt(MInv, M.nrow());
	// std::cout << "Expected sqrt(M) * V =" << m * V << std::endl;
	// std::cout << "Expected sqrt(M) * sqrt(M)_T=" << m * m.transpose() << std::endl;

	// // Test sqrt
	// SimTK::Matrix test;
	// matter.calcMInvSqrt(state, test);

	// std::cout << "Test SqrtMI=" << test * V << std::endl;
	// std::cout << "Test SqrtMI * SqrtMI=" << test * test.transpose() << std::endl;

	// // Get detM
	// SimTK::Real detM = 1.0;
	// SimTK::Vector DetV(nu);
	// Real* newDetM = new Real(1.0);

	// matter.calcDetM(state, V, DetV, newDetM);
	// std::cout << "newDetM: " << *newDetM << std::endl;
	// //==========

	// // ---- Verify with a correct algorithm (used to be Eigen) ----------
	// // Eigen M determinant - Eigen::MatrixXd
	// std::vector<SimTK::Real> EiM(nu * nu);
	// std::vector<SimTK::Real> EiD0(6 * 6);
	// SimTK::Real JainDetM;
	// for(int i=0; i<nu; i++){
	// 	for(int j=0; j<nu; j++){
	// 		EiM[i * nu + j] = M(i, j);
	// 	}
	// }
	// SimTK::Real EiDetM = det(&EiM[0], nu);
	// std::cout << "EiDetM= " << EiDetM << std::endl;

	// // Simulate for a bit
	// ts.initialize(state);
	// ts.stepTo(0.12); // 0.12ps

	return 0;
}