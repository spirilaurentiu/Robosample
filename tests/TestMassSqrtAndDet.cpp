/* -------------------------------------------------------------------------- *
 *                 SimTK Molmodel Example: Simple Protein                     *
 * -------------------------------------------------------------------------- *
 * This is the first example from the Molmodel User's Guide. It creates a     *
 * small protein (a five-residue peptide), simulates it and generates a live  *
 * animation while it is running.                                             *
 *                                                                            *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"
#include <iostream>
#include <exception>

// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>
// #include <Eigen/LU>

using namespace SimTK;

int trianglize(double **m, int n)
{
    int sign = 1;
    for (int i = 0; i < n; i++) {
        int max = 0;
 
        for (int row = i; row < n; row++)
            if (fabs(m[row][i]) > fabs(m[max][i]))
                max = row;
 
        if (max) {
            sign = -sign;
            double *tmp = m[i];
            m[i] = m[max], m[max] = tmp;
        }
 
        if (!m[i][i]) return 0;
 
        for (int row = i + 1; row < n; row++) {
            double r = m[row][i] / m[i][i];
            if (!r)    continue;
 
            for (int col = i; col < n; col ++)
                m[row][col] -= m[i][col] * r;
        }
    }
    return -sign;
}
 
double det(double *in, int n)
{
    double *m[n];
    m[0] = in;
 
    for (int i = 1; i < n; i++)
        m[i] = m[i - 1] + n;
 
    int sign = trianglize(m, n);
    if (!sign)
        return 0;
 
    double p = 1;
    for (int i = 0; i < n; i++)
        p *= m[i][i];
    return p * sign;
}

// Print matrix nicely
void printMat(SimTK::Matrix m, unsigned int precision)
{
    std::cout << std::fixed << std::setprecision(precision);
    for(int i = 0; i < m.nrow(); i++){
        for(int j = 0; j < m.nrow(); j++){
            std::cout << m[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

SimTK::Matrix matrixSqrt(SimTK::Matrix& DI, int dof)
{
    // DI Cholesky Decomposition: 
    int i, j, k;
    SimTK::Matrix sqrtDI(dof, dof);
    if(dof>1){
      double s1, s2;
      SimTK::Matrix L(dof, dof), De(dof, dof);
      L.setTo(0);
      De.setTo(0);

      for(j=0; j<dof; j++){
        L(j,j) = 1.0;
        // Diagonal
        s1 = 0.0;
        for(k=0; k<dof; k++){
          s1 += (L(j,k) * L(j,k)) * De(k,k);
        }
        De(j,j) = DI(j,j) - s1;
        // Off diagonal
        for(i=(j+1); i<dof; i++){
          s2 = 0.0;
          for(k=0; k<j; k++){
            s2 += L(i,k) * L(j,k) * De(k,k);
          }
          L(i,j) = (1 / De(j,j)) * (DI(i,j) - s2);
        }
      }

      // Square root the diagonal
      for(j=0; j<dof; j++){De(j,j) = sqrt(De(j,j));}

    //   std::cout << "choleksy de\n";
    //   for(int row = 0; row < dof; row++)
    //   {
    //       for(int col = 0; col < dof; col++)
    //       {
    //           std::cout << De(row, col) << ' ';
    //       }

    //       std::cout << '\n';
    //   }

    //   std::cout << '\n';

     // Get LsqrtDe
    sqrtDI = L *  De ;

    //   std::cout << "L = \n";
    //   printMat(L * De * L.transpose(), 3);
    //   std::cout << "\n\n";
    }
    else{
      sqrtDI(0,0) = sqrt(DI(0,0));
    }

    return sqrtDI;
}

SimTK::Matrix lapacc(SimTK::Matrix& DI, int dof)
{
  constexpr SimTK::Real K = 1;

        // dof = 7;

        // SimTK::Real A_rawData[dof * dof] = {1.00, -1.00, 0.00, -0.00, 0.00, 0.00, -0.00
        // , -1.00, 2.00, -1.00, -0.00, -0.00, 0.00, 0.00
        // , 0.00, -1.00, 2.00, -1.00, -0.00, -0.00, 0.00
        // , -0.00, -0.00, -1.00, 2.00, -1.00, 0.00, -0.00
        // , 0.00, -0.00, -0.00, -1.00, 2.00, -1.00, 0.00
        // , 0.00, 0.00, -0.00, 0.00, -1.00, 2.00, -1.00
        // , -0.00, 0.00, 0.00, -0.00, 0.00, -1.00, 1.00};

        // Matrix must be copied because LAPACK will overwrite input data.
        SimTK::Real A_rawData[dof * dof];
        for(int i = 0; i < dof; i++) {
            for(int j = 0; j < dof; j++) {
                A_rawData[(i * dof) + j] = DI[i][j];
            }
        }

        SimTK::Real wr[dof], wi[dof];
        SimTK::Real vl[dof * dof], vr[dof * dof];

        SimTK::Real Out[dof * dof];
        for(auto& i : Out) i = 0;

        SimTK::Real Out1[dof * dof];
        for(auto& i : Out1) i = 0;

        int lwork = -1;
        int info;

        SimTK::Real wkopt;
        dgeev_('V', 'V', dof, A_rawData, dof, &wr[0], &wi[0], vl, dof, vr, dof, &wkopt, lwork, info);

        //std::cout << "info 1 is " << info << std::endl;

        lwork = static_cast<int>(wkopt);
        std::vector<SimTK::Real> work_rawData(lwork);
        dgeev_('V', 'V', dof, A_rawData, dof, &wr[0], &wi[0], vl, dof, vr, dof, &work_rawData[0], lwork, info);

        std::cout << "D is\n";
        for(int i = 0; i < dof; i++)
        {
            for(int j = 0; j < dof; j++)
            {
                if(i == j) std::cout << wr[i] << ' ';
                else std::cout << "0 ";
            }

            std::cout << '\n';
        }
        std::cout << "\n\n";

        //std::cout << "info 2 is " << info << std::endl;

        // std::cout << "wi = ";
        // for(auto i : wi)
        // {
        //     std::cout << i << ' ';
        // }

        //std::cout << std::endl;

        // SimTK::Matrix sqrtD(dof, dof), sqrtDI(dof, dof);
        // sqrtD.setToZero();
        // for(int i = 0; i < dof; i++) {
        //     sqrtD[i][i] = std::sqrt(wr[i]);
        // }
        
        // dgemm_('N', 'T', dof, dof, dof, K, vl, dof, &sqrtD[0][0], dof, K, Out, dof);
        // dgemm_('T', 'T', dof, dof, dof, K, Out, dof, vr, dof, K, Out1, dof);

        // for(int i = 0; i < dof; i++)
        // {
        //     for(int j = 0; j < dof; j++)
        //     {
        //         sqrtDI[i][j] = Out1[(i*dof) + j];
        //     }   
        // }

        // Compute the square root
        SimTK::Matrix VL(dof, dof), VR(dof, dof), sqrtD(dof, dof);
        sqrtD.setToZero();
        VL.setToZero();
        VR.setToZero();

        for(int i = 0; i < dof; i++) {
            for(int j = 0; j < dof; j++) {
                VL[i][j] = vl[(i*3) + j];
                VR[i][j] = vr[(i*3) + j];
            }
        }
        for(int i = 0; i < dof; i++) {
            sqrtD[i][i] = std::sqrt(wr[i]);
        }

        // std::cout << "lapacc sqrtD\n";
        // printMat(sqrtD, 3);
        // std::cout << '\n';

        // std::cout << "lapacc VR\n";
        // printMat(VR, 3);
        // std::cout << '\n';

        //std::cout << "m should be " << std::endl;
        auto sqrt_m = VL.transpose() * sqrtD * VR;

    return sqrt_m;
}

template<typename T>
void printArrayAsMat(T* array, int n, std::string s)
{
    if(!s.empty())
    {
        std::cout << s << '\n';
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            std::cout << array[i * n + j] << ' ';
        }

        std::cout << '\n';
    }

    std::cout << "\n\n";
}

int main() {
try {
    CompoundSystem system; // Extends MolecularMechanicsSystem (Molmodel)
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();
    Protein protein("G");
    protein.assignBiotypes();
    system.adoptCompound(protein);

    for (unsigned int r=0 ; r<protein.getNumBonds(); r++){
        protein.setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
    }

    //system.modelCompounds(); 
    system.addEventReporter(new Visualizer::Reporter(system, 0.020));
    State state = system.realizeTopology();
    LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);

    // Get sqrt(MInv)
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    for (int i=0; i < nu; ++i){
       V[i] = i;
    }

    system.realize(state, SimTK::Stage::Position);
    matter.multiplyBySqrtMInv(state, V, SqrtMInvV);
    std::cout << "V= " << V << std::endl;
    std::cout << "SqrtMInvV=  " << SqrtMInvV << std::endl;
    //==========


    // Get M
    SimTK::Matrix M(nu, nu);
    matter.calcM(state, M);
    //==========

    // Get the inverse M
    SimTK::Matrix MInv;
    matter.calcMInv(state, MInv);

    // Compute expected sqrt
    auto m = matrixSqrt(MInv, M.nrow());
    std::cout << "Expected SqrtMI=" << m * V << std::endl;
    std::cout << "Expected SqrtMI * SqrtMI=" << m * m.transpose() << std::endl;

    // Test sqrt
    SimTK::Matrix test;
    matter.calcMInvSqrt(state, test);

    std::cout << "Test SqrtMI=" << test * V << std::endl;
    std::cout << "Test SqrtMI * SqrtMI=" << test * test.transpose() << std::endl;

    // Get detM
    SimTK::Real detM = 1.0;
    SimTK::Vector DetV(nu);
    Real* newDetM = new Real(1.0);

    matter.calcDetM(state, V, DetV, newDetM);
    std::cout << "newDetM: " << *newDetM << std::endl;
    //==========

    // ---- Verify with a correct algorithm (used to be Eigen) ----------
    // Eigen M determinant - Eigen::MatrixXd
    std::vector<SimTK::Real> EiM(nu * nu);
    std::vector<SimTK::Real> EiD0(6 * 6);
    SimTK::Real JainDetM;
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiM[i * nu + j] = M(i, j);
        }
    }
    SimTK::Real EiDetM = det(&EiM[0], nu);
    std::cout << "EiDetM= " << EiDetM << std::endl;

    // Simulate for a bit
    ts.initialize(state);
    ts.stepTo(0.12); // 0.12ps

    return 0;
} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" 
              << std::endl;
    return 1;
}

}
