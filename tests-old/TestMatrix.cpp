/* -------------------------------------------------------------------------- *
 *                      Simbody(tm) Example: Long Pendulum                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/*                      Simbody precision
This example shows how to build a linked chain of bodies programmatically,
simulate it, and produce a simple animation while it is simulating. */

#include "Simbody.h"
#include <iostream>

using namespace SimTK;

// Print matrix nicely
template <int M, class E, int CS, int RS> inline
void printMat(Mat<M,M,E,CS,RS>m, unsigned int precision)
//void printMat(Mat<3, 3, Real, 3, 1> m, unsigned int precision)
{
    std::cout << std::fixed << std::setprecision(precision);
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            std::cout << m[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Round matrix entries
void roundRealMatUseStr(Mat<3, 3, Real, 3, 1> &m, unsigned int precision)
{
    std::cout << std::fixed << std::setprecision(precision);
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            m[i][j] = std::stod(std::to_string(m[i][j]));
        }
    }
}

// Matrix square root
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;
Mat<3, 3, Complex, 3, 1> generalMatrixSqrt(Mat<3, 3, fComplex, 3, 1> &m) {

    Mat<3, 3, Complex, 3, 1> sqrt_m;
    int M = 3;

    // Perform Eigenvalue / Eigenvector decomposition
    char jobvl = 'V';    // Compute Eigenvectors too
    char jobvr = 'V';    // Store lower triangle of A

    int N = M;      // Order of the matrix A
    int lda = N;        // Leading dimension 0f A

    //Raw* A_rawData = reinterpret_cast<Raw*>(&m(0,0));
    fComplex *A_rawData = new fComplex[lda*N]; //

    for(int i = 0; i < M; i++) {
        for(int j = 0; j < M; j++) {
            A_rawData[(i*3) + j] = m[i][j];
        }
    }

    fComplex *w_eigenValues = new fComplex[N];
    int lwork = -1;
    fComplex* work_rawData; // = new fComplex[lwork]; // 2*N

    int ldvl = N;
    fComplex *vl = new fComplex[ldvl*N];

    int ldvr = N;
    fComplex *vr = new fComplex[ldvr*N];

    float *rwork = new float[2*N];
    int info;

    //cheev_(jobz, uplo, N, A_rawData, lda, w_eigenValues, work_rawData, lwork, rwork, info, 1, 1);
    fcomplex wkopt;
    cgeev_(jobvl, jobvr, N, A_rawData, lda, w_eigenValues, vl, ldvl, vr, ldvr,
           reinterpret_cast<complex<float> *>(&wkopt), lwork, rwork, info, 1, 1);
    //cgeev_(jobvl, jobvr, N, A_rawData, lda, w_eigenValues, vl, ldvl, vr, ldvr, &wkopt, lwork, rwork, info, 1, 1);
    lwork = (int)wkopt.re;
    work_rawData = new fComplex[lwork];
    cgeev_(jobvl, jobvr, N, A_rawData, lda, w_eigenValues, vl, ldvl, vr, ldvr, work_rawData, lwork, rwork, info, 1, 1);

    std::cout << "A_rawData cheev_" << std::endl;
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < M; j++) {
            std::cout << A_rawData[(i*3) + j] << " " ;
        }
        std::cout << std::endl;
    }

    std::cout << "info " << info << std::endl;
    std::cout << "w_eigenValues " << std::endl;
    for(int i = 0; i < N; i++){
        std::cout << w_eigenValues[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "vl" << std::endl;
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < M; j++) {
            std::cout << vl[(i*3) + j] << " " ;
        }
        std::cout << std::endl;
    }

    std::cout << "vr" << std::endl;
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < M; j++) {
            std::cout << vr[(i*3) + j] << " " ;
        }
        std::cout << std::endl;
    }

    // Compute the square root
    Mat<3, 3, Complex, 3, 1> VL(0), VR(0), sqrtD(0);
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < M; j++) {
            VL[i][j] = vl[(i*3) + j];
            VR[i][j] = vr[(i*3) + j];
        }
    }
    for(int i = 0; i < M; i++) {
        sqrtD[i][i] = std::sqrt(w_eigenValues[i]);
    }

    std::cout << "m should be " << std::endl;
    sqrt_m = VL.transpose() * sqrtD * VR;
    std::cout << sqrt_m << std::endl;

    Vec3 oneVector(1);
    std::cout << "sqrt_M * oneVector " << sqrt_m * oneVector << std::endl;


    return sqrt_m;
}


template <typename someType>
someType multiplyAnything(someType x, someType y)
{
    return x * y;
}

//
int main() {

    /***********************
     * Inverse of a real matrix
     ************************/
     /**/
    Mat<3, 3, Real, 3, 1> Mreal(1);

/*    std::cout << "Version 0" << std::endl;
    Mreal[0][0] = 5.109256853263803854758862144081; Mreal[0][1] =  6.631074837371017771658898709575; Mreal[0][2] = 3.410527705547373056305104910280;
    Mreal[1][0] = 6.631074837371017771658898709575; Mreal[1][1] = 11.951585331463572714483234449290; Mreal[1][2] = -2.078098713674341624368935299572;
    Mreal[2][0] = 3.410527705547373500394314760342; Mreal[2][1] = -2.078098713674342512547354999697; Mreal[2][2] = 14.923205855272648534537438536063;
    std::cout << "M" << std::endl; printMat(Mreal, 30);
    Mat<3, 3> invM = lapackInverse(Mreal);
    std::cout << "M * invM" << std::endl;
    printMat(Mreal * invM, 2);
    std::cout << "=============================================" << std::endl;

    std::cout << "Version 1" << std::endl;
    std::cout << "M rounded" << std::endl;
    roundRealMatUseStr(Mreal, 5);
    printMat(Mreal, 5);
    invM = lapackInverse(Mreal);
    std::cout << "M * invM" << std::endl; printMat(Mreal * invM, 2);
    std::cout << "=============================================" << std::endl;*/

/*
    std::cout << "Version 2" << std::endl;
    Mreal[0][0] = 5.109256; Mreal[0][1] =  6.631074; Mreal[0][2] = 3.410527;
    Mreal[1][0] = 6.631074; Mreal[1][1] = 11.951585; Mreal[1][2] = -2.078098;
    Mreal[2][0] = 3.410527; Mreal[2][1] = -2.078098; Mreal[2][2] = 14.923205;
    std::cout << "M" << std::endl; printMat(Mreal, 30);
    invM = lapackInverse(Mreal);

    //std::cout << "invM" << invM << std::endl;
    std::cout << "M * invM" << std::endl; printMat(Mreal * invM, 2);
    std::cout << "=============================================" << std::endl;
    */

    /*
     * Complex numbers in Simbody
     */
/*    fComplex x(1,2), y(3, 4);
    std::cout << "Conversion from Complex to std::complex "
    << std::complex<float>(x) << " " << std::complex<float>(y) << std::endl;

    Complex xy;
    xy = multiplyAnything(x, y);
    std::cout << "multiplication of arbitrary types: " << xy << std::endl;*/

    /***********************
     * Sqrt of any real matrix can be a complex matrix
     ************************/

    Mat<3, 3, fComplex, 3, 1> Mcomplex(0);
    Mcomplex[0][0] = -5078669.6653899600; Mcomplex[0][1] =  720530.3913155640; Mcomplex[0][2] = 1131653.7208911600;
    Mcomplex[1][0] =   720530.3913155640; Mcomplex[1][1] = -102222.1297841450; Mcomplex[1][2] = -160551.9930878970;
    Mcomplex[2][0] =  1131653.7208911600; Mcomplex[2][1] = -160551.9930878970; Mcomplex[2][2] = -252158.1972449670;
    std::cout << "Mcomplex" << std::endl; printMat(Mcomplex, 5);

    Mat<3, 3, Complex, 3, 1> sqrtM;
    sqrtM = generalMatrixSqrt(Mcomplex);
    //std::cout << "sqrtM" << std::endl; printMat(sqrtM, 5);

    SpatialMat SpatialM(1);
    std::cout << "Spatial matrix " << std::endl;
    std::cout << SpatialM << std::endl;
    std::cout << SpatialM[0] << std::endl;
    std::cout << SpatialM[0][0] << std::endl;
    std::cout << det(SpatialM[0][0]) << std::endl;

    Rotation R;
    R.setRotationFromAngleAboutX(SimTK_DEGREE_TO_RADIAN * 30);
    std::cout << "R: " << R << std::endl;

    // END
    return 0;
}
