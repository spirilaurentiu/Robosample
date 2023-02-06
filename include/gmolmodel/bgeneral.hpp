#ifndef BGENERAL_H_
#define BGENERAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <list>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <assert.h>
#include <cmath>
#include <functional>

// Molmodel specific headers

#include "Molmodel.h"
#include "mol.h"

#include "SimTKcommon.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/GrinPointer.h"
#include "molmodel/internal/units.h"

#include "Robo.hpp"

//#ifndef DEBUG
//#define DEBUG 1
//#endif


#ifndef TARGET_TYPE
#define TARGET_TYPE double
#endif

using namespace std;


// Smart trace
template <typename T>
void trace_impl(std::ostream &os, T &&t)
{
    os << t;
}

template <typename T, typename... Args>
void trace_impl(std::ostream &os, T &&t, Args &&... args)
{
    os << t << ' ';
    trace_impl(os, std::forward<Args>(args)...);
}

template <typename... Args>
void trace(Args &&... args)
{
    std::cout << "[TRACE] ";
    trace_impl(std::cout, std::forward<Args>(args)...);
    std::cout << std::endl;
}

// ----------

#define sq(x)		((x)*(x))

//#ifndef sqr
//#define sqr(x) ((x)*(x))
//#endif

#ifndef ANG_360_TO_180
#define ANG_360_TO_180(x) (((x)>180) ? ((x)-360) : (x))
#endif

// Check versus Velocity Verlet in cart coords
enum struct VELOCITY_VERLET : int {
  SERV_QX, SERV_QY, SERV_QZ,
  SERV_VX, SERV_VY, SERV_VZ,
  SERV_AX, SERV_AY, SERV_AZ,
  VASS_QX, VASS_QY, VASS_QZ,
  VASS_VX, VASS_VY, VASS_VZ,
  VASS_AX, VASS_AY, VASS_AZ,
  SAYS_QX, SAYS_QY, SAYS_QZ,
  SAYS_VX, SAYS_VY, SAYS_VZ,
  SAYS_AX, SAYS_AY, SAYS_AZ
};

/**************************************
 * 		Debugging Functions   *
 **************************************/

std::string exec(const char* cmd);

/**************************************
 * 		General Functions             *
 **************************************/

inline int RandomIntRange(int min, int max)
{
  assert(max > min);
  return rand()%(max-min+1) + min;
}

inline float RandomRealRange(float min, float max)
{
  assert(max > min);
  return min + (((float)rand()) / (float)RAND_MAX) * (max-min);
}

inline double RandomRealRange(double min, double max)
{
  assert(max > min);
  return min + (((double)rand()) / (double)RAND_MAX) * (max-min);
}

inline long double RandomRealRange(long double min, long double max)
{
  assert(max > min);
  return min + (((long double)rand()) / (long double)RAND_MAX) * (max-min);
}


//double round(double r);

// Log-sum-exp function
double LSE2(double, double);

/*
 *  Aminoacids 3 letter notation to 1 letter notation
 */
char aa321 (const char *aa);

/*
 * aa321 inverse
 * ?????????????????????????
 * ?????????????????????????
 */
char aa123 (char *dest, char aa);
char aa123 (std::string& dest, char aa);

/*
 *  If not, puts '/' at the end of RESULT
 */
void bCheck_path (char *result, const char *path);

/*
 * Puts START-STOP SRC substring into DEST
 * Character is counted from 0
 * Pointers have to be allocated
 */
int bExtractStr (char *dest, const char *src, int start, int stop);

/*
 * Tolower bExtractStr
 */
int bExtractTolowerStr ( char *dest, const char *src, int start, int stop );

/*
 * Modifies all characters in a string to NULL
 * Returns how many characters has modified
 */
int bZeroStr(char *dest);

/*
 * Modifies NO_ELEM elements in an array to NULL
 * Returns how many elements has modified
 */int bZeroCharArray(char *array, int no_elem);


/*
 * Awk Substr
 */
int bSubstr (char *dest, const char *src, int start, int no_chars);

std::string to_lower(std::string str);
std::string to_upper(std::string str);

/*
 * Left trim
 */
static inline std::string &ltrim(std::string &s) {
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) {return !std::isspace(c);}));
		return s;
}

/*
 * Right trim
 */
static inline std::string &rtrim(std::string &s) {
		s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {return !std::isspace(c);}).base(), s.end());
		return s;
}

/*
 * Trim from both ends
 */
static inline std::string &trim(std::string &s) {
		return ltrim(rtrim(s));
}

/*
 * Decimal prefix of zeros to limit
 */
string decimal_prefix(double inp_no, long int limit);

/*
 * Floating point numbers are aprox the same
 */
bool AreSame(double a, double b, double EPSILON);

/*
 * Check if pair, triple or quadruple is the same irrespective of order
 * from the chemical point of view (ex.: a-b-c-d != a-c-b-d)
 */
bool IsTheSameBond(const std::vector<SimTK::DuMM::AtomClassIndex>& tar, const std::vector<SimTK::DuMM::AtomClassIndex>& ref);
bool IsTheSameAngle(const std::vector<SimTK::DuMM::AtomClassIndex>& tar, const std::vector<SimTK::DuMM::AtomClassIndex>& ref);
bool IsTheSameTorsion(const std::vector<SimTK::DuMM::AtomClassIndex>& tar, const std::vector<SimTK::DuMM::AtomClassIndex>& ref);

/*
 * Given a frame F1 expressed in another frame G and a station v1 expressed 
 * in G return another frame F2 with origin in v1, aligne along F1 v1 vector
 * with the X axis and pointing towards F1
 */
SimTK::Transform alignFlipAndTranslateFrameAlongXAxis(SimTK::Transform G_X_F1, SimTK::Vec3 G_v1);


/*
 * Convert spatial maatrix (Mat< 2, 2, Mat33 >) to 6x6 matrix (Mat<6,6>)
 * Replaces inf and nan with zeros
 */
bool SpatialMat2Mat66(SimTK::SpatialMat SM, SimTK::Mat<6,6>& destination);

/*
 * Compute numerical matrix inverse with Eigen if possible
 */
bool NumericalInverse(SimTK::Matrix M, SimTK::Matrix& MInv, int nrows, int ncols);
bool NumericalLeftInverse(SimTK::Matrix M, SimTK::Matrix& MLeftInv, int nrows, int ncols);
bool NumericalRightInverse(SimTK::Matrix M, SimTK::Matrix& MRightInv, int nrows, int ncols);

/**
 * Compute determinant
*/
double cstyle_det(double *cstyle_matrix, int dim);

/**
 * Print C++ vector
*/
template <typename S>
 
void PrintCppVector(const vector<S>& vec,
	std::string sep = " ",
	std::string ending = "\n") 
{
    // Iterating over all elements of vector
    for (const auto &elem : vec) {
        std::cout << sep << elem;
    }
	std::cout << ending;

}

template <typename S>
void PrintCppVector(const vector < vector <S>>& vec,
	std::string sep = " ",
	std::string ending = "\n") 
{
    // Iterating over all elements of vector
    for (const auto &v : vec) {
		for(const auto &elem : v){
			std::cout << sep << elem;
		}
		std::cout << std::endl;
        
    }
	std::cout << ending;
}

/*
 * Get the block corresponding to a body from an H-like matrix
 * Body index "which" starts from 0, 0 being the Ground body.
 */
SimTK::Matrix& SOA_GetHLikeElement(SimTK::Matrix inMatrix, int which, SimTK::Matrix& outMatrix);
SimTK::Matrix& SOA_GetHstarLikeElement(SimTK::Matrix inMatrix, int which, SimTK::Matrix& outMatrix);

/*
 * Convert spatial vector to 6-dim vector
 */
SimTK::Vector& SOA_SpatialVec2Vector(SimTK::SpatialVec in, SimTK::Vector& out);

/*
 * Convert spatial matrix to Mat66
 */
void SOA_SpatialMat2Mat66(const SimTK::SpatialMat& in, SimTK::Mat66& out);

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(const SimTK::Matrix& M,
	int nrows, int ncols,
	int decimal_places,
	std::string header);
	
void PrintBigMat(SimTK::Mat33 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Mat44 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Mat55 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Mat66 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(
	const SimTK::Vector& M, int nrows,
	int decimal_places = 3, std::string header = "unknown matrix");

/*
 * Print Spatial Matrix
 */
void PrintSpatialVec(SimTK::SpatialVec V, int decimal_places, std::string header);
void PrintSpatialMat(SimTK::SpatialMat M, int decimal_places, std::string header);

/*
 * Print Transform
 */
void PrintTransform(SimTK::Transform T, int decimal_places);
/*
 * Angle
 */
SimTK::Real bAngle(SimTK::Vec3& pos0, SimTK::Vec3& pos1, SimTK::Vec3& pos2);

/*
 * Dihedral angle
 */
SimTK::Real bDihedral(SimTK::Vec3& pos0, SimTK::Vec3& pos1, SimTK::Vec3& pos2, SimTK::Vec3& pos3);

/**  Get a unique name based on number **/
std::string GetUniqueName(int key);

/** Magnitude (norm) of a vector of reals **/
SimTK::Real magnitude(std::vector<SimTK::Real>& V);
SimTK::Real magSq(std::vector<SimTK::Real>& V);

/** Normalize a std vector of double.
 * The result is stored in the input vector.
 * If the magnitude is 0, the vector is left untouched.
**/
void normalize(std::vector<SimTK::Real>& V);


/*
 * Thermodynamics
 */

enum struct ThermostatName : int { // Thermostats
	NONE,
	ANDERSEN,
	BERENDSEN,
	LANGEVIN,
	NOSE_HOOVER
};

/*
 * Simulation
 */

enum struct IntegratorName : int { // Integrators
	EMPTY,
	VERLET,
	EULER,
	EULER2,
	CPODES,
	RUNGEKUTTA,
	RUNGEKUTTA2,
	RUNGEKUTTA3,
	RUNGEKUTTAFELDBERG,
	OMMVV,
	NOF_INTEGRATORS
};

const std::unordered_map<std::string, IntegratorName>
IntegratorNameS{
	{"EMPTY", IntegratorName::EMPTY},
	{"VERLET", IntegratorName::VERLET},
	{"EULER", IntegratorName::EULER},
	{"EULER2", IntegratorName::EULER2},
	{"CPODES", IntegratorName::CPODES},
	{"RUNGEKUTTA", IntegratorName::RUNGEKUTTA},
	{"RUNGEKUTTA2", IntegratorName::RUNGEKUTTA2},
	{"RUNGEKUTTA3", IntegratorName::RUNGEKUTTA3},
	{"RUNGEKUTTAFELDBERG", IntegratorName::RUNGEKUTTAFELDBERG},
	{"OMMVV", IntegratorName::OMMVV}
};

/*
 * Statistics
 */
/** The type of distribution to draw a random number from **/
enum struct GmolRandDistributionType : int {
	UNIFORM,
	NORMAL
};

// Samplers
enum struct SamplerName : int {
	EMPTY,
	MC,
	HMC,
	LAHMC
};

enum struct JointType : int {
	LINEAR,
	ANGULAR180,
	ANGULAR360,
	QUATERNION_a,
	QUATERNION_b,
	QUATERNION_c,
	QUATERNION_d
};

/*
 * k-means clustering
 * http://www.goldsborough.me/c++/python/cuda/2017/09/10/20-32-46-exploring_k-means_in_python,_c++_and_cuda/
 */
struct Point {
  double x{0}, y{0};
};

using DataFrame = std::vector<Point>;

// Gmolmodel version of mean
SimTK::Real bMean(std::vector<SimTK::Real> v);
SimTK::Real bVariance(std::vector<SimTK::Real> v);
SimTK::Real bStdev(std::vector<SimTK::Real> v);

/** Gmolmodel version of correlation coefficient */
SimTK::Real bCorr(std::vector<SimTK::Real> V, std::vector<SimTK::Real> W);

/** Circular mean as in 2014 Fenwick **/
SimTK::Real circMean(std::vector<SimTK::Real> phi);

/** Circular correlation coefficient as in 2014 Fenwick **/
SimTK::Real circCorr(std::vector<SimTK::Real> phi, std::vector<SimTK::Real> psi);

/** 1 - circCorr **/
SimTK::Real circDist(std::vector<SimTK::Real> phi, std::vector<SimTK::Real> psi);

double squared_l2_distance(Point first, Point second);

DataFrame k_means(const DataFrame& data,
				  size_t k,
				  size_t number_of_iterations);

double update_variance(double current_value, 
	double last_variance, 
	double last_average,
	int n);


// STD linear algebra



// Print


void bPrintVec(std::vector<double> &src);

// Assign
std::vector<double>& bCopyVec(std::vector<double> &src,
			     std::vector<double> &dest);

// Magnitude
double bNorm(std::vector<double>& V);

// Normalize U and put it in V
std::vector<double>& bNormalize(std::vector<double> &U, std::vector<double> &V);
std::vector<double>& bNormalizeInPlace(std::vector<double> &U);

// Dot product
double bDot(std::vector<double> &u, std::vector<double> &v);

// Multiply by scalar in place
std::vector<double>& bMulByScalar(std::vector<double>& V, double scalar);

// Multiply by scalar and put in W
std::vector<double>& bMulByScalar(std::vector<double>& V, double scalar, std::vector<double>& W);

// Porject U on V and put in in p_UV
std::vector<double>& proj(std::vector<double>& u, std::vector<double>& v, std::vector<double>& p_uv);

std::vector<double>& bAddScalar(std::vector<double>& V, double scalar);

std::vector<double>& bAddVector(std::vector<double>& V, std::vector<double>& W);

// Substract V from W and put it in W
std::vector<double>& bSubstractVector(std::vector<double>& V, std::vector<double>& W);
// Matrix
typedef std::vector<std::vector<double>> bMatrix;
bMatrix& bCopyMat(bMatrix& src, bMatrix& dest);

// Print
void bPrintMat(bMatrix src);

// Transpose
bMatrix& bTransposeInPlace(bMatrix&);

// Multiply U by M and put it in V
std::vector<double>& bMulVecByMatrix(std::vector<double> &U,
				     bMatrix& M,
				     std::vector<double> &V);
// Gram–Schmidt
bMatrix& gram_schmidt(bMatrix& M, bMatrix& es);




#ifndef MONTECARLOSAMPLER
#define MONTECARLOSAMPLER MC
#endif

#ifndef HAMILTONIANMONTECARLOSAMPLER
#define HAMILTONIANMONTECARLOSAMPLER HMC
#endif

#ifndef LOOKAHEADHMCSAMPLER
#define LOOKAHEADHMCSAMPLER LAHMC
#endif

#endif /*BGENERAL_H_*/
