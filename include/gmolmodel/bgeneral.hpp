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

#define sq(x)		((x)*(x))

//#ifndef sqr
//#define sqr(x) ((x)*(x))
//#endif

#ifndef ANG_360_TO_180
#define ANG_360_TO_180(x) (((x)>180) ? ((x)-360) : (x))
#endif

// Check versus Velocity Verlet in cart coords
enum{
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


double round(double r);

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

/*
 * Left trim
 */
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

/*
 * Right trim
 */
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
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
SimTK::Mat66& SOA_SpatialMat2Mat66(SimTK::SpatialMat in, SimTK::Mat66& out);

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Matrix M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Mat33 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Mat44 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Mat66 M, int nrows, int ncols, int decimal_places, std::string header);
void PrintBigMat(SimTK::Vector M, int nrows, int decimal_places, std::string header);

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
 * Dihedral angle
 */
SimTK::Real bDihedral(SimTK::Vec3& pos0, SimTK::Vec3& pos1, SimTK::Vec3& pos2, SimTK::Vec3& pos3);

/**  Get a unique name based on number **/
std::string GetUniqueName(int key);


/*
 * Thermodynamics
 */

enum ThermostatName { // Thermostats
    NONE,
    ANDERSEN,
    BERENDSEN,
    LANGEVIN,
    NOSE_HOOVER
};

/*
 * Simulation
 */

enum IntegratorName { // Samplers
    VERLET,
    EULER,
    EULER2,
    CPODES,
    RUNGEKUTTA,
    RUNGEKUTTA2,
    RUNGEKUTTA3,
    RUNGEKUTTAFELDBERG
};

/*
 * Statistics
 */
/** The type of distribution to draw a random number from **/
enum GmolRandDistributionType {
    UNIFORM,
    NORMAL
};


enum SamplerName { // Samplers
    MC,
    HMC
};

#ifndef MONTECARLOSAMPLER
#define MONTECARLOSAMPLER MC
#endif

#ifndef HAMILTONIANMONTECARLOSAMPLER
#define HAMILTONIANMONTECARLOSAMPLER HMC
#endif


#endif /*BGENERAL_H_*/
