#include "bgeneral.hpp"

/**************************************
 * 		General Functions             *
 **************************************/


double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}
 


//
/*
 *  Aminoacids 3 letter notation to 1 letter notation
 */
char aa321 ( const char *aa )// converteste aa din cod de 3 litere in cod de o litera
{
if 	( (!(strcmp(aa, "ALA")) ) || (!(strcmp(aa, "ala")) ) )
	return 'A';
else if ( (!(strcmp(aa, "ARG")) ) || (!(strcmp(aa, "arg")) ) )
	return 'R';
else if ( (!(strcmp(aa, "ASN")) ) || (!(strcmp(aa, "asn")) ) )
	return 'N';
else if ( (!(strcmp(aa, "ASP")) ) || (!(strcmp(aa, "asp")) ) )
	return 'D';
else if ( (!(strcmp(aa, "CYS")) ) || (!(strcmp(aa, "cys")) ) )
	return 'C';
else if ( (!(strcmp(aa, "GLN")) ) || (!(strcmp(aa, "gln")) ) )
	return 'Q';
else if ( (!(strcmp(aa, "GLU")) ) || (!(strcmp(aa, "glu")) ) )
	return 'E';
else if ( (!(strcmp(aa, "GLY")) ) || (!(strcmp(aa, "gly")) ) )
	return 'G';
else if ( (!(strcmp(aa, "HIS")) ) || (!(strcmp(aa, "his")) ) )
	return 'H';
else if ( (!(strcmp(aa, "ILE")) ) || (!(strcmp(aa, "ile")) ) )
	return 'I';
else if ( (!(strcmp(aa, "LEU")) ) || (!(strcmp(aa, "leu")) ) )
	return 'L';
else if ( (!(strcmp(aa, "LYS")) ) || (!(strcmp(aa, "lys")) ) )
	return 'K';
else if ( (!(strcmp(aa, "MET")) ) || (!(strcmp(aa, "met")) ) )
	return 'M';
else if ( (!(strcmp(aa, "PHE")) ) || (!(strcmp(aa, "phe")) ) )
	return 'F';
else if ( (!(strcmp(aa, "PRO")) ) || (!(strcmp(aa, "pro")) ) )
	return 'P';
else if ( (!(strcmp(aa, "SER")) ) || (!(strcmp(aa, "ser")) ) )
	return 'S';
else if ( (!(strcmp(aa, "THR")) ) || (!(strcmp(aa, "thr")) ) )
	return 'T';
else if ( (!(strcmp(aa, "TRP")) ) || (!(strcmp(aa, "trp")) ) )
	return 'W';
else if ( (!(strcmp(aa, "TYR")) ) || (!(strcmp(aa, "tyr")) ) )
	return 'Y';
else if ( (!(strcmp(aa, "VAL")) ) || (!(strcmp(aa, "val")) ) )
	return 'V';
else return 'x';
}

/*
 * aa321 inverse
 */
char aa123 (char *dest, char aa)
{
	if(!dest) return 1;
	if((aa == 'A') || (aa == 'a')){
		strcpy(dest, "ALA");
	}
	else if((aa == 'R') || (aa == 'r')){
		strcpy(dest, "ARG");
	}
	else if((aa == 'N') || (aa == 'n')){
		strcpy(dest, "ASN");
	}
	else if((aa == 'D') || (aa == 'd')){
		strcpy(dest, "ASP");
	}
	else if((aa == 'C') || (aa == 'c')){
		strcpy(dest, "CYS");
	}
	else if((aa == 'Q') || (aa == 'q')){
		strcpy(dest, "GLN");
	}
	else if((aa == 'E') || (aa == 'e')){
		strcpy(dest, "GLU");
	}
	else if((aa == 'G') || (aa == 'g')){
		strcpy(dest, "GLY");
	}
	else if((aa == 'H') || (aa == 'h')){
		strcpy(dest, "HIS");
	}
	else if((aa == 'I') || (aa == 'i')){
		strcpy(dest, "ILE");
	}
	else if((aa == 'L') || (aa == 'l')){
		strcpy(dest, "LEU");
	}
	else if((aa == 'K') || (aa == 'k')){
		strcpy(dest, "LYS");
	}
	else if((aa == 'M') || (aa == 'm')){
		strcpy(dest, "MET");
	}
	else if((aa == 'F') || (aa == 'f')){
		strcpy(dest, "PHE");
	}
	else if((aa == 'P') || (aa == 'p')){
		strcpy(dest, "PRO");
	}
	else if((aa == 'S') || (aa == 's')){
		strcpy(dest, "SER");
	}
	else if((aa == 'T') || (aa == 't')){
		strcpy(dest, "THR");
	}
	else if((aa == 'W') || (aa == 'w')){
		strcpy(dest, "TRP");
	}
	else if((aa == 'Y') || (aa == 'y')){
		strcpy(dest, "TYR");
	}
	else if((aa == 'V') || (aa == 'v')){
		strcpy(dest, "VAL");
	}
	else strcpy(dest, "xxx");

	return 0;
}

/*
 * aa321 inverse
 */
char aa123 (std::string& dest, char aa)
{
	if((aa == 'A') || (aa == 'a')){
		dest = "ALA";
	}
	else if((aa == 'R') || (aa == 'r')){
		dest = "ARG";
	}
	else if((aa == 'N') || (aa == 'n')){
		dest = "ASN";
	}
	else if((aa == 'D') || (aa == 'd')){
		dest = "ASP";
	}
	else if((aa == 'C') || (aa == 'c')){
		dest = "CYS";
	}
	else if((aa == 'Q') || (aa == 'q')){
		dest = "GLN";
	}
	else if((aa == 'E') || (aa == 'e')){
		dest = "GLU";
	}
	else if((aa == 'G') || (aa == 'g')){
		dest = "GLY";
	}
	else if((aa == 'H') || (aa == 'h')){
		dest = "HIS";
	}
	else if((aa == 'I') || (aa == 'i')){
		dest = "ILE";
	}
	else if((aa == 'L') || (aa == 'l')){
		dest = "LEU";
	}
	else if((aa == 'K') || (aa == 'k')){
		dest = "LYS";
	}
	else if((aa == 'M') || (aa == 'm')){
		dest = "MET";
	}
	else if((aa == 'F') || (aa == 'f')){
		dest = "PHE";
	}
	else if((aa == 'P') || (aa == 'p')){
		dest = "PRO";
	}
	else if((aa == 'S') || (aa == 's')){
		dest = "SER";
	}
	else if((aa == 'T') || (aa == 't')){
		dest = "THR";
	}
	else if((aa == 'W') || (aa == 'w')){
		dest = "TRP";
	}
	else if((aa == 'Y') || (aa == 'y')){
		dest = "TYR";
	}
	else if((aa == 'V') || (aa == 'v')){
		dest = "VAL";
	}
	else dest = "xxx";

	return 0;
}


/*
 *  If not, puts '/' at the end of RESULT
 */
void bCheck_path ( char *result, const char *path )
{
	int i;
	strcpy ( result, path );
	i = strlen ( result );
	if ( result[i-1] != '/' )	result[i] = '/';
	result[i+1] = '\0';	
}

/*
 * Puts START-STOP SRC substring into DEST
 * Character is counted from 0
 * Pointers have to be allocated
 */
//indexarea e de la 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//inclusiv caracterul cu indexul stop => '<='
//ai grija: pointerii trebuie sa fie alocati!!!!!!!!!
int bExtractStr ( char *dest, const char *src, int start, int stop )
{
        int i,j;
        for ( i=0,j=start;j<=stop;i++,j++ )     dest[i] = src[j];
        dest[i] = '\0';
        return 0;
}

/*
 * Tolower bExtractStr
 */
int bExtractCaseStr ( char *dest, const char *src, int start, int stop )
{
        int i,j;
        for ( i=0,j=start;j<=stop;i++,j++ )     dest[i] = tolower(src[j]);
        dest[i] = '\0';
        return 0;
}

/*
 * Modifies all characters in a string to NULL
 * Returns how many characters has modified
 */
int bZeroStr(char *dest)
{
	if (dest){
		int len = strlen(dest);
		int i; 
		for(i=0; i<=len; i++){
			dest[i] = '\0';
		}
		return i;
	}
	return 0;
}

/*
 * Modifies NO_ELEM elements in an array to NULL
 * Returns how many elements has modified
 */
int bZeroCharArray(char *array, int no_elem)
{
	if(!array){
		printf("bZeroArray(): NULL pointer passed as argument\n");
		exit(1);	
	}
	int i;
	for(i=0; i<no_elem; i++){
		array[i] = '\0';
	}
	return i;
}

/*
 * Awk Substr
 */
int bSubstr (char *dest, const char *src, int start, int no_chars)
{
	int src_len=strlen(src);
	if((start >= 0) && (no_chars >= 0)){
		int stop=0, k=0, t=0;;
		stop = start + no_chars;
		for(k=start;k<stop;k++,t++){
			if(k>src_len) break;
			dest[t] =  src[k];
		}
		return 0;
	}
	return 1;
}


/*
 * SEQRES line to 1-letter code aa
 * ?????????????????????
 * ?????????????????????
 *
int seqres321(char *dest, const char *src){
	int i;
	char *buffer = new char[4];
	char ter[1];
	for(i=19;i<=70;i+=4){
		bSubstr(buffer,src,i,3);
		ter[1]=aa321(buffer);
		strcat(dest, ter);
	}
	delete buffer;
	return 0;
}
*/

/*
 * Decimal prefix of zeros to limit
 */
string decimal_prefix(double inp_no, long int limit)
{
  if(abs(limit) > 1E+20){
    return string("");
  }
  if((inp_no > (limit/10)) || (inp_no < 0)){
    return string("");
  }

  ostringstream ss;

  int digits = 0;
  long int num = limit;
  while(num){
    num /= 10;
    ++digits;
  }

  string prefix(digits, '\0');
  if(int(inp_no) == 0){
      ss << (limit/10);
      prefix = ss.str();
  }
  for(long int varlim = 10; varlim <= limit; varlim *= 10){
    if((int(inp_no) >= int(varlim/10)) && (inp_no < varlim)){
      ss << (limit/varlim);
      prefix = ss.str();
    } 
  }

  return prefix.substr(1, prefix.length());
}

bool AreSame(double a, double b, double EPSILON)
{
    return fabs(a - b) < EPSILON;
}

/*
 * Given a frame F1 expressed in another frame G and a station v1 expressed
 * in G return another frame F2 with origin in v1, aligne along F1 v1 vector
 * with the X axis and pointing towards F1
 */
SimTK::Transform alignFlipAndTranslateFrameAlongXAxis(SimTK::Transform G_X_F1, SimTK::Vec3 G_v1)
{

    // Get F1v1 vector
    SimTK::Vec3 G_F1v1 = (G_v1 - G_X_F1.p());

    // Re-express F1v1 in F1
    SimTK::Vec3 F1_F1v1 = ~(G_X_F1.R()) * G_F1v1;

    // Get rotation axis and angle with dot and cross products
    SimTK::Real cosAngle = SimTK::dot(SimTK::UnitVec3(F1_F1v1), SimTK::UnitVec3(1, 0, 0));
    //std::cout << "F1_F1v1 " << F1_F1v1 << " F1_F1v1 normed " << SimTK::UnitVec3(F1_F1v1) << " cos with (1,0,0) " << cosAngle << std::endl;
    assert(cosAngle < 1.1); 
    assert(cosAngle > -1.1);
    if (cosAngle > 1.0) cosAngle = 1.0;
    if (cosAngle < -1.0) cosAngle = -1.0;
    SimTK::Angle rotAngle = std::acos( cosAngle );
    SimTK::UnitVec3 F1_rotAxis(SimTK::cross(SimTK::UnitVec3(F1_F1v1), SimTK::UnitVec3(1, 0, 0)));

    // Put results into a new Transform
    SimTK::Transform F1_X_F2( SimTK::Rotation((-1.0 * rotAngle) + SimTK::Pi, F1_rotAxis), F1_F1v1) ;

    // ALign YAxis with XAxis
    // Express everything in F2
    SimTK::Transform F2_X_F1 = ~F1_X_F2;
    SimTK::Vec3 F2_F2YAxis = SimTK::Vec3(0, 1, 0);
    SimTK::Vec3 F2_F2Orig = SimTK::Vec3(0, 0, 0);
    SimTK::Vec3 F2_F1F2 = F2_X_F1.p();
    SimTK::Vec3 F2_F1XAxis = F2_X_F1.R() * SimTK::Vec3(1, 0, 0);
    rotAngle = bDihedral(F2_F2YAxis, F2_F2Orig, F2_F1F2, F2_F1XAxis);
    SimTK::Transform F2_X_F3 = SimTK::Rotation(rotAngle, SimTK::UnitVec3(1, 0, 0));
    // SimTK::Transform G_X_F3 = G_X_F2 * F2_X_F3;
    SimTK::Transform F1_X_F3 = F1_X_F2 * F2_X_F3;

    return F1_X_F3;

}

/*
 * Convert spatial maatrix (Mat< 2, 2, Mat33 >) to 6x6 matrix (Mat<6,6>)
 * Replaces inf and nan with zeros
 */
bool SpatialMat2Mat66(SimTK::SpatialMat SM, SimTK::Mat<6,6>& M66)
{
    int ti = -1; // 6x6
    int tj = -1; // 6x6
    for(int i = 0; i < 2; i++){ // i for SpatialMatrix
            for(int k = 0; k < 3; k++){ // k for 3x3
        for(int j = 0; j < 2; j++){ // j for SpatialMatrix
                for(int l = 0; l < 3; l++){ // l for 3x3
                    tj++;
                    tj = tj % 6;
                    if(!tj){ti++;}
                    if(std::isinf(SM[i][j][k][l]) || std::isnan(SM[i][j][k][l])){
                        M66[ti][tj] = 0.0;
                    }else{
                        M66[ti][tj] = SM[i][j][k][l];
                    }
                }
            }
        }
    }
    return true;
}

/*
 * Compute numerical matrix inverse with Eigen if possible
 */
bool NumericalInverse(SimTK::Matrix M, SimTK::Matrix& MInv, int nrows, int ncols)
{
	assert(!"Not implemented");
	return false;
	/*
    assert((nrows == ncols) && "Robo::NumericalInverse: need a square matrix");
    Eigen::MatrixXd EiM(nrows, ncols);
    Eigen::MatrixXd EiMInv(nrows, ncols);

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            EiM(i, j) = M[i][j];
        }
    }
    
    EiMInv = EiM.inverse();

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            MInv[i][j] = EiMInv(i, j);
        }
    }
    return true;
    */
}

/*
 * Compute numerical left inverse with Eigen if possible
 */
bool NumericalLeftInverse(SimTK::Matrix M, SimTK::Matrix& MLeftInv, int nrows, int ncols)
{
	assert(!"Not implemented");
	return false;
	/*
    Eigen::MatrixXd EiM(nrows, ncols);
    Eigen::MatrixXd EiMLeftInv(ncols, nrows);

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            EiM(i, j) = M[i][j];
        }
    }
    
    EiMLeftInv = (EiM.transpose() * EiM).inverse() * EiM.transpose();

    for(int i = 0; i < ncols; i++){
        for(int j = 0; j < nrows; j++){
            MLeftInv[i][j] = EiMLeftInv(i, j);
        }
    }
    return true;
   */ 
}

/*
 * Compute numerical right inverse with Eigen if possible
 */
bool NumericalRightInverse(SimTK::Matrix M, SimTK::Matrix& MRightInv, int nrows, int ncols)
{
	assert(!"Not implemented");
	return false;
	/*
    Eigen::MatrixXd EiM(nrows, ncols);
    Eigen::MatrixXd EiMRightInv(ncols, nrows);

    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            EiM(i, j) = M[i][j];
        }
    }
    
    EiMRightInv = EiM.transpose() * (EiM * EiM.transpose()).inverse();

    for(int i = 0; i < ncols; i++){
        for(int j = 0; j < nrows; j++){
            MRightInv[i][j] = EiMRightInv(i, j);
        }
    }
    return true;
   */ 
}

/*
 * Convert spatial vector to 6-dim vector
 */
SimTK::Vector& SOA_SpatialVec2Vector(SimTK::SpatialVec in, SimTK::Vector& out)
{
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            out[i*3 + j] = in[i][j];
        }
    }
    return out;
}

/*
 * Convert spatial matrix to Mat66
 */
SimTK::Mat66& SOA_SpatialMat2Mat66(SimTK::SpatialMat in, SimTK::Mat66& out)
{
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    out[i*3 + k][j*3 + l] = in[i][j][k][l];
                }
            }
        }
    }
}

/*
 * Get the block corresponding to a body from an H-like matrix.
 * Body index "which" starts from 0, 0 being the Ground body.
 */
SimTK::Matrix& SOA_GetHstarLikeElement(SimTK::Matrix inMatrix, int which, SimTK::Matrix& outMatrix)
{
    if(which == 0){ // Ground
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i][j];
            }
        }
    }else if(which == 1){ // First body
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i + 6][j];
            }
        }
    }else{
        for(int i = 0; i < 6; i++){
            outMatrix[i][0] = inMatrix[(which * 6) + i][4 + which];
        }
    }
    return outMatrix;
}

/*
 * Get the block corresponding to a body from an H-like matrix.
 * Body index "which" starts from 0, 0 being the Ground body.
 */
SimTK::Matrix& SOA_GetHLikeElement(SimTK::Matrix inMatrix, int which, SimTK::Matrix& outMatrix)
{
    if(which == 0){ // Ground
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i][j];
            }
        }
    }else if(which == 1){ // First body
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                outMatrix[i][j] = inMatrix[i][j + 6];
            }
        }
    }else{
        for(int i = 0; i < 6; i++){
            outMatrix[0][i] = inMatrix[4 + which][(which * 6) + i];
        }
    }
    return outMatrix;
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Matrix M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat33 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat44 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat55 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Big Matrices separated by spaces
 */
void PrintBigMat(SimTK::Mat66 M, int nrows, int ncols, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Print Spatial Matrix
 */
void PrintSpatialMat(SimTK::SpatialMat M, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    std::cout << M[i][j][k][l] << " ";
                }
                std::cout << "  ";
            }
        }
        std::cout << std::endl;
    }
}

/*
 * Print Spatial Vector
 */
void PrintSpatialVec(SimTK::SpatialVec M, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places) << std::fixed;
    for(int i = 0; i < 2; i++){
        for(int k = 0; k < 3; k++){
            std::cout << M[i][k] << " ";
        }
    }
    std::cout << std::endl;
}

/*
 * Print Big Vector separated by spaces
 */
void PrintBigMat(SimTK::Vector V, int nrows, int decimal_places, std::string header)
{
    std::cout << header << std::endl;
    std::cout << std::setprecision(decimal_places) << std::fixed;
    for(int i = 0; i < nrows; i++){
            std::cout << V[i] << " ";
    }
    std::cout << std::endl;
}

/*
 * Print Transform
 */
void PrintTransform(SimTK::Transform T, int decimal_places)
{
    const SimTK::Mat44 M = T.toMat44();
    std::cout << std::setprecision(decimal_places) << std::fixed;
    for(int i = 0; i < 4; i++){
        for(int k = 0; k < 4; k++){
            std::cout << M(i, k) << " ";
        }
        std::cout << std::endl;
    }

}

/*
 * Dihedral angle
 */
SimTK::Real bDihedral(SimTK::Vec3& pos0, SimTK::Vec3& pos1, SimTK::Vec3& pos2, SimTK::Vec3& pos3){
 
  //std::cout << "bDihedral pos " << pos0 << ' ' << pos1 << ' ' << pos2 << ' ' << pos3 << std::endl;
 
  SimTK::Vec3 diffs[3];
  SimTK::Vec3 normals[2];
  double dots[2];
  double psin, pcos;

  diffs[0] = pos1 - pos0;
  diffs[1] = pos2 - pos1;
  diffs[2] = pos3 - pos2;

  normals[0] = diffs[0] % diffs[1];
  normals[1] = diffs[1] % diffs[2];

  dots[0] = (normals[0][0] * diffs[2][0]) + (normals[0][1] * diffs[2][1]) + (normals[0][2] * diffs[2][2]);
  dots[1] = (diffs[1][0] * diffs[1][0]) + (diffs[1][1] * diffs[1][1]) + (diffs[1][2] * diffs[1][2]);

  psin = dots[0] * std::sqrt(dots[1]);
  pcos = (normals[0][0] * normals[1][0]) + (normals[0][1] * normals[1][1]) + (normals[0][2] * normals[1][2]);

  //std::cout << "bDihedral diffs " << diffs[0] << ' ' << diffs[1] << ' ' << diffs[2] << std::endl;
  //std::cout << "bDihedral normals " << normals[0] << ' ' << normals[1] << std::endl;
  //std::cout << "bDihedral dots " << dots[0] << ' ' << dots[1] << std::endl;
  //std::cout << "bDihedral psin pcos " << psin << ' ' << pcos << std::endl;
  //std::cout << " bDihedral return " << atan2(psin, pcos) << "  ==========" << std::endl; 

  return atan2(psin, pcos);
}

/**  Get a unique name based on number **/
// Assign a unique name specific to Gmolmodel. There are 60 available
// ASCII readble characters: 0-9, A-Z and a-z. This gives a 12.960.000
// of possible 4 character combinations in a number of the form
// a*60^3 + b*60^2 + c*60^1 + d. However the readble characters do not
// form a continuous interval in the ASCII table so they have to be
// spread.
std::string GetUniqueName(int nameCounter) {

    std::string string_name;
    std::string aStr, bStr, cStr, dStr;
    int a=65, b=65, c=65, d=65;
    int aRest=0, bRest=0, cRest=0;

    a = int(nameCounter / std::pow(25, 3));
    aStr = (char)(a + 65);
    aRest = nameCounter % int(std::pow(25, 3));

    b = int(aRest / std::pow(25, 2));
    bStr = (char)(b + 65);
    bRest = aRest % int(std::pow(25, 2));

    c = int(bRest / std::pow(25, 1));
    cStr = (char)(c + 65);
    cRest = bRest % int(std::pow(25, 1));

    d = int(cRest / std::pow(25, 0));
    dStr = (char)(d + 65);

    string_name = aStr + bStr + cStr + dStr;

    return string_name;
}


/** Magnitude (norm) of a vector of reals **/
SimTK::Real magnitude(std::vector<SimTK::Real>& V) {
	SimTK::Real mag = 0;
	for(int i = 0; i < V.size(); i++){
		mag += (V[i] * V[i]);
	}

	return std::sqrt(mag);
}

/** Magnitude (norm) of a vector of reals **/
SimTK::Real magSq(std::vector<SimTK::Real>& V) {
	SimTK::Real mag = 0;
	for(int i = 0; i < V.size(); i++){
		mag += (V[i] * V[i]);
	}

	return mag;
}

/** Normalize a std vector of double. If the magnitude is 0 it returns 
the vector back.**/
// TODO: optimize
std::vector<SimTK::Real>& normalize(std::vector<SimTK::Real>& V) {

	SimTK::Real mag = 0;
	for(int i = 0; i < V.size(); i++){
		mag += (V[i] * V[i]);
	}
	mag = std::sqrt(mag);

	if(mag != 0){
		for(int i = 0; i < V.size(); i++){
			V[i] /= mag;
		}
	}

	return V;
}


//#include <algorithm>
//#include <cstdlib>
//#include <limits>
//#include <random>
//#include <vector>

//struct Point {
//  double x{0}, y{0};
//};

//using DataFrame = std::vector<Point>;


//double square(double value) {
//  return value * value;
//}

using SimTK::square;

// Circular mean or mean direction as defined in Jammalamadaka and Sengupta
SimTK::Real circMean(std::vector<SimTK::Real> phi){
	SimTK::Real sum_sin = 0;
	SimTK::Real sum_cos = 0;
	for(int i = 0; i < phi.size(); i++){
		sum_sin += std::sin(phi[i]);
		sum_cos += std::cos(phi[i]);
	}
	return std::atan2(sum_sin, sum_cos);
}

SimTK::Real circCorr(std::vector<SimTK::Real> phi, std::vector<SimTK::Real> psi){
	assert(phi.size() == psi.size());
	SimTK::Real phi_mean, psi_mean;
	SimTK::Real num = 0;
	SimTK::Real sum_sinSq_phi = 0;
	SimTK::Real sum_sinSq_psi = 0;
	SimTK::Real denom = 0;

	phi_mean = circMean(phi);
	psi_mean = circMean(psi);
	
	for(int i = 0; i < phi.size(); i++){
		num += std::sin(phi[i] - phi_mean) * std::sin(psi[i] - psi_mean);
		sum_sinSq_phi += SimTK::square(std::sin(phi[i] - phi_mean));
		sum_sinSq_psi += SimTK::square(std::sin(psi[i] - psi_mean));
	}
	denom = std::sqrt(sum_sinSq_phi * sum_sinSq_psi);

	return num / denom;

}

SimTK::Real circDist(std::vector<SimTK::Real> phi, std::vector<SimTK::Real> psi){
	return 1 - circCorr(phi, psi);
}


double squared_l2_distance(Point first, Point second) {
  return square(first.x - second.x) + square(first.y - second.y);
}

DataFrame k_means(const DataFrame& data,
                  size_t k,
                  size_t number_of_iterations) {
  static std::random_device seed;
  static std::mt19937 random_number_generator(seed());
  std::uniform_int_distribution<size_t> indices(0, data.size() - 1);

  // Pick centroids as random points from the dataset.
  DataFrame means(k);
  for (auto& cluster : means) {
    cluster = data[indices(random_number_generator)];
  }

  std::vector<size_t> assignments(data.size());
  for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
    // Find assignments.
    for (size_t point = 0; point < data.size(); ++point) {
      double best_distance = std::numeric_limits<double>::max();
      size_t best_cluster = 0;
      for (size_t cluster = 0; cluster < k; ++cluster) {
        const double distance =
            squared_l2_distance(data[point], means[cluster]);
        if (distance < best_distance) {
          best_distance = distance;
          best_cluster = cluster;
        }
      }
      assignments[point] = best_cluster;
    }

    // Sum up and count points for each cluster.
    DataFrame new_means(k);
    std::vector<size_t> counts(k, 0);
    for (size_t point = 0; point < data.size(); ++point) {
      const auto cluster = assignments[point];
      new_means[cluster].x += data[point].x;
      new_means[cluster].y += data[point].y;
      counts[cluster] += 1;
    }

    // Divide sums by counts to get new centroids.
    for (size_t cluster = 0; cluster < k; ++cluster) {
      // Turn 0/0 into 0/1 to avoid zero division.
      const auto count = std::max<size_t>(1, counts[cluster]);
      means[cluster].x = new_means[cluster].x / count;
      means[cluster].y = new_means[cluster].y / count;
    }
  }

  return means;
}






