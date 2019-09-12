#ifndef __READAMBERINPUT_HPP__
#define __READAMBERINPUT_HPP__

#ifndef TARGET_TYPE
#define TARGET_TYPE double
#endif


#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <utility>

class readAmberInput{

  // TO DO: impropers ?!?
  // Errors to be continued...
  // Epot calc & units



    public:

      static int AtomNameSize;
      static TARGET_TYPE chargemMultiplier;
      static int AmberIndexMultiplier;
      static int AmberIndexDiff;

      void readAmberFiles(std::string inpcrdfile, std::string prmtopfile);

      int getNumberAtoms();
      int getNumberBonds();
      int getNumberAngles();
      int getNumberDihedrals();

      std::string getAtomsName(int p);
      std::string getAtomsNameAlias(int p);
      TARGET_TYPE getAtomsXcoord(int p);
      TARGET_TYPE getAtomsYcoord(int p);
      TARGET_TYPE getAtomsZcoord(int p);
      TARGET_TYPE getAtomsMass(int p);
      TARGET_TYPE getAtomsRadii(int p); // A
      TARGET_TYPE getAtomsRVdW(int p);  // A
      TARGET_TYPE getAtomsEpsilon(int p); // kcal
      // amber internal charge
      // q = amber int charge / chargemMultiplier;
      // chargemMultiplier = 18.2223;
      TARGET_TYPE getAtomsCharge(int p);

      // Natmos x Natomd matrix containing 0 for excluded pairs
      bool getNonBondedAtomsMatrix(int at1, int at2);

      // BONDS
      int getBondsAtomsIndex1(int bond);
      int getBondsAtomsIndex2(int bond);
      // Ebond = 1/2  * k (r - r0)^2
      TARGET_TYPE getBondsForceK(int bond);  // kcal / mol^-1 A^-2
      TARGET_TYPE getBondsEqval(int bond);   // A

      // ANGLE
      int getAnglesAtomsIndex1(int angle);
      int getAnglesAtomsIndex2(int angle);
      int getAnglesAtomsIndex3(int angle);
      // Eangle = 1/2  * k  (ang - ang0)^2
      TARGET_TYPE getAnglesForceK(int angle);   // kcal / mol^-1 rad^-2
      TARGET_TYPE getAnglesEqval(int angle);

      // DIHEDRALS
      // Edih = k * cos (n*phi + psi)
      int getDihedralsAtomsIndex1(int dih);
      int getDihedralsAtomsIndex2(int dih);
      int getDihedralsAtomsIndex3(int dih);
      int getDihedralsAtomsIndex4(int dih);

      //HOREA
      int getDihedralsAtomsIndex(int dihIndex, int atomIndx);
      void GeneratePairStartAndLen();

      std::vector < std::pair<int, int> > getPairStartAndLen();

      TARGET_TYPE getDihedralsForceK(int dih);  // kcal / mol^-1 rad^-2
      TARGET_TYPE getDihedralsEqval(int dih); // rad
      TARGET_TYPE getDihedralsPeriod(int dih);
      TARGET_TYPE getDihedralsPhase(int dih);

      // Other
      // TARGET_TYPE getDihedralsSCEE(int p);
      // TARGET_TYPE getDihedralsSCNB(int p);




  private:

  // General


    std::ifstream inpcrd;
    std::ifstream prmtop;

    int NumberAtoms;
    int NumberAtoms2; // for checking purposes
    int NumberBonds;
    int NumberAngles;
    int NumberDihedrals;
    int NumberExcludedAtoms;

    // ATOMS
    // 4 characters std::strings left aligned - "CA  "
    std::vector<std::string> AtomsName;
    std::vector<std::string> AtomsNameAlias;  //
    std::vector<TARGET_TYPE> AtomsXcoord;
    std::vector<TARGET_TYPE> AtomsYcoord;
    std::vector<TARGET_TYPE> AtomsZcoord;
    std::vector<TARGET_TYPE> AtomsMass;
    std::vector<TARGET_TYPE> AtomsRadii; // A
    std::vector<TARGET_TYPE> AtomsRVdW;  // A
    std::vector<TARGET_TYPE> AtomsEpsilon; // kcal
    // amber internal charge
    // q = amber int charge / chargemMultiplier;
    // chargemMultiplier = 18.2223;
    std::vector<TARGET_TYPE> AtomsCharge;


    // BONDS
    std::vector <std::vector <int> > BondsAtomsIndex;
    // Ebond = 1/2  * k (r - r0)^2
    std::vector<TARGET_TYPE> BondsForceK;  // kcal / mol^-1 A^-2
    std::vector<TARGET_TYPE> BondsEqval;   // A

    // ANGLE
    std::vector <std::vector <int> > AnglesAtomsIndex;

    // Eangle = 1/2  * k  (ang - ang0)^2
    std::vector<TARGET_TYPE> AnglesForceK;   // kcal / mol^-1 rad^-2
    std::vector<TARGET_TYPE> AnglesEqval;

    // DIHEDRALS
    std::vector <std::vector <int> > DihedralsAtomsIndex;
    std::vector < std::pair<int, int> > PairStartAndLen;
    // Edih = k * cos (n*phi + psi)
    std::vector<TARGET_TYPE> DihedralsForceK;
    std::vector<TARGET_TYPE> DihedralsEqval;
    std::vector<TARGET_TYPE> DihedralsPeriod;
    std::vector<TARGET_TYPE> DihedralsPhase;

    // Other
    // std::vector<TARGET_TYPE> SCEEScaleFactor;
    // std::vector<TARGET_TYPE> SCNBScaleFactor;






    // some temp data
    int i;
    TARGET_TYPE temp_val;
    std::string line;

    int FlagPointers[32];
    int NumberTypes;
    int NumberBondsH;
    int NumberBondsNONH;
    int NumberAnglesH;
    int NumberAnglesNONH;
    int NumberDihedralsH;
    int NumberDihedralsNONH;
    int NumberBondsTypes;
    int NumberAnglesTypes;
    int NumberDihedralsTypes;
    int LennardJonesTypes;

    std::vector<int> AtomsTypesIndex;
    std::vector<TARGET_TYPE> AtomsScreenGBIS;


    std::vector<TARGET_TYPE> tempBond_K;
    std::vector<TARGET_TYPE> tempBond_eq;
    std::vector<TARGET_TYPE> tempAngle_K;
    std::vector<TARGET_TYPE> tempAngle_eq;
    std::vector<TARGET_TYPE> tempDihedral_K;
    std::vector<TARGET_TYPE> tempDihedral_per;
    std::vector<TARGET_TYPE> tempDihedral_phase;
    std::vector<TARGET_TYPE> tempDihedral_SCEE;
    std::vector<TARGET_TYPE> tempDihedral_SCNB;
    TARGET_TYPE acoef;
    TARGET_TYPE bcoef;
    int t1;
    int t2;
    int t3;
    int t4;
    int t5;


    std::vector<TARGET_TYPE> LennardJonesNonbondParmIndex;
    std::vector<TARGET_TYPE> tempLJONES_ACOEFF;
    std::vector<TARGET_TYPE> tempLJONES_BCOEFF;


    std::vector<int> NumberExcludedAtomsList;
    std::vector <std::vector <bool> > NonBondedAtomsMatrix;


    // readAmberFiles internal field readers
    void readInpcrd();
    void readPointers();
    void readAtomsName();
    void readAtomsNameAlias();
    void readAtomsCharge();
    void readAtomsMass();
    void readAtomsRadii();
    void readAtomsTypesIndex();
    void readAtomsScreenGBIS();

    void readNonbondedParmIndex();
    void readTempBondsForceK();
    void readTempBondsEqval();
    void readTempAnglesForceK();
    void readTempAnglesEqval();
    void readTempDihedralsForceK();
    void readTempDihedralsPeriod();
    void readTempDihedralsPhase();
    void readTempSCEEScaleFactor();
    void readTempSCNBScaleFactor();

    void readBonds(int nr);
    void readAngles(int nr);
    void readDihedrals(int nr);

    void readLennardJonesACOEF();
    void readLennardJonesBCOEF();

    void readLennardJonesRVdWEpsilon();

    void readNumberExcludedAtomsList();
    void readExcludedAtomsList();




};


#endif // __READAMBERINPUT_HPP__
