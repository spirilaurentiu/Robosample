#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <utility>
#include "Molmodel.h"

template <typename T>
struct Interval {
    int start;
    int end;
    T data;

    // Overloading the less than operator for set comparison
    bool operator<(const Interval& other) const {
        return end < other.start;
    }
};

template <typename T>
class IntervalTree {
public:
    // Insert an interval into the tree
    void insert(const Interval<T>& interval) {
        intervals.insert(interval);
    }

    // Find the string for a given number
    const T& get(int number) const {
        Interval<T> dummyInterval;
        dummyInterval.start = number;
        dummyInterval.end = number;

        auto it = intervals.lower_bound(dummyInterval);

        if (it != intervals.end() && number >= it->start && number <= it->end) {
            return it->data;
        }

        // Return an empty string or an indication that the number was not found in any interval
        return notFound;
    }

private:
    std::set<Interval<T>> intervals;
    T notFound {};
};


// TO DO: impropers ?!?
// Errors to be continued...
// Epot calc & units
class readAmberInput{
public:
    void readAmberFiles(const std::string& inpcrdfile, const std::string& prmtopfile);



    int getNumberAtoms() const;
    int getNumberBonds() const;
    int getNumberAngles() const;
    int getNumberDihedrals() const;

    const std::string& getAtomsName(int p) const;
    const std::string& getAtomsType(int p) const;
    SimTK::Real getAtomsXcoord(int p) const;
    SimTK::Real getAtomsYcoord(int p) const;
    SimTK::Real getAtomsZcoord(int p) const;
    SimTK::Real getAtomsMass(int p) const;
    SimTK::Real getAtomsRadii(int p) const; // A
    SimTK::Real getAtomsRVdW(int p) const;  // A
    SimTK::Real getAtomsEpsilon(int p) const; // kcal

    SimTK::Real getAtomicNumber(int p) const;

    // amber internal charge
    // q = amber int charge / chargemMultiplier;
    // chargemMultiplier = 18.2223;
    SimTK::Real getAtomsCharge(int p) const;

    // Natmos x Natomd matrix containing 0 for excluded pairs
    bool getNonBondedAtomsMatrix(int at1, int at2) const;

    // BONDS
    int getBondsAtomsIndex1(int bond) const;
    int getBondsAtomsIndex2(int bond) const;
    // Ebond = 1/2  * k (r - r0)^2
    SimTK::Real getBondsForceK(int bond) const;  // kcal / mol^-1 A^-2
    SimTK::Real getBondsEqval(int bond) const;   // A

    // ANGLE
    int getAnglesAtomsIndex1(int angle) const;
    int getAnglesAtomsIndex2(int angle) const;
    int getAnglesAtomsIndex3(int angle) const;
    // Eangle = 1/2  * k  (ang - ang0)^2
    SimTK::Real getAnglesForceK(int angle) const;   // kcal / mol^-1 rad^-2
    SimTK::Real getAnglesEqval(int angle) const;

    // DIHEDRALS
    // Edih = k * cos (n*phi + psi)
    int getDihedralsAtomsIndex1(int dih) const;
    int getDihedralsAtomsIndex2(int dih) const;
    int getDihedralsAtomsIndex3(int dih) const;
    int getDihedralsAtomsIndex4(int dih) const;

    //HOREA
    int getDihedralsAtomsIndex(int dihIndex, int atomIndx);
    void GeneratePairStartAndLen();

    std::vector < std::pair<int, int> > getPairStartAndLen() const;

    SimTK::Real getDihedralsForceK(int dih) const;  // kcal / mol^-1 rad^-2
    SimTK::Real getDihedralsEqval(int dih) const; // rad
    SimTK::Real getDihedralsPeriod(int dih) const;
    SimTK::Real getDihedralsPhase(int dih) const;

    const std::string& getResidueLabel(int i) const;
    int getResidueIndex(int i) const;

    // Other
    // SimTK::Real getDihedralsSCEE(int p);
    // SimTK::Real getDihedralsSCNB(int p);




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
    std::vector<std::string> AtomsTypes;  //
    std::vector<SimTK::Real> AtomsXcoord;
    std::vector<SimTK::Real> AtomsYcoord;
    std::vector<SimTK::Real> AtomsZcoord;
    std::vector<SimTK::Real> AtomsMass;
    std::vector<SimTK::Real> AtomsRadii; // A
    std::vector<SimTK::Real> AtomsRVdW;  // A
    std::vector<SimTK::Real> AtomsEpsilon; // kcal
    // amber internal charge
    // q = amber int charge / chargemMultiplier;
    // chargemMultiplier = 18.2223;
    std::vector<SimTK::Real> AtomsCharge;


    // BONDS
    std::vector <std::vector <int> > BondsAtomsIndex;
    // Ebond = 1/2  * k (r - r0)^2
    std::vector<SimTK::Real> BondsForceK;  // kcal / mol^-1 A^-2
    std::vector<SimTK::Real> BondsEqval;   // A

    // ANGLE
    std::vector <std::vector <int> > AnglesAtomsIndex;

    // Eangle = 1/2  * k  (ang - ang0)^2
    std::vector<SimTK::Real> AnglesForceK;   // kcal / mol^-1 rad^-2
    std::vector<SimTK::Real> AnglesEqval;

    // DIHEDRALS
    std::vector <std::vector <int> > DihedralsAtomsIndex;
    std::vector < std::pair<int, int> > PairStartAndLen;
    // Edih = k * cos (n*phi + psi)
    std::vector<SimTK::Real> DihedralsForceK;
    std::vector<SimTK::Real> DihedralsEqval;
    std::vector<SimTK::Real> DihedralsPeriod;
    std::vector<SimTK::Real> DihedralsPhase;

    // Other
    // std::vector<SimTK::Real> SCEEScaleFactor;
    // std::vector<SimTK::Real> SCNBScaleFactor;

    std::vector<std::string> residueLabels;
    std::vector<int> residuePointers;
    IntervalTree<std::string> residueLabelTree;
    IntervalTree<int> residueIndexTree;

    // some temp data
    int global_i = 0;
    SimTK::Real temp_val;
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
    std::vector<SimTK::Real> AtomsScreenGBIS;

    std::vector<SimTK::Real> AtomicNumbers;

    std::vector<SimTK::Real> tempBond_K;
    std::vector<SimTK::Real> tempBond_eq;
    std::vector<SimTK::Real> tempAngle_K;
    std::vector<SimTK::Real> tempAngle_eq;
    std::vector<SimTK::Real> tempDihedral_K;
    std::vector<SimTK::Real> tempDihedral_per;
    std::vector<SimTK::Real> tempDihedral_phase;
    std::vector<SimTK::Real> tempDihedral_SCEE;
    std::vector<SimTK::Real> tempDihedral_SCNB;
    SimTK::Real acoef;
    SimTK::Real bcoef;
    int t1;
    int t2;
    int t3;
    int t4;
    int t5;


    std::vector<SimTK::Real> LennardJonesNonbondParmIndex;
    std::vector<SimTK::Real> tempLJONES_ACOEFF;
    std::vector<SimTK::Real> tempLJONES_BCOEFF;


    std::vector<int> NumberExcludedAtomsList;
    std::vector <std::vector <bool> > NonBondedAtomsMatrix;


    // readAmberFiles internal field readers
    void readInpcrd();
    void readPointers();
    void readAtomsName();
    void readAtomsTypes();
    void readAtomsCharge();
    void readAtomsMass();
    void readAtomsRadii();
    void readAtomsTypesIndex();
    void readAtomsScreenGBIS();

    void readResidueLabels();
    void readResiduePointers();

    void readAtomicNumber();

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
