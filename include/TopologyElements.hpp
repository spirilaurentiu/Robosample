#pragma once

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"

// struct AtomSpec {}

class Atom {
public:
    Atom() = default;
    
    Atom(int globalIndex, int molIx, int atomicNumber, SimTK::Real charge, SimTK::mdunits::Mass mass, SimTK::Real vdw, SimTK::Real lj, const std::string& resName, int resIx, SimTK::Real x, SimTK::Real y, SimTK::Real z, const std::string& name, bool root) :
        globalIndex(globalIndex),
        moleculeIndex(molIx),
        atomicNumber(atomicNumber),
        charge(charge),
        mass(mass),
        vdwRadius(vdw),
        LJWellDepth(lj),
        residueName(resName),
        residueIndex(resIx),
        name(name),
        root(root) {
        coords = SimTK::Vec3(x, y, z);
        dummAtomClassIndex = SimTK::DuMM::AtomClassIndex(globalIndex);
        chargedAtomTypeIndex = SimTK::DuMM::ChargedAtomTypeIndex(globalIndex);
    }

    int getGlobalIndex() const { return globalIndex; }
    SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex() const { return chargedAtomTypeIndex; }
    SimTK::DuMM::AtomClassIndex getDummAtomClassIndex() const { return dummAtomClassIndex; }

    int getMoleculeIndex() const { return moleculeIndex; }

    void addNeighborGlobalIndex(int neighbor) {
        neighborsGlobalIndices.push_back(neighbor);
        availableBonds++;
    }
    const std::vector<int>& getNeighborsGlobalIndices() const { return neighborsGlobalIndices; }
    int getNumBondsInvolved() const { return neighborsGlobalIndices.size(); }

    void addInvolvedBondGlobalIndex(int bondIdx) { bondsInvolvedGlobalIndex.push_back(bondIdx); }
    const std::vector<int>& getBondsInvolvedIndices() const { return bondsInvolvedGlobalIndex; }

    int getNumAvailableBonds() const { return availableBonds; }
    void decrementAvailableBonds() {
        SimTK_ASSERT_ALWAYS(availableBonds > -1, "No more available bonds to decrement.");
        --availableBonds;
    }

    int getResidueIndex() const { return residueIndex; }
    const std::string& getResidueName() const { return residueName; }

    // SimTK::BiotypeIndex getBiotypeIndex() const { return biotypeIndex; }
    // void setBiotypeIndex(SimTK::BiotypeIndex bIdx) { biotypeIndex = bIdx; }

    SimTK::DuMM::AtomIndex getDuMMAtomIndex() const { return dAIx; }
    void setDuMMAtomIndex(SimTK::DuMM::AtomIndex dIdx) { dAIx = dIdx; }

    int getAtomicNumber() const { return atomicNumber; }
    SimTK::Real getCharge() const { return charge; }
    SimTK::mdunits::Mass getMass() const { return mass; }
    SimTK::Real getVdwRadius() const { return vdwRadius; }
    SimTK::Real getLJWellDepth() const { return LJWellDepth; }

    std::string getElementName() const;
    std::string getElementSymbol() const;

    SimTK::Real getX() const { return coords[0]; }
    void setX(SimTK::Real x) { coords[0] = x; }

    SimTK::Real getY() const { return coords[1]; }
    void setY(SimTK::Real y) { coords[1] = y; }

    SimTK::Real getZ() const { return coords[2]; }
    void setZ(SimTK::Real z) { coords[2] = z; }

    const SimTK::Vec3& getCoords() const { return coords; }
    void setCoords(const SimTK::Vec3& c) { coords = c; }

    SimTK::Compound::AtomIndex getCompoundAtomIndex() const { return compoundAtomIndex; }
    void setCompoundAtomIndex(SimTK::Compound::AtomIndex cIdx) { compoundAtomIndex = cIdx; }

    const SimTK::Compound::SingleAtom& getSingleAtom() const { return *compoundSingleAtom; }
    void setSingleAtom(const SimTK::Element &element);

    SimTK::Compound::Name getName() const { return name; }
    bool isRoot() const { return root; }

private:

    int globalIndex = SimTK::InvalidIndex; // amber index
    int parentAtomGlobalIndex = SimTK::InvalidIndex; // amber index of parent atom
    int moleculeIndex = SimTK::InvalidIndex;
    
    std::vector<int> neighborsGlobalIndices;
    std::vector<int> bondsInvolvedGlobalIndex;
    int availableBonds = std::numeric_limits<int>::min();

    int residueIndex = std::numeric_limits<int>::min();
    std::string residueName; // Residue and chain

    SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    // SimTK::BiotypeIndex biotypeIndex;
    SimTK::DuMM::AtomIndex dAIx; // ??????????????????????????????????????????????????????????????????????????????????????????????????????????????

    int atomicNumber = std::numeric_limits<int>::min(); // atomic number
    SimTK::Real charge = SimTK::NaN;
    SimTK::mdunits::Mass mass = SimTK::NaN;
    SimTK::Real vdwRadius = SimTK::NaN;
    SimTK::Real LJWellDepth = SimTK::NaN;

    SimTK::Vec3 coords = SimTK::Vec3(SimTK::NaN);

    // wasted 6 hours trying to make this unique_ptr or allocated on the stack
    // after more hours wasted, i read this article: https://www.cppstories.com/2014/05/vector-of-objects-vs-vector-of-pointers/
    // as expected, it should be faster to have this allocated on the stack
    // however, taking into account that this is called only a few times, it is not worth the trouble
    // there are two use cases for this: addition of bond centers and building the molecule graph
    // in the first case, the calls are sequential, so we waste some time
    // in the second case, the calls are random and the article shows that it is actually better to have the objects allocated on the heap
    SimTK::Compound::SingleAtom* compoundSingleAtom = nullptr;
    SimTK::Compound::AtomIndex compoundAtomIndex; // this is the local index in the compound, not the global index. it's different from globalIndex

    SimTK::Compound::Name name;

    bool root = false; // specifies if this atom is a base atom
};


class BondLink {
public:
    BondLink() = default;

    BondLink(int parentAtomGlobalIndex, int childAtomGlobalIndex, int bondGlobalIndex, int moleculeIndex, bool ringClosing, SimTK::Real forceK, SimTK::Real forceEquil)
        : parentAtomGlobalIndex(parentAtomGlobalIndex),
          childAtomGlobalIndex(childAtomGlobalIndex),
          bondGlobalIndex(bondGlobalIndex),
          moleculeIndex(moleculeIndex),
          ringClosing(ringClosing),
          forceK(forceK),
          forceEquil(forceEquil) {
    }

    bool operator==(const BondLink& other) const {
        return (parentAtomGlobalIndex == other.parentAtomGlobalIndex && childAtomGlobalIndex == other.childAtomGlobalIndex) ||
               (parentAtomGlobalIndex == other.childAtomGlobalIndex && childAtomGlobalIndex == other.parentAtomGlobalIndex);
    }

    void addBondMobility(SimTK::BondMobility::Mobility someMobility) { mobilities.push_back(someMobility); }
	void setBondMobility(SimTK::BondMobility::Mobility someMobility, int world) { mobilities[world] = someMobility; }
    SimTK::BondMobility::Mobility getBondMobility(int world) const { return mobilities[world]; }

    SimTK::Real getUScaleFactor(int world) const { return uScaleFactors[world]; }
	void addUScaleFactor(SimTK::Real u) { uScaleFactors.push_back(u); }
	void setUScaleFactor(int world, SimTK::Real u) { uScaleFactors[world] = u; }

    int getParentAtomGlobalIndex() const { return parentAtomGlobalIndex; }
    int getChildAtomGlobalIndex() const { return childAtomGlobalIndex; }
    int getBondGlobalIndex() const { return bondGlobalIndex; }
    int getMoleculeIndex() const { return moleculeIndex; }
    bool isRingClosing() const { return ringClosing; }
    SimTK::Compound::BondIndex getCompoundBondIndex() const { return compoundBondIndex; }

    SimTK::Real getForceK() const { return forceK; }
    SimTK::Real getForceEquil() const { return forceEquil; }

private:
	std::vector<SimTK::BondMobility::Mobility> mobilities;
	std::vector<SimTK::Real> uScaleFactors = { 1.0f };

	// These will correspond to Atom.number
	int parentAtomGlobalIndex = std::numeric_limits<int>::min();
	int childAtomGlobalIndex = std::numeric_limits<int>::min();
    int bondGlobalIndex = 0; // amber index
    int moleculeIndex = -111111;
	bool ringClosing = false;
	SimTK::Compound::BondIndex compoundBondIndex = std::numeric_limits<SimTK::Compound::BondIndex>::min();

    SimTK::Real forceK = std::numeric_limits<SimTK::Real>::min();
	SimTK::Real forceEquil = std::numeric_limits<SimTK::Real>::min();
};

class BondAngle {
public:
    BondAngle() = default;

    BondAngle(int firstGlobalIndex, int secondGlobalIndex, int thirdGlobalIndex, SimTK::Real k, SimTK::Real equil)
        : firstGlobalIndex(firstGlobalIndex), secondGlobalIndex(secondGlobalIndex), thirdGlobalIndex(thirdGlobalIndex), k(k), equil(equil) {
    }

    int getFirstGlobalIndex() const { return firstGlobalIndex; }
    int getSecondGlobalIndex() const { return secondGlobalIndex; }
    int getThirdGlobalIndex() const { return thirdGlobalIndex; }
    SimTK::Real getK() const { return k; }
    SimTK::Real getEquil() const { return equil; }

private:
	int firstGlobalIndex = std::numeric_limits<int>::min();
	int secondGlobalIndex = std::numeric_limits<int>::min();
	int thirdGlobalIndex = std::numeric_limits<int>::min();

	SimTK::Real k = std::numeric_limits<SimTK::Real>::min();
	SimTK::Real equil = std::numeric_limits<SimTK::Real>::min();
};

class BondTorsion {
public:
    BondTorsion() = default;

    BondTorsion(int firstGlobalIndex, int secondGlobalIndex, int thirdGlobalIndex, int fourthGlobalIndex, bool improper, const std::array<SimTK::Real, 4>& k, const std::array<SimTK::Real, 4>& phase, const std::array<int, 4>& period)
        : firstGlobalIndex(firstGlobalIndex), secondGlobalIndex(secondGlobalIndex), thirdGlobalIndex(thirdGlobalIndex), fourthGlobalIndex(fourthGlobalIndex), improper(improper), k(k), phase(phase), period(period) {
    }

    int getFirstGlobalIndex() const { return firstGlobalIndex; }
    int getSecondGlobalIndex() const { return secondGlobalIndex; }
    int getThirdGlobalIndex() const { return thirdGlobalIndex; }
    int getFourthGlobalIndex() const { return fourthGlobalIndex; }

    const std::array<int, 4>& getPeriod() const { return period; }
    const std::array<SimTK::Real, 4>& getK() const { return k; }
    const std::array<SimTK::Real, 4>& getPhase() const { return phase; }
    int getNum() const { return num; }
    bool isImproper() const { return improper; }

private:
	// These values are filled according to num (see below)
	std::array<SimTK::Real, 4> k { 0, 0, 0, 0 };
	std::array<SimTK::Real, 4> phase { 0, 0, 0, 0 };
	std::array<int, 4> period { 0, 0, 0, 0 };

	// How many impropers with these four indices are present here
	int num = 0; 

	int firstGlobalIndex = std::numeric_limits<int>::min();
	int secondGlobalIndex = std::numeric_limits<int>::min();
	int thirdGlobalIndex = std::numeric_limits<int>::min();
	int fourthGlobalIndex = std::numeric_limits<int>::min();
	bool improper = false;
};
