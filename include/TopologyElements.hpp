#pragma once

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"

struct AtomSpec {
    // Indices
    int globalIndex = SimTK::InvalidIndex;
    int moleculeIndex = SimTK::InvalidIndex;
    int residueIndex = -1;
    int atomClassIndex;
    int chargedAtomTypeIndex;

    // Names
    std::string atomName, residueName, atomClassName, chargedAtomName;

    // Connectivity
    std::vector<int> neighborsGlobalIndices;
    int availableBonds = 0;
    bool root = false;

    // Physical properties
    int atomicNumber = 0;
    SimTK::Real chargeInE = 0.0;
    SimTK::mdunits::Mass massInDaltons {};
    SimTK::Real vdwRadiusInNm = 0.0;
    SimTK::Real vdwWellDepthInKJ = 0.0;
    SimTK::Real x = 0.0, y = 0.0, z = 0.0;
};

struct BondLinkSpec {
    int parentAtomGlobalIndex = std::numeric_limits<int>::min();
	int childAtomGlobalIndex = std::numeric_limits<int>::min();
    int bondGlobalIndex = 0; // amber index
    int moleculeIndex = -111111;
	bool ringClosing = false;

    SimTK::Real forceK = std::numeric_limits<SimTK::Real>::min();
	SimTK::Real forceEquil = std::numeric_limits<SimTK::Real>::min();
};

class Atom {
public:
    Atom() = default;
    
    Atom(const AtomSpec& spec) : atomSpec(spec) {
        coords = SimTK::Vec3(atomSpec.x, atomSpec.y, atomSpec.z);

        element = SimTK::Element(getAtomicNumber(), getElementName(), getElementSymbol(), getMassInDaltons());
    }

    int getGlobalIndex() const { return atomSpec.globalIndex; }

    SimTK::DuMM::AtomClassIndex getAtomClassIndex() const { return SimTK::DuMM::AtomClassIndex(atomSpec.atomClassIndex); }
    const std::string& getAtomClassName() const { return atomSpec.atomClassName; }

    SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex() const { return SimTK::DuMM::ChargedAtomTypeIndex(atomSpec.chargedAtomTypeIndex); }
    const std::string& getChargedAtomName() const { return atomSpec.chargedAtomName; }

    int getMoleculeIndex() const { return atomSpec.moleculeIndex; }

    const std::vector<int>& getNeighborsGlobalIndices() const { return atomSpec.neighborsGlobalIndices; }
    int getNumBondsInvolved() const { return atomSpec.neighborsGlobalIndices.size(); }

    int getNumAvailableBonds() const { return atomSpec.availableBonds; }
    void decrementAvailableBonds() {
        SimTK_ASSERT_ALWAYS(atomSpec.availableBonds > -1, "No more available bonds to decrement.");
        --atomSpec.availableBonds;
    }

    int getResidueIndex() const { return atomSpec.residueIndex; }
    const std::string& getResidueName() const { return atomSpec.residueName; }

    SimTK::BiotypeIndex getBiotypeIndex() const { return biotypeIndex; }
    void setBiotypeIndex(SimTK::BiotypeIndex bIdx) { biotypeIndex = bIdx; }

    SimTK::DuMM::AtomIndex getDuMMAtomIndex() const { return dAIx; }
    void setDuMMAtomIndex(SimTK::DuMM::AtomIndex dIdx) { dAIx = dIdx; }

    int getAtomicNumber() const { return atomSpec.atomicNumber; }
    SimTK::Real getChargeInE() const { return atomSpec.chargeInE; }
    SimTK::mdunits::Mass getMassInDaltons() const { return atomSpec.massInDaltons; }
    SimTK::Real getVdwRadiusInNm() const { return atomSpec.vdwRadiusInNm; }
    SimTK::Real getVdwWellDepthInKJ() const { return atomSpec.vdwWellDepthInKJ; }

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
    void setSingleAtom();

    SimTK::Compound::AtomPathName getAtomName() const { return atomSpec.atomName; }
    bool isRoot() const { return atomSpec.root; }

    const SimTK::Element& getElement() const { return element; }

private:

    AtomSpec atomSpec;
    SimTK::Element element;

    SimTK::BiotypeIndex biotypeIndex;
    SimTK::DuMM::AtomIndex dAIx; // ??????????????????????????????????????????????????????????????????????????????????????????????????????????????

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
};


class BondLink {
public:
    BondLink() = default;

    BondLink(const BondLinkSpec& spec) : bondSpec(spec) {
    }

    bool operator==(const BondLink& other) const {
        return (bondSpec.parentAtomGlobalIndex == other.bondSpec.parentAtomGlobalIndex && bondSpec.childAtomGlobalIndex == other.bondSpec.childAtomGlobalIndex) ||
               (bondSpec.parentAtomGlobalIndex == other.bondSpec.childAtomGlobalIndex && bondSpec.childAtomGlobalIndex == other.bondSpec.parentAtomGlobalIndex);
    }

    void addBondMobility(SimTK::BondMobility::Mobility someMobility) { mobilities.push_back(someMobility); }
	void setBondMobility(SimTK::BondMobility::Mobility someMobility, int world) { mobilities[world] = someMobility; }
    SimTK::BondMobility::Mobility getBondMobility(int world) const { return mobilities[world]; }

    // SimTK::Real getUScaleFactor(int world) const { return uScaleFactors[world]; }
	// void addUScaleFactor(SimTK::Real u) { uScaleFactors.push_back(u); }
	// void setUScaleFactor(int world, SimTK::Real u) { uScaleFactors[world] = u; }

    int getParentAtomGlobalIndex() const { return bondSpec.parentAtomGlobalIndex; }
    int getChildAtomGlobalIndex() const { return bondSpec.childAtomGlobalIndex; }
    int getBondGlobalIndex() const { return bondSpec.bondGlobalIndex; }
    int getMoleculeIndex() const { return bondSpec.moleculeIndex; }
    bool isRingClosing() const { return bondSpec.ringClosing; }

    SimTK::Real getForceK() const { return bondSpec.forceK; }
    SimTK::Real getForceEquil() const { return bondSpec.forceEquil; }

private:
    BondLinkSpec bondSpec;

	std::vector<SimTK::BondMobility::Mobility> mobilities;
	// std::vector<SimTK::Real> uScaleFactors = { 1.0f };
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
