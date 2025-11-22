#pragma once

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"

// struct AtomClassDefinition {
//     std::string atomTypeName; // Equivalent to AMBER atom type name (eg CT C CA CM CC CV CW CR etc), not AMBER atom name (eg N, CA, C, O, C1, C2, H1 etc)
//     SimTK::Real vdwRadiusInNm = 0.0;
//     SimTK::Real vdwWellDepthInKJ = 0.0;
//     int atomClassIndex = 0;
//     int atomicNumber = 0;
//     int expectedValence = 0;
// };

// struct ChargedAtomTypeDefinition {
//     std::string biotypeAtomName, biotypeResidueName;
//     SimTK::Real partialChargeInE = 0.0;
//     int chargedAtomTypeIndex = 0;
//     int atomClassIndex = 0;
// };

// struct BondStretchDefinition {
//     int atomClassIndex1 = 0, atomClassIndex2 = 0;
//     SimTK::Real stiffnessInKJperNmSq = 0.0;
// 	SimTK::Real nominalLengthInNm = 0.0;
// };

// struct BondBendDefinition {
//     
//     SimTK::Real stiffnessInKJPerRadSq = 0.0;
//     SimTK::Real nominalAngleInDeg = 0.0;
// };

// struct BondTorsionDefinition {
//     int atomClassIndex1 = 0, atomClassIndex2 = 0, atomClassIndex3 = 0, atomClassIndex4 = 0;
//     SimTK::Real ampInKJ = 0.0;
//     SimTK::Real phaseInDegrees = 0.0;
//     int periodicity = 0;
//     bool improper = false;
// };



struct AtomDefinition {
    // Indices
    int globalIndex = 0;
    int moleculeIndex = 0;
    int residueIndex = 0;
    int atomClassIndex = 0;
    int chargedAtomTypeIndex = 0;

    // Names
    std::string biotypeAtomName; // AMBER atom name (eg N, CA, C, O, C1, C2, H1 etc) + ':' + valence (number of actual bonds, not typical valence)
    std::string atomClassName; // AMBER atom type name (eg CT C CA CM CC CV CW CR etc)
    std::string chargedAtomName; // Biotype name for charged atom type
    std::string residueName; // AMBER residue name (eg ALA, GLY, SER, THR etc)
    std::string uniqueAtomName; // LYS2_NZ_23:3 (23 is the atom index in the entire molecule as specified by the prmtop file, :3 is the valence)

    // Connectivity
    std::vector<int> neighborsGlobalIndices;
    
    bool root = false;

    // Physical properties
    int atomicNumber = 0;
    SimTK::Real chargeInE = 0.0;
    SimTK::Real massInDaltons = 0.0;
    SimTK::Real vdwRadiusInNm = 0.0, sigmaInNm = 0.0;
    SimTK::Real vdwWellDepthInKJ = 0.0;
    SimTK::Real x = 0.0, y = 0.0, z = 0.0;
};

struct BondStretchDefinition {
    int parentAtomGlobalIndex = 0;
	int childAtomGlobalIndex = 0;
    int bondGlobalIndex = 0;
    int moleculeIndex = 0;
	bool ringClosing = false;

    SimTK::Real stiffnessInKJPerNmSq = 0.0;
	SimTK::Real nominalLengthInNm = 0.0;
};

struct BondBendDefinition {
    int globalIndex1 = 0, globalIndex2 = 0, globalIndex3 = 0;

	SimTK::Real stiffnessInKJPerRadSq = 0.0;
    SimTK::Real nominalAngleInDeg = 0.0;
};

struct BondTorsionDefinition {
    int globalIndex1 = 0, globalIndex2 = 0, globalIndex3 = 0, globalIndex4 = 0;
    SimTK::Real ampInKJ = 0.0;
    SimTK::Real phaseInDegrees = 0.0;
    int periodicity = 0;
    bool improper = false;
};

class Atom {
public:
    Atom() = default;
    
    Atom(const AtomDefinition& spec) : atomSpec(spec) {
        coords = SimTK::Vec3(atomSpec.x, atomSpec.y, atomSpec.z);
        availableBonds = atomSpec.neighborsGlobalIndices.size();
        element = SimTK::Element(getAtomicNumber(), getElementName(), getElementSymbol(), getMassInDaltons());
    }

    int getGlobalIndex() const { return atomSpec.globalIndex; }

    SimTK::DuMM::AtomClassIndex getAtomClassIndex() const { return SimTK::DuMM::AtomClassIndex(atomSpec.atomClassIndex); }
    SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex() const { return SimTK::DuMM::ChargedAtomTypeIndex(atomSpec.chargedAtomTypeIndex); }

    SimTK::Compound::AtomPathName getAtomName() const { return atomSpec.biotypeAtomName; }
    const std::string& getAtomClassName() const { return atomSpec.atomClassName; }
    const std::string& getChargedAtomName() const { return atomSpec.chargedAtomName; }
    const std::string& getResidueName() const { return atomSpec.residueName; }
    const std::string& getUniqueAtomName() const { return atomSpec.uniqueAtomName; }

    int getMoleculeIndex() const { return atomSpec.moleculeIndex; }

    const std::vector<int>& getNeighborsGlobalIndices() const { return atomSpec.neighborsGlobalIndices; }
    int getNumBondsInvolved() const { return atomSpec.neighborsGlobalIndices.size(); }

    int getNumAvailableBonds() const { return availableBonds; }
    void decrementAvailableBonds() {
        SimTK_ASSERT_ALWAYS(availableBonds > -1, "No more available bonds to decrement.");
        --availableBonds;
    }

    int getResidueIndex() const { return atomSpec.residueIndex; }

    SimTK::BiotypeIndex getBiotypeIndex() const { return biotypeIndex; }
    void setBiotypeIndex(SimTK::BiotypeIndex bIdx) { biotypeIndex = bIdx; }

    SimTK::DuMM::AtomIndex getDuMMAtomIndex() const { return dAIx; }
    void setDuMMAtomIndex(SimTK::DuMM::AtomIndex dIdx) { dAIx = dIdx; }

    int getAtomicNumber() const { return atomSpec.atomicNumber; }
    SimTK::Real getChargeInE() const { return atomSpec.chargeInE; }
    SimTK::mdunits::Mass getMassInDaltons() const { return atomSpec.massInDaltons; }
    SimTK::Real getVdwRadiusInNm() const { return atomSpec.vdwRadiusInNm; }
    SimTK::Real getSigmaInNm() const { return atomSpec.sigmaInNm; }
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
    void createSingleAtom();

    bool isRoot() const { return atomSpec.root; }

    const SimTK::Element& getElement() const { return element; }

private:

    AtomDefinition atomSpec;
    SimTK::Element element;

    int availableBonds = 0;

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


class BondStretch {
public:
    BondStretch() = default;

    BondStretch(const BondStretchDefinition& spec) : bondSpec(spec) {}

    bool operator==(const BondStretch& other) const {
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

    SimTK::Real getStiffnessInKJPerNmSq() const { return bondSpec.stiffnessInKJPerNmSq; }
    SimTK::Real getNominalLengthInNm() const { return bondSpec.nominalLengthInNm; }

private:
    BondStretchDefinition bondSpec;

	std::vector<SimTK::BondMobility::Mobility> mobilities;
	// std::vector<SimTK::Real> uScaleFactors = { 1.0f };
};

class BondBend {
public:
    BondBend() = default;

    BondBend(const BondBendDefinition & spec) : angleSpec(spec) {}

    int getGlobalIndex1() const { return angleSpec.globalIndex1; }
    int getGlobalIndex2() const { return angleSpec.globalIndex2; }
    int getGlobalIndex3() const { return angleSpec.globalIndex3; }
    SimTK::Real getStiffnessInKJPerRadSq() const { return angleSpec.stiffnessInKJPerRadSq; }
    SimTK::Real getNominalAngleInDeg() const { return angleSpec.nominalAngleInDeg; }

private:
	BondBendDefinition angleSpec;
};

class BondTorsion {
public:
    BondTorsion() = default;

    BondTorsion(const BondTorsionDefinition& spec) : torsionSpec(spec) {}

    int getGlobalIndex1() const { return torsionSpec.globalIndex1; }
    int getGlobalIndex2() const { return torsionSpec.globalIndex2; }
    int getGlobalIndex3() const { return torsionSpec.globalIndex3; }
    int getGlobalIndex4() const { return torsionSpec.globalIndex4; }

    int getPeriodicity() const { return torsionSpec.periodicity; }
    SimTK::Real getAmpInKJ() const { return torsionSpec.ampInKJ; }
    SimTK::Real getPhaseInDegrees() const { return torsionSpec.phaseInDegrees; }
    bool isImproper() const { return torsionSpec.improper; }

private:
    BondTorsionDefinition torsionSpec;
};
