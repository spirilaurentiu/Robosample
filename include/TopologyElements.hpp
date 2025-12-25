#pragma once

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"


struct AtomDefinition {
    // Indices
    int globalIndex = 0;
    int prmtopIndex = 0;
    int moleculeIndex = 0;
    int residueIndex = 0;

    // Atom class identifiers
    std::string atomClassName;
    int atomClassIndex = 0;

    // Charged atom class identifiers
    std::string chargedAtomTypeName;
    int chargedAtomTypeIndex = 0;

    // Names
    std::string residueName;
    std::string uniqueAtomName;

    // Connectivity
    std::vector<int> neighborsGlobalIndices;
    bool root = false;

    // Physical properties
    int atomicNumber = 0;
    SimTK::Real chargeInE = 0.0;
    SimTK::Real massInDaltons = 0.0;
    SimTK::Real vdwRadiusInNm = 0.0, sigmaInNm = 0.0;
    SimTK::Real vdwWellDepthInKJ = 0.0;
    SimTK::Real x_nm = 0.0, y_nm = 0.0, z_nm = 0.0;
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
    int globalIndex1 = -1, globalIndex2 = -1, globalIndex3 = -1, globalIndex4 = -1;
    bool improper = false;

    SimTK::Real ampInKJ_1 = -1.0, phaseInDegrees_1 = -1.0, periodicity_1 = -1;
    SimTK::Real ampInKJ_2 = -1.0, phaseInDegrees_2 = -1.0, periodicity_2 = -1;
    SimTK::Real ampInKJ_3 = -1.0, phaseInDegrees_3 = -1.0, periodicity_3 = -1;
    SimTK::Real ampInKJ_4 = -1.0, phaseInDegrees_4 = -1.0, periodicity_4 = -1;
    SimTK::Real ampInKJ_5 = -1.0, phaseInDegrees_5 = -1.0, periodicity_5 = -1;
};

class Atom {
public:
    Atom() = default;
    
    Atom(const AtomDefinition& spec) : atomSpec(spec) {
        CoordsInNm = SimTK::Vec3(atomSpec.x_nm, atomSpec.y_nm, atomSpec.z_nm);
        availableBonds = atomSpec.neighborsGlobalIndices.size();
        element = SimTK::Element(getAtomicNumber(), getElementName(), getElementSymbol(), getMassInDaltons());
    }

    int getGlobalIndex() const { return atomSpec.globalIndex; }
    int getPrmtopIndex() const { return atomSpec.prmtopIndex; }

    const std::string& getAtomClassName() const { return atomSpec.atomClassName; }
    SimTK::DuMM::AtomClassIndex getAtomClassIndex() const { return SimTK::DuMM::AtomClassIndex(atomSpec.atomClassIndex); }

    const std::string& getChargedAtomTypeName() const { return atomSpec.chargedAtomTypeName; }
    SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex() const { return SimTK::DuMM::ChargedAtomTypeIndex(atomSpec.chargedAtomTypeIndex); }

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

    SimTK::Real getXInNm() const { return CoordsInNm[0]; }
    void setXInNm(SimTK::Real x_nm) { CoordsInNm[0] = x_nm; }

    SimTK::Real getYInNm() const { return CoordsInNm[1]; }
    void setYInNm(SimTK::Real y_nm) { CoordsInNm[1] = y_nm; }

    SimTK::Real getZInNm() const { return CoordsInNm[2]; }
    void setZInNm(SimTK::Real z_nm) { CoordsInNm[2] = z_nm; }

    const SimTK::Vec3& getCoordsInNm() const { return CoordsInNm; }
    void setCoordsInNm(const SimTK::Vec3& c) { CoordsInNm = c; }

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

    SimTK::Vec3 CoordsInNm = SimTK::Vec3(SimTK::NaN);

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
    bool isImproper() const { return torsionSpec.improper; }

    SimTK::Real getAmpInKJ_1() const { return torsionSpec.ampInKJ_1; }
    SimTK::Real getPhaseInDegrees_1() const { return torsionSpec.phaseInDegrees_1; }
    int getPeriodicity_1() const { return torsionSpec.periodicity_1; }

    SimTK::Real getAmpInKJ_2() const { return torsionSpec.ampInKJ_2; }
    SimTK::Real getPhaseInDegrees_2() const { return torsionSpec.phaseInDegrees_2; }
    int getPeriodicity_2() const { return torsionSpec.periodicity_2; }

    SimTK::Real getAmpInKJ_3() const { return torsionSpec.ampInKJ_3; }
    SimTK::Real getPhaseInDegrees_3() const { return torsionSpec.phaseInDegrees_3; }
    int getPeriodicity_3() const { return torsionSpec.periodicity_3; }

    SimTK::Real getAmpInKJ_4() const { return torsionSpec.ampInKJ_4; }
    SimTK::Real getPhaseInDegrees_4() const { return torsionSpec.phaseInDegrees_4; }
    int getPeriodicity_4() const { return torsionSpec.periodicity_4; }

    SimTK::Real getAmpInKJ_5() const { return torsionSpec.ampInKJ_5; }
    SimTK::Real getPhaseInDegrees_5() const { return torsionSpec.phaseInDegrees_5; }
    int getPeriodicity_5() const { return torsionSpec.periodicity_5; }

private:
    BondTorsionDefinition torsionSpec;
};
