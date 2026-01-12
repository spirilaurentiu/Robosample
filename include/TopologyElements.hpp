#pragma once

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"


struct RoboAtomDefinition {
    // Indices
    int globalIndex = 0;
    int prmtopIndex = 0;
    int compoundAtomIndex = 0;
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
    SimTK::Real solventRadiusInNm = 0.0;
    SimTK::Real screen = 0.0;
    SimTK::Real x_nm = 0.0, y_nm = 0.0, z_nm = 0.0;
};

struct RoboBondStretchDefinition {
    int parentAtomGlobalIndex = 0, childAtomGlobalIndex = 0;
    int parentCompoundAtomIndex = 0, childCompoundAtomIndex = 0;

    int bondGlobalIndex = 0;
    int moleculeIndex = 0;
	bool ringClosing = false;

    SimTK::Real stiffnessInKJPerNmSq = 0.0;
	SimTK::Real nominalLengthInNm = 0.0;
};

struct RoboBondBendDefinition {
    int globalIndex1 = 0, globalIndex2 = 0, globalIndex3 = 0;
    int compoundAtomIndex1 = 0, compoundAtomIndex2 = 0, compoundAtomIndex3 = 0;
    int moleculeIndex = 0;
	SimTK::Real stiffnessInKJPerRadSq = 0.0;
    SimTK::Real nominalAngleInDeg = 0.0;
};

struct RoboBondTorsionDefinition {
    int globalIndex1 = -1, globalIndex2 = -1, globalIndex3 = -1, globalIndex4 = -1;
    int compoundAtomIndex1 = -1, compoundAtomIndex2 = -1, compoundAtomIndex3 = -1, compoundAtomIndex4 = -1;
    int moleculeIndex = 0;
    bool improper = false;
    SimTK::Real ampInKJ_1 = -1.0, phaseInDegrees_1 = -1.0, periodicity_1 = -1;
    SimTK::Real ampInKJ_2 = -1.0, phaseInDegrees_2 = -1.0, periodicity_2 = -1;
    SimTK::Real ampInKJ_3 = -1.0, phaseInDegrees_3 = -1.0, periodicity_3 = -1;
    SimTK::Real ampInKJ_4 = -1.0, phaseInDegrees_4 = -1.0, periodicity_4 = -1;
    SimTK::Real ampInKJ_5 = -1.0, phaseInDegrees_5 = -1.0, periodicity_5 = -1;
};

class RoboAtom {
public:
    RoboAtom() = default;
    
    RoboAtom(const RoboAtomDefinition& spec) : atomSpec(spec) {
        CoordsInNm = SimTK::Vec3(atomSpec.x_nm, atomSpec.y_nm, atomSpec.z_nm);
        availableBonds = atomSpec.neighborsGlobalIndices.size();
        element = SimTK::Element(getAtomicNumber(), getElementName(), getElementSymbol(), getMassInDaltons());
    }

    int getGlobalIndex() const { return atomSpec.globalIndex; }
    int getPrmtopIndex() const { return atomSpec.prmtopIndex; }
    SimTK::Compound::AtomIndex getCompoundAtomIndex() const { return SimTK::Compound::AtomIndex(atomSpec.compoundAtomIndex); }

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

    int getAtomicNumber() const { return atomSpec.atomicNumber; }
    SimTK::Real getChargeInE() const { return atomSpec.chargeInE; }
    SimTK::mdunits::Mass getMassInDaltons() const { return atomSpec.massInDaltons; }
    SimTK::Real getVdwRadiusInNm() const { return atomSpec.vdwRadiusInNm; }
    SimTK::Real getSigmaInNm() const { return atomSpec.sigmaInNm; }
    SimTK::Real getVdwWellDepthInKJ() const { return atomSpec.vdwWellDepthInKJ; }
    SimTK::Real getSolventRadiusInNm() const { return atomSpec.solventRadiusInNm; }
    SimTK::Real getScreen() const { return atomSpec.screen; }

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

    const SimTK::Compound::SingleAtom& getSingleAtom() const { return *compoundSingleAtom; }
    void createSingleAtom();

    bool isRoot() const { return atomSpec.root; }

    const SimTK::Element& getElement() const { return element; }

private:

    RoboAtomDefinition atomSpec;
    SimTK::Element element;

    int availableBonds = 0;

    SimTK::BiotypeIndex biotypeIndex;

    SimTK::Vec3 CoordsInNm = SimTK::Vec3(SimTK::NaN);

    // wasted 6 hours trying to make this unique_ptr or allocated on the stack
    // after more hours wasted, i read this article: https://www.cppstories.com/2014/05/vector-of-objects-vs-vector-of-pointers/
    // as expected, it should be faster to have this allocated on the stack
    // however, taking into account that this is called only a few times, it is not worth the trouble
    // there are two use cases for this: addition of bond centers and building the molecule graph
    // in the first case, the calls are sequential, so we waste some time
    // in the second case, the calls are random and the article shows that it is actually better to have the objects allocated on the heap
    SimTK::Compound::SingleAtom* compoundSingleAtom = nullptr;
};


class RoboBondStretch {
public:
    RoboBondStretch() = default;

    RoboBondStretch(const RoboBondStretchDefinition& spec) : bondSpec(spec) {}

    bool operator==(const RoboBondStretch& other) const {
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

    SimTK::Compound::AtomIndex getParentCompoundAtomIndex() const { return SimTK::Compound::AtomIndex(bondSpec.parentCompoundAtomIndex); }
    SimTK::Compound::AtomIndex getChildCompoundAtomIndex() const { return SimTK::Compound::AtomIndex(bondSpec.childCompoundAtomIndex); }

    int getBondGlobalIndex() const { return bondSpec.bondGlobalIndex; }
    int getMoleculeIndex() const { return bondSpec.moleculeIndex; }
    bool isRingClosing() const { return bondSpec.ringClosing; }

    SimTK::Real getStiffnessInKJPerNmSq() const { return bondSpec.stiffnessInKJPerNmSq; }
    SimTK::Real getNominalLengthInNm() const { return bondSpec.nominalLengthInNm; }

private:
    RoboBondStretchDefinition bondSpec;

	std::vector<SimTK::BondMobility::Mobility> mobilities;
	// std::vector<SimTK::Real> uScaleFactors = { 1.0f };
};

class RoboBondBend {
public:
    RoboBondBend() = default;

    RoboBondBend(const RoboBondBendDefinition & spec) : angleSpec(spec) {}

    int getGlobalIndex1() const { return angleSpec.globalIndex1; }
    int getGlobalIndex2() const { return angleSpec.globalIndex2; }
    int getGlobalIndex3() const { return angleSpec.globalIndex3; }

    SimTK::Compound::AtomIndex getCompoundAtomIndex1() const { return SimTK::Compound::AtomIndex(angleSpec.compoundAtomIndex1); }
    SimTK::Compound::AtomIndex getCompoundAtomIndex2() const { return SimTK::Compound::AtomIndex(angleSpec.compoundAtomIndex2); }
    SimTK::Compound::AtomIndex getCompoundAtomIndex3() const { return SimTK::Compound::AtomIndex(angleSpec.compoundAtomIndex3); }

    int getMoleculeIndex() const { return angleSpec.moleculeIndex; }

    SimTK::Real getStiffnessInKJPerRadSq() const { return angleSpec.stiffnessInKJPerRadSq; }
    SimTK::Real getNominalAngleInDeg() const { return angleSpec.nominalAngleInDeg; }

private:
	RoboBondBendDefinition angleSpec;
};

class RoboBondTorsion {
public:
    RoboBondTorsion() = default;

    RoboBondTorsion(const RoboBondTorsionDefinition& spec) : torsionSpec(spec) {}

    int getGlobalIndex1() const { return torsionSpec.globalIndex1; }
    int getGlobalIndex2() const { return torsionSpec.globalIndex2; }
    int getGlobalIndex3() const { return torsionSpec.globalIndex3; }
    int getGlobalIndex4() const { return torsionSpec.globalIndex4; }

    SimTK::Compound::AtomIndex getCompoundAtomIndex1() const { return SimTK::Compound::AtomIndex(torsionSpec.compoundAtomIndex1); }
    SimTK::Compound::AtomIndex getCompoundAtomIndex2() const { return SimTK::Compound::AtomIndex(torsionSpec.compoundAtomIndex2); }
    SimTK::Compound::AtomIndex getCompoundAtomIndex3() const { return SimTK::Compound::AtomIndex(torsionSpec.compoundAtomIndex3); }
    SimTK::Compound::AtomIndex getCompoundAtomIndex4() const { return SimTK::Compound::AtomIndex(torsionSpec.compoundAtomIndex4); }

    int getMoleculeIndex() const { return torsionSpec.moleculeIndex; }

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
    RoboBondTorsionDefinition torsionSpec;
};
