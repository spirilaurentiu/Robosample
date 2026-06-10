#pragma once

#include <array>

#include "CompoundSystem.h"
#include "Molmodel.h"
#include "Simbody.h"
#include "bgeneral.hpp"

struct RoboAtomPhysics {
    SimTK::Real chargeInE = SimTK::NaN;
    SimTK::Real massInDaltons = SimTK::NaN;
    SimTK::Real vdwRadiusInNm = SimTK::NaN;
    SimTK::Real vdwWellDepthInKJ = SimTK::NaN;
    SimTK::Real sigmaInNm = SimTK::NaN;
    SimTK::Real solventRadiusInNm = SimTK::NaN;
    SimTK::Real screen = SimTK::NaN;

    RoboAtomPhysics() = default;
    RoboAtomPhysics(SimTK::Real chargeInE_,
                    SimTK::Real massInDaltons_,
                    SimTK::Real vdwRadiusInNm_,
                    SimTK::Real vdwWellDepthInKJ_,
                    SimTK::Real sigmaInNm_,
                    SimTK::Real solventRadiusInNm_,
                    SimTK::Real screen_)
        : chargeInE(chargeInE_)
        , massInDaltons(massInDaltons_)
        , vdwRadiusInNm(vdwRadiusInNm_)
        , vdwWellDepthInKJ(vdwWellDepthInKJ_)
        , sigmaInNm(sigmaInNm_)
        , solventRadiusInNm(solventRadiusInNm_)
        , screen(screen_) {
    }
};

struct RoboAtomIdentity {
    std::string uniqueAtomName;
    std::string residueName;
    std::string atomClassName;
    std::string chargedAtomTypeName;

    int globalIndex = std::numeric_limits<int>::min();
    int prmtopIndex = std::numeric_limits<int>::min();
    int moleculeIndex = std::numeric_limits<int>::min();
    int residueIndex = std::numeric_limits<int>::min();

    int nonbondedIndex = std::numeric_limits<int>::min();

    SimTK::Compound::AtomIndex compoundAtomIndex;
    SimTK::DuMM::AtomClassIndex atomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;

    RoboAtomIdentity() = default;
    RoboAtomIdentity(std::string uniqueAtomName_,
                     std::string residueName_,
                     std::string atomClassName_,
                     std::string chargedAtomTypeName_,
                     int globalIndex_,
                     int prmtopIndex_,
                     int moleculeIndex_,
                     int residueIndex_,
                     int nonbondedIndex_,
                     int compoundAtomIndex_,
                     int atomClassIndex_,
                     int chargedAtomTypeIndex_)
        : uniqueAtomName(uniqueAtomName_)
        , residueName(residueName_)
        , atomClassName(atomClassName_)
        , chargedAtomTypeName(chargedAtomTypeName_)
        , globalIndex(globalIndex_)
        , prmtopIndex(prmtopIndex_)
        , moleculeIndex(moleculeIndex_)
        , residueIndex(residueIndex_)
        , nonbondedIndex(nonbondedIndex_)
        , compoundAtomIndex(SimTK::Compound::AtomIndex(compoundAtomIndex_))
        , atomClassIndex(SimTK::DuMM::AtomClassIndex(atomClassIndex_))
        , chargedAtomTypeIndex(SimTK::DuMM::ChargedAtomTypeIndex(chargedAtomTypeIndex_)) {
    }
};

struct RoboAtomElement {
    SimTK::Element element;
    std::string elementName;
    std::string elementSymbol;
    int atomicNumber = std::numeric_limits<int>::min();

    RoboAtomElement() = default;
    RoboAtomElement(std::string elementName_, std::string elementSymbol_, int atomicNumber_)
        : elementName(elementName_)
        , elementSymbol(elementSymbol_)
        , atomicNumber(atomicNumber_) {
    }
};

struct RoboAtomConnectivity {
    std::vector<int> neighborsGlobalIndices;
    int numAvailableBonds = std::numeric_limits<int>::min();
    int numBondsInvolved = std::numeric_limits<int>::min();
    bool root = false;

    RoboAtomConnectivity() = default;
    RoboAtomConnectivity(std::vector<int> neighborsGlobalIndices_, bool root_)
        : neighborsGlobalIndices(neighborsGlobalIndices_)
        , root(root_) {
    }
};

struct RoboAtom {
    RoboAtomIdentity identity;
    RoboAtomElement elementInfo;
    RoboAtomPhysics physics;
    RoboAtomConnectivity connectivity;

    SimTK::Vec3 position{SimTK::NaN};

    SimTK::BiotypeIndex biotypeIndex;
    SimTK::Compound::SingleAtom* compoundSingleAtom = nullptr;

    RoboAtom() = default;
    RoboAtom(RoboAtomIdentity identity_,
             RoboAtomElement elementInfo_,
             RoboAtomPhysics physics_,
             RoboAtomConnectivity connectivity_,
             std::array<SimTK::Real, 3> position_)
        : identity(identity_)
        , elementInfo(elementInfo_)
        , physics(physics_)
        , connectivity(connectivity_) {
        connectivity.numAvailableBonds = connectivity.neighborsGlobalIndices.size();
        connectivity.numBondsInvolved = connectivity.neighborsGlobalIndices.size();
        elementInfo.element = SimTK::Element(elementInfo.atomicNumber,
                                             elementInfo.elementName,
                                             elementInfo.elementSymbol,
                                             physics.massInDaltons);
        position = SimTK::Vec3(position_[0], position_[1], position_[2]);
    }

    void createSingleAtom();
};

struct RoboBond {
    std::vector<SimTK::BondMobility::Mobility> mobilities;
    std::array<SimTK::Compound::AtomIndex, 2> compoundAtomIndices;
    std::string dihedralType;

    std::array<int, 2> globalIndices = {-1, -1};
    std::array<int, 2> prmtopIndices = {-1, -1};
    SimTK::Real stiffnessInKJPerNmSq = 0.0;
    SimTK::Real nominalLengthInNm = 0.0;
    int moleculeIndex = 0;
    bool ringClosing = false;

    RoboBond() = default;
    RoboBond(std::array<int, 2> globalIndices_,
             std::array<int, 2> prmtopIndices_,
             std::array<int, 2> compoundAtomIndices_,
             SimTK::Real stiffnessInKJPerNmSq_,
             SimTK::Real nominalLengthInNm_,
             int moleculeIndex_,
             bool ringClosing_,
             const std::string& dihedralType_)
        : globalIndices(globalIndices_)
        , prmtopIndices(prmtopIndices_)
        , stiffnessInKJPerNmSq(stiffnessInKJPerNmSq_)
        , nominalLengthInNm(nominalLengthInNm_)
        , moleculeIndex(moleculeIndex_)
        , ringClosing(ringClosing_)
        , dihedralType(dihedralType_) {
        compoundAtomIndices[0] = SimTK::Compound::AtomIndex(compoundAtomIndices_[0]);
        compoundAtomIndices[1] = SimTK::Compound::AtomIndex(compoundAtomIndices_[1]);
    }

    void addBondMobility(SimTK::BondMobility::Mobility someMobility) {
        mobilities.push_back(someMobility);
    }
    void setBondMobility(SimTK::BondMobility::Mobility someMobility, int world) {
        mobilities[world] = someMobility;
    }
    SimTK::BondMobility::Mobility getBondMobility(int world) const {
        return mobilities[world];
    }
};

struct RoboAngle {
    std::array<SimTK::Compound::AtomIndex, 3> compoundAtomIndices;

    std::array<int, 3> globalIndices{-1, -1, -1};
    std::array<int, 3> prmtopIndices{-1, -1, -1};
    int moleculeIndex = 0;
    SimTK::Real stiffnessInKJPerRadSq = 0.0;
    SimTK::Real nominalAngleInDeg = 0.0;

    RoboAngle() = default;
    RoboAngle(std::array<int, 3> globalIndices_,
              std::array<int, 3> prmtopIndices_,
              std::array<int, 3> compoundAtomIndices_,
              int moleculeIndex_,
              SimTK::Real stiffnessInKJPerRadSq_,
              SimTK::Real nominalAngleInDeg_)
        : compoundAtomIndices{SimTK::Compound::AtomIndex(compoundAtomIndices_[0]),
                              SimTK::Compound::AtomIndex(compoundAtomIndices_[1]),
                              SimTK::Compound::AtomIndex(compoundAtomIndices_[2])}
        , globalIndices(globalIndices_)
        , prmtopIndices(prmtopIndices_)
        , moleculeIndex(moleculeIndex_)
        , stiffnessInKJPerRadSq(stiffnessInKJPerRadSq_)
        , nominalAngleInDeg(nominalAngleInDeg_) {
    }
};

struct RoboPeriodicTorsionTerm {
    SimTK::Real amplitudeKJ = 0.0;
    SimTK::Real phaseDeg = 0.0;
    int periodicity = -1;

    RoboPeriodicTorsionTerm() = default;
    RoboPeriodicTorsionTerm(SimTK::Real amplitudeKJ_, SimTK::Real phaseDeg_, int periodicity_)
        : amplitudeKJ(amplitudeKJ_)
        , phaseDeg(phaseDeg_)
        , periodicity(periodicity_) {
    }
};

struct RoboPeriodicTorsion {
    std::array<RoboPeriodicTorsionTerm, 5> terms;
    std::array<int, 4> globalIndices{-1, -1, -1, -1};
    std::array<int, 4> prmtopIndices{-1, -1, -1, -1};
    std::array<SimTK::Compound::AtomIndex, 4> compoundAtomIndices;
    int moleculeIndex = -1;
    int numTerms = 0;
    bool improper = false;

    RoboPeriodicTorsion(std::array<int, 4> globalIndices_,
                        std::array<int, 4> prmtopIndices_,
                        std::array<int, 4> compoundAtomIndices_,
                        int moleculeIndex_,
                        bool isImproper,
                        std::vector<RoboPeriodicTorsionTerm> inputTerms)
        : globalIndices(globalIndices_)
        , prmtopIndices(prmtopIndices_)
        , moleculeIndex(moleculeIndex_)
        , improper(isImproper) {
        if (inputTerms.size() > 5) {
            throw std::runtime_error("Periodic torsion supports at most 5 terms");
        }

        numTerms = static_cast<int>(inputTerms.size());

        for (std::size_t i = 0; i < inputTerms.size(); ++i) {
            terms[i] = inputTerms[i];
        }

        for (int i = 0; i < 4; ++i) {
            compoundAtomIndices[i] = SimTK::Compound::AtomIndex(compoundAtomIndices_[i]);
        }
    }
};

struct RoboHarmonicImproperTorsion {
    std::array<int, 4> globalIndices{-1, -1, -1, -1};
    std::array<int, 4> prmtopIndices{-1, -1, -1, -1};
    std::array<SimTK::Compound::AtomIndex, 4> compoundAtomIndices;
    int moleculeIndex = -1;
    SimTK::Real stiffnessInKJPerRadSq = 0.0;
    SimTK::Real nominalAngleInRad = 0.0;

    RoboHarmonicImproperTorsion(std::array<int, 4> globalIndices_,
                                std::array<int, 4> prmtopIndices_,
                                std::array<int, 4> compoundAtomIndices_,
                                int moleculeIndex_,
                                SimTK::Real stiffnessInKJPerRadSq_,
                                SimTK::Real nominalAngleInRad_)
        : globalIndices(globalIndices_)
        , prmtopIndices(prmtopIndices_)
        , moleculeIndex(moleculeIndex_)
        , stiffnessInKJPerRadSq(stiffnessInKJPerRadSq_)
        , nominalAngleInRad(nominalAngleInRad_) {
        for (int i = 0; i < 4; ++i) {
            compoundAtomIndices[i] = SimTK::Compound::AtomIndex(compoundAtomIndices_[i]);
        }
    }
};

class TopologyRange {
    // Stores [begin, end) pairs for each type
    std::array<std::pair<int, int>, (int)TopologyRangeType::NofTopologyRangeTypes> ranges;

    public:
    TopologyRange(std::vector<int> startCounts) {
        ranges[(int)TopologyRangeType::Atom] = {startCounts[0], startCounts[0]};
        ranges[(int)TopologyRangeType::Bond] = {startCounts[1], startCounts[1]};
        ranges[(int)TopologyRangeType::Angle] = {startCounts[2], startCounts[2]};
        ranges[(int)TopologyRangeType::PeriodicTorsion] = {startCounts[3], startCounts[3]};
        ranges[(int)TopologyRangeType::ImproperHarmonicTorsion] = {startCounts[4], startCounts[4]};
    }

    void close(std::vector<int> endCounts) {
        ranges[(int)TopologyRangeType::Atom].second = endCounts[0];
        ranges[(int)TopologyRangeType::Bond].second = endCounts[1];
        ranges[(int)TopologyRangeType::Angle].second = endCounts[2];
        ranges[(int)TopologyRangeType::PeriodicTorsion].second = endCounts[3];
        ranges[(int)TopologyRangeType::ImproperHarmonicTorsion].second = endCounts[4];
    }

    [[nodiscard]] auto getRange(TopologyRangeType type) const -> const std::pair<int, int>& {
        return ranges[(int)type];
    }
};

struct ZMatrixRow {
    std::array<int, 4> globalIndices{-1, -1, -1, -1};
    std::array<SimTK::Compound::AtomIndex, 4> compoundAtomIndices;
    int moleculeIndex = -1;

    ZMatrixRow(std::array<int, 4> globalIndices_, std::array<int, 4> compoundAtomIndices_, int moleculeIndex_)
        : globalIndices(globalIndices_)
        , moleculeIndex(moleculeIndex_) {
        compoundAtomIndices[0] = SimTK::Compound::AtomIndex(compoundAtomIndices_[0]);

        if (compoundAtomIndices_[1] != -1) {
            compoundAtomIndices[1] = SimTK::Compound::AtomIndex(compoundAtomIndices_[1]);
        }
        if (compoundAtomIndices_[2] != -1) {
            compoundAtomIndices[2] = SimTK::Compound::AtomIndex(compoundAtomIndices_[2]);
        }
        if (compoundAtomIndices_[3] != -1) {
            compoundAtomIndices[3] = SimTK::Compound::AtomIndex(compoundAtomIndices_[3]);
        }
    }
};

using ZMatrix = std::vector<ZMatrixRow>;

enum NonbondedMethod : std::uint8_t {
    NoCutoff = 0,
    CutoffNonPeriodic,
};

struct CMAPGrid {
    std::vector<SimTK::Real> energy;
    int size = std::numeric_limits<int>::max();
};

struct CMAPTorsion {
    int mapIndex = std::numeric_limits<int>::max();

    // Torsion A atoms
    int torsionAAtom1GlobalIndex = std::numeric_limits<int>::max();
    int torsionAAtom2GlobalIndex = std::numeric_limits<int>::max();
    int torsionAAtom3GlobalIndex = std::numeric_limits<int>::max();
    int torsionAAtom4GlobalIndex = std::numeric_limits<int>::max();

    // Torsion B atoms
    int torsionBAtom1GlobalIndex = std::numeric_limits<int>::max();
    int torsionBAtom2GlobalIndex = std::numeric_limits<int>::max();
    int torsionBAtom3GlobalIndex = std::numeric_limits<int>::max();
    int torsionBAtom4GlobalIndex = std::numeric_limits<int>::max();
};

struct UreyBradley {
    int atom1GlobalIndex = std::numeric_limits<int>::max();
    int atom3GlobalIndex = std::numeric_limits<int>::max();
    SimTK::Real stiffnessInKJPerNmSq = SimTK::NaN;
    SimTK::Real nominalLengthInNm = SimTK::NaN;
};

struct Scaling14 {
    int atom1GlobalIndex = std::numeric_limits<int>::max();
    int atom4GlobalIndex = std::numeric_limits<int>::max();
    SimTK::Real chargeProduct = SimTK::NaN;
    SimTK::Real epsilon = SimTK::NaN;
    SimTK::Real sigma = SimTK::NaN;
};

struct Exclusion {
    int atom1GlobalIndex = std::numeric_limits<int>::max();
    int atom2GlobalIndex = std::numeric_limits<int>::max();
};

struct SystemTopology {
    std::vector<int> rootAtomGlobalIndices;
    std::vector<TopologyRange> topologyRanges;
    std::vector<RoboAtom> atoms;
    std::vector<RoboBond> bonds;
    std::vector<RoboAngle> angles;
    std::vector<RoboPeriodicTorsion> periodicTorsions;
    std::vector<RoboHarmonicImproperTorsion> harmonicImproperTorsions;
    std::vector<CMAPGrid> cmapGrids;
    std::vector<CMAPTorsion> cmapTorsions;
    std::vector<UreyBradley> ureyBradleys;
    std::vector<Scaling14> scaling14s;
    std::vector<Exclusion> exclusions;
    std::vector<SimTK::RootMobility> rootMobilities;
};

struct ForceFieldParams {
    bool hasNBfix = false;
    int numTypes = 0;
    std::vector<SimTK::Real> aCoef;
    std::vector<SimTK::Real> bCoef;

    bool useGBSAOBC2 = true;
    SimTK::Real gbsaSolventDielectric = 78.5;
    SimTK::Real gbsaSoluteDielectric = 1.0;

    NonbondedMethod nonbondedMethod = NonbondedMethod::NoCutoff;
    SimTK::Real nonbondedCutoffInNm = 1.2;
};

struct SimulationSettings {
    SimTK::Real thermostatTemperatureInK = 300.0;
    SimTK::Real collisionFrequency = 1.0;
    int seed = 0;
};

using CanonicalBond = std::pair<std::size_t, std::size_t>;
using CanonicalAngle = std::array<std::size_t, 3>;
using CanonicalTorsion = std::array<std::size_t, 4>;

[[nodiscard]] static inline std::pair<std::size_t, std::size_t> canonicalizeBond(std::size_t i,
                                                                                 std::size_t j) noexcept {
    return {std::min(i, j), std::max(i, j)};
}

[[nodiscard]] static inline auto canonicalizeAngle(std::size_t i, std::size_t j, std::size_t k) noexcept
    -> std::array<std::size_t, 3> {
    return {std::min(i, k), j, std::max(i, k)};
}

[[nodiscard]] static inline auto
canonicalizeTorsion(std::size_t i, std::size_t j, std::size_t k, std::size_t l) noexcept
    -> std::array<std::size_t, 4> {
    const std::array<std::size_t, 4> forward{i, j, k, l};
    const std::array<std::size_t, 4> reverse{l, k, j, i};
    return (forward < reverse) ? forward : reverse;
}
