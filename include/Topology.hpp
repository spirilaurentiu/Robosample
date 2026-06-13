#pragma once

#include <cstddef>
#include <unordered_map>

#include "CompoundSystem.h"
#include "TopologyElements.hpp"

using CompoundAtomIndexPair = std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>;

/** Topological information (bonds graph) for one molecule.
It maps to one compound in Molmodel thus it is derived
from Molmodel Compound class.
It does the following things:
   - loads information from input files such as Amber input prmtop / inpcrd
   - adds parameters to a DuMM force field which belongs to the World class
     because one DuMM class should be used for multiple molecules
   - constructs the graph based on a list of Atom objects each of
     which already contains bonding information from the input files
   - defines the rigid bodies based on input files provided by the users.
Contains a list of atoms bAtomList which consists of Atom
objects **/
class Topology : public SimTK::Compound {
    public:
    Topology(const SimTK::Compound::Name& name,
             SimTK::CompoundSystem::CompoundIndex compoundIndex,
             int rootGlobalAtomIx,
             SimTK::RootMobility rootMobility);

    ~Topology() override = default;

    void setAtoms(Span<RoboAtom> atoms);
    [[nodiscard]] auto getAtoms() const -> Span<RoboAtom> {
        return subAtomList;
    }
    [[nodiscard]] auto updAtoms() -> Span<RoboAtom> {
        return subAtomList;
    }

    void setBonds(Span<RoboBond> bonds);
    [[nodiscard]] auto getBonds() const -> Span<RoboBond> {
        return subBondList;
    }
    [[nodiscard]] auto updBonds() -> Span<RoboBond> {
        return subBondList;
    }
    auto getBondByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0, SimTK::Compound::AtomIndex cAIx1) const
        -> const RoboBond&;

    void setAngles(Span<RoboAngle> angles);
    [[nodiscard]] auto getAngles() const -> Span<RoboAngle> {
        return subAngleList;
    }
    [[nodiscard]] auto updAngles() -> Span<RoboAngle> {
        return subAngleList;
    }
    auto getAngleByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0,
                                     SimTK::Compound::AtomIndex cAIx1,
                                     SimTK::Compound::AtomIndex cAIx2) const -> const RoboAngle&;

    void setPeriodicTorsions(Span<RoboPeriodicTorsion> periodicTorsions);
    [[nodiscard]] auto getPeriodicTorsions() const -> Span<RoboPeriodicTorsion> {
        return subPeriodicTorsions;
    }
    [[nodiscard]] auto updPeriodicTorsions() -> Span<RoboPeriodicTorsion> {
        return subPeriodicTorsions;
    }

    void setImproperHarmonicTorsions(Span<RoboHarmonicImproperTorsion> improperHarmonicTorsions);
    [[nodiscard]] auto getImproperHarmonicTorsions() const -> Span<RoboHarmonicImproperTorsion> {
        return subImproperHarmonicTorsions;
    }
    [[nodiscard]] auto updImproperHarmonicTorsions() -> Span<RoboHarmonicImproperTorsion> {
        return subImproperHarmonicTorsions;
    }

    [[nodiscard]] auto getRootMobility() const -> SimTK::RootMobility {
        return rootMobility;
    }

    /**
     * @brief Get the name of this molecule
     * @return name of the molecule
     */
    [[nodiscard]] auto getName() const -> const std::string& {
        return this->name;
    }

    /**
     * @brief Get own CompoundIndex in CompoundSystem
     *
     * This is equivalent to the index of the molecule in the CompoundSystem's vector of Compounds.
     *
     * @return CompoundIndex
     */
    [[nodiscard]] auto getCompoundIndex() const -> SimTK::CompoundSystem::CompoundIndex {
        return compoundIndex;
    }

    /**
     * @brief Computes log(sin^2(pitch)) for the root atom's orientation in ground frame.
     *
     * Extracts the pitch angle from the atom’s quaternion orientation.
     * The result is numerically stabilized near pitch ≈ 0 or ±π,
     * where sin(pitch) → 0 and log(sin²(pitch)) would diverge.
     *
     * @note This quantity may represent an orientation regularization term.
     * @return log(sin²(pitch)), computed safely with analytic limits near singularities.
     */
    [[nodiscard]] auto calcLogSineSqrGamma2(const SimTK::State& quatState,
                                            const SimTK::SimbodyMatterSubsystem& matter) const
        -> SimTK::Real {
        // Get atom transform and convert to quaternion (w,x,y,z)
        const SimTK::Transform X =
            calcAtomFrameInGroundFrameThroughSimbody(quatState, rootCompoundAtomIx, matter);
        const SimTK::Quaternion quat = X.R().convertRotationToQuaternion();

        const SimTK::Real w = quat[0];
        const SimTK::Real x = quat[1];
        const SimTK::Real y = quat[2];
        const SimTK::Real z = quat[3];

        // Compute sin(pitch) = 2(wy - zx)
        SimTK::Real sinPitch = 2.0 * ((w * y) - (z * x));

        // Clamp to account for floating-point drift outside [-1, 1]
        sinPitch = std::clamp(sinPitch, SimTK::Real(-1.0), SimTK::Real(1.0));

        // Compute pitch and evaluate the stable log(sin²)
        const SimTK::Real pitch = std::asin(sinPitch);
        return safeLogSineSqr(pitch);
    }

    [[nodiscard]] auto calcLogDetMBATGamma2Contribution(const SimTK::State& quatState,
                                                        const SimTK::SimbodyMatterSubsystem& matter) const
        -> SimTK::Real {
        SimTK::Transform X = calcAtomFrameInGroundFrameThroughSimbody(quatState, rootCompoundAtomIx, matter);
        SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();

        SimTK::Real w = quat[0];
        SimTK::Real x = quat[1];
        SimTK::Real y = quat[2];
        SimTK::Real z = quat[3];
        SimTK::Real sinPitch = 2 * (w * y - z * x);

        std::cout << std::setprecision(20) << std::fixed;
        std::cout << "sinpitch " << sinPitch << std::endl;
        std::cout << "sinpitchsq" << sinPitch * sinPitch << std::endl;

        SimTK::Real pitch = std::asin(sinPitch);
        std::cout << "pitch " << pitch << std::endl;
        // if(pitch < 0){
        //	pitch = pitch + (2*SimTK_PI);
        //	std::cout << "sin converted pitch " << std::sin(pitch) << std::endl;
        // }

        if (sinPitch < SimTK::Eps) { // consider using SimTK::Eps
            return -SimTK::Infinity;
        }
        SimTK::Real result = std::log(sinPitch * sinPitch);
        return result;
    }

    /**
     * @brief Get a reference to the atom object in the atom list of this Compound.
     *
     * @param cAIx Compound Atom Index. This is in range [0, num_atoms-1] for this Compound. Not to confuse
     * with the global atom index.
     * @return Reference to the Atom object.
     */
    [[nodiscard]] const RoboAtom& getAtom(SimTK::Compound::AtomIndex cAIx) const;

    /**
     * @brief Get the bonded neighbor atom in the parent mobilized body.
     * @param aIx Compound Atom Index
     * @return Compound atom index of the root
     */
    auto getChemicalParentOfMobodRootAtom(SimTK::Compound::AtomIndex aIx,
                                          const SimTK::SimbodyMatterSubsystem& matter,
                                          const SimTK::DuMMForceFieldSubsystem& dumm) const
        -> SimTK::Compound::AtomIndex;

    /**
     * @brief Get atom Top level transform from the existing Topology map
     * @param cAIx: atom Compound AtomIndex
     * @return Atom's Top level transform
     */
    auto getTopTransform(SimTK::Compound::AtomIndex cAIx) const -> const SimTK::Transform& {
        return aIx2TopTransform[cAIx];
    }

    // Return mbx by calling DuMM functions
    auto getAtomMobilizedBodyIndexThroughDumm(SimTK::Compound::AtomIndex aIx,
                                              const SimTK::DuMMForceFieldSubsystem& dumm) const
        -> SimTK::MobilizedBodyIndex;

    // Get atom location on mobod through DuMM functions
    auto getAtomLocationInMobilizedBodyFrameThroughDumm(SimTK::Compound::AtomIndex aIx,
                                                        const SimTK::DuMMForceFieldSubsystem& dumm) const
        -> SimTK::Vec3;

    auto calcAtomLocationInGroundFrameThroughSimbody(SimTK::Compound::AtomIndex aIx,
                                                     const SimTK::DuMMForceFieldSubsystem& dumm,
                                                     const SimTK::SimbodyMatterSubsystem& matter,
                                                     const SimTK::State& someState) const -> SimTK::Vec3;

    /**
     * @brief Create a mapping between the local compound atom indices and the global atom indices.
     */
    void loadIndicesMaps(const SimTK::Compound::AtomTargetLocations& atomTargets);

    auto getAtomFrameCache() const -> const std::vector<SimTK::Transform>& {
        return atomFrameCache;
    }
    auto updAtomFrameCache() -> std::vector<SimTK::Transform>& {
        return atomFrameCache;
    }

    private:
    Span<RoboAtom> subAtomList;
    Span<RoboBond> subBondList;
    Span<RoboAngle> subAngleList;
    Span<RoboPeriodicTorsion> subPeriodicTorsions;
    Span<RoboHarmonicImproperTorsion> subImproperHarmonicTorsions;

    std::unordered_map<std::string, std::size_t> atomName2bond, atomName2angle;


    // Map Atom number to aIx
    std::vector<int> compound2GlobalAtomIndex;
    std::vector<std::pair<int, SimTK::Compound::AtomIndex>> global2CompoundAtomIndex;


    std::vector<std::pair<CompoundAtomIndexPair, int>> aIxPair2Bonds;

    // Atom frames in Top frame
    std::vector<SimTK::Transform> atomFrameCache;

    std::string name;

    /** Every Compound has an index which is the position in the dvector
     * of Compounds in CompoundSystem
     */
    SimTK::CompoundSystem::CompoundIndex compoundIndex;

    std::size_t rootGlobalAtomIx = 0;              // in global atom index
    SimTK::Compound::AtomIndex rootCompoundAtomIx; // in subAtomList index

    bool flipAllChirality = false;
    SimTK::RootMobility rootMobility = SimTK::RootMobility::Weld;

    // Map aIx to its Transform Default top transform
    std::vector<SimTK::Transform> aIx2TopTransform;

    auto calcAtomFrameInGroundFrameThroughSimbody(const SimTK::State& state,
                                                  SimTK::Compound::AtomIndex cAIx,
                                                  const SimTK::SimbodyMatterSubsystem& matter) const
        -> SimTK::Transform {
        // Body this atom is welded to (the value World wrote via setAtomMobilizedBodyIndex).
        const SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(cAIx);
        const SimTK::MobilizedBody& body = matter.getMobilizedBody(mbx);

        // G_X_B from Simbody, B_X_atom from the stored per-atom frame-in-body.
        const SimTK::Transform& G_X_B = body.getBodyTransform(state);
        const SimTK::Transform& B_X_atom = getFrameInMobilizedBodyFrame(cAIx);

        return G_X_B * B_X_atom;
    }
};
