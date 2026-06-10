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
    Topology(const SystemTopology& systemTopology,
             SimTK::CompoundSystem::CompoundIndex compoundIndex,
             SimTK::RootMobility rootMobility);

    ~Topology() override = default;

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
    [[nodiscard]] auto calcLogSineSqrGamma2(const SimTK::State& quatState) const -> SimTK::Real;

    [[nodiscard]] auto calcLogDetMBATGamma2Contribution(const SimTK::State& quatState) const -> SimTK::Real;

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
     * @brief Calculate all atom frames in top frame. It avoids calling
     * calcDefaultAtomFrameInCompoundFrame multiple times. This has to be called
     * every time the coordinates change though.
     * @param :
     * @return
     */
    void calcAtomsTopTransforms();

    /**
     * @brief
     * @return
     */
    void printTopTransforms();

    /**
     * @brief Get atom Top level transform from the existing Topology map
     * @param cAIx: atom Compound AtomIndex
     * @return Atom's Top level transform
     */
    auto getTopTransform(SimTK::Compound::AtomIndex cAIx) const -> const SimTK::Transform&;

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

    auto matchAtomTargetLocations(const SimTK::Compound::AtomTargetLocations& atomTargets)
        -> SimTK::Transform;
    auto getMatchError(const SimTK::Compound::AtomTargetLocations& atomTargets) -> SimTK::Real;

    void writeAtomListPdb(std::string dirname,
                          std::string prefix,
                          std::string suffix,
                          int maxNofDigits,
                          int index) const;

    /**
     * @brief Create a mapping between the local compound atom indices and the global atom indices.
     */
    void loadIndicesMaps(const SimTK::Compound::AtomTargetLocations& atomTargets);

    /**
     * @brief Get the global atom index from the local compound atom index.
     * @param cAIx Compound Atom Index
     * @return Global Atom Index
     */
    auto getGlobalAtomIndex(SimTK::Compound::AtomIndex cAIx) -> int;

    /** Print atom to MobilizedBodyIndex and bond to Compound::Bond index
     * maps **/
    void printMaps();

    auto getAtomFrameCache() const -> const std::vector<SimTK::Transform>& {
        return atomFrameCache;
    }
    auto updAtomFrameCache() -> std::vector<SimTK::Transform>& {
        return atomFrameCache;
    }

    private:
    std::unordered_map<std::string, std::size_t> atomName2bond, atomName2angle;

    // Map aIx to its Transform Default top transform
    std::vector<SimTK::Transform> aIx2TopTransform;

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
};
