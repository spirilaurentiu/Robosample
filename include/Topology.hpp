#pragma once

#include "TrivalentAtomTetra.hpp"
#include "TopologyElements.hpp"
#include "server.hpp"

using CompoundAtomIndexPair = std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>;

inline CompoundAtomIndexPair canonical(SimTK::Compound::AtomIndex a, SimTK::Compound::AtomIndex b) {
	return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

/** Topological information (bonds graph) for one molecule.
It maps to one compound in Molmodel thus it is derived 
from Molmodel Compound class.
It does the following things:
   - loads information from input files such as Amber input prmtop / inpcrd
   - adds parameters to a DuMM force field which belongs to the World class
	 because one DuMM class should be used for multiple molecules
   - contructs the graph based on a list of Atom objects each of 
	 which already contains bonding information from the input files
   - defines the rigid bodies based on imput files provided by the users.
Contains a list of atoms bAtomList which consists of Atom 
objects **/
class Topology : public SimTK::Compound {
public:

	Topology(const SimTK::Compound::Name& name, SimTK::CompoundSystem::CompoundIndex compoundIndex, int rootGlobalAtomIx);

	~Topology() override = default;

	void setAtoms(Span<RoboAtom> atoms) { subAtomList = atoms; }
	const Span<RoboAtom> getAtoms() const { return subAtomList; }
	Span<RoboAtom> updAtoms() { return subAtomList; }

	void setBonds(Span<RoboBondStretch> bonds) { subBondList = bonds; }
	const Span<RoboBondStretch> getBonds() const { return subBondList; }
	Span<RoboBondStretch> updBonds() { return subBondList; }

	void setAngles(Span<RoboBondBend> angles) { subAngleList = angles; }
	const Span<RoboBondBend> getAngles() const { return subAngleList; }
	Span<RoboBondBend> updAngles() { return subAngleList; }

	void setTorsions(Span<RoboBondTorsion> torsions) { subTorsionList = torsions; }
	const Span<RoboBondTorsion> getTorsions() const { return subTorsionList; }
	Span<RoboBondTorsion> updTorsions() { return subTorsionList; }

	/**	
	* @brief Get the name of this molecule
	* @return name of the molecule
	*/
	const std::string getName() const {return this->name;}

	/**	
	* @brief Get own CompoundIndex in CompoundSystem
	* 
	* This is equivalent to the index of the molecule in the CompoundSystem's vector of Compounds.
	* 
	* @return CompoundIndex
	*/
	inline SimTK::CompoundSystem::CompoundIndex getCompoundIndex() const {
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
	SimTK::Real calcLogSineSqrGamma2(const SimTK::State &quatState) const;

	SimTK::Real calcLogDetMBATGamma2Contribution(const SimTK::State& quatState) const;

	/**
	 * @brief Get a reference to the atom object in the atom list of this Compound.
	 * 
	 * @param cAIx Compound Atom Index. This is in range [0, num_atoms-1] for this Compound. Not to confuse with the global atom index.
	 * @return Reference to the Atom object.
	 */
	const RoboAtom& getAtom(SimTK::Compound::AtomIndex cAIx) const;

	/**
	 * @brief Get a reference to the bond object in the bond list of this Compound.
	 * 
	 * This is not the Compound Atom Index but the global atom index.
	 * 
	 * @param aIx0 Global Atom Index of one atom in the bond.
	 * @param aIx1 Global Atom Index of the other atom in the bond.
	 * @return Reference to the BondLink object.
	 */
	const RoboBondStretch& getBondByGlobalAtomIndex(int aIx0, int aIx1) const;

	/**
	 * @brief Get a reference to the bond object in the bond list of this compound using Compound Atom Indices.
	 * 
	 * This is not the global atom index but the local Compound Atom Index.
	 * 
	 * @param cAIx0 Compound Atom Index of one atom in the bond.
	 * @param cAIx1 Compound Atom Index of the other atom in the bond.
	 * 
	 * @return Reference to the BondLink object.
	 */
	const RoboBondStretch& getBondByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0, SimTK::Compound::AtomIndex cAIx1) const;

	/**	
	* @brief Get the bonded neighbor atom in the parent mobilized body.
	* @param aIx Compound Atom Index
	* @return Compound atom index of the root
	*/
	SimTK::Compound::AtomIndex getChemicalParentOfMobodRootAtom(SimTK::Compound::AtomIndex aIx, const SimTK::SimbodyMatterSubsystem& matter, const SimTK::DuMMForceFieldSubsystem& dumm) const;

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
	const SimTK::Transform& getTopTransform(SimTK::Compound::AtomIndex cAIx) const;

	// Interface to access the maps

	// Retunr mbx by calling DuMM functions
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexThroughDumm(SimTK::Compound::AtomIndex aIx, const SimTK::DuMMForceFieldSubsystem& dumm) const;

	// Get atom location on mobod through DuMM functions
	SimTK::Vec3 getAtomLocationInMobilizedBodyFrameThroughDumm(SimTK::Compound::AtomIndex aIx, const SimTK::DuMMForceFieldSubsystem& dumm) const;

	SimTK::Vec3 calcAtomLocationInGroundFrameThroughSimbody(SimTK::Compound::AtomIndex aIx, const SimTK::DuMMForceFieldSubsystem& dumm, const SimTK::SimbodyMatterSubsystem& matter, const SimTK::State& someState) const;

	SimTK::Transform matchAtomTargetLocations(const SimTK::Compound::AtomTargetLocations& atomTargets);
	SimTK::Real getMatchError(const SimTK::Compound::AtomTargetLocations& atomTargets);

	void writeAtomListPdb(std::string dirname,
			              std::string prefix,
			              std::string sufix,
			              int maxNofDigits,
			              int index) const;

    /**
	 * @brief Create a mapping between the local compound atom indices and the global atom indices.
	 */
	void loadIndicesMaps(void);
	
	/**
	 * @brief Get the global atom index from the local compound atom index.
	 * @param cAIx Compound Atom Index
	 * @return Global Atom Index
	 */
	int getGlobalAtomIndex(SimTK::Compound::AtomIndex cAIx);
        
	/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index
	 * maps **/
	void printMaps();

	const std::vector<SimTK::Transform>& getAtomFrameCache() const { return atomFrameCache; }
	std::vector<SimTK::Transform>& updAtomFrameCache() { return atomFrameCache; }

private:
	Span<RoboAtom> subAtomList;
	Span<RoboBondStretch> subBondList;
	Span<RoboBondBend> subAngleList;
	Span<RoboBondTorsion> subTorsionList;

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

	std::size_t rootGlobalAtomIx = 0; // in global atom index
	SimTK::Compound::AtomIndex rootCompoundAtomIx; // in subAtomList index

	bool flipAllChirality = false;
};
