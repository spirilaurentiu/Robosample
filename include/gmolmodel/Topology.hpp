#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

/* -------------------------------------------------------------------------- *
 *			         Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling	                                              *
 */

#include "TrivalentAtomTetra.hpp"
#include "bSpecificAtom.hpp"
#include "bBond.hpp"
#include "server.hpp"

/** Topological information (bonds graph) for one molecule.
It maps to one compound in Molmodel thus it is derived 
from Molmodel Compound class.
It does the following things:
   - loads information from input files such as Amber input prmtop / inpcrd
   - adds parameters to a DuMM force field which belongs to the World class
	 because one DuMM class should be used for multiple molecules
   - contructs the graph based on a list of bSpecificAtom objects each of 
	 which already contains bonding information from the input files
   - defines the rigid bodies based on imput files provided by the users.
Contains a list of atoms bAtomList which consists of bSpecificAtom 
objects **/
class Topology : public SimTK::Compound{
public:

	/** Default Constructor. Sets the name of this molecule to 'no_name '.
	The name has no particular function and is not guaranteed to be unique.**/
	Topology();

	/** Constructor that sets the name of the molecule. The name has no 
	particular function and is not guaranteed to be unique. **/
	explicit Topology(std::string nameOfThisMolecule);

	/** Default Destructor. **/
	virtual ~Topology();

	/** Set atoms properties from a reader: number, name, element, initial
	 * name, force field type, charge, coordinates, mass, LJ parameters **/
	void SetGmolAtomPropertiesFromReader(readAmberInput *amberReader);

	/** Set bonds properties from reader: bond indeces, atom neighbours **/
	void SetGmolBondingPropertiesFromReader(readAmberInput *amberReader);

	/** Set atoms Molmodel types (Compound::SingleAtom derived) based on
	 * their valence **/
	void SetGmolAtomsMolmodelTypes();
	void SetGmolAtomsMolmodelTypesTrial();


	/** Reads data from a specific reader (readAmberInput for now) object **/
	void loadAtomAndBondInfoFromReader(readAmberInput *amberReader);

	/** Print atom list **/
	void PrintAtomList(int whichWorld);


	/** Biotype is a Molmodel hook that is usually used to look up molecular
	 force field specific parameters for an atom type. Gmolmodel defines a
	 new Biotype for each atom. The only thing that is specified is the element
	 with info about name, atomic number, valence and mass. **/
	void bAddBiotypes(
			//std::string resName, 
			readAmberInput *amberReader
			//, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
	setBiotypeChargedAtomType for every atom. These Molmodel functions contain
	information regarding the force field parameters. **/
	void bAddDummAtomClasses(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** Calls DuMM defineBondStretch. **/
	void bAddDummBondParams(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** Calls DuMM defineBondBend. **/
	void bAddDummAngleParams(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** Calls DuMM defineBondTorsion for 1, 2 and 3 periodicities **/
	void bAddDummTorsionParams(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** Adds force field parameters read by the inputReader to DuMM **/
	void bAddDummParams(
		readAmberInput *amberReader
		, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** Print Molmodel specific types as introduced in Gmolmodel **/
	const void PrintMolmodelAndDuMMTypes(SimTK::DuMMForceFieldSubsystem& dumm) const;

	/** Build the molecular tree without the cycle closing bonds **/
	void buildAcyclicGraph(bSpecificAtom *node, bSpecificAtom *previousNode);

	/** After building the acyclic molecular tree close the remaining bonds **/
	void addRingClosingBonds();

	/** Match Default configuration with the coordinates loaded from
	 * the input reader **/
	void matchDefaultConfigurationWithAtomList(
			SimTK::Compound::MatchStratagem matchStratagem);

	/** Builds the Compound's tree, closes the rings, matches the configuration
	on the graph using using Molmodels matchDefaultConfiguration and sets the
	general flexibility of the molecule. **/
	void buildGraphAndMatchCoords(int argRoot);

	// Set flexibility according to flexibility file
	void setFlexibility(std::string argRegimen, std::string flexFN, int whichWorld);

	// Set scale factors for U entries according to flexibility file
	void setUScaleFactorsToBonds(std::string flexFN);

	// Get regimen keyword
	std::string getRegimen();

	// Interface:
	/** Get the name of this molecule **/
	const std::string getName() const {return this->name;}

	/** Set the name of this molecule **/
	void setName(std::string nameOfThisMolecule){
		this->name = nameOfThisMolecule;
	}

	/** Get own CompoundIndex in CompoundSystem **/
	const SimTK::CompoundSystem::CompoundIndex &getCompoundIndex() const;

	/** Set the compoundIndex which is the position in the vector of Compounds
 * of the CompoundSystem **/
	void setCompoundIndex(const SimTK::CompoundSystem::CompoundIndex &compoundIndex);

	/** Compute BAT determinant
	**/
	bool checkIfTripleUnorderedAreEqual(
			std::vector<Compound::AtomIndex> &first,
			std::vector<Compound::AtomIndex> &second);

	void loadTriples(void);
	SimTK::Real calcLogSineSqrGamma2(const SimTK::State &quatState);
	SimTK::Real calcLogDetMBATGamma2Contribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATDistsContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATDistsMassesContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATAnglesContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATMassesContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATInternal(const SimTK::State& someState);
	SimTK::Real calcLogDetMBAT(const SimTK::State&);

	/** Get the number of atoms. **/
	int getNAtoms() const;

	/** Get the number of bonds. **/
	int getNBonds() const;

	/** Get a pointer to an atom object in the atom list inquiring
	by number **/
	bSpecificAtom * getAtomByNumber(int number) const;

	/** Get a pointer to an atom object in the atom list inquiring
	by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
	bSpecificAtom * updAtomByAtomIx(int aIx);

	/** Get a pointer to an atom object in the atom list inquiring
	by atom name **/
	bSpecificAtom * getAtomByName(std::string name) const;

	/** Get the neighbours in the graph **/
	std::vector<bSpecificAtom *> getNeighbours(int) const;

	/** Get the bonded neighbor atom in the parent mobilized body **/
	SimTK::Compound::AtomIndex
	getChemicalParent(
		SimTK::SimbodyMatterSubsystem *matter,
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	/**  **/
	std::vector<SimTK::Transform>
	calcMobodTransforms(
		SimTK::SimbodyMatterSubsystem *matter,
		SimTK::Compound::AtomIndex aIx,
		const SimTK::State& someState,
		SimTK::DuMMForceFieldSubsystem& dumm,
		int whichWorld);

	/** Calculate all atom frames in top frame. It avoids calling 
	calcDefaultAtomFrameInCompoundFrame multiple times. This has 
	to be called every time the coordinates change though. **/
	void calcTopTransforms(void);
	
	/**  **/
	void printTopTransforms(void);

	/**  **/
	SimTK::Transform getTopTransform(SimTK::Compound::AtomIndex);

	/** **/
	const bBond& getBond(int, int) const;

	/** Get bond order. **/
	int getBondOrder(int, int) const;

	// Interface to access the maps
	/** Get MobilizedBody to AtomIndex map **/
	std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex >
	getMbx2aIx(){
		return mbx2aIx;
	}

	// Retunr mbx by calling DuMM functions
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexThroughDumm(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	// Retunr mbx from an olresdy saved map inside Topology
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexFromMap(
		SimTK::Compound::AtomIndex aIx);

	// Get atom location on mobod through DuMM functions
	SimTK::Vec3 getAtomLocationInMobilizedBodyFrameThroughDumm(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	// 
	SimTK::Vec3 calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm,
		SimTK::SimbodyMatterSubsystem& matter,
		const SimTK::State& someState);

	/** Get AtomIndex to MobilizedBodyIndex map **/
	std::map< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex >
	getAIx2mbx(){
		return aIx2mbx;
	}

	/** Get the number of MobilizedBodies in this Compound **/
	std::size_t getNofMobilizedBodies() const {
		return mbx2aIx.size();
	}

	void writeAtomListPdb(std::string dirname,
			              std::string prefix,
			              std::string sufix,
			              int maxNofDigits,
			              int index) const;

	/** To be removed. *Create MobilizedBodyIndex vs Compound::AtomIndex
	 * maps. **/
	void loadMobodsRelatedMaps();

        /** Compound AtomIndex to bAtomList number **/
	void loadCompoundAtomIx2GmolAtomIx(void);
	
	/**  **/
	int getNumber(SimTK::Compound::AtomIndex);
        
	/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index
	 * maps **/
	void printMaps();

	/** Get coordinates **/
	void getCoordinates(
			std::vector<SimTK::Real> Xs,
			std::vector<SimTK::Real> Ys,
			std::vector<SimTK::Real> Zs);

public:

	// TODO move to Context
	int natoms;
	std::vector<bSpecificAtom> bAtomList;

	int nbonds;
	std::vector<bBond> bonds;

	int nTriples;
	std::vector< std::vector<Compound::AtomIndex> > triples;

	// Map mbx2aIx contains only atoms at the origin of mobods
	std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > mbx2aIx;

	// Map aIx is redundant in MobilizedBodyIndeces
	std::map< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex > aIx2mbx;

	// Map aIx to its Transform Default top transform
	std::map< SimTK::Compound::AtomIndex, SimTK::Transform > aIx2TopTransform;

	// Map bSpecificAtom number to aIx
	std::map< SimTK::Compound::AtomIndex, int > CompoundAtomIx2GmolAtomIx;

private:

	int bSpecificAtomRootIndex;

	std::string regimen;
	std::string name;

	/** Every Compound has an index which is the position in the vector
	 * of Compounds in CompoundSystem
	 */
	SimTK::CompoundSystem::CompoundIndex compoundIndex;

	int nofProcesses;
	int baseSetFlag;
	int baseAtomNumber;

	// Gmolmodel to Molmodel (and inverse) bond mappings
	std::map< SimTK::Compound::BondIndex, int > bondIx2GmolBond;
	std::map< int,  SimTK::Compound::BondIndex> GmolBond2bondIx;

};




#endif //TOPOLOGY_H_
