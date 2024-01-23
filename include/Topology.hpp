#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

/* -------------------------------------------------------------------------- *
 *			         Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling	                                              *
 */

#include "ElementCache.hpp"
#include "InternalCoordinates.hpp"
#include "TrivalentAtomTetra.hpp"
#include "bSpecificAtom.hpp"
#include "bBond.hpp"
#include "server.hpp"


/** Helper class for Topology class, used as a key in an AtomClass related
map **/
class AtomClassParams {

  public:
	// Parameters
	int atomicNumber = 0;
	int valence = 0;
	SimTK::Real vdwRadius = 0;
	SimTK::Real LJWellDepth = 0;

	// Constructor
	AtomClassParams(int a, int v, SimTK::Real vdw, SimTK::Real l) :
		atomicNumber(a), valence(v), vdwRadius(vdw), LJWellDepth(l) {}

	// Dump function
	const void dump(void) const { std::cout 
		<< " atomicNumber " << atomicNumber 
		<< " valence " << valence 
		<< " vdwRadius " << vdwRadius 
		<< " LJWellDepth " << LJWellDepth << std::endl;
	}

	// Equal operator
	bool operator==(const AtomClassParams& other) const {

		if ( 	( atomicNumber == other.atomicNumber ) &&
			( valence == other.valence ) &&
			( std::abs(vdwRadius - other.vdwRadius) < 0.0000001 ) &&
			( std::abs(LJWellDepth - other.LJWellDepth) < 0.0000001) ){
			
			return true;
		}else{
			return false;
		}
	}

	// Sort operator
	bool operator<(const AtomClassParams& other) const {
		if ( atomicNumber < other.atomicNumber ){
			return true;
		}else if ( atomicNumber == other.atomicNumber ){
		
			if( valence < other.valence){
				return true;
			}else if( valence == other.valence ){


				if( vdwRadius < other.vdwRadius ){
					return true;
				}else if( vdwRadius == other.vdwRadius ){


					if( LJWellDepth < other.LJWellDepth ){
						return true;
					}else if( LJWellDepth == other.LJWellDepth ){
						// Should throw an error ?
						return true;
					}else{
						return false; // LJWellDepth
					}			

				}else{
					return false; // vdwRadius
				}			

			}else{
				return false; // valence
			}			

		}else{
			return false; // atomicNumber
		}
	}

};

/** Helper class for Topology class, used as a value in an AtomClass related
map **/
class AtomClassId {

  public:
	SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
	std::string name = "noName";

	AtomClassId(int i, std::string s) :
		dummAtomClassIndex(i), name(s) {}
};

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

	void setAtomList(
		std::vector<bSpecificAtom>& argAtomList)
	{		
		(this->bAtomList) = (argAtomList);
		natoms = (this->bAtomList).size();
	}

	/** Default Constructor. Sets the name of this molecule to 'no_name '.
	The name has no particular function and is not guaranteed to be unique.**/
	Topology();

	/** Constructor that sets the name of the molecule. The name has no 
	particular function and is not guaranteed to be unique. **/
	explicit Topology(std::string nameOfThisMolecule);

	/** Default Destructor. **/
	virtual ~Topology();

	/** Set atoms properties from a reader: number, name, element, initial
	 * name, force field type, charge, coordinates, mass, LJ parameters.
	 * This does not set anything in Compund or DuMM. **/
	void SetGmolAtomPropertiesFromReader(readAmberInput *amberReader);

	/** Set bonds properties from reader: bond indeces, atom neighbours **/
	void SetGmolBondingPropertiesFromReader(readAmberInput *amberReader);

	/** Set atoms Molmodel types (Compound::SingleAtom derived) based on
	 * their valence **/
	void SetGmolAtomsCompoundTypes();


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
			//readAmberInput *amberReader
			//, SimTK::DuMMForceFieldSubsystem& dumm
	);

	/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
	setBiotypeChargedAtomType for every atom. These Molmodel functions contain
	information regarding the force field parameters. **/
	void generateDummAtomClasses(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
			, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	);

	/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
	setBiotypeChargedAtomType for every atom. These Molmodel functions contain
	information regarding the force field parameters. **/
	void transferDummAtomClasses(
		std::string resName
		, readAmberInput *amberReader
		, SimTK::DuMMForceFieldSubsystem& dumm
		, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	);

	/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
	setBiotypeChargedAtomType for every atom. These Molmodel functions contain
	information regarding the force field parameters. **/
	void transferDummChargedAtomClasses(
		std::string resName
		, readAmberInput *amberReader
		, SimTK::DuMMForceFieldSubsystem& dumm
		, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	);

	/** Calls DuMM defineBondStretch. **/
	void bAddDummBondParams(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
			, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
	);

	/** Calls DuMM defineBondBend. **/
	void bAddDummAngleParams(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
			, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
	);

	/** Calls DuMM defineBondTorsion for 1, 2 and 3 periodicities **/
	void bAddDummTorsionParams(
			std::string resName
			, readAmberInput *amberReader
			, SimTK::DuMMForceFieldSubsystem& dumm
			, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
			, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs
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
	void generateDummParams(
		readAmberInput *amberReader
		, SimTK::DuMMForceFieldSubsystem& dumm
		, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs

	);

	/** Transfer already generated force field parameters to DuMM **/
	void transferDummParams(
		readAmberInput *amberReader
		, SimTK::DuMMForceFieldSubsystem& dumm
		, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
		, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs

	);

	/** Print Molmodel specific types as introduced in Gmolmodel **/
	const void PrintMolmodelAndDuMMTypes(SimTK::DuMMForceFieldSubsystem& dumm) const;

	/** Build the molecular tree without the cycle closing bonds **/
	void buildAcyclicGraph(bSpecificAtom *node, bSpecificAtom *previousNode);

	/** */
	void buildAcyclicGraphWrap(bSpecificAtom* root);

	/** After building the acyclic molecular tree close the remaining bonds **/
	void addRingClosingBonds();

	/** Match Default configuration with the coordinates loaded from
	 * the input reader **/
	void matchDefaultConfigurationWithAtomList(
			SimTK::Compound::MatchStratagem matchStratagem);

	/**
	 * Check if the provided atom is a possible root and if not find one
	*/
	bSpecificAtom* findARoot(int argRoot);

	std::vector<bSpecificAtom>::iterator findARoot(
		std::vector<bSpecificAtom>::iterator bAtomListBeg,
		std::vector<bSpecificAtom>::iterator bAtomListEnd);

	void findARoot(std::vector<bSpecificAtom>& bAtomListArg);

	/**
	 * Generate an AtomIndex to Top Transforms map
	*/
	void generateAIx2TopXMaps( void );
	void generateAIx2TopXMaps_SP_NEW( void );

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

	// Helper function for calcLogDetMBATAnglesContribution
	// Finds all triple runs - TODO VERY INEFFICIENT
	void loadTriples(void);
	void loadTriples_SP_NEW(void);
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

	/** Get the bonded upstream neighbor atom **/
	SimTK::Compound::AtomIndex
	getNeighbourWithSmallerAIx(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	/** Get the bonded neighbor atom in the parent mobilized body **/
	SimTK::Compound::AtomIndex
	getChemicalParent_IfIAmRoot(
		SimTK::SimbodyMatterSubsystem *matter,
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	/**  **/
	//std::vector<SimTK::Transform>
	//calcMobodToMobodTransforms(
	//	SimTK::SimbodyMatterSubsystem *matter,
	//	SimTK::Compound::AtomIndex aIx,
	//	const SimTK::State& someState,
	//	SimTK::DuMMForceFieldSubsystem& dumm,
	//	int whichWorld);

	/** Calculate all atom frames in top frame. It avoids calling 
	calcDefaultAtomFrameInCompoundFrame multiple times. This has 
	to be called every time the coordinates change though. **/
	void calcTopTransforms(void);
	
	/**  **/
	void printTopTransforms(void);

	/**  **/
	SimTK::Transform getTopTransform(SimTK::Compound::AtomIndex);

	/** **/
	bool checkBond(int, int);

	/** **/
	const bBond& getBond(int, int) const;

	/** Get bond order. **/
	int getBondOrder(int, int) const;

	// Interface to access the maps
	/** Get MobilizedBody to AtomIndex map **/
	//std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex >
	//getMbx2aIx(){
	//	return mbx2aIx;
	//}

	// Alternatives
	SimTK::Transform 
		calcDefaultAtomFrameInCompoundFrameThroughDuMM(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm,
		SimTK::SimbodyMatterSubsystem& matter,
		const SimTK::State& someState);

	// Retunr mbx by calling DuMM functions
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexThroughDumm(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	// Retunr mbx from an olresdy saved map inside Topology
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexFromMap(
		SimTK::Compound::AtomIndex aIx,
		int whichWorld);

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
	std::map< SimTK::Compound::AtomIndex, std::vector<SimTK::MobilizedBodyIndex> >
	getAIx2mbx(){
		return aIx2mbx;
	}

	void writeAtomListPdb(std::string dirname,
			              std::string prefix,
			              std::string sufix,
			              int maxNofDigits,
			              int index) const;

	/** To be removed. *Create MobilizedBodyIndex vs Compound::AtomIndex
	 * maps. **/
	void loadAIx2MbxMap();
	void loadAIx2MbxMap_SP_NEW();
	//void loadMbx2AIxMap();

        /** Compound AtomIndex to bAtomList number **/
	void loadCompoundAtomIx2GmolAtomIx(void);
	void loadCompoundAtomIx2GmolAtomIx_SP_NEW(void);
	
	/**  **/
	int getNumber(SimTK::Compound::AtomIndex);
        
	/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index
	 * maps **/
	void printMaps();

	/** Get coordinates **/
	void getCoordinates(
			std::vector<SimTK::Real>& Xs,
			std::vector<SimTK::Real>& Ys,
			std::vector<SimTK::Real>& Zs);

	void setSubAtomList(
		std::vector<bSpecificAtom>::iterator beginArg,
		std::vector<bSpecificAtom>::iterator endArg,
		ELEMENT_CACHE& elementCacheArg);

	void setAtomList(void);		

	void setSubBondList(
		std::vector<bBond>::iterator beginArg,
		std::vector<bBond>::iterator endArg);

	void setBondList(void);
	
public:

	void BAT();

	// Atoms
	int natoms;
	std::vector<bSpecificAtom> bAtomList;
	std::vector<bSpecificAtom>::iterator atomsBeg_It;
	std::vector<bSpecificAtom>::iterator atomsEnd_It;
	size_t atomsBeg_Ix;
	size_t atomsEnd_Ix;
	array_view<std::vector<bSpecificAtom>::iterator> subAtomList;

	// Bonds
	int nbonds;
	std::vector<bBond> bonds;
	std::vector<bBond>::iterator bondsBeg_It;
	std::vector<bBond>::iterator bondsEnd_It;
	size_t bondsBeg_Ix;
	size_t bondsEnd_Ix;
	array_view<std::vector<bBond>::iterator> subBondList;


	// Triples
	int nTriples;
	std::vector< std::vector<Compound::AtomIndex> > triples;

	// Map mbx2aIx contains only atoms at the origin of mobods
	//std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > mbx2aIx;

	// Map aIx is redundant in MobilizedBodyIndeces // TODO remove 
	std::map< SimTK::Compound::AtomIndex, std::vector<SimTK::MobilizedBodyIndex> > aIx2mbx;

	// Map aIx to its Transform Default top transform
	std::map< SimTK::Compound::AtomIndex, SimTK::Transform > aIx2TopTransform;

	// Map bSpecificAtom number to aIx
	std::map< SimTK::Compound::AtomIndex, int > CompoundAtomIx2GmolAtomIx;

	// Gmolmodel to Molmodel (and inverse) bond mappings
	std::map< SimTK::Compound::BondIndex, int > bondIx2GmolBond;
	std::map< int,  SimTK::Compound::BondIndex> GmolBond2bondIx;


	int nofProcesses;
	int baseSetFlag;
	int baseAtomNumber;

	int bSpecificAtomRootIndex;

	// Atom frames in Top frame
	std::vector<SimTK::Transform> atomFrameCache;


private:


	std::string regimen;
	std::string name;

	/** Every Compound has an index which is the position in the vector
	 * of Compounds in CompoundSystem
	 */
	SimTK::CompoundSystem::CompoundIndex compoundIndex;



	//std::map<AtomClassParams, AtomClassId> aClassParams2aClassId;

	ELEMENT_CACHE elementCache;
};




#endif //TOPOLOGY_H_
