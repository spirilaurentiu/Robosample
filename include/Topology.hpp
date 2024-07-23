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

	/** Default Constructor. Sets the name of this molecule to 'no_name '.
	The name has no particular function and is not guaranteed to be unique.**/
	Topology();

	/** Constructor that sets the name of the molecule. The name has no 
	particular function and is not guaranteed to be unique. **/
	explicit Topology(std::string nameOfThisMolecule);

	/** Default Destructor. **/
	virtual ~Topology();

	/** Print atom list **/
	void PrintAtomList(int whichWorld);

	/** Print Molmodel specific types as introduced in Gmolmodel **/
	void PrintMolmodelAndDuMMTypes(SimTK::DuMMForceFieldSubsystem& dumm) const;

	/**
	 * Generate an AtomIndex to Top Transforms map
	*/
	void generateAIx2TopXMaps( void );

	/**	
	* @brief Get the name of this molecule
	* @param 
	* @return name of the molecule
	*/
	const std::string getName() const {return this->name;}

	/** Set the name of this molecule **/
	void setName(std::string nameOfThisMolecule){
		this->name = nameOfThisMolecule;
	}

	/**	
	* @brief Get own CompoundIndex in CompoundSystem
	* @param 
	* @return CompoundIndex
	*/
	/**  **/
	const SimTK::CompoundSystem::CompoundIndex &getCompoundIndex() const;

	/**	
	* @brief Set the compoundIndex which is the position in the vector of
	* Compounds of the CompoundSystem
	* @param compoundIndex 
	* @return
	*/
	/**  **/
	void setCompoundIndex(const SimTK::CompoundSystem::CompoundIndex &compoundIndex);

	/** Compute BAT determinant
	**/
	bool checkIfTripleUnorderedAreEqual(
			std::vector<Compound::AtomIndex> &first,
			std::vector<Compound::AtomIndex> &second);

	// Helper function for calcLogDetMBATAnglesContribution
	// Finds all triple runs - TODO VERY INEFFICIENT
	void loadTriples_SP_NEW(void);

	SimTK::Real calcLogSineSqrGamma2(const SimTK::State &quatState);
	SimTK::Real calcLogDetMBATGamma2Contribution(const SimTK::State&);

	SimTK::Real calcLogDetMBATDistsContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATAnglesContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATMassesContribution(const SimTK::State&);
	SimTK::Real calcLogDetMBATInternal(const SimTK::State& someState);

	/** Get the number of atoms. **/
	int getNAtoms() const;

	/** Get the number of bonds. **/
	int getNBonds() const;

	/** Get a pointer to an atom object in the atom list inquiring
	by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
	bSpecificAtom * updAtomByAtomIx(int aIx);

	/** Get a pointer to an atom object in the atom list inquiring
	by atom name **/
	bSpecificAtom * getAtomByName(std::string name) const;

	/** Get the neighbours in the graph **/
	std::vector<bSpecificAtom *> getNeighbours(int) const;

	/**	
	* @brief Get the bonded neighbor atom in the parent mobilized body.
	* @param aIx Compound Atom Index
	* @return Compound atom index of the root
	*/
	/**  **/
	SimTK::Compound::AtomIndex
	getChemicalParent_IfIAmRoot(
		SimTK::SimbodyMatterSubsystem *matter,
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

	/**	
	* @brief Calculate all atom frames in top frame. It avoids calling 
	* calcDefaultAtomFrameInCompoundFrame multiple times. This has to be called
	* every time the coordinates change though.
	* @param : 
	* @return
	*/
	void calcAtomsTopTransforms(void);
	
	/**	
	* @brief 
	* @return
	*/
	void printTopTransforms(void);

	/**	
	* @brief Get atom Top level transform from the existing Topology map
	* @param cAIx: atom Compound AtomIndex
	* @return Atom's Top level transform
	*/
	SimTK::Transform getTopTransform_FromMap(SimTK::Compound::AtomIndex cAIx);

	/**	
	* @brief 
	* @return
	*/
	bool checkBond(int, int);

	/**	
	* @brief 
	* @return
	*/
	const bBond& getBond(int, int) const;

	// Interface to access the maps

	// Retunr mbx by calling DuMM functions
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexThroughDumm(
		SimTK::Compound::AtomIndex aIx,
		SimTK::DuMMForceFieldSubsystem& dumm);

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

    /** Compound AtomIndex to bAtomList number **/
	void loadCompoundAtomIx2GmolAtomIx(void);
	
	/**  **/
	int getNumber(SimTK::Compound::AtomIndex cAIx);
        
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

	//void setAtomList(void);		

	void setSubBondList(
		std::vector<bBond>::iterator beginArg,
		std::vector<bBond>::iterator endArg);

	//void setBondList(void);

    void setBondMappings(std::unordered_map<int, SimTK::Compound::BondIndex>& argBondMapping) {
		bondMapping = &argBondMapping;
	}


public:

	//void BAT();

	// Atoms
	int natoms;
	//std::vector<bSpecificAtom> bAtomList;
	std::vector<bSpecificAtom>::iterator atomsBeg_It;
	std::vector<bSpecificAtom>::iterator atomsEnd_It;
	size_t atomsBeg_Ix;
	size_t atomsEnd_Ix;
	array_view<std::vector<bSpecificAtom>::iterator> subAtomList;

	// Bonds
	int nbonds;
	//std::vector<bBond> bonds;
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
	std::unordered_map<int, SimTK::Compound::BondIndex>* bondMapping;


	int nofProcesses;
	int baseSetFlag;
	int baseAtomNumber;

	int bSpecificAtomRootIndex;

	// Atom frames in Top frame
	std::vector<SimTK::Transform> atomFrameCache;


private:

	std::string name;

	/** Every Compound has an index which is the position in the vector
	 * of Compounds in CompoundSystem
	 */
	SimTK::CompoundSystem::CompoundIndex compoundIndex;



	//std::map<AtomClassParams, AtomClassId> aClassParams2aClassId;

	ELEMENT_CACHE elementCache;

	std::size_t rootAtomIx = 0;
};




#endif //TOPOLOGY_H_
