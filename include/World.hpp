#ifndef WORLD_H_
#define WORLD_H_

/* -------------------------------------------------------------------------- *
 *		                       Robosampling                           *
 * -------------------------------------------------------------------------- *
 * This is part of Robosample		                                      *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include <thread>
#include <array>
//#include <random>
#include <math.h>

//#include <Eigen/Dense>
//#include <unsupported/Eigen/MatrixFunctions>
//#include <Eigen/Eigenvalues>
//#include <Eigen/LU>
//using Eigen::MatrixXd;

#include "Simbody.h"
#include "Molmodel.h"

#include "readAmberInput.hpp"

#include "ParaMolecularDecorator.hpp"
#include "FixmanTorque.hpp"

#ifndef BaseSampler
#define BaseSampler HMCSampler
#endif

#include "server.hpp"
#include "Topology.hpp"
//#include "LAHMCSampler.hpp"
#include "HMCSampler.hpp"
#include "ConformationalSearch.hpp"

// TODO write a pdb writer for all the Compounds
// TODO move this in Topology since they work only for one Compound
// The following use PdbStructure for a Compound's Default Configuration
void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix);

void writePdb(	  SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix);

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

void writePdb(	  SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

void writePdb(SimTK::PdbStructure pdb, const char *FN);

//==============================================================================
//                   CLASS TaskSpace
//==============================================================================
/**
 *  Contains a Symbody task space and additional data
 **/
class StationTaskLaurentiu{
    friend class Context;

public:

    StationTaskLaurentiu();

private:

};



enum class ROOT_MOBILITY : int {
	FREE = 0,
	CARTESIAN,
	WELD,
	FREE_LINE,
	BALL,
	PIN
};

struct BOND_FLEXIBILITY {
	BOND_FLEXIBILITY() = default;
	BOND_FLEXIBILITY(int i, int j, BondMobility::Mobility mobility)
		: i(i), j(j), mobility(mobility) {}
		
	int i = -1;
	int j = -1;
	BondMobility::Mobility mobility = BondMobility::Default;
};

//==============================================================================
//                   CLASS World
//==============================================================================
/**
 *  Contains a Symbody system and additional data that define a regimen
 **/

class Context;

class World {
public:
	// --- Structural functions ---
	/** Constructor **/
	explicit World(	int worldIndex,
					int requestedNofMols,
					bool isVisual=true,
					SimTK::Real visualizerFrequency = 0.0015);

	void setFlexibilities(const std::vector<BOND_FLEXIBILITY>& flexibilities);
	const std::vector<BOND_FLEXIBILITY>& getFlexibilities() const;

	void generateDummParams(const std::vector<bSpecificAtom>& atoms,
		const std::vector<bBond>& bonds,
		const std::vector<DUMM_ANGLE>& dummAngles,
		const std::vector<DUMM_TORSION>& dummTorsions,
		const ELEMENT_CACHE& elementCache);

	/** Creates a topology object and based on amberReader forcefield
	 parameters - defines Biotypes; - adds BAT parameters to DuMM **/
	void AddMolecule(
		readAmberInput *amberReader,
		//std::string rbFN,
		//std::string flexFN,
		//std::string regimenSpec,
		std::string argRoot
		//, std::string argRootMobility
		);

	//
	void AddBiotypes(int which, readAmberInput *amberReader);
	
	void BuildTopologyGraph(int which, std::string argRoot);
	void AllocateCoordBuffers(int natoms);

	/** Adopts a topology **/
	void adoptTopology(int which);

	/** Calls CompoundSystem.modelCompounds and realizes Topology
	To be called after loading all Compounds. **/
	void modelTopologies(std::string GroundToCompoundMobilizerType);

	SimTK::Real getRecommendedTimesteps(void);

	//=========================================================================
	//                   CONSTRAINTS
	//=========================================================================

	/** Add contact constraints to specific bodies **/
	void addRodConstraint(State& someState);

	/** Add contact constraints to specific bodies **/
	const SimTK::State& addSpeedConstraint(int prmtopIndex);

	//=========================================================================
	//                   TaskSpace Functions
	//=========================================================================

	/** Allocate memory for a task space consisting of a set of body indeces
 	* station on the bodies expresed in both guest (target) and host 
 	* and the difference between them
	*/
	void addTaskSpaceLS(void);

	/** Update target task space */
	void updateTaskSpace(const State& someState);

	/** Get delta stationP */
	SimTK::Array_<SimTK::Vec3>& 
	getTaskSpaceStationPInGuest(void);

	/** Get delta stationP */
	SimTK::Array_<SimTK::Vec3>& 
	getTaskSpaceStationPInHost(void);

	/** Get delta stationP */
	SimTK::Array_<SimTK::Vec3>& 
	getTaskSpaceDeltaStationP(void);

	/** Calc station Jacobian */
	void calcStationJacobian(const State& someState,
        SimTK::Matrix_<SimTK::Vec3>& JS) const;


	/** Add contact surfaces to bodies **/
	const SimTK::State& addContacts(int prmtopIndex);

	//=========================================================================
	//                   Membrane-Related Functions
	//=========================================================================
	
	//-------------------------------------------------------------------------
	/** @name Contacts - Contacts  **/
	/**@{**/

	/**	Add contact surfaces to bodies 
		By (my own) convention, the atoms that are in Z>0 are in constant 
		contact	with the half-space with clique1, so set those in Clique 0 to
		avoid large energies/movements.
		If a prmtopIndex of (-1) is encountered, that means that 
		particlular topology needs to be skipped.**/
	void addMembrane(const SimTK::Real halfThickness);
	/**@}**/
	//-------------------------------------------------------------------------

	//=========================================================================
	//                   CONTACTS Functions
	//=========================================================================
	/** @name Membrane - Membrane  **/
	/**@{**/

	/** Add a membrane represented by a contact surface 
	We can approximate a membrane-like environment (mechanically speaking,
 	i.e. no electrostatic interactions (!yet!) ) via an elastic environment.
	In practice, we realize this by way of 4 overlapping half-spaces. 
	2 of them occupy the Z>0 space and are translated on the Z axis by the
	value of halfThickness. The other 2 occupy the Z<0 spaceand are translated
	by the same value.
	The way this works is we select 4 subsets of atoms: 
		*) 1 that is "below" the membrane, and can't go up on the Z axis:
			these can be charged atoms that can't cross the hydrophobic core of the membrane
		*) 1 that is "above" the membrane, and can't go higher on the Z axis:
			these can be hydrophobic patches that can't "escape" the hydrophobic patch
			Think of a hydrophobic helix that's lodged in the membrane
		*) The other 2 are the same, but in reverse ("above" the membrane, but can't cross down
			and "below" but can't go lower") 
	Each of these subsets are only affected by one half-space. The contact cliques are:
		0) Z>0, translated by +halfThickness on the Z axis
		1) Z>0, translated by -halfThickness on the Z axis
		2) Z<0, translated by +halfThickness on the Z axis
		3) Z<0, translated by -halfThickness on the Z axis **/
	void addContacts(const std::vector<int>& prmtopIndex, const int topologyIx, 
		const SimTK::ContactCliqueId cliqueId);	/**@}**/
	//-------------------------------------------------------------------------



	/** Realize Topology for this World **/
	const SimTK::State& realizeTopology();

	/** Assign a scale factor for generalized velocities to every mobilized
	body **/
	void setUScaleFactorsToMobods(void);

	/** Load CompoundAtomIndex to Gmolmodel atom index map **/
	void loadCompoundRelatedMaps();

	/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
	void loadMbx2AIxMap();

	/** Get the number of molecules **/
	int getNofMolecules() const;

	// These are no longer needed TODO: delete
	/** Get MobilizedBody to AtomIndex map **/
	std::map< SimTK::MobilizedBodyIndex, std::pair<int, SimTK::Compound::AtomIndex>>&
	getMbx2aIx();

	/** Get the number of MobilizedBodies in this Compound **/
	std::size_t getNofMobilizedBodies() const ;

	/** Get U scale factor for the mobilized body **/
	SimTK::Real getMobodUScaleFactor(SimTK::MobilizedBodyIndex& ) const;
	//...............

	/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index
	 * maps **/
	void printMaps();

	//...............


	// --- Inter-world functions: Pass configurations among Worlds
	// ELIZA
	SimTK::Vec3 calcAtomLocationInGroundFrameThroughOMM(const SimTK::DuMM::AtomIndex&);

	// Get geometric center of a subset of atoms
	// TEODOR
	SimTK::Vec3 getGeometricCenterOfSelection(
		const SimTK::State & state);



	float setSphereRadius (float argRadius);

	/**
	 * @brief Get the current Compound Cartesian coords.
	 * @param state state.
	 * @details Return a 2D vector representing all the coordinates of this World.
 	* The first dimension represents the molecules (topologies) and the second
 	* dimension (inner) represents the coordinates. The second inner dimension
 	* type is pair of bSpecificAtom* and a Vec3. Thus, besides coordinates, it
 	* contains all the information in bSpecificAtom as well. The bottleneck here
 	* is the calcAtomLocationInGroundFrame from Compound.
	*/
	std::vector<std::vector<
	std::pair<bSpecificAtom *, SimTK::Vec3> > >
		getAtomsLocationsInGround(SimTK::State&);

	/** Get the current Compound Cartesian coordinates using Simbody **/
	std::vector<std::vector<
	std::pair<bSpecificAtom *, SimTK::Vec3> > >
		getCurrentAtomsLocationsInGround(void);

	/** Nice print helper for get/setAtomsLocations */
	void PrintAtomsLocations(const std::vector<std::vector<
		std::pair<bSpecificAtom *, SimTK::Vec3> > >& someAtomsLocations);
	void WriteRst7FromTopology(std::string FN);


	//=========================================================================
	//                   RECONSTRUCTION-Related Functions
	//=========================================================================
	
	//-------------------------------------------------------------------------
	/** @name Contacts - Contacts  **/
	/**@{**/

	/**	Description **/

	void PrintFullTransformationGeometry(std::string indS, const SimTK::State&,
		bool x_pf_r = true, bool x_fm_r = true, bool x_bm_r = true,
		bool x_pf_p = true, bool x_fm_p = true, bool x_bm_p = true);
	
	/**
	 * RMSD function
	*/
	SimTK::Real RMSD(
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations	
	) const ;

	/**
	 * Maximum distance between two corresponding atoms
	*/
	std::pair<int, SimTK::Real> maxAtomDeviation(
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations	
	) const ;

	/**
	 * Helper for setAtoms Locations This function is only intended for root atoms!!
	*/
	std::vector<SimTK::Transform>calcMobodToMobodTransforms(
		Topology& topology,
		SimTK::Compound::AtomIndex aIx,
		const SimTK::State& someState);

	/**
	 * Calc X_FM transforms for reconstruction
	*/
	SimTK::Transform calcX_FMTransforms(
		Topology& topology,
		SimTK::Compound::AtomIndex aIx,
		const SimTK::State& someState);


	// REFAC ----------------------------------------------------------------------

	/**	
	* @brief Takes coordinates from molecule topoIx and puts them into atomTargets
	* @param otherWorldsAtomsLocations: Pairs of (atom, and its position) within
	* 		 a vector of Topologies
	* @param atomTargets: map of atoms' CompoundAtomIndex to their positions
	* @return
	*/
	void
	extractAtomTargets(
		int topoIx,
		const std::vector<std::vector<
		std::pair<bSpecificAtom *, SimTK::Vec3> > >& otherWorldsAtomsLocations,
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets);

	/*!
	* <!-- Compound matchDefaultConfiguration for molecule topoIx -->
	*/
	SimTK::Transform
	setAtoms_Compound_Match(
		int topoIx,
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets);

	/*!
	* <!-- Set atoms' frames in mobods. Also get locations in mobods for 
	* further use -->
	*/
	void
	setAtoms_Compound_FramesAndLocsInMobods(
		int topoIx,
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets,
		SimTK::Vec3* locs);

	void
	setAtoms_SetDuMMStations(
		int topoIx,
		SimTK::Vec3* locs
	);

	SimTK::State&
	setAtoms_XPF_XBM(
		SimTK::State& someState,
		int topoIx
	);

	SimTK::State&
	setAtoms_MassProperties(
		SimTK::State& someState,
		int topoIx
	);

	SimTK::State&
	setAtoms_XFM(
		SimTK::State& someState,
		int topoIx
	);

	/** Set Compound, MultibodySystem and DuMM configurations according to
	some other World's atoms **/
	SimTK::State&
	setAtomsLocationsInGround_REFAC(SimTK::State&,
		const std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		otherWorldsAtomsLocations);

	// REFAC ----------------------------------------------------------------------

	/**@}**/


	/** Return own CompoundSystem **/
	CompoundSystem *getCompoundSystem() const;

	SimTK::GeneralForceSubsystem* getGeneralForceSubsystem() const {
		return forces.get();
	}
	SimTK::SimbodyMatterSubsystem* getSimbodyMatterSubsystem() const {
		return matter.get();
	}

	/** Set own Compound system **/
	// TODO find a solution for the old one
	void setCompoundSystem(CompoundSystem *compoundSystem);

	/** Update Gmolmodel bSpecificAtom Cartesian coordinates according to
	Molmodel Compound which in turn relizes Position and uses matter
	 to calculate locations. **/
	void updateAtomListsFromCompound(const SimTK::State &state);

	/** Access to molecule (Topology) objects
	Get a readble reference of one of the molecules **/
	const Topology& getTopology(std::size_t moleculeNumber) const;

	/** Get a writeble reference of one of the molecules **/
	Topology& updTopology(std::size_t moleculeNumber);
	//...............

	//.......................
	// --- Thermodynamics ---
	//.......................
	/** Get the World (macro) temperature **/
	SimTK::Real getTemperature();

	/** Set the World (macro) temperature **/
	void setTemperature(SimTK::Real);
	//...............

	/** Set the World (macro) temperature **/
	void setBoostTemperature(SimTK::Real);
	//...............

	// --- Simulation ---
	/** Get/Set seed for reproducibility. **/
	void setSeed(uint32_t argSeed);

	/** Use the Fixman torque as an additional force subsystem.
	Careful not have different temperatures for World and Fixman Torque. **/
	void addFixmanTorque();

	/** Amber like scale factors. **/
	void setAmberForceFieldScaleFactors();

	/** Set a global scaling factor for all the terms the forcefield **/
	void setGlobalForceFieldScaleFactor(SimTK::Real);

	/** Set GBSA implicit solvent scale factor **/
	void setGbsaGlobalScaleFactor(SimTK::Real);

	SimTK::DuMMForceFieldSubsystem& getForceField();

	/** Get a writeble pointer to the DuMM force field **/
	SimTK::DuMMForceFieldSubsystem * updForceField();

	/** Return true if the Fixman torque flag is set **/
	bool isUsingFixmanTorque() const;
	//...............

	// DOESN'T WORK WITH OPENMM
	SimTK::Real CalcFullPotentialEnergyIncludingRigidBodies(void);
	SimTK::Real CalcPotentialEnergy(void);

	// Calculate Fixman potential
	SimTK::Real calcFixman();

	/**
	 *  Generate a proposal
	 * */
	bool generateProposal(void);

	/** Generate a number of samples **/
	bool generateSamples(int howMany, std::stringstream& worldOutStream);
	//...............

	//...................
	// --- Statistics ---
	//...................
	/** How many samples did we have so far **/
	std::size_t getNofSamples() const;

	/** Sampler manipulation functions **/
	std::size_t getNofSamplers() const;

	/** Add a sampler to the World **/
	bool addSampler(SamplerName samplerName,
		SampleGenerator generator,
		IntegratorName integratorName,
		ThermostatName thermostatName,
		SimTK::Real timestep,
		int mdStepsPerSample,
		int mdStepsPerSampleStd,
		SimTK::Real boostTemperature,
		int boostMDSteps,
		int distort,
		int work,
		int flow,
		bool useFixmanPotential);

	void useOpenMM(bool ommvv, SimTK::Real boostTemp, SimTK::Real timestep);

	// TODO Use Sampler polymorphism
	/** Get a sampler based on its position in the samplers vector **/
	BaseSampler * getSampler(std::size_t which) const;

	/** Get a writable sampler based on its position in the samplers vector **/
	BaseSampler * updSampler(std::size_t which);

	/** Get writble pointer to FixmanTorque implementation **/
	FixmanTorque * updFixmanTorque();

	/** Get pointer to FixmanTorque implementation **/
	FixmanTorque * getFixmanTorque() const;

	/** Allocate space for containers that keep statistics if we're doing any **/
	void allocateStatsContainers(void);

	// Get the (potential) energy transfer
	// If any of the Q, U or tau is actively modifyied by the sampler
	// the Jacobian of that transformation will be included too
	SimTK::Real getWorkOrHeat(void);

	// Get the (potential) energy transfer in the form of work
	// If any of the Q, U or tau is actively modifyied by the sampler
	// the Jacobian of that transformation will be included too
	SimTK::Real getWork(void);

	// Set initial values of X_PF or X_BM
	void setTransformsMeansToIni(void);

	// Set initial values of X_PF or X_BM
	void setTransformsMeansToCurrent(SimTK::State& someState);

	// Set initial values of X_PF or X_BM
	void setTransformsMeans(const std::vector<SimTK::Real>& givenX_PF,
		const std::vector<SimTK::Real>& givenX_BM);

	// Set X_PF and X_BM related values
	void setTransformsMeansToMin(readAmberInput &amberReader);

	// Get X_PF and X_BM related values
	void getTransformsStatistics(SimTK::State& someState);

	/** Update transforms means given a previous value */
	void updateTransformsMeans(SimTK::State& someState);

	/** Get X_PF means */
	std::vector<SimTK::Real>& getX_PFMeans(void);

	/** Get X_BM means */
	std::vector<SimTK::Real>& getX_BMMeans(void);

	/**
	 * Calculate bond length and angle deviations from their means
	*/ 
	void calcBendStretchDeviations(
		SimTK::State& someState,
		std::vector<SimTK::Real>& X_PFdiffs,
		std::vector<SimTK::Real>& X_BMdiffs
	);
	
	// Print bond lengths and angle bends
	void traceBendStretch(SimTK::State& someState);

	// Print X_PF
	void PrintAcosX_PFs(void);

	// Print X_PF
	void PrintNormX_BMs(void);

	// Print X_PF means
	void PrintAcosX_PFMeans(void);

	// Print X_PF means
	void PrintNormX_BMMeans(void);

	//...............

	// REORIENT

	SimTK::Transform& getReorientTransformInAnotherBody(
		const State &someState,
		const MobilizedBody &inBodyA,
		const MobilizedBody &ofBodyB,
		const SimTK::Transform &reorientAB,
		SimTK::Transform& X_FMprim);

	//...............

	// -- Debugging / helper functions ---
	/** Print information about Simbody systems **/
	void PrintSimbodyStateCache(SimTK::State&);

	/** Print a Compound Cartesian coordinates as given by
	 * Compound::calcAtomLocationInGroundFrame **/
	void printPoss(const SimTK::Compound& c, SimTK::State& someState);

	/** Print a Compound Cartesian velocities as given by
 * Compound::calcAtomVelocityInGroundFrame **/
	void printVels(const SimTK::Compound& c, SimTK::State& someState);

	/** Print a Compound Cartesian coordinates and velocities
	 * as given by Compound::calcAtomLocationInGroundFrame and
	 * Compound::calcAtomVelocityInGroundFrame**/
	void printPossVels(const SimTK::Compound& c, SimTK::State& someState);
	//...............


	// RANDOM_WALK related functions; we don't need getter, since we only
	// use these values inside the scope of World.
	void setTopologyIXs(std::vector<int> topologyIXs);
	void setAmberAtomIXs(std::vector<std::vector<int>> AmberAtomIXs);

	// const std::vector<SimTK::Real&>& getCppQs(void){
	// 	return CppQs;
	// }
	// std::vector<SimTK::Real&>& updCppQs(void){
	// 	return CppQs;
	// }

	/** Get references to Qs **/
	void initializeCppQs(void){

		// // 
		// SimTK::State& currentState = integ->updAdvancedState();
		// const SimTK::Vector & allQs = matter->getQ(currentState);
		// int NQ = matter->getNQ(currentState);
		// // 
		// for(int qIx = 0; qIx < NQ; qIx++){
		// 	CppQs[qIx] = allQs[qIx];
		// }

	}


public:

	// --- The three S: Study, System and State related ---
	/** System->MultibodySystem->MolecularMechanicsSystems->CompoundSystem **/
	std::unique_ptr<SimTK::CompoundSystem> compoundSystem;

	/** Subsystem->SimbodyMatterSubsystem **/
	std::unique_ptr<SimTK::SimbodyMatterSubsystem> matter;

	/** Subsystem->ForceSubsystem->GeneralForceSubsystem **/
	std::unique_ptr<SimTK::GeneralForceSubsystem> forces;

	/** Subsystem->ForceSubsystem->DuMMForceFieldSubsystem **/
	std::unique_ptr<SimTK::DuMMForceFieldSubsystem> forceField;

	// std::vector<std::unique_ptr<SimTK::ConformationalController>> controller;
	// std::vector<std::unique_ptr<SimTK::Force::Custom>> controlForce;

	/** Nof molecules **/
	int moleculeCount;

	/** Molecules (topologies<-Compounds) objects **/
	//std::vector<bMoleculeReader *> moleculeReaders;
	std::vector<Topology>* topologies = nullptr;
	std::vector<std::string> roots;
	std::vector<std::string> rootMobilities;

	/** Joint types **/
	//std::map< SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility> mbx2mobility;

	/** Vectors of Cartesian coordinates **/
	std::vector<SimTK::Real> Xs;
	std::vector<SimTK::Real> Ys;
	std::vector<SimTK::Real> Zs;

	/** This vector stores a configuration if is needed for later use **/
	SimTK::Transform *TVector;

	/** Topologies graphs as tables - to be removed **/
	int **mbxTreeMat;    // tree representing the bonding
	SimTK::Real *branchMassVec; // branch masses self body included
	//...............

	// --- Thermodynamics ---
	SimTK::Real temperature;
	//...............

	// --- Simulation ---
	std::unique_ptr<SimTK::VerletIntegrator> integ;
	std::unique_ptr<SimTK::TimeStepper> ts;
	std::vector<std::unique_ptr<BaseSampler>> samplers;

	// Contact related
	std::unique_ptr<ContactTrackerSubsystem> tracker;
	std::unique_ptr<CompliantContactSubsystem> contactForces;
	ContactCliqueId clique1;
	std::unique_ptr<MobilizedBody::Weld> membrane;
	std::unique_ptr<Body::Rigid> memBody;
	//...............

	// --- Statistics ---
	std::vector<SimTK::Real> acosX_PF00;
	std::vector<SimTK::Real> normX_BMp;
	std::vector<SimTK::Real> acosX_PF00_means;
	std::vector<SimTK::Real> normX_BMp_means;

	//std::vector<SimTK::Real> CppQs;

	//...............

	// // --- Graphics ---
	bool visual;

	// Our decorations
	std::unique_ptr<ParaMolecularDecorator> paraMolecularDecorator;

	// Decoration subsystem
	std::unique_ptr<SimTK::DecorationSubsystem> decorations;

	// Visualizer
	std::unique_ptr<SimTK::Visualizer> visualizer;

	// Visualizer reporter
	std::unique_ptr<SimTK::Visualizer::Reporter> visualizerReporter;
	//...............

	// --- Mixing data ---
	int ownWorldIndex;
	//...............

	/** Get writble pointer to Fixman Torque and other forces**/
	// std::unique_ptr<FixmanTorque> FixmanTorqueImpl;
	// std::unique_ptr<SimTK::Force::Custom> FixmanTorqueForce;

	std::unique_ptr<SimTK::Force::Custom> FixmanTorqueForce;
	// std::unique_ptr<FixmanTorque> FixmanTorqueImpl;
	std::unique_ptr<SimTK::Force::Custom> FixmanTorqueExtForce;
	// std::unique_ptr<FixmanTorqueExt> FixmanTorqueExtImpl;

	FixmanTorque* FixmanTorqueImpl = nullptr;
	// SimTK::Force::Custom* FixmanTorqueForce = nullptr;
	FixmanTorqueExt* FixmanTorqueExtImpl = nullptr;
	// SimTK::Force::Custom* FixmanTorqueExtForce = nullptr;

	//Task Space
	SimTK::Array_<SimTK::MobilizedBodyIndex> onBodyB;
	SimTK::Array_<SimTK::Vec3> taskStationPInGuest;
	SimTK::Array_<SimTK::Vec3> taskStationPInHost;
	SimTK::Array_<SimTK::Vec3> taskDeltaStationP;

	// Constraints
	std::vector<std::pair <SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex> > rodBodies;
	SimTK::Array_<SimTK::Vec3> conStationPInGuest;
	SimTK::Array_<SimTK::Vec3> conStationPInHost;
	SimTK::Array_<SimTK::Vec3> conDeltaStationP;
	SimTK::Array_<SimTK::Constraint::Rod> rodConstraints;


	/**
	 * Define some convenient transforms
	*/

	// SimTK::UnitVec3 constXAxis(1, 0, 0);
	// SimTK::UnitVec3 constYAxis(0, 1, 0);
	// SimTK::UnitVec3 constZAxis(0, 0, 1);
	// SimTK::UnitVec3 constOriginVec(0, 0, 0);

	// X axis to Z axis switch
	const SimTK::Transform X_to_Z 
		=  SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);
	const SimTK::Transform Z_to_X = ~X_to_Z;

	// Y axis to Z axis switch
	const SimTK::Transform Y_to_Z =
		SimTK::Transform(SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::XAxis));
	const SimTK::Transform Z_to_Y = ~Y_to_Z;

	// X axis to X axis switch
	const SimTK::Transform Y_to_X =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::ZAxis);
	const SimTK::Transform X_to_Y = ~Y_to_X;

	//
	void setSamplesPerRound(int samples);
	int getSamplesPerRound() const;

	void setDistortOption(int distort);
	int getDistortOption() const;

	void setRootMobility(ROOT_MOBILITY rootMobility);
	const SimTK::String& getRootMobility() const;


	// BAT ====================================================================

    // Getter for myContext
    const Context* getMyContext() const {
        return myContext;
    }

    // Setter for myContext
    void setMyContext(Context* context) {
        myContext = context;
    }

    // Updater for myContext
    Context* updMyContext(void) {
        return myContext;
    }	


	const int getOwnIndex(void) const{
		return ownWorldIndex;
	}

	/*!
	* <!--	 -->
	*/
    const std::pair<int, SimTK::Compound::AtomIndex>&
	getMobodRootAtomIndex(SimTK::MobilizedBodyIndex mbIndex) const
	{

        auto it = mbx2aIx.find(mbIndex);

        if (it != mbx2aIx.end()) {

            return it->second;

        } else { 
            
            return errorTopoAtomPair;

        }
    }

	/*!
	* <!--	 -->
	*/
    std::pair<int, SimTK::Compound::AtomIndex>&
	updMobodRootAtomIndex(SimTK::MobilizedBodyIndex mbIndex)
	{

        auto it = mbx2aIx.find(mbIndex);

        if (it != mbx2aIx.end()) {

            return it->second;

        } else {
			
        	return errorTopoAtomPair;
        }
    }

    // Assign the atom index to the mobilized body index 
    void setAtomIndex(
		SimTK::MobilizedBodyIndex mbIndex,
		int topoIx,
		SimTK::Compound::AtomIndex aIndex) 
	{
        mbx2aIx[mbIndex] = std::pair<int, SimTK::Compound::AtomIndex> {topoIx, aIndex}; 
    }

	// BAT --------------------------------------------------------------------

	//...........................
	// --- Drilling functions ---
	//...........................
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<std::vector<double>>& getEnergies_drl_bon();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<std::vector<double>>& getEnergies_drl_ang();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<std::vector<double>>& getEnergies_drl_tor();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<std::vector<double>>& getEnergies_drl_n14();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<std::vector<double>>& getEnergies_drl_vdw();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<std::vector<double>>& getEnergies_drl_cou();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<OpenMM::Vec3>& getForces_drl_bon();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<OpenMM::Vec3>& getForces_drl_ang();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<OpenMM::Vec3>& getForces_drl_tor();
	
	/**	
	* @brief Drill
	* @param
	*/
	const std::vector<OpenMM::Vec3>& getForces_drl_n14();

	
	/**	
	* @brief Drill
	* @param
	*/
	void printDrilling(void);



	void setDuMMAtomIndexes(void);

	SimTK::Compound::AtomIndex getCompoundAtomIndex(SimTK::DuMM::AtomIndex);


private:

	// Map mbx2aIx contains only atoms at the origin of mobods
	// topology index and atom index
	std::map< SimTK::MobilizedBodyIndex, std::pair<int, SimTK::Compound::AtomIndex>> mbx2aIx;

	// Map mbx2aIx contains only atoms at the origin of mobods
	//std::map<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex> mbx2aIx;

	// Maps a generalized velocity scale factor for every mobod
	std::map< SimTK::MobilizedBodyIndex, SimTK::Real > mbx2uScale;

	// Binding Site Data: Topologies, AtomIx
	std::vector<int> topologyIXs;
	std::vector<std::vector<int>> amberAtomIXs;

	// Track the Stage of the system
	SimTK::Stage currStage;

	bool useFixmanTorque = false;
	int samplesPerRound = 0;

	Random32 randomEngine;
	SimTK::String rootMobilizer;

	std::vector<BOND_FLEXIBILITY> flexibilities;


	// Context
	Context *myContext;

	// Default return value for non-existing topology atom, pair
	std::pair<int, SimTK::Compound::AtomIndex> errorTopoAtomPair{-1, SimTK::Compound::AtomIndex(SimTK::InvalidIndex)};

};

#endif /*WORLD_H_*/
