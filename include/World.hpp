#pragma once

/* -------------------------------------------------------------------------- *
 *		                       Robosampling                           *
 * -------------------------------------------------------------------------- *
 * This is part of Robosample		                                      *
 */

#include "TopologyElements.hpp"
#include <functional>
#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include <thread>
#include <array>
#include <math.h>

#include "Simbody.h"
#include "Molmodel.h"
#include "ParaMolecularDecorator.hpp"
#include "FixmanTorque.hpp"

#ifndef BaseSampler
#define BaseSampler HMCSampler
#endif

#include "server.hpp"
#include "Topology.hpp"
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

class Context;

template <typename T>
std::string vecToString(const std::vector<T>& v) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        oss << std::scientific
            << std::showpos
            << std::setprecision(2)
            << v[i];
        if (i + 1 < v.size()) oss << ", ";
    }
    oss << "]";
    return oss.str();
}

inline std::string atomsToString(const std::vector<std::string>& atoms) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < atoms.size(); ++i) {
        oss << atoms[i];
        if (i + 1 < atoms.size()) oss << ", ";
    }
    oss << "]";
    return oss.str();
}

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

// Describes bonds involving a root atom
struct RootAtomBond {
	std::size_t topologyIndex = 0;

	SimTK::Compound::AtomIndex childCAIx;
	SimTK::Compound::AtomIndex parentCAIx;

	SimTK::MobilizedBodyIndex childMBIx;
	SimTK::MobilizedBodyIndex parentMBIx;
};

// Describes atom bonds and angles between two linked rigid bodies
struct RigidBodyAtomBond {
	std::size_t topologyIndex = 0;
	SimTK::BondMobility::Mobility mobility = SimTK::BondMobility::Default;

	std::size_t childAtomGlobalIndex = -1;
	std::size_t parentAtomGlobalIndex = -1;
	std::size_t grandParentAtomGlobalIndex = -1;

	SimTK::Compound::AtomIndex childCAIx;
	SimTK::Compound::AtomIndex parentCAIx;
	SimTK::Compound::AtomIndex grandParentCAIx;

	// Compound atom index of the root atom in the parent rigid body
	SimTK::Compound::AtomIndex parentMobodRootCAIx;

	SimTK::MobilizedBodyIndex childMBIx;
	SimTK::MobilizedBodyIndex parentMBIx;
	SimTK::MobilizedBodyIndex grandParentMBIx;
};

struct RigidBond {
	std::size_t topologyIndex = 0;
	SimTK::Compound::AtomIndex childCAIx;
	SimTK::Compound::AtomIndex parentCAIx;
};

struct RigidAngle {
	std::size_t topologyIndex = 0;
	SimTK::Compound::AtomIndex cAIx1, cAIx2, cAIx3;
};

struct RigidTorsion {
	std::size_t topologyIndex = 0;
	SimTK::Compound::AtomIndex cAIx1, cAIx2, cAIx3, cAIx4;
};

enum class ROOT_MOBILITY : int {
	FREE = 0,
	CARTESIAN,
	WELD,
	FREE_LINE,
	BALL,
	PIN
};

//==============================================================================
//                   CLASS World
//==============================================================================
/**
 *  Contains a Symbody system and additional data that define a regimen
 **/
class World {
public:

	void setAtomTargetLocationsToState(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);
	void updateFramesFromTopologies();

	explicit World(int worldIndex, Span<Topology> topo, bool testing, const ZMatrix& _zMatrix, bool isVisual=true, SimTK::Real visualizerFrequency = 0.0015);

	const std::vector<SimTK::Compound::AtomTargetLocations>& getAtomTargetLocationsCache() const {
		return atomTargetLocaltionsCache;
	}

	void setMobodLocks(const std::vector<std::vector<SimTK::MobilizedBodyIndex>>& mobodLocks) {
		this->mobodLocks = mobodLocks;
	}

	void generateDummParams(
		const std::vector<RoboAtom>& atoms,
		const std::vector<RoboBond>& bonds,
		const std::vector<RoboAngle>& angles,
		const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions,
		const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions
	);
	
	void modelTopologies();

	SimTK::Real getRecommendedTimesteps(void);

	//=========================================================================
	//                   CONSTRAINTS
	//=========================================================================

	/** Add contact constraints to specific bodies **/
	void addRodConstraint(SimTK::State& someState);

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
	void updateTaskSpace(const SimTK::State& someState);

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
	void calcStationJacobian(const SimTK::State& someState,
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

	/** Assign a scale factor for generalized velocities to every mobilized
	body **/
	void setUScaleFactorsToMobods(void);

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
	SimTK::Vec3 getGeometricCenterOfSelection(const SimTK::State & state);

	/** Nice print helper for get/setAtomsLocations */
	void PrintAtomsLocations(const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >& someAtomsLocations);
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
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations	
	) const ;

	/**
	 * Maximum distance between two corresponding atoms
	*/
	std::pair<int, SimTK::Real> maxAtomDeviation(
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
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
	 *
	*/
	SimTK::BondMobility::Mobility determineMobilityFrom_H(SimTK::MobilizedBodyIndex mbx, SimTK::State& someState);

	/**
	 *
	*/
	SimTK::Real getRootAngle(
		Topology& topology,
		SimTK::Compound::AtomIndex rootAIx,
		const SimTK::State& someState
	);

	/**
	 * Calc X_FM transforms for reconstruction
	*/
	SimTK::Transform calcX_FMTransforms(
		Topology& topology,
		SimTK::Compound::AtomIndex aIx,
		const SimTK::State& someState);

	/** Update Gmolmodel Atom Cartesian coordinates according to
	Molmodel Compound which in turn relizes Position and uses matter
	 to calculate locations. **/
	void updateAtomListsFromSimbody(const SimTK::State &state);

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

	const SimTK::CompoundSystem& getCompoundSystem() const { return *compoundSystem; }
	SimTK::CompoundSystem& updCompoundSystem() { return *compoundSystem; }

	const SimTK::SimbodyMatterSubsystem& getMatterSubsystem() const { return *matter; }
	SimTK::SimbodyMatterSubsystem& updMatterSubsystem() { return *matter; }

	const SimTK::GeneralForceSubsystem& getForces() const { return *forces; }
	SimTK::GeneralForceSubsystem& updForces() { return *forces; }

	const SimTK::DuMMForceFieldSubsystem& getForceField() const { return *forceField; }
	SimTK::DuMMForceFieldSubsystem& updForceField() { return *forceField; }

	const SimTK::VerletIntegrator& getIntegrator() const { return *integrator; }
	SimTK::VerletIntegrator& updIntegrator() { return *integrator; }

	const SimTK::TimeStepper& getTimeStepper() const { return *timeStepper; }
	SimTK::TimeStepper& updTimeStepper() { return *timeStepper; }

	/** Return true if the Fixman torque flag is set **/
	bool isUsingFixmanTorque() const;
	//...............

	// Calculate Fixman potential
	SimTK::Real calcFixman();

	/** Generate a number of samples **/
	bool generateSamples(int howMany, std::stringstream& worldOutStream, const std::string& header, bool verbose);

	//...................
	// --- Statistics ---
	//...................
	/** How many samples did we have so far **/
	std::size_t getNofSamples() const;

	/** Sampler manipulation functions **/
	std::size_t getNofSamplers() const;

	/** Add a sampler to the World **/
	bool addSampler(SamplerName samplerName,
		IntegratorType integratorType,
		ThermostatName thermostatName,
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

	// Get the (potential) energy transfer
	// If any of the Q, U or tau is actively modifyied by the sampler
	// the Jacobian of that transformation will be included too
	SimTK::Real getWorkOrHeat(void);

	// Get the (potential) energy transfer in the form of work
	// If any of the Q, U or tau is actively modifyied by the sampler
	// the Jacobian of that transformation will be included too
	SimTK::Real getWork(void) const;

	// Set initial values of X_PF or X_BM
	void setTransformsMeansToIni(void);

	// Set initial values of X_PF or X_BM
	void setTransformsMeansToCurrent(SimTK::State& someState);

	// Set initial values of X_PF or X_BM
	void setTransformsMeans(const std::vector<SimTK::Real>& givenX_PF,
		const std::vector<SimTK::Real>& givenX_BM);

	// Set X_PF and X_BM related values
	void setTransformsMeansToMin(AmberReader &amberReader);

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
		const SimTK::State &someState,
		const SimTK::MobilizedBody &inBodyA,
		const SimTK::MobilizedBody &ofBodyB,
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

	void PrintDefaultTransforms() const;
	void PrintAllTransforms() const;
	void PrintXFMs() const;

	void PrintXBMps() const;
	const SimTK::Vector & getBMps();
	const SimTK::Vector & getPFrs();

	// Get Qs
	int getNQs(void);
	int getNUs(void);
	const SimTK::Vector & getAdvancedQs();
	const void PrintAdvancedQs() const;

	const SimTK::Vector & getAdvancedUs();

	void PrintBATFromSimbody() const;
	void calcSimbodyBAT_TODEL(std::vector<std::vector<int>>& ZMatrix, std::vector<SimTK::Real>& BONDLengths, std::vector<SimTK::Real>& ANGLEBends, std::vector<SimTK::Real>& TORSIONAngles);
	void calcSimbodyBAT(std::vector<std::vector<int>>& ZMatrix, std::vector<SimTK::Real>& BONDLengths, std::vector<SimTK::Real>& ANGLEBends, std::vector<SimTK::Real>& TORSIONAngles);

	void setFlexibilites(const std::vector<std::vector<SimTK::MobilizedBodyIndex>>& flexibilities_UNCHAINED) {
		this->flexibilities_UNCHAINED = flexibilities_UNCHAINED;
	}

public:

	// The three S: Study, System and State related.
	
	// System->MultibodySystem->MolecularMechanicsSystems->CompoundSystem
	// This is non-copyable and non-movable, so we use a unique_ptr because World needs to be copyable/movable.
	std::unique_ptr<SimTK::CompoundSystem> compoundSystem;

	// Subsystem->SimbodyMatterSubsystem
	// This is non-copyable and non-movable, so we use a unique_ptr because World needs to be copyable/movable.
	std::unique_ptr<SimTK::SimbodyMatterSubsystem> matter;

	// Subsystem->ForceSubsystem->GeneralForceSubsystem
	// This is non-copyable and non-movable, so we use a unique_ptr because World needs to be copyable/movable.
	std::unique_ptr<SimTK::GeneralForceSubsystem> forces;

	// Subsystem->ForceSubsystem->DuMMForceFieldSubsystem
	// This is non-copyable and non-movable, so we use a unique_ptr because World needs to be copyable/movable.
	std::unique_ptr<SimTK::DuMMForceFieldSubsystem> forceField;

	// --- Simulation ---
	std::unique_ptr<SimTK::VerletIntegrator> integrator;
	std::unique_ptr<SimTK::TimeStepper> timeStepper;
	// TODO they belong to Sampler, not World


	std::vector<std::unique_ptr<BaseSampler>> samplers;

	std::vector<std::vector<SimTK::MobilizedBodyIndex>> flexibilities_UNCHAINED;

	SimTK::Vector BMps;
	SimTK::Vector PFrs;


	// std::vector<std::unique_ptr<SimTK::ConformationalController>> controller;
	// std::vector<std::unique_ptr<SimTK::Force::Custom>> controlForce;

	/** Nof molecules **/
	std::size_t numMolecules = 0;
	std::size_t numAtoms = 0;


	/** Molecules (topologies<-Compounds) objects **/
	Span<Topology> topologies;
	std::vector<std::string> roots;
	std::vector<std::string> rootMobilitiesStr;

	/** Joint types **/
	//std::map< SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility> mbx2mobility;

	// 
	std::vector<std::vector<int>> zMatrixTable;
	std::vector<std::vector<SimTK::Real>> zMatrixBAT;

	// --- Thermodynamics ---
	SimTK::Real temperature;

	// // Contact related
	// std::unique_ptr<ContactTrackerSubsystem> tracker;
	// std::unique_ptr<CompliantContactSubsystem> contactForces;
	// ContactCliqueId clique1;
	// std::unique_ptr<MobilizedBody::Weld> membrane;
	// std::unique_ptr<Body::Rigid> memBody;
	// //...............

	// --- Statistics ---
	std::vector<SimTK::Real> acosX_PF00;
	std::vector<SimTK::Real> normX_BMp;
	std::vector<SimTK::Real> acosX_PF00_means;
	std::vector<SimTK::Real> normX_BMp_means;

	//std::vector<SimTK::Real> CppQs;

	//...............

	// // --- Graphics ---
	bool visual;

	// // Our decorations
	// std::unique_ptr<ParaMolecularDecorator> paraMolecularDecorator;

	// Decoration subsystem
	std::unique_ptr<SimTK::DecorationSubsystem> decorations;

	// Visualizer
#ifdef BUILD_VISUALIZER
	std::unique_ptr<SimTK::Visualizer> visualizer;

	// Visualizer reporter
	std::unique_ptr<SimTK::Visualizer::Reporter> visualizerReporter;
#endif

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

	const int getOwnIndex(void) const{
		return ownWorldIndex;
	}

	const std::vector<RootAtomBond>& getRootAtomBonds() const {
		return rootAtomBonds;
	}

	const std::vector<RigidBodyAtomBond>& getRigidBodyAtomBonds() const {
		return rigidBodyAtomBonds;
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

	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	//std::vector<std::vector<int>>& zMatrixTable;
	//std::vector<std::vector<SimTK::Real>>& zMatrixBAT;

	//setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value);
	//void calcZMatrixBAT(SimTK::State& someState);

	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	void addTestRigidBond(int globalIndex1, int globalIndex2) {
		testRigidBonds.push_back({globalIndex1, globalIndex2});
	}

	void addTestNonRigidBond(int globalIndex1, int globalIndex2) {
		testNonRigidBonds.push_back({globalIndex1, globalIndex2});
	}

	void addTestRigidAngle(int globalIndex1, int globalIndex2, int globalIndex3) {
		testRigidAngles.push_back({globalIndex1, globalIndex2, globalIndex3});
	}

	void addTestNonRigidAngle(int globalIndex1, int globalIndex2, int globalIndex3) {
		testNonRigidAngles.push_back({globalIndex1, globalIndex2, globalIndex3});
	}

	void addTestRigidPeriodicTorsion(int globalIndex1, int globalIndex2, int globalIndex3, int globalIndex4) {
		testRigidPeriodicTorsions.push_back({globalIndex1, globalIndex2, globalIndex3, globalIndex4});
	}

	void addTestNonRigidPeriodicTorsion(int globalIndex1, int globalIndex2, int globalIndex3, int globalIndex4) {
		testNonRigidPeriodicTorsions.push_back({globalIndex1, globalIndex2, globalIndex3, globalIndex4});
	}

	void addTestRigidHarmonicTorsion(int globalIndex1, int globalIndex2, int globalIndex3, int globalIndex4) {
		testRigidHarmonicTorsions.push_back({globalIndex1, globalIndex2, globalIndex3, globalIndex4});
	}

	void addTestNonRigidHarmonicTorsion(int globalIndex1, int globalIndex2, int globalIndex3, int globalIndex4) {
		testNonRigidHarmonicTorsions.push_back({globalIndex1, globalIndex2, globalIndex3, globalIndex4});
	}

	const std::vector<std::array<int, 2>>& getTestRigidBonds() const {
		return testRigidBonds;
	}

	const std::vector<std::array<int, 2>>& getTestNonRigidBonds() const {
		return testNonRigidBonds;
	}

	const std::vector<std::array<int, 3>>& getTestRigidAngles() const {
		return testRigidAngles;
	}

	const std::vector<std::array<int, 3>>& getTestNonRigidAngles() const {
		return testNonRigidAngles;
	}

	const std::vector<std::array<int, 4>>& getTestRigidPeriodicTorsions() const {
		return testRigidPeriodicTorsions;
	}

	const std::vector<std::array<int, 4>>& getTestNonRigidPeriodicTorsions() const {
		return testNonRigidPeriodicTorsions;
	}

	const std::vector<std::array<int, 4>>& getTestRigidHarmonicTorsions() const {
		return testRigidHarmonicTorsions;
	}

	const std::vector<std::array<int, 4>>& getTestNonRigidHarmonicTorsions() const {
		return testNonRigidHarmonicTorsions;
	}
	

	const std::vector<std::vector<SimTK::Real>>& getMatchAtomTargetLocationsResiduals() const {
		SimTK_ASSERT_ALWAYS(testing, "getMatchAtomTargetLocationsResiduals called outside testing mode");
		return matchAtomTargetLocationsResiduals;
	}

	const std::vector<SimTK::Real>& getCumulativeCartesianDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getCumulativeCartesianDisplacements called outside testing mode");
		return cumulativeCartesianDisplacements;
	}

	const std::vector<SimTK::Real>& getCumulativeBondDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getCumulativeBondDisplacements called outside testing mode");
		return cumulativeBondDisplacements;
	}

	const std::vector<SimTK::Real>& getCumulativeAngleDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getCumulativeAngleDisplacements called outside testing mode");
		return cumulativeAngleDisplacements;
	}

	const std::vector<SimTK::Real>& getCumulativeTorsionDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getCumulativeTorsionDisplacements called outside testing mode");
		return cumulativeTorsionDisplacements;
	}

	const std::vector<SimTK::Real>& getOpenMMCumulativeCartesianDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getOpenMMCumulativeCartesianDisplacements called outside testing mode");
		return openmmCumulativeCartesianDisplacements;
	}

	const std::vector<SimTK::Real>& getOpenMMCumulativeBondDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getOpenMMCumulativeBondDisplacements called outside testing mode");
		return openmmCumulativeBondDisplacements;
	}

	const std::vector<SimTK::Real>& getOpenMMCumulativeAngleDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getOpenMMCumulativeAngleDisplacements called outside testing mode");
		return openmmCumulativeAngleDisplacements;
	}

	const std::vector<SimTK::Real>& getOpenMMCumulativeTorsionDisplacements() const {
		SimTK_ASSERT_ALWAYS(testing, "getOpenMMCumulativeTorsionDisplacements called outside testing mode");
		return openmmCumulativeTorsionDisplacements;
	}

	const std::vector<std::pair<bool, SimTK::Real>>& getAcceptanceRMSD() const {
		SimTK_ASSERT_ALWAYS(testing, "getRMSD called outside testing mode");
		return acceptanceRMSD;
	}

	const std::vector<SimTK::Real>& getRigidBodyBondRMSDInNm() const {
		SimTK_ASSERT_ALWAYS(testing, "getRigidBodyBondRMSDInNm called outside testing mode");
		return rigidBodyBondRMSDInNm;
	}

	const std::vector<SimTK::Real>& getRigidBodyAngleDriftInRad() const {
		SimTK_ASSERT_ALWAYS(testing, "getRigidBodyAngleDriftInRad called outside testing mode");
		return rigidBodyAngleDriftInRad;
	}

	const std::vector<SimTK::Real>& getRigidBodyProperTorsionDriftInRad() const {
		SimTK_ASSERT_ALWAYS(testing, "getRigidBodyProperTorsionDriftInRad called outside testing mode");
		return rigidBodyProperTorsionDriftInRad;
	}

	const std::vector<SimTK::Real>& getRigidBodyImproperTorsionDriftInRad() const {
		SimTK_ASSERT_ALWAYS(testing, "getRigidBodyImproperTorsionDriftInRad called outside testing mode");
		return rigidBodyImproperTorsionDriftInRad;
	}

	bool isOverconstrained(const SimTK::State& state, std::ostream& out) const;

private:
	SimTK::Real findDecorrelationTime(const SimTK::State& state, int equilSteps, int tuneSteps, SimTK::Real timestep);
	void optimizeCoordinates(std::vector<std::vector<double>>& coordinates, const std::vector<double>& scale) const;
	SimTK::Real integratedAutocorrelation(const std::vector<SimTK::Real>& x, int maxLag = -1) const;
	bool tuned = false;

	bool testing = false;

	void checkCoordinateTransfer(
		const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets,
		std::vector<SimTK::Real>& matchAtomTargetLocationsResiduals,
		SimTK::Real& cumulDiffCartesian,
		SimTK::Real& cumulDiffBonds,
		SimTK::Real& cumulDiffAngles,
		SimTK::Real& cumulDiffDihedrals
	);

	std::vector<std::vector<SimTK::Real>> matchAtomTargetLocationsResiduals;
	std::vector<SimTK::Real> cumulativeCartesianDisplacements;
	std::vector<SimTK::Real> cumulativeBondDisplacements;
	std::vector<SimTK::Real> cumulativeAngleDisplacements;
	std::vector<SimTK::Real> cumulativeTorsionDisplacements;

	std::vector<SimTK::Real> openmmCumulativeCartesianDisplacements;
	std::vector<SimTK::Real> openmmCumulativeBondDisplacements;
	std::vector<SimTK::Real> openmmCumulativeAngleDisplacements;
	std::vector<SimTK::Real> openmmCumulativeTorsionDisplacements;

	std::vector<SimTK::Compound::AtomTargetLocations> atomTargetLocaltionsCache;
	std::vector<SimTK::Compound::AtomTargetLocations> atomTargetLocaltionsCacheOld;
	std::vector<std::pair<bool, SimTK::Real>> acceptanceRMSD;

	std::vector<RigidBond> rigidBonds;
	std::vector<SimTK::Real> rigidBodyBondRMSDInNm;

	std::vector<RigidAngle> rigidAngles;
	std::vector<SimTK::Real> rigidBodyAngleDriftInRad;

	std::vector<RigidTorsion> rigidProperTorsions, rigidImproperTorsions;
	std::vector<SimTK::Real> rigidBodyProperTorsionDriftInRad, rigidBodyImproperTorsionDriftInRad;

	// Map mbx2aIx contains only atoms at the origin of mobods
	// topology index and atom index
	std::map< SimTK::MobilizedBodyIndex, std::pair<int, SimTK::Compound::AtomIndex>> mbx2aIx;
	std::vector<RootAtomBond> rootAtomBonds;
	std::vector<RigidBodyAtomBond> rigidBodyAtomBonds;

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

	std::vector<std::vector<SimTK::MobilizedBodyIndex>> mobodLocks;

	// Default return value for non-existing topology atom, pair
	std::pair<int, SimTK::Compound::AtomIndex> errorTopoAtomPair{-1, SimTK::Compound::AtomIndex(SimTK::InvalidIndex)};

	std::reference_wrapper<const ZMatrix> zMatrix;

	std::vector<std::array<int, 2>> testRigidBonds, testNonRigidBonds;
	std::vector<std::array<int, 3>> testRigidAngles, testNonRigidAngles;
	std::vector<std::array<int, 4>> testRigidPeriodicTorsions, testNonRigidPeriodicTorsions, testRigidHarmonicTorsions, testNonRigidHarmonicTorsions;
};
