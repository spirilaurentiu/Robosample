#pragma once

#include <vector>
#include "Simbody.h"
#include "bSpecificAtom.hpp"
#include "Topology.hpp"
#include "InternalCoordinates.hpp"

class Replica{
public:
	Replica(int index,
			std::vector<bSpecificAtom>& atoms_,
			std::vector<int>& roots_,
			std::vector<Topology>& topologies_,
            InternalCoordinates& internCoords_,
			std::vector<std::vector<int>>& zMatrixTable_
			) :
			  atoms(atoms_)
			, roots(roots_)
			, topologies(topologies_)
			, internCoords(internCoords_)
			, zMatrixTable(zMatrixTable_)
	{}

	const std::vector<std::vector<std::pair <bSpecificAtom*, SimTK::Vec3>>>& getAtomsLocationsInGround() const;
	const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>& get_WORK_AtomsLocationsInGround() const;

	// Reserve memory and set values
	void setAtomsLocationsInGround(const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>& otherAtomsLocations);

	// Reserve memory and set values
	void set_WORK_AtomsLocationsInGround(const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>& otherAtomsLocations);

	// This assumes allocation has been done already
	void updAtomsLocationsInGround(const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>& otherAtomsLocations);

	// Transfers work coordinates into regular cooordinates
	void updAtomsLocationsInGround_FromWORK();

	void upd_WORK_AtomsLocationsInGround(const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>& otherAtomsLocations);

	void set_WORK_PotentialEnergy_New(SimTK::Real somePotential);

	SimTK::Real get_WORK_Jacobian() const;
	void set_WORK_Jacobian(SimTK::Real);

	void setPotentialEnergy_FromWORK();

	// Load atomLocations coordinates into the front world
	void restoreCoordinates();

	// Stores coordinates from front world into atomsLocations
	void storeCoordinates();

	SimTK::Real getPotentialEnergy() const;
	void setPotentialEnergy(SimTK::Real somePotential);

	SimTK::Real getTransferedEnergy() const;
	void setTransferedEnergy(SimTK::Real workArg);

	SimTK::Real get_WORK_PotentialEnergy_New() const;
    // void set_WORK_LastPotentialEnergy(SimTK::Real wpArg);

    SimTK::Real getFixman() const;
	void setFixman(SimTK::Real somePotential);

	void Print() const;
	void PrintCoordinates() const;
	void Print_WORK_Coordinates() const;
	
	void PrintRst7(void) const;
	void WriteRst7(std::string FN) const;

    /** @name Z Matrix and BAT functions
	*/

    /**@{**/
	
	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

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

	// Function to find and return the value for a given AtomIndex
	SimTK::Vec3
	findAtomTarget(
		const std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets,
		SimTK::Compound::AtomIndex searchIndex);

    void setZMatrixTable(const std::vector<std::vector<int>>& newZMatrixTable);

    // zmatrixbat_ Setter for a specific entry
    void setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value) ;

    // zmatrixbat_ Function to get a given row
    const std::vector<SimTK::Real>& getZMatrixBATRow(size_t rowIndex) const;

    // zmatrixbat_ Function to get a given row
    std::vector<SimTK::Real>& updZMatrixBATRow(size_t rowIndex) ;

	// Allocate Z Matrix BAT
	void reallocZMatrixBAT(void);

	// zmatrixbat_
	void
	calcZMatrixBAT(const std::vector< std::vector<
	std::pair <bSpecificAtom *, SimTK::Vec3 > > >&
		otherWorldsAtomsLocations);

	void calcZMatrixBAT_WORK(void);

    // zmatrixbat_ Function to get the value for a given row and column in zMatrixBAT
    SimTK::Real getZMatrixBATValue(size_t rowIndex, size_t colIndex) const ;

    // zmatrixbat_ Function to print the zMatrixBAT
    void PrintZMatrixBAT() const ;

    // zmatrixbat_  Function to add a new row to the zMatrixBAT
    void addZMatrixBATRow(const std::vector<SimTK::Real>& newRow);

	/**
	* @brief zmatrixbat_ Get log of the Cartesian->BAT Jacobian
	* @param
	*/
	SimTK::Real
	calcInternalBATJacobianLog(void);

    // Incrementer function for nofSamples
    void incrementNofSamples();

   std::vector<std::vector<SimTK::Real>>& getZMatrixBATPointer() {
		scout("Address of zMatrixBAT ") << &zMatrixBAT << eolf;
        return (zMatrixBAT);
    }


	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	/**@}**/

private:

	int myIndex = 0;

	// Replica configurations
	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>> atomsLocations;

	// Replica configurations
	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>> WORK_atomsLocations;

	// Replica potential energy
	SimTK::Real potential; // TODO: turn into a vector for worlds
	SimTK::Real WORK_potential; // TODO: turn into a vector for worlds

	SimTK::Real transferedEnergy; // TODO: turn into a vector for worlds
	SimTK::Real workJacobiansContributions; // TODO: turn into a vector for worlds

	SimTK::Real FixmanPotential; // TODO: turn into a vector for worlds


	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////
	std::vector<bSpecificAtom>& atoms;

	std::vector<int>& roots;
	std::vector<Topology>& topologies;
	InternalCoordinates& internCoords;
	std::vector<std::vector<int>>& zMatrixTable;

	std::vector<std::vector<SimTK::Real>> zMatrixBAT;
	

	// BAT
	int nofSamples = 0;	

	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

};
