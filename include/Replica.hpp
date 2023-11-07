#pragma once

#include <vector>
#include "Simbody.h"
#include "bSpecificAtom.hpp"

class Replica{
public:
	Replica(int index);

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
};