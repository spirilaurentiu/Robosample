#include "Replica.hpp"

Replica::Replica(int argIndex) : myIndex(argIndex)
{
}

// Get coordinates from this replica
const std::vector<std::vector<std::pair <bSpecificAtom*, SimTK::Vec3>>>&
	Replica::getAtomsLocationsInGround() const
{

	return atomsLocations;
}

// Get coordinates from this replica
const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
	Replica::get_WORK_AtomsLocationsInGround() const
{
	return WORK_atomsLocations;
}

// Set the coordinates of this replica
// Also allocate memory
void Replica::setAtomsLocationsInGround(
	const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		otherAtomsLocations)
{
	for (auto& topology : otherAtomsLocations){
		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		currentTopologyInfo.reserve(topology.size());

		for(auto& otherAtom : topology){
			bSpecificAtom* batom = otherAtom.first;
			SimTK::Vec3 location = otherAtom.second;
			std::pair<bSpecificAtom*, SimTK::Vec3> atomLocationPair(batom, location);
			currentTopologyInfo.push_back(atomLocationPair);
		}

		atomsLocations.emplace_back(currentTopologyInfo);
	}

	//atomsLocations = otherAtomsLocations;
}

// Set the coordinates of this replica
// Also allocate memory
void Replica::set_WORK_AtomsLocationsInGround(
	const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		otherAtomsLocations)
{
	for (auto& topology : otherAtomsLocations){
		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		currentTopologyInfo.reserve(topology.size());

		for(auto& otherAtom : topology){
			bSpecificAtom* batom = otherAtom.first;
			SimTK::Vec3 location = otherAtom.second;
			std::pair<bSpecificAtom*, SimTK::Vec3> atomLocationPair(batom, location);
			currentTopologyInfo.push_back(atomLocationPair);
		}

		WORK_atomsLocations.emplace_back(currentTopologyInfo);
	}
}

// Update the coordinates of this replica
void Replica::updAtomsLocationsInGround(
	const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		otherAtomsLocations)
{
	int i = -1;
	for (auto& topology : otherAtomsLocations){
		i += 1;

		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>& myTopology =
			atomsLocations[i];

		int j = -1;
		for(auto& otherAtom : topology){
			j += 1;

			std::pair<bSpecificAtom *, SimTK::Vec3>& myAtom = myTopology[j];

			myAtom.first  = otherAtom.first;
			myAtom.second = otherAtom.second;

		}
	}
}

void Replica::updAtomsLocationsInGround_FromWORK()
{
	/* atomsLocations.insert(WORK_atomsLocations.end(),
		WORK_atomsLocations.begin(),
		WORK_atomsLocations.end()); */

		std::cout << "before" << std::endl;
		PrintCoordinates();
		Print_WORK_Coordinates();
		std::cout << std::flush;

	atomsLocations = WORK_atomsLocations;

		std::cout << "after" << std::endl;
		PrintCoordinates();
		Print_WORK_Coordinates();
		std::cout << std::flush;


}

// Update the coordinates of this replica
void Replica::upd_WORK_AtomsLocationsInGround(
	const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		otherAtomsLocations)
{
	int i = -1;
	for (auto& topology : otherAtomsLocations){
		i += 1;

		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>& myTopology = WORK_atomsLocations[i];

		int j = -1;
		for(auto& otherAtom : topology){
			j += 1;

			std::pair<bSpecificAtom *, SimTK::Vec3>&
			myAtom = myTopology[j];

			myAtom.first  = otherAtom.first;
			myAtom.second = otherAtom.second;
		}
	}
}

void Replica::set_WORK_PotentialEnergy_New(SimTK::Real somePotential) {
    WORK_potential = somePotential;
}

SimTK::Real Replica::get_WORK_Jacobian() const
{
	return this->workJacobiansContributions;
}

void Replica::set_WORK_Jacobian(SimTK::Real inpJac)
{
	this->workJacobiansContributions = inpJac;
}

void Replica::setPotentialEnergy_FromWORK() {
    this->potential = this->WORK_potential;
}

// Load atomLocations coordinates into the front world
void Replica::restoreCoordinates() {
}

// Stores coordinates from front world into atomsLocations
void Replica::storeCoordinates() {
}

SimTK::Real Replica::getPotentialEnergy() const {
	return potential;
}

void Replica::setPotentialEnergy(SimTK::Real somePotential) {
	potential = somePotential;
}

SimTK::Real Replica::getTransferedEnergy() const {
	return transferedEnergy;
}

void Replica::setTransferedEnergy(SimTK::Real workArg) {
	this->transferedEnergy = workArg;
}

SimTK::Real Replica::get_WORK_PotentialEnergy_New() const {
	return this->WORK_potential;
}

// void Replica::set_WORK_LastPotentialEnergy(SimTK::Real wpArg) {
//     this->WORK_potential = wpArg;
// }

SimTK::Real Replica::getFixman() const {
	return FixmanPotential;
}

void Replica::setFixman(SimTK::Real somePotential) {
	FixmanPotential = somePotential;
}

void Replica::Print() const {
}

void Replica::PrintCoordinates() const
{
	for(auto& topology : atomsLocations) {
		for(auto& atomCoordinates : topology) {
			std::cout
			//<< (atomCoordinates.first)->inName
			<< (atomCoordinates.first)->getCompoundAtomIndex() << " "
			<< atomCoordinates.second[0] << " "
			<< atomCoordinates.second[1] << " "
			<< atomCoordinates.second[2] << std::endl;
		}
	}
}

void Replica::Print_WORK_Coordinates() const
{
	for(auto& topology : WORK_atomsLocations) {
		for(auto& atomCoordinates : topology) {
			std::cout
			//<< (atomCoordinates.first)->inName
			<< (atomCoordinates.first)->getCompoundAtomIndex() << " "
			<< atomCoordinates.second[0] << " "
			<< atomCoordinates.second[1] << " "
			<< atomCoordinates.second[2] << std::endl;
		}
	}
}

/**
 * Write coordinates to a rst7 file
*/
void Replica::WriteRst7(std::string FN) const

{
	FILE *File = fopen(FN.c_str(), "w+");

	int Natoms = 0;
	for(auto& topology : atomsLocations){
		Natoms += topology.size();
	}

	fprintf(File, "TITLE: Created by Robosample with %d atoms\n", Natoms);
	fprintf(File, "%6d\n", Natoms);

	int atomCnt = -1;
	for(auto& topology : atomsLocations){
		for(auto& atomCoordinates : topology){
			++atomCnt;
			fprintf(File, "%12.7f%12.7f%12.7f", 
				atomCoordinates.second[0] * 10.0,
				atomCoordinates.second[1] * 10.0,
				atomCoordinates.second[2] * 10.0);

			if(atomCnt % 2 == 1){
				fprintf(File, "\n");
			}

		}
	}

	if(atomCnt % 2 == 0){
		fprintf(File, "\n");
	}

	fflush(File);
	fclose(File);
}

/**
 * Print coordinates in Amber rst7 format
 */
void Replica::PrintRst7(void) const
{
	int Natoms = 0;
	for(auto& topology : atomsLocations){
		Natoms += topology.size();
	}

	printf("TITLE: Created by Robosample with %d atoms\n", Natoms);
	printf("%6d\n", Natoms);

	int atomCnt = -1;
	for(auto& topology : atomsLocations){
		for(auto& atomCoordinates : topology){
			++atomCnt;
			printf("%12.7f%12.7f%12.7f", 
				atomCoordinates.second[0],
				atomCoordinates.second[1],
				atomCoordinates.second[2]);

			if(atomCnt % 2 == 1){
				printf("\n");
			}

		}
	}

	if(atomCnt % 2 == 0){
		printf("\n");
	}

}
