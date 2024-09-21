#include "Replica.hpp"

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

	atomsLocations = WORK_atomsLocations;

}

// Update the coordinates of this replica
void Replica::upd_WORK_AtomsLocationsInGround(
	const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		otherAtomsLocations)
{
	int i = -1;
	for (auto& topology : otherAtomsLocations){
		i += 1;

		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>& myTopology =
			WORK_atomsLocations[i];

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

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> Replica::getCoordinates() const {
	std::vector<double> x, y, z;

	for (const auto& topology : atomsLocations) {
		for(auto& atomCoordinates : topology) {
			x.push_back(atomCoordinates.second[0]);
			y.push_back(atomCoordinates.second[1]);
			z.push_back(atomCoordinates.second[2]);
		}
	}

	return std::make_tuple(x, y, z);
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




// ===========================================================================
// ===========================================================================
// ZMatrix BAT
// ===========================================================================
// ===========================================================================

/*!
 * <!-- Takes coordinates from molecule topoIx and puts them into atomTargets
 * otherWorldsAtomsLocations: Pairs of (atom, and its position) within
 * a vector of Topologies
 * atomTargets: a map of atoms' Comopund atom index to positions-->
*/
void
Replica::extractAtomTargets(
	int topoIx,
	const std::vector<std::vector<std::pair
		<bSpecificAtom *, SimTK::Vec3> > >& otherWorldsAtomsLocations,
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets
)
{
	for(std::size_t j = 0; j < otherWorldsAtomsLocations[topoIx].size(); j++){
		auto atomIndex = otherWorldsAtomsLocations[topoIx][j].first->getCompoundAtomIndex();
		auto location = otherWorldsAtomsLocations[topoIx][j].second;
		atomTargets.insert(std::make_pair(atomIndex, location));
	}
}

/*!
 * <!--	zmatrixbat_ Function to find and return the value for a given 
 * AtomIndex -->
*/
SimTK::Vec3
Replica::findAtomTarget(
	const std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets,
	SimTK::Compound::AtomIndex searchIndex)
{
	auto it = atomTargets.find(searchIndex);

	if (it != atomTargets.end()) {
		return it->second;
	} else {
		return SimTK::Vec3(SimTK::NaN);
	}
}

void Replica::setZMatrixTable(const std::vector<std::vector<int>>& newZMatrixTable) {
    zMatrixTable = newZMatrixTable;
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Replica::setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value) {

	// Set the value at the specified position
	zMatrixBAT[rowIndex][colIndex] = value;
}

/*!
 * <!--	zmatrixbat_ -->
*/
const std::vector<SimTK::Real>& Replica::getZMatrixBATRow(size_t rowIndex) const {

	// assert(rowIndex >= 0);
	// assert(rowIndex < static_cast<int>(zMatrixBAT.size()));

	// Check if the indices are within bounds
	return zMatrixBAT[rowIndex];

}

/*!
 * <!--	zmatrixbat_ -->
*/
std::vector<SimTK::Real>& Replica::updZMatrixBATRow(size_t rowIndex) {

	// assert(rowIndex >= 0);
	// assert(rowIndex < static_cast<int>(zMatrixBAT.size()));

	// Check if the indices are within bounds
	return zMatrixBAT[rowIndex];

}

/*!
 * <!--	zmatrixbat_ -->
*/


/*!
 * <!-- zmatrixbat_ Allocate Z Matrix BAT -->
*/
void Replica::reallocZMatrixBAT(void){

	zMatrixBAT.resize(zMatrixTable.size());
	for (auto& row : zMatrixBAT) {
		row.resize(3, SimTK::NaN);
	}

}




/*!
 * <!-- zmatrixbat_ Calculate Z-matrix -->
*/
void
Replica::calcZMatrixBAT_WORK(void)
{

	// Iterate molecules
	int allCnt = 0;

	int topoIx = 0;

	// Get locations of this molecule
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
	for (int i = 0; i < topologies.size(); i++) {
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
		extractAtomTargets(i, get_WORK_AtomsLocationsInGround(), temp);
		atomTargets.insert(temp.begin(), temp.end());
	}

	int rowCnt = 0;
	SimTK::Real bondLength, bondBend, bondTorsion;
	for (const auto& row : zMatrixTable) {
		
		bondLength = SimTK::NaN;
		bondBend = SimTK::NaN;
		bondTorsion = SimTK::NaN;

		SimTK::Compound::AtomIndex a0_cAIx, a1_cAIx;
		
		// Calculate bond length
		a0_cAIx = atoms[row[0]].getCompoundAtomIndex();
		a1_cAIx = atoms[row[1]].getCompoundAtomIndex();
		SimTK::Vec3 a0loc = findAtomTarget(atomTargets, a0_cAIx);
		SimTK::Vec3 a1loc = findAtomTarget(atomTargets, a1_cAIx);

		SimTK::Vec3 v_a0a1 = a0loc - a1loc;
		SimTK::Real bondLength = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

		if(row[2] >= 0){

			SimTK::Compound::AtomIndex a2_cAIx;
			a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
			SimTK::Vec3 a2loc = findAtomTarget(atomTargets, a2_cAIx);

			// Calculate angle
			UnitVec3 v1(v_a0a1);
			UnitVec3 v2(a2loc - a1loc);

			Real dotProduct = SimTK::dot(v1, v2);
			assert(dotProduct < 1.1);
			assert(dotProduct > -1.1);
			if (dotProduct > 1.0) dotProduct = 1.0;
			if (dotProduct < -1.0) dotProduct = -1.0;
			bondBend = std::acos(dotProduct);

			if(row[3] >= 0){
				SimTK::Compound::AtomIndex
					a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
				SimTK::Vec3 a3loc = findAtomTarget(atomTargets, a3_cAIx);

				bondTorsion = bDihedral(a0loc, a1loc, a2loc, a3loc);
			}

		} // angle
		
		setZMatrixBATValue(rowCnt, 0, bondLength);
		setZMatrixBATValue(rowCnt, 1, bondBend);
		setZMatrixBATValue(rowCnt, 2, bondTorsion);

		if(row[3] == -2){
			topoIx++;
		}

		rowCnt++;

	} // every zMatrix row		

}



/*!
 * <!-- zmatrixbat_ Calculate Z-matrix -->
*/
void
Replica::calcZMatrixBAT(
	const std::vector< std::vector<
		std::pair <bSpecificAtom *, SimTK::Vec3 > > >&
		otherWorldsAtomsLocations)
{

	// Iterate molecules
	int allCnt = 0;

	int topoIx = 0;

	// Get locations of this molecule
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
	for (int i = 0; i < topologies.size(); i++) {
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
		extractAtomTargets(i, otherWorldsAtomsLocations, temp);
		atomTargets.insert(temp.begin(), temp.end());
	}

	int rowCnt = 0;
	SimTK::Real bondLength, bondBend, bondTorsion;
	for (const auto& row : zMatrixTable) {
		
		bondLength = SimTK::NaN;
		bondBend = SimTK::NaN;
		bondTorsion = SimTK::NaN;

		SimTK::Compound::AtomIndex a0_cAIx, a1_cAIx;
		
		// Calculate bond length
		a0_cAIx = atoms[row[0]].getCompoundAtomIndex();
		a1_cAIx = atoms[row[1]].getCompoundAtomIndex();
		SimTK::Vec3 a0loc = findAtomTarget(atomTargets, a0_cAIx);
		SimTK::Vec3 a1loc = findAtomTarget(atomTargets, a1_cAIx);

		SimTK::Vec3 v_a0a1 = a0loc - a1loc;
		SimTK::Real bondLength = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

		if(row[2] >= 0){

			SimTK::Compound::AtomIndex a2_cAIx;
			a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
			SimTK::Vec3 a2loc = findAtomTarget(atomTargets, a2_cAIx);

			// Calculate angle
			UnitVec3 v1(v_a0a1);
			UnitVec3 v2(a2loc - a1loc);

			Real dotProduct = SimTK::dot(v1, v2);
			assert(dotProduct < 1.1);
			assert(dotProduct > -1.1);
			if (dotProduct > 1.0) dotProduct = 1.0;
			if (dotProduct < -1.0) dotProduct = -1.0;
			bondBend = std::acos(dotProduct);

			if(row[3] >= 0){
				SimTK::Compound::AtomIndex
					a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
				SimTK::Vec3 a3loc = findAtomTarget(atomTargets, a3_cAIx);

				bondTorsion = bDihedral(a0loc, a1loc, a2loc, a3loc);
			}

		} // angle
		
		setZMatrixBATValue(rowCnt, 0, bondLength);
		setZMatrixBATValue(rowCnt, 1, bondBend);
		setZMatrixBATValue(rowCnt, 2, bondTorsion);

		if(row[3] == -2){
			topoIx++;
		}

		rowCnt++;

	} // every zMatrix row		

}


/*!
 * <!--	zmatrixbat_ -->
*/
SimTK::Real Replica::getZMatrixBATValue(size_t rowIndex, size_t colIndex) const {
	// Check if the indices are within bounds
	if (rowIndex < zMatrixBAT.size() && colIndex < zMatrixBAT[0].size()) {
		// Return the value at the specified position
		return zMatrixBAT[rowIndex][colIndex];
	} else {
		// Indices are out of bounds, handle this case accordingly
		return SimTK::NaN;
	}
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Replica::PrintZMatrixBAT() const {

	int bati = 0;
	for (const auto& row : zMatrixBAT) {

		scout("zm ") ; // indicator

		for(const auto tabValue : zMatrixTable[bati]){
			std::cout << tabValue << " ";
		}

		for (SimTK::Real value : row) {
			std::cout << std::setw(9) << value << " ";
		}
		std::cout << std::endl;

		bati++;
	}

}

/*!
 * <!--	zmatrixbat_ -->
*/
void Replica::addZMatrixBATRow(const std::vector<SimTK::Real>& newRow) {
	zMatrixBAT.push_back(newRow);
}

/*!
 * <!-- zmatrixbat_ BAT JAcobian -->
*/
SimTK::Real
Replica::calcInternalBATJacobianLog(void)
	{

		// Get log of the Cartesian->BAT Jacobian
		SimTK::Real logJacBAT = 0.0;

		for(size_t zCnt = 0; zCnt = zMatrixBAT.size(); zCnt++){

				// Get bond term
				SimTK::Real currBond = zMatrixBAT[zCnt][0];
				
				if(currBond != SimTK::NaN){
				
					logJacBAT += 4.0 * std::log(currBond);
				}

				// Get the angle term
				SimTK::Real currAngle = zMatrixBAT[zCnt][1];

				if(currAngle != SimTK::NaN){

					logJacBAT += 2.0 * std::log(std::sin(currAngle));
					
				}

		}

		return logJacBAT;

	}

/*!
 * <!--	zmatrixbat_ Incrementer function implementation -->
*/
void Replica::incrementWorldsNofSamples() {
    ++allWorldsNofSamples;
}

/*!
 * <!--	zmatrixbat_ Incrementer function implementation -->
*/
void Replica::incrementWorldsNofSamples(int howMany) {
    allWorldsNofSamples += howMany;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ZMatrix BAT
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------




