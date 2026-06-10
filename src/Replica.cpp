#include "Replica.hpp"

#include "Compound.h"

// Set the coordinates of this replica
void Replica::setAtomsLocationsInGround(
    const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) {
    atomsLocations = atomTargets;

    // Allocate memory for coordinate buffers
    if (x.size() == 0) {
        for (const auto& topology : atomsLocations) {
            for (const auto& locations : topology) {
                x.push_back(0.0);
                y.push_back(0.0);
                z.push_back(0.0);
            }
        }
    }

    // Update coordinate buffers
    int numAtomsInMols = 0;
    for (const auto& topology : atomsLocations) {
        for (SimTK::Compound::AtomIndex atomIndex(0); atomIndex < topology.size(); ++atomIndex) {
            x[atomIndex + numAtomsInMols] = topology[atomIndex][0];
            y[atomIndex + numAtomsInMols] = topology[atomIndex][1];
            z[atomIndex + numAtomsInMols] = topology[atomIndex][2];
        }
        numAtomsInMols += topology.size();
    }
}

void Replica::updAtomsLocationsInGround_FromWORK() {
    atomsLocations = WORK_atomsLocations;

    // Allocate memory for coordinate buffers
    if (x.size() == 0) {
        for (const auto& topology : atomsLocations) {
            for (const auto& locations : topology) {
                x.push_back(0.0);
                y.push_back(0.0);
                z.push_back(0.0);
            }
        }
    }

    // Update coordinate buffers
    int numAtomsInMols = 0;
    for (const auto& topology : atomsLocations) {
        for (SimTK::Compound::AtomIndex atomIndex(0); atomIndex < topology.size(); ++atomIndex) {
            x[atomIndex + numAtomsInMols] = topology[atomIndex][0];
            y[atomIndex + numAtomsInMols] = topology[atomIndex][1];
            z[atomIndex + numAtomsInMols] = topology[atomIndex][2];
        }
        numAtomsInMols += topology.size();
    }
}


// Load atomLocations coordinates into the front world
void Replica::restoreCoordinates() {
}

// Stores coordinates from front world into atomsLocations
void Replica::storeCoordinates() {
}

void Replica::Print() const {
    // SimTK_ASSERT_ALWAYS(false, "Replica::Print not implemented yet.");
}

void Replica::PrintCoordinates() const {
    SimTK_ASSERT_ALWAYS(false, "Replica::PrintCoordinates not implemented yet.");

    // for(auto& topology : atomsLocations) {
    // 	for(auto& atomCoordinates : topology) {
    // 		std::cout
    // 		//<< (atomCoordinates.first)->inName
    // 		<< (atomCoordinates.first)->identity.compoundAtomIndex << " "
    // 		<< atomCoordinates.second[0] << " "
    // 		<< atomCoordinates.second[1] << " "
    // 		<< atomCoordinates.second[2] << std::endl;
    // 	}
    // }
}

void Replica::Print_WORK_Coordinates() const {
    SimTK_ASSERT_ALWAYS(false, "Replica::Print_WORK_Coordinates not implemented yet.");

    // for(auto& topology : WORK_atomsLocations) {
    // 	for(auto& atomCoordinates : topology) {
    // 		std::cout
    // 		//<< (atomCoordinates.first)->inName
    // 		<< (atomCoordinates.first)->identity.compoundAtomIndex << " "
    // 		<< atomCoordinates.second[0] << " "
    // 		<< atomCoordinates.second[1] << " "
    // 		<< atomCoordinates.second[2] << std::endl;
    // 	}
    // }
}

/**
 * Write coordinates to a rst7 file
 */
void Replica::WriteRst7(const std::string& fileName) const {
    SimTK_ASSERT_ALWAYS(false, "Replica::WriteRst7 not implemented yet.");

    // FILE *File = fopen(fileName.c_str(), "w+");

    // int Natoms = 0;
    // for(auto& topology : atomsLocations){
    // 	Natoms += topology.size();
    // }

    // fprintf(File, "TITLE: Created by Robosample with %d atoms\n", Natoms);
    // fprintf(File, "%6d\n", Natoms);

    // int atomCnt = -1;
    // for(auto& topology : atomsLocations){
    // 	for(auto& atomCoordinates : topology){
    // 		++atomCnt;
    // 		fprintf(File, "%12.7f%12.7f%12.7f",
    // 			atomCoordinates.second[0] * 10.0,
    // 			atomCoordinates.second[1] * 10.0,
    // 			atomCoordinates.second[2] * 10.0);

    // 		if(atomCnt % 2 == 1){
    // 			fprintf(File, "\n");
    // 		}

    // 	}
    // }

    // if(atomCnt % 2 == 0){
    // 	fprintf(File, "\n");
    // }

    // fflush(File);
    // fclose(File);
}

/**
 * Print coordinates in Amber rst7 format
 */
void Replica::PrintRst7() const {
    SimTK_ASSERT_ALWAYS(false, "Replica::PrintRst7 not implemented yet.");

    // int Natoms = 0;
    // for(auto& topology : atomsLocations){
    // 	Natoms += topology.size();
    // }

    // printf("TITLE: Created by Robosample with %d atoms\n", Natoms);
    // printf("%6d\n", Natoms);

    // int atomCnt = -1;
    // for(auto& topology : atomsLocations){
    // 	for(auto& atomCoordinates : topology){
    // 		++atomCnt;
    // 		printf("%12.7f%12.7f%12.7f",
    // 			atomCoordinates.second[0],
    // 			atomCoordinates.second[1],
    // 			atomCoordinates.second[2]);

    // 		if(atomCnt % 2 == 1){
    // 			printf("\n");
    // 		}

    // 	}
    // }

    // if(atomCnt % 2 == 0){
    // 	printf("\n");
    // }
}

/*!
 * <!-- zmatrixbat_ BAT JAcobian -->
 */
auto Replica::calcInternalBATJacobianLog() -> SimTK::Real {
    // Get log of the Cartesian->BAT Jacobian
    SimTK::Real logJacBAT = 0.0;

    for (size_t zCnt = 0; zCnt < zMatrixBAT.size(); zCnt++) {
        // Get bond term
        SimTK::Real currBond = zMatrixBAT[zCnt][0];

        if (currBond != SimTK::NaN) {
            logJacBAT += 4.0 * std::log(currBond);
        }

        // Get the angle term
        SimTK::Real currAngle = zMatrixBAT[zCnt][1];

        if (currAngle != SimTK::NaN) {
            logJacBAT += 2.0 * std::log(std::sin(currAngle));
        }
    }

    return logJacBAT;
}
