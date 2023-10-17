#include "Context.hpp"

constexpr auto DECIMAL_PLACES = 5;
bool areEqual(SimTK::Real a, SimTK::Real b) {
    SimTK::Real epsilon = std::pow(10, -DECIMAL_PLACES);
    return std::abs(a - b) < epsilon;
}

int main(int argc, char **argv){
    const std::string inpcrdfile = "2but/ligand.inpcrd";
    const std::string prmtopfile = "2but/ligand.prmtop";

    std::vector<std::string> atomName = { "O1", "C1", "C2", "C3", "C4", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10" };
    std::vector<std::string> atomType = { "oh", "c3", "c3", "c3", "c3", "h1", "hc", "hc", "hc", "hc", "hc", "hc", "hc", "hc", "ho" };
    std::vector<SimTK::Real> atomX = { -0.518, -0.481, 0.735, -1.776, 2.041, -0.395, 0.773, 0.648, -1.913, -1.791, -2.638, 2.173, 2.076, 2.888, -1.287 };
    std::vector<SimTK::Real> atomY = { 1.356, 0.02, -0.7, -0.683, 0.006, 0.067, -1.73, -0.746, -0.686, -1.716, -0.155, 0.089, 1.012, -0.558, 1.796 };
    std::vector<SimTK::Real> atomZ = { -0.152, 0.343, -0.244, -0.041, 0.094, 1.435, 0.128, -1.336, -1.128, 0.319, 0.381, 1.177, -0.336, -0.308, 0.248 };
    std::vector<SimTK::Real> atomMass = { 16., 12.01, 12.01, 12.01, 12.01, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008 };
    std::vector<SimTK::Real> atomCharge = { -10.9479578, 2.51649963, -1.35573912, -2.38894353, -1.67827383, 0.55942461, 0.87831486, 0.87831486, 0.765937936, 0.765937936, 0.765937936, 0.66875841, 0.66875841, 0.66875841, 7.23425310 };
    std::vector<SimTK::Real> atomRadii = { 1.500, 1.700, 1.700, 1.700, 1.700, 1.300, 1.300, 1.300, 1.300, 1.300, 1.300, 1.300, 1.300, 1.300, 0.800 };
    std::vector<SimTK::Real> atomRVdw = { 1.7210000632, 1.9080000762, 1.9080000762, 1.9080000762, 1.9080000762, 1.3870000413, 1.487000048, 1.487000048, 1.487000048, 1.487000048, 1.487000048, 1.487000048, 1.487000048, 1.487000048, 0.5 };
    std::vector<SimTK::Real> atomEpsilon = { 0.2104000002, 0.1093999999, 0.1093999999, 0.1093999999, 0.1093999999, 0.0157000001, 0.0157, 0.0157, 0.0157, 0.0157, 0.0157, 0.0157, 0.0157, 0.0157, 0.0 };

    std::vector<std::size_t> bondInd1 = { 0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 0, 1, 1, 2 };
    std::vector<std::size_t> bondInd2 = { 14, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 2, 3, 4 };
    std::vector<SimTK::Real> bondForce = { 369.6, 335.9, 337.3, 337.3, 337.3, 337.3, 337.3, 337.3, 337.3, 337.3, 314.1, 303.1, 303.1, 303.1 };
    std::vector<SimTK::Real> bondEq = { 0.974, 1.093, 1.092, 1.092, 1.092, 1.092, 1.092, 1.092, 1.092, 1.092, 1.426, 1.535, 1.535, 1.535 };

    std::vector<std::size_t> angleInd1 = { 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 6, 8, 8, 9, 11, 11, 12, 0, 0, 1, 2 };
    std::vector<std::size_t> angleInd2 = { 1, 0, 2, 2, 3, 3, 3, 1, 4, 4, 4, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 2, 1 };
    std::vector<std::size_t> angleInd3 = { 5, 14, 6, 7, 8, 9, 10, 5, 11, 12, 13, 5, 6, 7, 7, 9, 10, 10, 12, 13, 13, 2, 3, 4, 3 };
    std::vector<SimTK::Real> angleForce = { 50.97, 47.09, 46.37, 46.37, 46.37, 46.37, 46.37, 46.36, 46.37, 46.37, 46.37, 46.36, 46.37, 46.37, 39.43, 39.43, 39.43, 39.43, 39.43, 39.43, 39.43, 67.72, 67.72, 63.21, 63.21 };
    std::vector<SimTK::Real> angleEq = { 1.917769, 1.887749, 1.920736, 1.920736, 1.920736, 1.920736, 1.920736, 1.921085, 1.920736, 1.920736, 1.920736, 1.921085, 1.920736, 1.920736, 1.891065, 1.891065, 1.891065, 1.891065, 1.891065, 1.891065, 1.891065, 1.909915, 1.909915, 1.930859, 1.930859 };

    std::vector<std::size_t> dihedralInd1 = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 14, 14, 2, 2, 2, 14, 14, 3, 3, 4, 14, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 0, 3, 3, 3 };
    std::vector<std::size_t> dihedralInd2 = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0, 1, 1, 1, 0, 0, 1, 1, 2, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1 };
    std::vector<std::size_t> dihedralInd3 = { 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 1, 1, 3, 3, 3, 1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2 };
    std::vector<std::size_t> dihedralInd4 = { 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 12, 13, 2, 2, 8, 9, 10, 3, 3, 6, 7, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 11, 12, 13, 4, 4, 4, 4 };
    std::vector<SimTK::Real> dihedralForce = { 0.25, 0, 0.25, 0, 0.25, 0, 0.25, 0, 0.25, 0, 0.16, 0.16, 0.16, 0.25, 0.16, 0.16, 0.16, 0.16, 0.25, 0.16, 0.16, 0.16, 0.155556, 0.166667, 0.155556, 0.155556, 0.155556, 0.155556, 0.155556, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.155556, 0.2, 0.25, 0.18 };
    std::vector<SimTK::Real> dihedralPeriod = { 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 3, 3, 3, 1, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 2, 3 };
    std::vector<SimTK::Real> dihedralPhase = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.141594, 3.141594, 0 };

    std::vector<std::vector<bool>> nonBondedMatrix = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1 },
        { 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1 },
        { 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1 },
        { 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1 },
        { 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1 },
        { 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0 }
    };

    readAmberInput reader;
    reader.readAmberFiles(inpcrdfile, prmtopfile);

    // general info
    if (reader.getNumberAtoms() != 15) {
        std::cout << "incorrect number of atoms read" << std::endl;
        return 1;
    }

    if (reader.getNumberBonds() != 14) {
        std::cout << "incorrect number of bonds read" << std::endl;
        return 1;
    }

    if (reader.getNumberAngles() != 25) {
        std::cout << "incorrect number of angles read" << std::endl;
        return 1;
    }

    if (reader.getNumberDihedrals() != 39) {
        std::cout << "incorrect number of dihedrals read" << std::endl;
        return 1;
    }

    // test atoms
    for (std::size_t i = 0; i < reader.getNumberAtoms(); i++) {
        if (reader.getAtomsName(i) != atomName[i]) {
            std::cout << "Atom name mismatch. Expected '" << atomName[i] << "', got '" << reader.getAtomsName(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getAtomsType(i) != atomType[i]) {
            std::cout << "Atom type mismatch. Expected '" << atomType[i] << "', got '" << reader.getAtomsType(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsXcoord(i), atomX[i])) {
            std::cout << "Atom X coord mismatch. Expected '" << atomX[i] << "', got '" << reader.getAtomsXcoord(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsYcoord(i), atomY[i])) {
            std::cout << "Atom Y coord mismatch. Expected '" << atomY[i] << "', got '" << reader.getAtomsYcoord(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsZcoord(i), atomZ[i])) {
            std::cout << "Atom Z coord mismatch. Expected '" << atomZ[i] << "', got '" << reader.getAtomsZcoord(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsMass(i), atomMass[i])) {
            std::cout << "Atom mass mismatch. Expected '" << atomMass[i] << "', got '" << reader.getAtomsMass(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsCharge(i), atomCharge[i])) {
            std::cout << "Atom charge mismatch. Expected '" << atomCharge[i] << "', got '" << reader.getAtomsCharge(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsRadii(i), atomRadii[i])) {
            std::cout << "Atom radius mismatch. Expected '" << atomRadii[i] << "', got '" << reader.getAtomsRadii(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsRVdW(i), atomRVdw[i])) {
            std::cout << "Atom VdW radius mismatch. Expected '" << atomRVdw[i] << "', got '" << reader.getAtomsRVdW(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAtomsEpsilon(i), atomEpsilon[i])) {
            std::cout << "Atom epsilon mismatch. Expected '" << atomEpsilon[i] << "', got '" << reader.getAtomsEpsilon(i) << "'" << std::endl;
            return 1;
        }
    }

    // test bonds
    for (std::size_t i = 0; i < reader.getNumberBonds(); i++) {
        if (reader.getBondsAtomsIndex1(i) != bondInd1[i]) {
            std::cout << "Bond atom 1 index mismatch. Expected '" << bondInd1[i] << "', got '" << reader.getBondsAtomsIndex1(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getBondsAtomsIndex2(i) != bondInd2[i]) {
            std::cout << "Bond atom 2 index mismatch. Expected '" << bondInd2[i] << "', got '" << reader.getBondsAtomsIndex2(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getBondsForceK(i), bondForce[i])) {
            std::cout << "Bond force constant mismatch. Expected '" << bondForce[i] << "', got '" << reader.getBondsForceK(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getBondsEqval(i), bondEq[i])) {
            std::cout << "Bond force equilibrium mismatch. Expected '" << bondEq[i] << "', got '" << reader.getBondsEqval(i) << "'" << std::endl;
            return 1;
        }
    }

    // test angles
    for (std::size_t i = 0; i < reader.getNumberAngles(); i++) {
        if (reader.getAnglesAtomsIndex1(i) != angleInd1[i]) {
            std::cout << "Angle atom 1 index mismatch. Expected '" << angleInd1[i] << "', got '" << reader.getAnglesAtomsIndex1(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getAnglesAtomsIndex2(i) != angleInd2[i]) {
            std::cout << "Angle atom 2 index mismatch. Expected '" << angleInd2[i] << "', got '" << reader.getAnglesAtomsIndex2(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getAnglesAtomsIndex3(i) != angleInd3[i]) {
            std::cout << "Angle atom 3 index mismatch. Expected '" << angleInd3[i] << "', got '" << reader.getAnglesAtomsIndex3(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAnglesForceK(i), angleForce[i])) {
            std::cout << "Angle force mismatch. Expected '" << angleForce[i] << "', got '" << reader.getAnglesForceK(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getAnglesEqval(i), angleEq[i])) {
            std::cout << "Angle equilibrium mismatch. Expected '" << angleEq[i] << "', got '" << reader.getAnglesEqval(i) << "'" << std::endl;
            return 1;
        }
    }

    // test dihedrals
    for (std::size_t i = 0; i < reader.getNumberDihedrals(); i++) {
        if (reader.getDihedralsAtomsIndex1(i) != dihedralInd1[i]) {
            std::cout << "Dihedral atom 1 index mismatch. Expected '" << dihedralInd1[i] << "', got '" << reader.getDihedralsAtomsIndex1(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getDihedralsAtomsIndex2(i) != dihedralInd2[i]) {
            std::cout << "Dihedral atom 2 index mismatch. Expected '" << dihedralInd2[i] << "', got '" << reader.getDihedralsAtomsIndex2(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getDihedralsAtomsIndex3(i) != dihedralInd3[i]) {
            std::cout << "Dihedral atom 3 index mismatch. Expected '" << dihedralInd3[i] << "', got '" << reader.getDihedralsAtomsIndex3(i) << "'" << std::endl;
            return 1;
        }

        if (reader.getDihedralsAtomsIndex4(i) != dihedralInd4[i]) {
            std::cout << "Dihedral atom 4 index mismatch. Expected '" << dihedralInd4[i] << "', got '" << reader.getDihedralsAtomsIndex4(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getDihedralsForceK(i), dihedralForce[i])) {
            std::cout << "Dihedral force mismatch. Expected '" << dihedralForce[i] << "', got '" << reader.getDihedralsForceK(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getDihedralsPeriod(i), dihedralPeriod[i])) {
            std::cout << "Dihedral atom 1 index mismatch. Expected '" << dihedralPeriod[i] << "', got '" << reader.getDihedralsPeriod(i) << "'" << std::endl;
            return 1;
        }

        if (!areEqual(reader.getDihedralsPhase(i), dihedralPhase[i])) {
            std::cout << "Dihedral atom 1 index mismatch. Expected '" << dihedralPhase[i] << "', got '" << reader.getDihedralsPhase(i) << "'" << std::endl;
            return 1;
        }
    }

    // test non-bonded matrix
    for(std::size_t at0 = 0; at0 < reader.getNumberAtoms(); at0++)
    {
        for(std::size_t at1 = 0; at1 < reader.getNumberAtoms(); at1++) {
            if (nonBondedMatrix[at0][at1] != reader.getNonBondedAtomsMatrix(at0, at1)) {
                std::cout << "Nonbonded atom matrix mismatch between atoms " << at0 << " and " << at1 << std::endl;
            }
        }
    }

    return 0;
}
