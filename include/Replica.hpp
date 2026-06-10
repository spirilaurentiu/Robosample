#pragma once

#include <vector>

#include "Simbody.h"
#include "Topology.hpp"
#include "TopologyElements.hpp"

class Replica {
    public:
    Replica(int index,
            std::vector<RoboAtom>& atoms_,
            std::vector<int>& roots_,
            Span<Topology> topologies_,
            std::vector<std::vector<int>>& zMatrixTable_)
        : atoms(atoms_)
        , roots(roots_)
        , topologies(topologies_)
        // , internCoords(internCoords_)
        , zMatrixTable(zMatrixTable_)
        , zMatrixBAT() {
    }

    const std::vector<SimTK::Compound::AtomTargetLocations>& getAtomsLocationsInGround() const;
    const std::vector<SimTK::Compound::AtomTargetLocations>& get_WORK_AtomsLocationsInGround() const;

    // Reserve memory and set values
    void setAtomsLocationsInGround(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);

    // Reserve memory and set values
    void
    set_WORK_AtomsLocationsInGround(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);

    // Transfers work coordinates into regular coordinates
    void updAtomsLocationsInGround_FromWORK();

    void
    upd_WORK_AtomsLocationsInGround(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);


    SimTK::Real get_WORK_Jacobian() const;
    SimTK::Real& upd_WORK_Jacobian();

    void set_WORK_Jacobian(SimTK::Real);

    void setPotentialEnergy_FromWORK();

    // Load atomLocations coordinates into the front world
    void restoreCoordinates();

    // Stores coordinates from front world into atomsLocations
    void storeCoordinates();

    SimTK::Real getPotentialEnergy() const;
    void setPotentialEnergy(SimTK::Real somePotential);

    SimTK::Real getReferencePotentialEnergy() const;
    void setReferencePotentialEnergy(SimTK::Real somePotential);

    SimTK::Real getWORK() const;
    void setWORK(SimTK::Real workArg);
    SimTK::Real& updWORK();

    SimTK::Real get_WORK_PotentialEnergy_New() const;
    void set_WORK_PotentialEnergy_New(SimTK::Real somePotential);

    SimTK::Real get_WORK_ReferencePotentialEnergy_New() const;
    void set_WORK_ReferencePotentialEnergy_New(SimTK::Real somePotential);

    // void set_WORK_LastPotentialEnergy(SimTK::Real wpArg);

    SimTK::Real getFixman() const;
    void setFixman(SimTK::Real somePotential);

    SimTK::Real get_WORK_Fixman() const;
    void set_WORK_Fixman(SimTK::Real somePotential);

    void Print() const;
    void PrintCoordinates() const;
    void Print_WORK_Coordinates() const;

    const std::vector<SimTK::Real>& getX() const {
        return x;
    }
    const std::vector<SimTK::Real>& getY() const {
        return y;
    }
    const std::vector<SimTK::Real>& getZ() const {
        return z;
    }

    void PrintRst7() const;
    void WriteRst7(std::string FN) const;

    /** @name Z Matrix and BAT functions
     */

    /**@{**/

    //////////////////////////////////
    /////      Z Matrix BAT      /////
    //////////////////////////////////

    // Function to find and return the value for a given AtomIndex
    SimTK::Vec3 findAtomTarget(const SimTK::Compound::AtomTargetLocations& atomTargets,
                               SimTK::Compound::AtomIndex searchIndex) const;

    void setZMatrixTable(const std::vector<std::vector<int>>& newZMatrixTable);

    // zmatrixbat_ Setter for a specific entry
    void setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value);

    // zmatrixbat_ Function to get a given row
    const std::vector<SimTK::Real>& getZMatrixBATRow(size_t rowIndex) const;

    // zmatrixbat_ Function to get a given row
    std::vector<SimTK::Real>& updZMatrixBATRow(size_t rowIndex);

    // Allocate Z Matrix BAT
    void reallocZMatrixBAT();

    // zmatrixbat_
    void calcZMatrixBAT(
        const std::vector<std::vector<std::pair<RoboAtom*, SimTK::Vec3>>>& otherWorldsAtomsLocations);

    void calcZMatrixBAT_WORK();

    // zmatrixbat_ Function to get the value for a given row and column in zMatrixBAT
    SimTK::Real getZMatrixBATValue(size_t rowIndex, size_t colIndex) const;

    // zmatrixbat_ Function to print the zMatrixBAT
    void PrintZMatrixBAT() const;

    // zmatrixbat_  Function to add a new row to the zMatrixBAT
    void addZMatrixBATRow(const std::vector<SimTK::Real>& newRow);

    /**
     * @brief zmatrixbat_ Get log of the Cartesian->BAT Jacobian
     * @param
     */
    SimTK::Real calcInternalBATJacobianLog();

    // Incrementer function for nofSamples
    void incrementWorldsNofSamples();
    void incrementWorldsNofSamples(int howMany);
    void incrementNofSamples() {
        nofSamples++;
    }
    void incrementNofSamples(int howMany) {
        nofSamples += howMany;
    }

    std::vector<std::vector<SimTK::Real>>& getZMatrixBATPointer() {
        return (zMatrixBAT);
    }

    //////////////////////////////////
    /////      Z Matrix BAT      /////
    //////////////////////////////////

    /**@}**/

    private:
    int myIndex = 0;

    // Replica configurations
    std::vector<SimTK::Compound::AtomTargetLocations> atomsLocations;
    std::vector<SimTK::Compound::AtomTargetLocations> WORK_atomsLocations;

    // Replica potential energy
    SimTK::Real potential;
    SimTK::Real referencePotential;

    SimTK::Real WORK_potential;
    SimTK::Real referenceWORK_potential;

    SimTK::Real WORK;                       // TODO: turn into a vector for worlds
    SimTK::Real workJacobiansContributions; // TODO: turn into a vector for worlds

    SimTK::Real FixmanPotential;      // TODO: turn into a vector for worlds
    SimTK::Real WORK_FixmanPotential; // TODO: turn into a vector for worlds

    std::vector<SimTK::Real> x, y, z;

    //////////////////////////////////
    /////      Z Matrix BAT      /////
    //////////////////////////////////
    std::vector<RoboAtom>& atoms;

    std::vector<int>& roots;
    Span<Topology> topologies;
    // InternalCoordinates& internCoords;

    std::vector<std::vector<int>>& zMatrixTable;
    std::vector<std::vector<SimTK::Real>> zMatrixBAT;

    // std::vector<SimTK::QIndex> QIxs;

    // BAT
    int allWorldsNofSamples = 0;
    int nofSamples = 0;
};
