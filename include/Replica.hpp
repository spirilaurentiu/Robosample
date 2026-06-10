#pragma once

#include <vector>

#include "TopologyElements.hpp"

class Replica {
    public:
    Replica(const SystemTopology& systemTopology, int replicaIndex) {
        myIndex = replicaIndex;
    }

    [[nodiscard]] auto getAtomsLocationsInGround() const
        -> const std::vector<SimTK::Compound::AtomTargetLocations>& {
        return atomsLocations;
    }

    [[nodiscard]] auto get_WORK_AtomsLocationsInGround() const
        -> const std::vector<SimTK::Compound::AtomTargetLocations>& {
        return WORK_atomsLocations;
    }

    // Reserve memory and set values
    void setAtomsLocationsInGround(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);

    // Reserve memory and set values
    void
    set_WORK_AtomsLocationsInGround(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) {
        WORK_atomsLocations = atomTargets;
    }

    // Transfers work coordinates into regular coordinates
    void updAtomsLocationsInGround_FromWORK();

    void
    upd_WORK_AtomsLocationsInGround(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) {
        WORK_atomsLocations = atomTargets;
    }

    [[nodiscard]] auto get_WORK_PotentialEnergy_New() const -> SimTK::Real {
        return this->WORK_potential;
    }

    void set_WORK_PotentialEnergy_New(SimTK::Real somePotential) {
        WORK_potential = somePotential;
    }

    [[nodiscard]] auto get_WORK_ReferencePotentialEnergy_New() const -> SimTK::Real {
        return this->referenceWORK_potential;
    }

    void set_WORK_ReferencePotentialEnergy_New(SimTK::Real somePotential) {
        referenceWORK_potential = somePotential;
    }

    [[nodiscard]] auto get_WORK_Jacobian() const -> SimTK::Real {
        return this->workJacobiansContributions;
    }

    auto upd_WORK_Jacobian() -> SimTK::Real& {
        return this->workJacobiansContributions;
    }

    void set_WORK_Jacobian(SimTK::Real inpJac) {
        this->workJacobiansContributions = inpJac;
    }

    void setPotentialEnergy_FromWORK() {
        this->potential = this->WORK_potential;
    }

    // Load atomLocations coordinates into the front world
    void restoreCoordinates();

    // Stores coordinates from front world into atomsLocations
    void storeCoordinates();

    [[nodiscard]] auto getPotentialEnergy() const -> SimTK::Real {
        return potential;
    }

    void setPotentialEnergy(SimTK::Real somePotential) {
        potential = somePotential;
    }

    [[nodiscard]] auto getReferencePotentialEnergy() const -> SimTK::Real {
        return referencePotential;
    }

    void setReferencePotentialEnergy(SimTK::Real somePotential) {
        referencePotential = somePotential;
    }

    [[nodiscard]] auto getWORK() const -> SimTK::Real {
        return WORK;
    }

    auto updWORK() -> SimTK::Real& {
        return WORK;
    }

    void setWORK(SimTK::Real workArg) {
        this->WORK = workArg;
    }

    [[nodiscard]] auto getFixman() const -> SimTK::Real {
        return FixmanPotential;
    }

    void setFixman(SimTK::Real somePotential) {
        FixmanPotential = somePotential;
    }

    [[nodiscard]] auto get_WORK_Fixman() const -> SimTK::Real {
        return this->WORK_FixmanPotential;
    }

    void set_WORK_Fixman(SimTK::Real somePotential) {
        this->WORK_FixmanPotential = somePotential;
    }

    void Print() const;
    void PrintCoordinates() const;
    void Print_WORK_Coordinates() const;

    [[nodiscard]] auto getX() const -> const std::vector<SimTK::Real>& {
        return x;
    }
    [[nodiscard]] auto getY() const -> const std::vector<SimTK::Real>& {
        return y;
    }
    [[nodiscard]] auto getZ() const -> const std::vector<SimTK::Real>& {
        return z;
    }

    void PrintRst7() const;
    void WriteRst7(const std::string& fileName) const;

    /**
     * @brief zmatrixbat_ Get log of the Cartesian->BAT Jacobian
     * @param
     */
    auto calcInternalBATJacobianLog() -> SimTK::Real;

    // Incrementer function for nofSamples
    void incrementWorldsNofSamples() {
        ++allWorldsNofSamples;
    }
    void incrementWorldsNofSamples(int howMany) {
        allWorldsNofSamples += howMany;
    }
    void incrementNofSamples() {
        nofSamples++;
    }
    void incrementNofSamples(int howMany) {
        nofSamples += howMany;
    }

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

    // //////////////////////////////////
    // /////      Z Matrix BAT      /////
    // //////////////////////////////////

    std::vector<std::vector<int>>& zMatrixTable;
    std::vector<std::vector<SimTK::Real>> zMatrixBAT;

    int allWorldsNofSamples = 0;
    int nofSamples = 0;
};
