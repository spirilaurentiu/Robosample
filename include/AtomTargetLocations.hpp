#pragma once

#include "Compound.h"
#include "TopologyElements.hpp"
#include "bgeneral.hpp"

class AtomLocations {
    public:
    AtomLocations() = default;

    void initialize(const SystemTopology& systemTopology);

    [[nodiscard]] auto getLocation(int topoIx, SimTK::Compound::AtomIndex cIAx) const -> const SimTK::Vec3& {
        return locations[atomOffsets[topoIx] + cIAx];
    }

    [[nodiscard]] auto updLocation(int topoIx, SimTK::Compound::AtomIndex cIAx) -> SimTK::Vec3& {
        return locations[atomOffsets[topoIx] + cIAx];
    }

    [[nodiscard]] auto getLocationsForTopology(int topoIx) const -> Span<const SimTK::Vec3> {
        return {locations.data() + atomOffsets[topoIx], atomCounts[topoIx]};
    }

    [[nodiscard]] auto updLocationsForTopology(int topoIx) -> Span<SimTK::Vec3> {
        return {locations.data() + atomOffsets[topoIx], atomCounts[topoIx]};
    }

    private:
    std::vector<SimTK::Vec3> locations;

    std::vector<std::size_t> atomOffsets;
    std::vector<std::size_t> atomCounts;

    std::vector<std::size_t> bondAtom1CAIx;
    std::vector<std::size_t> bondAtom2CAIx;
    std::vector<int> bondCIx;
    std::vector<std::size_t> bondOffsets;
    std::vector<std::size_t> bondCounts;

    std::vector<std::size_t> angleOffsets;
    std::vector<std::size_t> angleCounts;

    std::vector<std::size_t> dihedralOffsets;
    std::vector<std::size_t> dihedralCounts;
};