#include "AtomTargetLocations.hpp"

void AtomLocations::initialize(const SystemTopology& systemTopology) {
    // Allocate enough memory for all atom locations
    const auto numAtoms = static_cast<int>(systemTopology.atoms.size());
    locations.resize(numAtoms);
}