#include "GibbsSweep.hpp"

#include <algorithm>

#include "World.hpp"

namespace robo::gibbs {

bool stepWorldInGround(World& w, std::vector<robo::Vec3>& coords, int numAtoms) {
    w.setAtomsLocationsInGround(coords);
    const bool accepted = w.generateSample();
    const robo::Vec3* p = w.getAtomsLocationsInGround();
    std::copy(p, p + numAtoms, coords.begin());
    return accepted;
}

} // namespace robo::gibbs
