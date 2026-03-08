#include "bgeneral.hpp"

#pragma once

struct EnergySnapshot {
    SimTK::Real potential = SimTK::NaN;
    SimTK::Real kinetic = SimTK::NaN;
    SimTK::Real fixman = SimTK::NaN;
    SimTK::Real logSineSqrGamma2 = SimTK::NaN;
    SimTK::Real total = SimTK::NaN;

    static constexpr SimTK::Real EXPLOSION_THRESHOLD = 10.0; 
    static constexpr SimTK::Real GEOMETRIC_LIMIT = 1e3;
	
    bool finite() const;
    bool checkTotal(SimTK::Real RT) const;
    bool checkKinetic() const;
    bool termExplosion(const EnergySnapshot& ref) const;
    bool geometryExplosion() const;
    bool hamiltonianDrift(const EnergySnapshot& ref, SimTK::Real beta, std::size_t ndofs) const;
    bool validate(const EnergySnapshot& ref, SimTK::Real RT, std::size_t ndofs) const;
};
