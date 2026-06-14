#pragma once
// ============================================================================
//  units_md.hpp -- boundary unit safety for robosample, built on the
//  nholthaus compile-time units library (third_party/units/units.h, MIT).
//
//  POLICY: typed quantities live ONLY at I/O boundaries (prmtop parsing, the
//  OpenMM transfer, DCD writing, timestep configuration) where unit-system
//  mistakes actually happen -- Angstrom vs nm, kcal vs kJ. The dynamics core
//  stays raw `double` in the consistent MD system (nm, dalton, ps, kJ/mol);
//  there is nothing to mix up once past the boundary, so paying for dimensional
//  analysis inside the ABA recursion would be cost with no benefit. nholthaus
//  unit_t is zero runtime overhead, but the ergonomic + compile-time cost is
//  only worth it at the edges. The conversions below are the exact spots that
//  have bitten us (e.g. the DCD nm->Angstrom x10).
//
//  Add the submodule and point the include path at it:
//      git submodule add https://github.com/nholthaus/units third_party/units
//      target_include_directories(robosample PRIVATE third_party/units/include)
// ============================================================================

#include "units.h" // nholthaus

namespace robo::mdunits {

// boundary unit aliases (use these on function signatures that cross an edge)
using Nanometer = units::length::nanometer_t;
using Angstrom = units::length::angstrom_t;
using Picosecond = units::time::picosecond_t;
using Femtosecond = units::time::femtosecond_t;
using Kilojoule = units::energy::kilojoule_t;     // per-mole is implicit in MD
using Kilocalorie = units::energy::kilocalorie_t; // ratio kcal:kJ is system-independent

// ---- boundary conversions: typed in, raw double out (for the core) ----------
inline auto angstromToNm(double valAngstrom) -> double {
    return Nanometer(Angstrom(valAngstrom)).value();
}
inline auto nmToAngstrom(double valNanometer) -> double {
    return Angstrom(Nanometer(valNanometer)).value(); // the DCD writer's x10
}
inline auto kcalToKj(double valKilocalorie) -> double {
    return Kilojoule(Kilocalorie(valKilocalorie)).value();
}
inline auto kjToKcal(double valKilojoule) -> double {
    return Kilocalorie(Kilojoule(valKilojoule)).value();
}
inline auto fsToPs(double valFemtosecond) -> double {
    return Picosecond(Femtosecond(valFemtosecond)).value();
}

// Mass (dalton == g/mol numerically) and energy-per-mole stay raw double: the
// prmtop already supplies them in the MD system, and they never cross into a
// different unit system, so a strong type would add friction with no payoff.

} // namespace robo::mdunits