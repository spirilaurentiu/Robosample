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

/**
 * @brief Boundary unit-safety helpers built on the nholthaus compile-time units
 *        library. Typed quantities live only at I/O edges (prmtop parse, OpenMM
 *        transfer, DCD write, timestep config); the dynamics core stays raw
 *        double in the consistent MD system (nm, dalton, ps, kJ/mol).
 *
 * @note The engine's inter-world coordinate currency is nm (INV-3). These
 *       conversions exist for the exact spots where a foreign unit system enters
 *       (Angstrom, kcal, fs); mass (dalton) and energy-per-mole never leave the
 *       MD system and so are not wrapped.
 */
namespace robo::mdunits {

// boundary unit aliases (use these on function signatures that cross an edge)
/** @brief Length in nanometers (the engine's internal length unit, INV-3). */
using Nanometer = units::length::nanometer_t;
/** @brief Length in Angstrom (prmtop / DCD-file length unit). */
using Angstrom = units::length::angstrom_t;
/** @brief Time in picoseconds (the engine's internal time unit). */
using Picosecond = units::time::picosecond_t;
/** @brief Time in femtoseconds (a common timestep-input unit). */
using Femtosecond = units::time::femtosecond_t;
/** @brief Energy in kilojoules; per-mole is implicit in the MD system. */
using Kilojoule = units::energy::kilojoule_t;
/** @brief Energy in kilocalories (Amber/CHARMM force-field input unit). */
using Kilocalorie = units::energy::kilocalorie_t;

// ---- boundary conversions: typed in, raw double out (for the core) ----------
/** @brief Convert an Angstrom length to nm (engine unit). */
inline auto angstromToNm(double valAngstrom) -> double {
    return Nanometer(Angstrom(valAngstrom)).value();
}
/** @brief Convert an nm length to Angstrom (e.g. the DCD writer's x10). */
inline auto nmToAngstrom(double valNanometer) -> double {
    return Angstrom(Nanometer(valNanometer)).value(); // the DCD writer's x10
}
/** @brief Convert a kcal/mol energy to kJ/mol (engine unit). */
inline auto kcalToKj(double valKilocalorie) -> double {
    return Kilojoule(Kilocalorie(valKilocalorie)).value();
}
/** @brief Convert a kJ/mol energy to kcal/mol. */
inline auto kjToKcal(double valKilojoule) -> double {
    return Kilocalorie(Kilojoule(valKilojoule)).value();
}
/** @brief Convert a femtosecond time to ps (engine unit). */
inline auto fsToPs(double valFemtosecond) -> double {
    return Picosecond(Femtosecond(valFemtosecond)).value();
}

// Mass (dalton == g/mol numerically) and energy-per-mole stay raw double: the
// prmtop already supplies them in the MD system, and they never cross into a
// different unit system, so a strong type would add friction with no payoff.

} // namespace robo::mdunits