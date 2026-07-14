#pragma once
// ============================================================================
//  TestPhysConstants.hpp -- the shared kB/kT300 literals every statistical
//  test built against a 300 K reference temperature was duplicating locally
//  (TESTS.md section 4). Values are byte-identical to the literals they
//  replace (0.0083144626 kJ/mol/K, the gas constant R in kJ/mol/K, times
//  300 K), so converting a `constexpr double kT300 = ...;` call site to this
//  header changes nothing about the drawn sample stream.
//
//  Scope: this header holds ONLY the temperature-independent gas constant and
//  the one reference-temperature product every duplicate call site used
//  (kT300). Per-test two-temperature sweeps (e.g. T1/T2 ratio checks) use
//  DIFFERENT numeric values across files (TestEnsembleValidation.cpp's
//  T1=280/T2=320 vs TestEquipartition.cpp's T1=250/T2=400) -- those are not
//  literal duplicates of one shared constant, so they stay local to each test
//  file rather than being forced into a single named pair here.
// ============================================================================

namespace rtest::phys {

// Gas constant R, kJ/(mol*K) -- AMBER/OpenMM unit convention used throughout
// the engine (RobotEngine/HmcDriver take RT in kJ/mol).
inline constexpr double kB = 0.0083144626;

// RT at the 300 K reference temperature every converted test draws at.
inline constexpr double kT300 = kB * 300.0;

} // namespace rtest::phys
