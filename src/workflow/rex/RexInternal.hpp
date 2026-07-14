#pragma once

// Shared REX-internal constant, split out of Context.cpp's former anonymous
// namespace (SPLIT-C4) so it is reachable from both ReplicaExchangeDriver.cpp
// and SwapAcceptance.cpp. Anonymous namespace => internal linkage per TU that
// includes this header (each gets its own copy), matching the original
// file-local constant's visibility.
namespace {
constexpr double kBoltzmann_kJ = 0.0083144626; // kJ/mol/K
} // namespace
