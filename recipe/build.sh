#!/usr/bin/env bash
set -xeuo pipefail

# Clear for new CMake configs
rm build/cuda-pgo-* -rf
rm profile-data -rf

# Build with PGO
cmake --preset cuda-pgo-train
cmake --build --preset cuda-pgo-train

# Train PGO
# python3 python/robosample/roborun.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6000 0 5 1
python3 python/robosample/roborun.py 2ala examples/2ala.prmtop examples/2ala.rst7 6000 0 100 1

# Use PGO data to build optimized binary
cmake --preset cuda-pgo-use
cmake --build --preset cuda-pgo-use

# Resolve the single pybind11 module path
SO_PATH=$(ls python/robosample/robo_bindings*.so)
BOLT_SO="robo_bindings.bolt.so"

# Store perf data along with PGO data for later analysis
PERF_DATA="profile-data/perf.data"
FDATA="profile-data/perf.fdata"

# Profile
perf record \
    -o "$PERF_DATA" \
    -e cycles:u -j any,u \
    python3 python/robosample/roborun.py \
    2ala examples/2ala.prmtop examples/2ala.rst7 6000 0 100 1

# Convert perf data
perf2bolt "$SO_PATH" \
    -p "$PERF_DATA" \
    -o "$FDATA"

# Run BOLT
llvm-bolt "$SO_PATH" \
    -o "$BOLT_SO" \
    -data="$FDATA" \
    -reorder-blocks=ext-tsp \
    -reorder-functions=hfsort \
    -dyno-stats

# Atomically replace original
mv -f "$BOLT_SO" "$SO_PATH"
