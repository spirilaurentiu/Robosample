#!/usr/bin/env bash
set -xeuo pipefail

mkdir -p build
cd build

cmake .. \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DPython3_EXECUTABLE=$PYTHON \
    -DCUDAToolkit_ROOT=/usr/local/cuda-12.8 \
    -DOPENMM_PLATFORM=CUDA \
    -DBUILD_VISUALIZER=OFF

ninja
ninja install

## Now install the pure Python wrapper
#cd ../python
#$PYTHON -m pip install . --no-deps --ignore-installed -vv
