#!/bin/bash'

cd openmm
mkdir build-release
cd build-release/
cmake -Wno-dev ..
make -j$(nproc)
sudo make install
