#!/bin/bash

# remove previous compilation results
sudo rm /usr/local/openmm -rf
sudo find /usr/local/lib -iname '*openmm*' -delete
rm openmm/build-release -rf

# compile
cd openmm
mkdir build-release
cd build-release/
cmake -Wno-dev ..
make -j$(nproc)
sudo make install
