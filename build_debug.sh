#!/bin/bash

# used to have
# cmake -Wno-dev -DCMAKE_BUILD_TYPE=Debug ..

# remove openmm previous compilation results
sudo rm /usr/local/openmm -rf
sudo find /usr/local/lib -iname '*openmm*' -delete
rm openmm/build-release -rf
rm openmm/build-debug -rf

# compile openmm
cd openmm
mkdir build-debug
cd build-debug/
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2))
sudo make install
cd ../../

# delete previous compilation from project
cd Molmodel/Simbody01
rm build-debug -rf
rm build-release -rf

cd ../
rm build-debug -rf
rm build-release -rf

cd ../
rm build-debug -rf
rm build-release -rf

# delete previous compilation from /usr/
sudo rm /usr/local/lib/cmake/simbody -rf
sudo rm /usr/local/lib/simbody -rf
sudo rm /usr/local/lib/pkgconfig/simbody.pc
sudo find /usr/local/lib -name 'libSimTK*' -delete

# compile project
cd Molmodel/Simbody01
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2))
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2))
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2))
sudo /sbin/ldconfig

# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs
