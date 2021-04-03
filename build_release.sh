#!/bin/bash

# notes for the future
# make -j$((`nproc`*2)+1) -> c++: fatal error: Killed signal terminated program cc1plus (i.e too much memory consumed)

# remove openmm previous openmm compilation
sudo rm /usr/local/openmm -rf
sudo find /usr/local/lib -iname '*openmm*' -delete
rm openmm/build-debug -rf
rm openmm/build-release -rf

# delete previous simbody compilation from /usr/
sudo rm /usr/local/lib/cmake/simbody -rf
sudo rm /usr/local/lib/simbody -rf
sudo rm /usr/local/lib/pkgconfig/simbody.pc
sudo find /usr/local/lib -name 'libSimTK*' -delete

# compile openmm
cd openmm
mkdir build-release
cd build-release/
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install
cd ../../

# delete previous compilation from project (Simbody)
cd Molmodel/Simbody01
rm build-debug -rf
rm build-release -rf

# remove Molmodel
cd ../
rm build-debug -rf
rm build-release -rf

# remove Robosample
cd ../
rm build-debug -rf
rm build-release -rf

# compile project (Simbody)
cd Molmodel/Simbody01
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install

# compile Molmodel
cd ../../
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install

# ensure this is correctly installed
sudo mkdir -p /usr/local/lib/plugins/ && sudo cp -f libOpenMMPlugin.so $_
rm /usr/local/lib/plugins/libOpenMMPlugin_d.so

# compile Robosample
cd ../../
mkdir build-release
cd build-release
mkdir pgo
cmake -DCMAKE_BUILD_TYPE=Release -DROBO_PGO=Generate -DROBO_PGO_PATH="$(pwd)/pgo/" ..
make -j$((`nproc`*2))
sudo /sbin/ldconfig

# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# run PGO on Robosample and recompile
./tests/Robosample inp.dock.test
cmake -DCMAKE_BUILD_TYPE=Release -DROBO_PGO=Use -DROBO_PGO_PATH="$(pwd)/pgo/" ..
make -j$((`nproc`*2))
sudo /sbin/ldconfig

# add tests
cp ../tests/test-memory.sh test-memory.sh
