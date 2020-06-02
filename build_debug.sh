#!/bin/bash

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
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
sudo /sbin/ldconfig

# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs
