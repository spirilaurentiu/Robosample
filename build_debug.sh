#!/bin/bash

# # remove openmm previous openmm compilation
# sudo rm /usr/local/openmm -rf
# sudo find /usr/local/lib -iname '*openmm*' -delete
# rm openmm/build-debug -rf
# rm openmm/build-release -rf

# # delete previous simbody compilation from /usr/
# sudo rm /usr/local/lib/cmake/simbody -rf
# sudo rm /usr/local/lib/simbody -rf
# sudo rm /usr/local/lib/pkgconfig/simbody.pc
# sudo find /usr/local/lib -name 'libSimTK*' -delete

<<<<<<< Updated upstream
# compile openmm
cd openmm
git checkout master && git pull
mkdir build-debug
cd build-debug/
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2)) 2> build-openmm.txt
sudo make install
cd ../../
=======
# # compile openmm
# cd openmm
# mkdir build-debug
# cd build-debug/
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# make -j$((`nproc`*2)) 2> build-openmm.txt
# sudo make install
# cd ../../
>>>>>>> Stashed changes

# # delete previous compilation from project
# cd Molmodel/Simbody01
# rm build-debug -rf
# rm build-release -rf

# cd ../
# rm build-debug -rf
# rm build-release -rf

# cd ../
# rm build-debug -rf
# rm build-release -rf

<<<<<<< Updated upstream
# compile project
cd Molmodel/Simbody01
git checkout master && git pull
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2)) 2> build-simbody.txt
sudo make install

cd ../../
git checkout master && git pull
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2)) 2> build-molmodel.txt
sudo make install
=======
# # compile project
# cd Molmodel/Simbody01
# mkdir build-debug
# cd build-debug
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# make -j$((`nproc`*2)) 2> build-simbody.txt
# sudo make install

# cd ../../
# mkdir build-debug
# cd build-debug
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# make -j$((`nproc`*2)) 2> build-molmodel.txt
# sudo make install
>>>>>>> Stashed changes

# # ensure this is correctly installed
# sudo mkdir -p /usr/local/lib/plugins/ && sudo cp -f libOpenMMPlugin_d.so $_
# rm /usr/local/lib/plugins/libOpenMMPlugin.so

<<<<<<< Updated upstream
cd ../../
git checkout master && git pull
=======
# cd ../../
rm build-debug -rf
>>>>>>> Stashed changes
mkdir build-debug
cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$((`nproc`*2)) 2> build-robosample.txt
sudo /sbin/ldconfig

# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# add tests
cp ../tests/test-memory.sh test-memory.sh

echo "Build done. See build logs in build-debug/"
