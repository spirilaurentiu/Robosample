#!/bin/bash

cd ../molmodel_legacy/simbody
mkdir build-debug
cd build-debug

cmake ..
make -j12
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake ..
make -j12
sudo make install

#<<<<<<< HEAD
#cd ../../openmm
#mkdir build-debug
#cd build-debug
#cmake ..
#make -j12
#sudo make install
#=======
#cd ../../openmm
#mkdir build-debug
#cd build-debug
#cmake ..
#make -j4
#sudo make install
#>>>>>>> 2d76ac7c3a99b1aed273e5f79b4f6d219034bda9

cd ../../build-debug/


cmake ..
make -j12
