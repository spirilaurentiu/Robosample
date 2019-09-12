#!/bin/bash

cd ../Molmodel/Simbody01
mkdir build-debug
cd build-debug
cmake ..
make -j4
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake ..
make -j4
sudo make install

cd ../../build-debug/
cmake ..
make -j4
