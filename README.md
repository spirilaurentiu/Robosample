# Prerequisites
Install prerequisites.
```
sudo apt-get update
sudo apt-get install cmake gfortran libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig ocl-icd-opencl-dev fftw2
```

Make sure you added the following to LD_LIBRARY_PATH:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib/plugins
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/plugins/libOpenMMPlugin.so

export CUDA_INC_DIR=/usr/local/cuda
export CUDA_ROOT=/usr/local/cuda
```

Also make sure that:
1. `cmake --version` is greater than 3.1.
1. `gcc` version can compile C++11.

#  Cloning the project
```
git clone --recurse-submodules https://github.com/spirilaurentiu/Robosample.git
cd Robosample
```

# Building OpenMM
```
cd openmm
mkdir build_debug
cd build_debug/
cmake ..
make -j4
sudo make install
cd ../../
```

# Building Robosample
```
mkdir build_debug
cd build_debug
bash ../cmake_regenerateAll.sh - build_debug will remain empty
cmake ..
make -j4
sudo /sbin/ldconfig
```

# Running the tests
```
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs
```
To run the test called `Robosample`:  `./tests/Robosample inp`
