# Installation

## Installing Xming
Download `Xming` from https://sourceforge.net/projects/xming/ and install it using the default options.

Open `.bashrc` with `vi ~/.bashrc` and add `export DISPLAY=:0` at the end of the file.

## Prerequisites
Install the dependencies:
```
sudo apt-get update
sudo apt-get install cmake gfortran libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig ocl-icd-opencl-dev fftw2
```

Add the following to `LD_LIBRARY_PATH`:
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

##  Cloning the project
```
git clone --recurse-submodules https://github.com/spirilaurentiu/Robosample.git
cd Robosample
```

## Building OpenMM
Run `build_openmm.sh`. It will ask for password after compilation.
```
bash build_openmm.sh
```

## Building Robosample
Run `build_debug.sh` or `bulild_release.sh`. Password will be required as a lot of files will be deleted from `/usr/`.
```
<<<<<<< HEAD
bash build_debug.sh
=======
mkdir build_debug
cd build_debug
bash ../cmake_regenerateAll.sh
cmake ..
make -j4
sudo /sbin/ldconfig
```
Keep in mind that after running `bash ../cmake_regenerateAll.sh`, `build_debug` will remain empty.

## Set up tests' environment
```
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs
cd ../
>>>>>>> origin/master
```

# Open the project using any IDE (e.g. Visual Studio Code)
Install [Visual Studio Code](https://code.visualstudio.com/) on Windows. Run `code .` in `Robosample`.

# Working on Robosample
After working on Robosample, it must be compiled as `Debug` or `Release` (the same flags as in section [Building Robosample](#building-robosample)). To compile as a different configuration, full recompilation is needed (see [Building Robosample](#building-robosample)).
```
make -j$(nproc)
```

<<<<<<< HEAD
# Running Robosample
`Robosample` is located in `build-debug/tests` or `build-release/tests`.
```
cd build-debug
./tests/Robosample inp
```
To change different parameters (use visualizer, use OpenMM etc) edit `inp` which is located in `build-debug` or `build-release`.
=======
# Running the tests
The tests are located in `Robosample/build_debug/tests`.

From `build_debug` (this is where you should be if you have just compiled the project), type `./tests/Robosample inp` to run the test called `Robosample`.

# Troubleshooting
`freeglut (simbody-visualizer_d): failed to open display ':0'`: make sure Xming is running.
>>>>>>> origin/master
