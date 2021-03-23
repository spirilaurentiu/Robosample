# Installation

## Installing dependencies
Install the dependencies:
```
sudo apt-get update
sudo apt-get install g++ cmake graphviz gfortran libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig ocl-icd-opencl-dev fftw2 libxmu-dev libxi-dev
```

Add the following to `LD_LIBRARY_PATH`:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib/plugins
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/plugins/libOpenMMPlugin.so

export CUDA_INC_DIR=/usr/local/cuda
export CUDA_ROOT=/usr/local/cuda
```
**WARNING:** It is recommended to run these exports in the WSL or native Linux terminal. When running in the Visual Studio Code terminal on WSL, these exports are not saved (not even in `.bashrc`) and must be restated before executing Robosample.

Also make sure that:
1. `cmake --version` is greater than 3.1.
1. `gcc` version can compile C++17.

##  Cloning the project
**WARNING:** Maximum performance is obtained if the folder project folder is installed in `/home/<user name>/` (see [this](https://docs.microsoft.com/en-us/windows/wsl/compare-versions#performance-across-os-file-systems) for more details). An antivirus program might interfere with WSL and cause performace drops. When running WSL, the user may opt to shut the running antivirus.
```
git clone --recurse-submodules https://github.com/spirilaurentiu/Robosample.git
cd Robosample
```

## Robosample branches
* `master` Stable version. Installed by default.
* `cpp_latest` Stable version with latest C++ standard. Faster than `master`.
* `experimental_cpp_latest` Unstable version with latest C++ version. This is the fastest branch.

To switch branches, go to `Robosample` and execute:
```
cd Molmodel/Simbody01/
git checkout experimental_cpp_latest

cd ../
git checkout experimental_cpp_latest

cd ../
git checkout experimental_cpp_latest
```

## Building Robosample
Run `build_debug.sh` or `build_release.sh`. Password will be required as a lot of files will be deleted from `/usr/`. **Do not run as sudo.**
```
bash build_release.sh
```

# Open the project using any IDE (e.g. Visual Studio Code)
Install [Visual Studio Code](https://code.visualstudio.com/) on Windows. Run `code .` in `Robosample`.

# Working on Robosample
After working on Robosample, it must be compiled as `Debug` or `Release` (the same flags as in section [Building Robosample](#building-robosample)). To compile as a different configuration, full recompilation is needed (see [Building Robosample](#building-robosample)).
```
make -j$((`nproc`*2))
```

# Running Robosample
`Robosample` is located in `build-debug/tests` or `build-release/tests`.
```
cd build-debug
./tests/Robosample inp
```
**WARNING**: To avoid dynamic library lookup failures, it is recommended to run Robosample in the WSL or native Linux terminal and not in Visual Studio Code terminal. See more details at [Installing dependencies](#installing-dependencies).  
  
To change different parameters (use visualizer, use OpenMM etc) edit `inp` which is located in `build-debug` or `build-release`.

# Using the visualizer
Robosample can show in real time the progress of the simulation. To do so, change `VISUAL FALSE` to `VISUAL TRUE` in `inp`.

When running in WSL, an X server such as Xming is needed in order for the program to draw to screen. In Windows, download and install [Xming](https://sourceforge.net/projects/xming/) using the default options. Open `.bashrc` with `vi ~/.bashrc` and add `export DISPLAY=:0` at the end of the file. Xming must be running when Robosample is started with a visualizer\.

## Troubleshooting the visualizer
`freeglut (simbody-visualizer_d): failed to open display ':0'`: make sure [Xming](#installing-xming) is running. If it does not work, change `VISUAL TRUE` to `VISUAL FALSE`.
