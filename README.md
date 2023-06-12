# Robosample: Generalized Coordinates Molecular Simulation Coupled with Gibbs Sampling (GCHMC)

Robosample is a C++ library based on Simbody and Molmodel, which uses high-speed robotics algorithms imlemented in Simbody and molecular modelling facilities in Molmodel to generate Markov Chain Monte Carlo moves coupled with Gibbs sampling able to reproduce the atomistic level detailed distribution of molecular systems.

![](drug.gif)

[More about the method.](https://pubmed.ncbi.nlm.nih.gov/28892630/)

# Installation
## Prerequisites
* Turn off all `conda` environments.
* Disable all antivirus programs (especially for WSL).
* Execute all commands in the native terminal if running under WSL.
* Install everything in `/home/<username>/` if running under WSL.

## Installing dependencies
Install the dependencies:
```
sudo apt-get update
sudo apt-get install git cmake graphviz gfortran libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig ocl-icd-opencl-dev fftw2 libxmu-dev libxi-dev clang ninja-build linux-tools-common linux-tools-generic linux-tools-`uname -r`
```

### Fill in exports exports:
```
export CUDA_INC_DIR=/usr/local/cuda
export CUDA_ROOT=/usr/local/cuda
```

##  Cloning the project
```
git clone --recurse-submodules https://github.com/spirilaurentiu/Robosample.git
cd Robosample
```

## Robosample branches
* `master` Stable version. Install by executing (from Robosample directory):
```
cd openmm/
git checkout master

cd ../Simbody01/
git checkout master

cd ../Molmodel/
git checkout master

cd ../
git checkout master
```
* `cpp_latest` - Stable version with latest C++ standard.
* `experimental_cpp_latest` - Unstable version with latest C++ version.

## Building Robosample
### Compiler
We have compiled Robosample with `gcc` and `clang`:
* `clang` compiles faster and produces marginally faster code.
* `GCC 11` produces ICE (internal compiler error) for OpenMM when using IPO. This is not the case with earlier versions (`GCC 7.5` works).

### OpenMM platform
OpenMM can use hardware acceleration. Robosample defaults with OpenCL. To set the platform, you can set it via the `cmake` command in the next step:
* `OPENMM_PLATFORM=CPU` for CPU.
* `OPENMM_PLATFORM=CUDA` for CUDA.
* `OPENMM_PLATFORM=OPENCL` for OpenCL.

### Building Robosample
```
mkdir build
cd build
cmake -G Ninja ../ -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D OPENMM_PLATFORM=OPENCL
ninja robosample
```

# Running Robosample
The executable is compiled in `build/`. Examples are available in the same folder. To see them, type:
```
ll inp.*
```

To run any of them, execute:
```
./robosample inp.2but
```

# BOLT
We have applied [LLVM-BOLT](https://github.com/llvm/llvm-project/tree/main/bolt) developed by Facebook.

## Installing BOLT
Downlad BOLT and compile it. A docker file is also available.
```
git clone https://github.com/llvm/llvm-project.git
mkdir build
cd build
cmake -G Ninja ../llvm-project/llvm -DLLVM_TARGETS_TO_BUILD="X86;AArch64" -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_ASSERTIONS=ON -DLLVM_ENABLE_PROJECTS="bolt"
ninja bolt
```

Add BOLT to `PATH`:
```
echo "PATH=$(pwd)/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
```

Allow intrumentation:
```
sudo echo "-1" > /proc/sys/kernel/perf_event_paranoid
```
This does not seem to fix the value forever. You will most likely need to change it every time you want to run BOLT.

## Instrumentation
We have discovered that running **only one** simulation round yields the best result. Also, using a larger system seems to be optimal. Perform the necessary changes in the input file and execute:
```
perf record -e cycles:u -j any,u -a -o perf.data ./robosample inp.aper
```

Convert the data into something that can be used by BOLT:
```
perf2bolt -p perf.data robosample -o perf.fdata
```

Optimize the binary:
```
llvm-bolt robosample -o robosample.bolt -data=perf.fdata -reorder-blocks=ext-tsp -reorder-functions=hfsort -split-functions -split-all-cold -split-eh -dyno-stats
```

Finally, the optimized binary can be run as:
```
./robosample.bolt inp.2but
```


# Development
Robosample is being developed using [Visual Studio Code](https://code.visualstudio.com/) on Windows. To start it, run `code .` in `Robosample`.

## Sanitizers
We use address and undefined behaviour sanitizers in our debug builds. To get the correct output, run:
```
echo "export ASAN_OPTIONS=detect_odr_violation=0:detect_leaks=0" >> ~/.bashrc
echo "export UBSAN_OPTIONS=print_stacktrace=1" >> ~/.bashrc
source ~/.bashrc
```

Set `detect_leaks=1` if you want to see memory leaks.

## Fun facts
To get the total number of lines in header and source files, execute this from the root directory:
```
find . -name '*.h' -o -name '*.cpp' | xargs wc -l
```

To see all exported symbols, use:
```
nm -an build/robosample | c++filt
```
