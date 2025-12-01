# changelog
* full machine setup (drivers and prerequisites). local conda package. debug in python
* python wrappers with type safety
* c++ unit safety
* single openmm with low vram footprint. openmm plugin and robosample openmm are now gone
* non-redundant atom/bond/angle/torsion types for low molmodel ram footprint
* internal coordinates moved to python

# TODO
* getTransformAndResidual() in Context.cpp - see nan when we have a perfect match



you are a software engineer and architect with over 20 years in scientific computing. you will be brutally honest with me and question whether every decision we make is the best or not since we are working with critical software components. you are going to rate how sure you are and list what you don't know exactly 

so we are developing robosample, a new molecular simulation software based on Generalized Coordinates Molecular Simulation Coupled with Gibbs Sampling (GCHMC) that uses robot mechanics and openmm (cuda) for energy calculations and integration.
this program is written in c++ and builds the sources of pybind11 (in cmakelists it's pybind11_add_module), openmm, simbody and molmodel into one single binary using llvm-clang, applies lto, then pgo, then llvm-bolt for maximum efficiency.
the build pipeline works (i tested it on several machines) but it forces the user (me or my colleagues from the lab) to compile robosample on their machines (this is messy because we're only two devs in the team who have the knowlegde to set up everything correctly)
it also comes with a handful of python modules that offer type safety (wrappers) and molecular functionality (like atom selection, user interface etc). 
for countless times we have faced deployment problems regarding dependencies compatibilities. we have thought that distributing this package (the .so and some .py files) as conda package might be helpful, but now i'm doubting my strategies.
i also understand that we also have the option of using wheels - not really sure why we would ever use them, i understand that they only work with pure python projects.
it is borderline impossible to compile robosample as with the tools provided by conda-build (gcc - not clang, no support for the other tools).
the .so links only cuda (system wide installation - from my understanding, cuda 12 cannot be installed as a conda package), blas and some other minor libraries. as a side questions: should we depend on system-wide or conda versions of cuda? what's the difference between them?
i don't believe the abi will ever change (we are distributing only the .so and only some minor functions from it are going to be called anyway).
what can i do? i already wasted like a week trying to make this shit work (and stopped dead in the face of conda-tools). how did openmm or pytorch do it?




let's go with the external route. i want to be able to compile my .so outside conda-build. i want to target python 3.12 and cuda 12.8
for this, i installed system-wide cuda toolkit 12.8.
we build using clang, cmake and ninja. do they have to be system-wide or can they be provided by conda? i believe system-wide, but tell me
robosample also depends on some other external libraries: linux-vdso.so.1 libblas.so.3 liblapack.so.3 libm.so.6 libcuda.so.1 libcufft.so.11 libstdc++.so.6 libgcc_s.so.1 libc.so.6 ld-linux-x86-64.so.2 libgfortran.so.5 libdl.so.2 libpthread.so.0 librt.so.1
right now, it links only the system wide libraries. is this ok since the .so will be distributed as package? do i need to build the .so linking the conda versions? i'm really confused about these


(what tasks for pgo and llvm)?




# Robosample: Generalized Coordinates Molecular Simulation Coupled with Gibbs Sampling (GCHMC)

Robosample is a C++ library based on Simbody and Molmodel, which uses high-speed robotics algorithms imlemented in Simbody and molecular modelling facilities in Molmodel to generate Markov Chain Monte Carlo moves coupled with Gibbs sampling able to reproduce the atomistic level detailed distribution of molecular systems.

![Docking with Robosample](drug.gif)

[More about the method.](https://pubmed.ncbi.nlm.nih.gov/28892630/)

## Installing dependencies (for developers)

### Installing the Nvidia driver
The only working driver is `proprietary`, not `open`. Remove the `open` driver and all CUDA Toolkit installations:
```bash
sudo apt --fix-broken install
sudo apt-get purge 'nvidia-*' 'cuda*'
sudo apt-get autoremove
sudo rm -rf /usr/local/cuda*
```

First, check what the recommended version is for your machine. To my knowledge, the last supported kernel is 6.14 (check with `uname -r`).
```
ubuntu-drivers devices
```

Install the recommended one (`nvidia-smi` will not work before you `sudo reboot`):
```
sudo apt install nvidia-driver-580
sudo reboot
```

### Installing CUDA Toolkit 12.8
Only this exact version can be used and is hard-coded inside the build files. It is forward compatible with any future CUDA toolkit versions. Download and follow instructions from [here](https://developer.nvidia.com/cuda-12-8-0-download-archive):
```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-8
```

At this point, you should have working `nvidia-smi`. However, `nvcc` is not available at the terminal since nothing from the CUDA Toolkit is in `$PATH`.
It can still be accessed via the absolute path: `/usr/local/cuda/bin/nvcc --version`.
This is no problem since Robosample references absolute pathways to this installation.

### Installing OpenGL (visualizer) and OpenCL (hardware acceleration)
Straightforward installation that does not interfere with CUDA:
```bash
sudo update
sudo apt-get install libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libxmu-dev libxi-dev ocl-icd-opencl-dev
```

### Miniforge

We will use `miniforge`, a variant of `miniconda` that comes with `mamba` installed (`conda` but with a faster solver).
To my knowledge, `miniconda` and `miniforge` are theoretically compatible and can run concurrently on the same machine.
However, uncertainty still looms over this, so we prefer so remove any `miniconda` installations before doing anything else.
```bash
conda deactivate
~/miniconda3/uninstall.sh
```

Download and install `miniforge` from their [GitHub page](https://github.com/conda-forge/miniforge) using:
```bash
cd ~
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
source ~/.bashrc
```

### Other dependencies
Install the dependencies:
```bash
sudo apt-get update
sudo apt-get install git cmake graphviz gfortran libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig fftw2 clang ninja-build linux-tools-common linux-tools-generic linux-tools-`uname -r`
```

### CMake
Minimum `CMake` version is 3.17. It can be tested with:
```bash
cmake --version
```

To install the correct version, head to the [CMake website](https://cmake.org/download/) and find the `.tar.gz` version for your operating system. The code below is an example for CMake 3.27 under Linux:
```bash
cd ~
wget https://github.com/Kitware/CMake/releases/download/v3.27.7/cmake-3.27.7-linux-x86_64.tar.gz
tar -xf cmake-3.27.7-linux-x86_64.tar.gz
rm cmake-3.27.7-linux-x86_64.tar.gz
```

The executable is located in the `bin` folder:
```bash
~/cmake-3.27.7-linux-x86_64/cmake-3.27.7-linux-x86_64/bin/cmake
```

### Ninja
There is no required Ninja version. It can also be replaced with Unix Makefiles. If installation from `apt-get` fails, downloading the binaries is recommended. Go to [Ninja website](https://ninja-build.org/) and find the [binary downloads](https://ninja-build.org/). The following is an example for version 1.11.1:
```bash
cd ~
wget https://github.com/ninja-build/ninja/releases/download/v1.11.1/ninja-linux.zip
unzip ninja-linux.zip
rm ninja-linux.zip
```

If used as intended further into the README, the executable must be run from the full path:

```bash
/home/myuser/ninja
```

## Download and compile Robosample
```bash
git clone --recurse-submodules https://github.com/spirilaurentiu/Robosample.git
cd Robosample

cd openmm
git checkout master
cd ../Molmodel
git checkout singularity
cd ../Simbody01
git checkout master
cd ../
git checkout refactor
```

Create the build environment:
```bash
mamba env create -f tools/robo_dev_py312.yaml
conda activate robo_dev_py312
```

OpenMM can use hardware acceleration. Robosample defaults with OpenCL. To set the platform, you can set it via the `cmake` command in the next step:
* `OPENMM_PLATFORM=CPU` for CPU.
* `OPENMM_PLATFORM=CUDA` for CUDA.
* `OPENMM_PLATFORM=OPENCL` for OpenCL.

```bash
mkdir -p build
cd build
cmake -G Ninja ../ -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D OPENMM_PLATFORM=CUDA
ninja
```

Assuming that CMake and Ninja have been installed as binaries and not from `apt-get`:
```bash
~/cmake-3.27.7-linux-x86_64/bin/cmake -G Ninja -DCMAKE_MAKE_PROGRAM=/home/myuser/ninja -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -D OPENMM_PLATFORM=CUDA
~/ninja robosample
```

If you want to use Unix Makefiles (please don't):
```bash
cmake -G "Unix Makefiles" ../ -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D OPENMM_PLATFORM=CUDA
make -j$(nproc)
```

## Running from the Python interface
```bash
conda install -c conda-forge mamba
mamba env create -f tools/robo_dev_py312.yaml
```

In order to run the program under VSCode debugger, place breakpoints in `simulate.py` and run with `Python: Current File (Debug Robosample Libraries)`. After it stops, run start `(gdb) Attach to Python`, enter the PID of the `python` instance (usually, it has `-X frozen_modules=OFF` as arguments) and press the `Continue (F5)` button of the VSCode debugger. Look in the `CALL STACK` section on the left side of the screen and press `(gdb) Attach to Python` and it should take you to the breakpoint placed inside a `.cpp` file. For more details, see [this](https://nadiah.org/2020/03/01/example-debug-mixed-python-c-in-visual-studio-code/).

## LLV-BOLT (Binary Optimization and Layout Tool)
We have applied [LLVM-BOLT](https://github.com/llvm/llvm-project/tree/main/bolt), improving the execution speed by rearranging code layout based on execution profiles from sampling profilers like `perf`.
Downlad BOLT and compile it:
```bash
cd ~
git clone https://github.com/llvm/llvm-project.git
cd llvm-project
mkdir build
cd build
cmake -G Ninja ../llvm -DLLVM_TARGETS_TO_BUILD="X86;AArch64" -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_ASSERTIONS=ON -DLLVM_ENABLE_PROJECTS="bolt"
ninja bolt
```

Add BOLT to `PATH`:
```bash
echo "PATH=$(pwd)/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
```

Allow instrumentation:
```bash
sudo sysctl kernel.perf_event_paranoid=-1
```

### Instrumentation
We have discovered that running **only one** simulation round yields the best result. Also, using a larger system seems to be optimal. Perform the necessary changes in the input file and execute:
```bash
perf record -e cycles:u -j any,u -a -o perf.data ./robosample.pgo.use inp.aper
```

Convert the data into something that can be used by BOLT:
```bash
perf2bolt -p perf.data robosample.pgo.use -o perf.fdata
```

Optimize the binary:
```bash
llvm-bolt robosample.pgo.use -o robosample.pgo.use.bolt -data=perf.fdata -reorder-blocks=ext-tsp -reorder-functions=hfsort -split-functions -split-all-cold -split-eh -dyno-stats
```

Compare the binaries:
```bash
time ./robosample inp.ala10
time ./robosample.bolt inp.ala10
```

## PGO (Profile Guided Optimization)
PGO requires us to compile to compile Robosample once, run it a few times and compile it again taking into account the hot code paths.
First compilation:
```bash
cmake -G Ninja ../ -D CMAKE_BUILD_TYPE=PGO_Train -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D OPENMM_PLATFORM=OPENCL
ninja robosample
```

Clear output of previous runs:
```bash
find . -name "*.gcda" -delete
```

Run the examples:
```bash
bash pgo.sh
```







END OF TIME HERE DON'T GO BELOW








## Sanitizers (**mandatory**)

We use address and undefined behaviour sanitizers in our debug builds. To get the correct output, run:

```bash
echo "export ASAN_OPTIONS=detect_odr_violation=0:detect_leaks=0:protect_shadow_gap=0" >> ~/.bashrc
echo "export UBSAN_OPTIONS=print_stacktrace=1" >> ~/.bashrc
source ~/.bashrc
```

Explaination:

* `detect_odr_violation`
* `detect_leaks=0` - OpenMM has some memory leaks. Set `detect_leaks=1` if you want to see memory all leaks.
* `protect_shadow_gap=0` - OpenCL and CUDA (which both use the NVIDIA driver) conflict with ASAN, as stated by [here](https://stackoverflow.com/a/68027496/3740613).
* `print_stacktrace=1`: show which lines trigger the undefined behaviour sanitizer (UBSAN).

## Fun facts

To get the total number of lines in header and source files, execute this from the root directory:

```bash
find . -name '*.h' -o -name '*.cpp' | xargs wc -l
```

To see all exported symbols, use:

```bash
nm -an build/robosample | c++filt
```
 

bash
```
cmake -G Ninja ../ -D CMAKE_BUILD_TYPE=PGO_Train -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D OPENMM_PLATFORM=OPENCL
ninja robosample
```
