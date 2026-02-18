# Robosample: Generalized Coordinates Molecular Simulation Coupled with Gibbs Sampling (GCHMC)

Robosample is a C++ library based on Simbody and Molmodel, which uses high-speed robotics algorithms imlemented in Simbody and molecular modelling facilities in Molmodel to generate Markov Chain Monte Carlo moves coupled with Gibbs sampling able to reproduce the atomistic level detailed distribution of molecular systems.

![Docking with Robosample](drug.gif)

[More about the method.](https://pubmed.ncbi.nlm.nih.gov/28892630/)

## Installing dependencies

### Installing the Nvidia driver for native Linux
The only working driver is `proprietary`, not `open`. Remove the `open` driver and all CUDA Toolkit installations:
```bash
sudo apt --fix-broken install
sudo apt-get purge 'nvidia-*' 'cuda*'
sudo apt-get autoremove
sudo rm -rf /usr/local/cuda*
```

First, check what the recommended version is for your machine. To my knowledge, the last supported kernel is `6.14` (check with `uname -r`).
```bash
ubuntu-drivers devices
```

Install the recommended one (`nvidia-smi` will not work before you `sudo reboot`):
```bash
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

At this point, you should have working `nvidia-smi`.

However, `nvcc` is not available at the terminal since nothing from the CUDA Toolkit is in `$PATH`. It can still be accessed via the absolute path e.g. `/usr/local/cuda/bin/nvcc --version`. This is no problem since Robosample CMake configuration references absolute pathways to `/usr/local/cuda/bin/`.

### [DEPRECATED] Installing OpenGL (visualizer)
Straightforward installation that does not interfere with CUDA:
```bash
sudo update
sudo apt-get install libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libxmu-dev libxi-dev
```

### Installing OpenCL for hardware acceleration
```bash
sudo update
sudo apt-get install ocl-icd-opencl-dev clinfo
```

### Other dependencies
```bash
sudo apt-get update
sudo apt-get install git cmake graphviz gfortran libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig fftw2 clang ninja-build linux-tools-common linux-tools-generic linux-tools-`uname -r`
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

Install `mamba`:
```bash
conda install conda-forge::mamba
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

## Installing Robosample

### Clone Robosample
```bash
git clone --recurse-submodules https://github.com/spirilaurentiu/Robosample.git
cd Robosample

cd openmm
git checkout master
cd ../Simbody01
git checkout master
cd ../Molmodel
git checkout refactor
cd ../
git checkout refactor
```

### Create a `mamba` environment:
```bash
mamba env create -f tools/robo_py312.yaml
conda activate robo_py312
```

If using `CUDA`, update the environment:
```bash
mamba env update -f tools/robo_py312_cuda.yaml
```

Test that `OpenMM` is installed correctly:
```bash
python -m openmm.testInstallation
```

If something goes wrong, delete this environment using:
```bash
conda deactivate
mamba env remove -n robo_py312
mamba clean --all
```

### Configuring Robosample
You can set the OpenMM hardware acceleration platform using `USE_CPU=ON`, `USE_OPENCL=ON` or `USE_CUDA=ON`.
```bash
mkdir -p build
cd build
cmake -G Ninja ../ -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D USE_CUDA=ON
```

If you want to use Unix Makefiles:
```bash
cmake -G "Unix Makefiles" ../ -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D OPENMM_PLATFORM=CUDA
make -j$(nproc)
```

When examining the output of the CMake configuration run, you should see:
- `Found Python3:` pointing to the `conda` package (inside `~/miniforge3/bin`), not the system-wide Python package. This is needed to ensure that we compile for a certain version of Python specified in the development environment.
- `Found OpenCL` pointing to the system-wide version (inside `/usr/lib/x86_64-linux-gnu/`), not the `conda` version.
- `Found CUDAToolkit` and `Check for working CUDA compiler` both pointing to the system-wide version (inside `/usr/local/cuda`), not the `conda` version.

### Compiling Robosample
Note that we call `ninja`, not `ninja robosample`.
```bash
ninja
```

### Running the program
```bash
cd build/
python3 roborun.py 2ala ../examples/2ala.prmtop ../examples/2ala.inpcrd 6000 1 1 1
```

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
