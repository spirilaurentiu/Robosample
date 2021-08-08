# Installation

## Installing dependencies
Install the dependencies:
```
sudo apt-get update
sudo apt-get install cmake graphviz gfortran libglfw3-dev freeglut3-dev libglew-dev libxmu-dev libeigen3-dev doxygen subversion libblas-dev liblapack-dev libboost-all-dev swig ocl-icd-opencl-dev fftw2 libxmu-dev libxi-dev
```

## Install a compiler
### For GCC, execute:
```
sudo apt-get install g++
```

### For Clang (recommended)
```
sudo apt-get install clang-12
```
or whatever version suits you and is C++17 capable.

For older machines (<20.04), follow [this](https://solarianprogrammer.com/2017/12/13/linux-wsl-install-clang-libcpp-compile-cpp-17-programs/) guide. A truncated version follows. We prefer this method because verison 9.0 is capable of compiling C++17 and can be easilty installed on many machines (14.04 to 18.04).
```
sudo apt install build-essential xz-utils curl
cd ~
```

Now download the compiler. Make sure the version suitable for yout Linux version is selected (14.04, 16.04, 18.04). You can check this verison with `cat /etc/issue`. Assuming 18.04, execute:
```
curl -SL http://releases.llvm.org/9.0.0/clang+llvm-9.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz | tar -xJC .
```

Continue with the installation.
```
sudo mv clang+llvm-9.0.0-x86_64-linux-gnu-ubuntu-18.04 /usr/local/clang_9.0.0
export PATH=/usr/local/clang_9.0.0/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/clang_9.0.0/lib:$LD_LIBRARY_PATH
```

By now, you should have these two commands available on the terminal:
```
clang
clang++
```

### Fill in `LD_LIBRARY_PATH` exports:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmm/lib/plugins
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/plugins/libOpenMMPlugin.so

export CUDA_INC_DIR=/usr/local/cuda
export CUDA_ROOT=/usr/local/cuda
```
**WARNING:** It is recommended to run these exports in the WSL or native Linux terminal. When running in the Visual Studio Code terminal on WSL, these exports are not saved (not even in `.bashrc`) and must be restated before executing Robosample.

Also make sure that `cmake --version` is greater than 3.1.

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
Run `build.sh`. A mandatory parameter is `-b` (`--build`) which specifies whether to build as debug or release. Other two mandatory parameters are `--cc` (C compiler) and `--cpp` (C++ compiler). You can specify full path to the compiler (`/usr/local/clang-12/bin/clang++`) or just pass what you would normally write to the command line (`clang++`). To set the number of CPU cores involved in compilation, use `-n` (`--nproc`, default is `nproc * 2`). On release builds, if you want to have an address sanitizer, compile with `-u` (`--ubsan`). **Do not run as sudo.**

Building for release with Clang:
```
bash build.sh -b release --cc clang --cpp clang++
```

Building for release with Clang:
```
bash build.sh -b release --cc /usr/local/clang-12/bin/clang --cpp /usr/local/clang-12/bin/clang++
```

Building for debug with GCC:
```
bash build.sh -b debug --cc gcc --cpp g++
```

Building with fewer cores (machines with limited resources, especially RAM, can fail if too many cores are used):
```
bash build.sh -b release --cc gcc --cpp g++ -n 2
```

Building with address sanitizer (already on if compiling for debug):
```
bash build.sh -b release --cc gcc --cpp g++ -u
```

To read more information, execute with the help flag (`-h` or `--help`):
```
bash build.sh -h
```

# Open the project using any IDE (e.g. Visual Studio Code)
Install [Visual Studio Code](https://code.visualstudio.com/) on Windows. Run `code .` in `Robosample`.

# Working on Robosample
After working on Robosample, it must be compiled as debug or release (the same flags as in section [Building Robosample](#building-robosample)). To compile as a different configuration, full recompilation is needed (see [Building Robosample](#building-robosample)).
```
cd build/debug/robosample
make -j$((`nproc`*2))
```

# Running Robosample
`Robosample` is located in `build/debug/robosample` or `build/release/robosample`.
```
cd build/release/robosample
./src/GMOLMODEL_robo inp.2but
```
**WARNING**: To avoid dynamic library lookup failures, it is recommended to run Robosample in the WSL or native Linux terminal and not in Visual Studio Code terminal. See more details at [Installing dependencies](#installing-dependencies).  
  
To change different parameters (use visualizer, use OpenMM etc) edit `inp.2but` which is located in `build-debug` or `build-release`.

# Using the visualizer
Robosample can show in real time the progress of the simulation. To do so, change `VISUAL FALSE` to `VISUAL TRUE` in `inp`.

When running in WSL, an X server such as Xming is needed in order for the program to draw to screen. In Windows, download and install [Xming](https://sourceforge.net/projects/xming/) using the default options. Open `.bashrc` with `vi ~/.bashrc` and add `export DISPLAY=:0` at the end of the file. Xming must be running when Robosample is started with a visualizer\.

## Troubleshooting the visualizer
`freeglut (simbody-visualizer_d): failed to open display ':0'`: make sure [Xming](#installing-xming) is running. If it does not work, change `VISUAL TRUE` to `VISUAL FALSE`.
