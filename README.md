# Robosample: Generalized Coordinates Molecular Simulation Coupled with Gibbs Sampling (GCHMC)

Robosample is a C++ library based on Simbody and Molmodel, which uses high-speed robotics algorithms imlemented in Simbody and molecular modelling facilities in Molmodel to generate Markov Chain Monte Carlo moves coupled with Gibbs sampling able to reproduce the atomistic level detailed distribution of molecular systems.

![Docking with Robosample](drug.gif)

[More about the method.](https://pubmed.ncbi.nlm.nih.gov/28892630/)

## Check hardware acceleration capabilities

On native Linux, output should not be empty:

```bash
lspci | grep -iE "vga|3d"
```

On Windows, open `Device Manager` and read the GPUs from `Display adapters`. Also, check WSL version (should be 2) in `Powershell`:

```bash
wsl --version
```

## CUDA

### Native Linux

Check that the driver (if already installed) works correctly:

- Userspace tool works, driver communicates with kernel and GPU is accessible: `nvidia-smi` should output a table showing your GPU name (e.g., RTX 4090) and `Driver Version: 5xx.xx`.

- Userspace tool exists: `command -v nvidia-smi` should point to `/usr/bin/nvidia-smi`. If empty, there is a `$PATH` or package issue.

- Confirm kernel module is active: `cat /proc/driver/nvidia/version` and, if `nvidia-smi` works, matches `nvidia-smi --query-gpu=driver_version --format=csv,noheader`.

- Kernel module is loaded: `lsmod | grep nvidia` which should output `nvidia, nvidia_drm, nvidia_modeset, nvidia_uvm`. If empty, the driver is not loaded.

- Verify module actually exists: `modinfo nvidia | grep filename` and can be inserted manually: `sudo modprobe nvidia`.

- `nouveau`, an open-source reverse-engineered Linux graphics driver for NVIDIA cards, is fully disabled. The following commands should return nothing: `lsmod | grep nouveau`, `cat /etc/modprobe.d/* | grep nouveau` and `lsinitramfs /boot/initrd.img-$(uname -r) | grep nouveau`.

- Kernel module loads successfully: `dmesg | grep -i nvidia`.

- Driver installed for current kernel: `uname -r` and `dkms status | grep nvidia` should match.

- Driver bound to GPU: `lspci -k | grep -A 3 -i nvidia` which should output `Kernel driver in use: nvidia`. If it says `nouveau`, the wrong driver bound to this GPU.

- Secure Boot is disabled: `mokutil --sb-state`. The DKMS-built nvidia.ko module may be blocked from loading unless it is signed and enrolled via MOK. This is one of the most common post-installation failure modes on UEFI systems.

- Driver is actually visible: `lspci | grep -i nvidia`.

Please note that `nvcc` is part of CUDA toolkit, not Nvidia driver, so we do not test for it here.

If any of these tests fail, it means that the driver is incorrectly installed. We begin by checking what drivers we currently have installed:

```bash
dpkg -l | grep -i nvidia
```

Remove the old driver:

```bash
sudo apt purge "*nvidia*"
sudo apt autoremove
```

Confirm that we uninstalled everything:

- `lsmod | grep nvidia` should be empty.

- `lspci -k | grep -A3 -i nvidia` should show `nouveau`.

Drivers can be downloaded either manually from the Nnvidia website or automatically by Ubuntu. However, manual downloads don't support DKMS (Dynamic Kernel Module Support), meaning that kernel updates can potentially break the driver installation. `ubuntu-drivers` is safer because it considers Ubuntu version GLIBC and kernel stability, whereas the website just looks at the GPU model. `ubuntu-drivers` uses DKMS, meaning that the driver automatically rebuilds itself when the kernel updates. Note that both methods respect the compatibility matrix.

First, check what the recommended version is for your machine:

```bash
sudo ubuntu-drivers devices
```

Install the recommended one (`nvidia-smi` will not work before you `sudo reboot`):

```bash
sudo ubuntu-drivers autoinstall
sudo reboot
```

After reboot, check for module load errors:

```bash
dmesg | grep -i nvidia
```

Check if the driver communicates with kernel and GPU is accessible:

```bash
nvidia-smi
```

The `CUDA Version` in the `nvidia-smi` is the maximum CUDA runtime API version supported by the driver. The driver does not install any CUDA toolkit.

### WSL2

In `Powershell`, confirm that `nvidia-smi` works:

```powershell
nvidia-smi
```

If not already done, download and install [normal Windows driver (Game Ready or Studio)](https://www.nvidia.com/Download/index.aspx). **Do not install Linux drivers**. After installation, **reboot Windows**.

Open Linux and check that:

- Windows driver is exposed to WSL: `nvidia-smi` should work and `command -v nvidia-smi` should point to Windows Nvidia driver reflection `/usr/lib/wsl/lib/nvidia-smi`.

- Linux Nvidia driver is not installed: `dpkg -l | grep nvidia` should be empty.

- Confirm kernel module is not activate: `/proc/driver/nvidia/version` and `dkms status | grep nvidia` should be empty.

- No kernel module attempts exist: `lsmod | grep nvidia` should be empty.

- GPU bridge device exists: `ls -l /dev/dxg`.

If any of these tests fail, it means that the driver is incorrectly installed. We begin by checking what drivers we currently have installed:

```bash
dpkg -l | grep -i nvidia
```

Remove the old driver:

```bash
sudo apt purge "*nvidia*"
sudo apt autoremove
```

Shutdown WSL from `Powershell` and restart (WSL does not have a real reboot). Upon restart, WSL will re-link the `/usr/lib/wsl/lib` directory from the Windows host.

```powershell
wsl --shutdown
```

## OpenCL - vendor based

OpenCL installation depends on the udnerlying hardware:

- Nvidia GPU: already installed with the divers.

- AMD GPU: TODO

- Intel GPU/CPU: install `sudo apt-get install intel-opencl-icd`

- AMD CPU: TODO

Other dependencies will be installed via `conda` (see below).

### Installing CUDA Toolkit 12.8

We inted to make this the oldest supported version. It is forward compatible with any future CUDA toolkit versions. Download and follow instructions [from here](https://developer.nvidia.com/cuda-12-8-0-download-archive):

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

### Create a `mamba` environment

Update `conda` and install in base `mamba` (a faster reimplementation of the `conda` package manager) and `conda-merge` (a tool for merging `conda` environment files into one file).

```bash
conda update -n base -c defaults conda
pip install --upgrade pip

conda install conda-forge::mamba
pip install conda-merge zstandard
```

#### CUDA Toolkit

CUDA Toolkit version is dependent on the major GCC version:

| CUDA Toolkit Version | Max Supported GCC Version | Min Supported GCC Version |
| -------------------- | ------------------------- | ------------------------- |
| **13.0, 13.1**       | **15**                    | 7.x                       |
| **12.8, 12.9**       | **14**                    | 6.x                       |
| **12.4, 12.5, 12.6** | **13.2**                  | 6.x                       |
| **12.1, 12.2, 12.3** | **12.2**                  | 6.x                       |
| **12**               | **12.1**                  | 6.x                       |
| **11.4.1 – 11.8**    | **11**                    | 6.x                       |
| **11.1 – 11.4.0**    | **10**                    | 6.x                       |
| **11**               | **9**                     | 6.x                       |
| **10.1, 10.2**       | **8**                     | 4.8.5                     |
| **9.2, 10.0**        | **7**                     | 4.8.5                     |
| **9.0, 9.1**         | **6**                     | 4.8.5                     |
| **8**                | **5.3**                   | 4.7                       |

#### Generating the right development environment `.yaml` file

The development environment is a combination of `.yaml` files from `envs/`:

- `envs/robo_py312.yaml`

- `envs/cuda*.yaml` with the correct `nvcc` - `gcc` version pair

- `envs/opencl.yaml`

- `envs/util.yaml` for non-essential, but useful packages.

We combine into an environment that contains all tools needed to configure and build the project. The name of the combined environment is the name of the last `.yaml` file in sequence (in our case, `robo_py312.yaml`):

```bash
conda-merge envs/cuda12.8.yaml envs/util.yaml envs/robo_py312.yaml > robo_py312.yaml
```

Now we create the environment:

```bash
mamba env create -f robo_py312.yaml
conda activate robo_py312
```

If something goes wrong, delete this environment using:

```bash
conda deactivate
mamba env remove -n robo_py312
mamba clean --all
```

#### Testing CUDA

```bash
nvcc --version && nvidia-smi && nvidia-smi topo -m && nvidia-smi -L
```

#### Testing OpenCL installation

```bash
python -c "import pyopencl as cl; print(cl.get_platforms())"
```

### Building Robosample

If using CUDA, the default architecture is `native` which is not portable. If you build for distribution, go to `CMakePresets.json` and replace `CMAKE_CUDA_ARCHITECTURES` to `70;75;80;86` (Volta, Turing, Ampere (DC) and Ampere (Consumer)).

The build system is organized into a matrix of Platforms and Build Types. You can combine them using the format `--preset <platform>-<type>`.

Platforms: `cpu`, `opencl`, `cuda`, `reference`.

Build types: `debug`, `release`, `relwithdebinfo`, `pgo-train`, `pgo-use`.

While in `robosample/`:

```bash
cmake --preset cuda-release
cmake --build --preset cuda-release
```

This will automatically install the Python bindings (`.so` file) into `robosample/build/robosample`.

### Running the program

We provide a series of examples in `robosample/examples`. To rung the program:

```bash
cd build/robosample
python3 robosample/roborun.py 2ala ../examples/2ala.prmtop ../examples/2ala.inpcrd 6000 0 1 1
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
