#!/bin/bash

PASS=0; FAIL=0; WARNS=()

check() {
    local label=$1; shift
    local output
    output=$("$@" 2>&1)
    if [ $? -eq 0 ]; then
        echo "  PASS  $label"
        ((PASS++))
    else
        echo "  FAIL  $label"
        WARNS+=("FAIL: $label\n$output")
        ((FAIL++))
    fi
}

# Verify we're inside a conda env
if [ -z "$CONDA_PREFIX" ]; then
    echo "ERROR: not inside a conda environment. Activate it first."
    exit 1
fi
echo "Environment : $CONDA_DEFAULT_ENV"
echo "Prefix      : $CONDA_PREFIX"
echo ""

# CUDA is conda-local, not system
echo "[ CUDA / GPU ]"
check "nvcc is conda-local"   bash -c "which nvcc | grep -q '$CONDA_PREFIX'"
check "nvcc version"          nvcc --version
check "libcudart loadable"    python -c "import ctypes; ctypes.CDLL('libcudart.so')"
check "clinfo"                clinfo

echo ""
echo "[ PyTorch CUDA ]"
check "torch.cuda available"  python -c "
import torch
assert torch.cuda.is_available(), 'no CUDA'
print(f'torch {torch.__version__}, cuda {torch.version.cuda}, device: {torch.cuda.get_device_name(0)}')
"
check "torch CUDA compute"    python -c "
import torch
assert torch.cuda.is_available(), 'no CUDA'
a = torch.randn(1000, 1000, device='cuda')
b = torch.randn(1000, 1000, device='cuda')
c = torch.mm(a, b)
torch.cuda.synchronize()
print(f'matmul OK — tensor on: {c.device}')
"
check "nccl"                  python -c "
import torch
v = torch.cuda.nccl.version()
print(f'NCCL version: {v}')
"

echo ""
echo "[ OpenMM CUDA ]"
check "openmm testInstallation" python -m openmm.testInstallation
check "openmm CUDA platform"    python -c "
import openmm as mm
import openmm.app as app
import openmm.unit as unit

system = mm.System()
system.addParticle(1.0)
force = mm.CustomExternalForce('0.5*k*(x^2+y^2+z^2)')
force.addGlobalParameter('k', 1.0)
force.addParticle(0, [])
system.addForce(force)

platform = mm.Platform.getPlatformByName('CUDA')
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
context = mm.Context(system, integrator, platform)
context.setPositions([[0.1, 0.0, 0.0]])
integrator.step(100)
state = context.getState(getEnergy=True)
print(f'platform: {context.getPlatform().getName()}, energy: {state.getPotentialEnergy()}')
"
check "openmm version+platforms" python -c "
import openmm
print(openmm.__version__)
print(openmm.Platform.getNumPlatforms(), 'platforms')
for i in range(openmm.Platform.getNumPlatforms()):
    print(' ', openmm.Platform.getPlatform(i).getName())
"

echo ""
echo "[ AmberTools ]"
check "tleap"               tleap -h
check "cpptraj"             cpptraj --version
check "antechamber"         antechamber -h
check "sqm"                 sqm -h
check "pytraj"              python -c "import pytraj; print(pytraj.__version__)"
check "parmed"              python -c "import parmed; print(parmed.__version__)"
check "netcdf4"             python -c "import netCDF4; print(netCDF4.__version__)"

echo ""
echo "[ MD / Analysis ]"
check "mdtraj"              python -c "import mdtraj; print(mdtraj.__version__)"
check "MDAnalysis"          python -c "import MDAnalysis; print(MDAnalysis.__version__)"
check "openmmtools"         python -c "import openmmtools; print(openmmtools.__version__)"

echo ""
echo "[ Cheminformatics ]"
check "rdkit"               python -c "import rdkit; print(rdkit.__version__)"
check "openbabel"           python -c "import openbabel; print(openbabel.__version__)"
check "pdb2pqr"             python -c "import pdb2pqr; print(pdb2pqr.__version__)"
check "pymol headless"      pymol -c

echo ""
echo "[ Build Toolchain ]"
check "OPENMM_CUDA_COMPILER set"      bash -c "[ -n '$OPENMM_CUDA_COMPILER' ]"
check "OPENMM_CUDA_COMPILER expanded" bash -c "echo '$OPENMM_CUDA_COMPILER' | grep -qv '\$CONDA_PREFIX'"
check "OPENMM_CUDA_COMPILER works"    bash -c "$OPENMM_CUDA_COMPILER --version"
check "CUDA_HOST_COMPILER set"        bash -c "[ -n '$CUDA_HOST_COMPILER' ]"
check "CUDA_HOST_COMPILER expanded"   bash -c "echo '$CUDA_HOST_COMPILER' | grep -qv '\$CONDA_PREFIX'"
check "CUDA_HOST_COMPILER works"      bash -c "$CUDA_HOST_COMPILER --version"
check "cmake"                         cmake --version
check "ninja"                         ninja --version
check "pybind11"                      python -c "import pybind11; print(pybind11.get_include())"

echo ""
echo "[ Science Stack ]"
check "numpy"           python -c "import numpy; print(numpy.__version__)" 
check "scipy"           python -c "import scipy; print(scipy.__version__)"
check "pandas"          python -c "import pandas; print(pandas.__version__)"
check "matplotlib"      python -c "import matplotlib; print(matplotlib.__version__)"
check "astropy"         python -c "import astropy; print(astropy.__version__)"
check "networkx"        python -c "import networkx; print(networkx.__version__)"
check "seaborn"         python -c "import seaborn; print(seaborn.__version__)"
check "sklearn"         python -c "import sklearn; print(sklearn.__version__)"
check "erfa"            python -c "import erfa; print(erfa.__version__)"

echo ""
echo "[ Modeller ]"
check "modeller"        python -c "import modeller; print(modeller.__version__)"

# Summary
echo ""
echo "===================================="
echo "  $PASS passed, $FAIL failed"
echo "===================================="

if [ ${#WARNS[@]} -gt 0 ]; then
    echo ""
    echo "Failure details:"
    for w in "${WARNS[@]}"; do
        echo -e "  $w"
        echo ""
    done
fi

[ $FAIL -eq 0 ]
