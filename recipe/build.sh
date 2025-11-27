#!/usr/bin/env bash
set -xeuo pipefail

## conda-build provides $PYTHON (host python)
#echo "Using build Python: $PYTHON"
#INC=$($PYTHON - <<'PY'
#import sysconfig
#print(sysconfig.get_paths()["include"])
#PY
#)
#
#LIBDIR=$($PYTHON - <<'PY'
#import sysconfig
#print(sysconfig.get_config_var("LIBDIR") or "")
#PY
#)

export CC=clang
export CXX=clang++

cmake ${SRC_DIR} \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DPython3_EXECUTABLE=$PYTHON \
    -DCUDAToolkit_ROOT=/usr/local/cuda-12.8 \
    -DOPENMM_PLATFORM=CUDA

ninja