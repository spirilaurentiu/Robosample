#!/bin/bash

# TODO: 
# update gitignore


######################################################################
# CHOOSE THE FOLLOWING VARIABLES ACCORDING TO YOUR SYSTEM & DESIRES
######################################################################

# Robosample HOME DIR
export ROBOSAMPLE_HOME=$(pwd);

# Build type = "Release" / "Debug" (argv 1)
export BUILD_TYPE=$1

export BUILD_DIRNAME='build-'${BUILD_TYPE,,}
export INSTALL_DIRNAME='install-'${BUILD_TYPE,,}


# 0/1 if you want to :
#	- checkout to a specific branch/commit 
#	- pull updates 

export SIMBODY_GIT_PULL=0
export SIMBODY_GIT_BRANCH="master"

export MOLMODEL_GIT_PULL=0
export MOLMODEL_GIT_BRANCH="calc_energy_openmm"

export OPENMM_GIT_PULL=0
export OPENMM_GIT_BRANCH="master"

export ROBOSAMPLE_GIT_PULL=0
export ROBOSAMPLE_GIT_BRANCH="calc_energy_openmm"


#	- regenerate cmake (for the first time usage OR if you modified CMake files)
#	- add to PATH variable to bashrc (for the first time usage ONLY !!!!!!)

export SIMBODY_REGENERATE=1
export MOLMODEL_REGENERATE=1
export OPENMM_REGENERATE=1
export ROBOSAMPLE_REGENERATE=1

export SIMBODY_EXPORT_PATH=1
export MOLMODEL_EXPORT_PATH=1
export OPENMM_EXPORT_PATH=1
export ROBOSAMPLE_EXPORT_PATH=1

######################################################################


################
#   SIMBODY    #
################

cd ${ROBOSAMPLE_HOME}/Molmodel/Simbody01

git checkout ${SIMBODY_GIT_BRANCH};
if [ ${SIMBODY_GIT_PULL} -eq 1 ];
	then git pull;
fi

mkdir -p ${BUILD_DIRNAME}
mkdir -p ${INSTALL_DIRNAME}

cd ${BUILD_DIRNAME}

if [ ${SIMBODY_REGENERATE} -eq 1 ]
  then
    rm -rf *
    cmake \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_INSTALL_PREFIX=../${INSTALL_DIRNAME} ..
fi

make -j$(nproc)
make install


if [ ${SIMBODY_EXPORT_PATH} -eq 1 ]
  then
    cd ../${INSTALL_DIRNAME}/
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$(pwd)/lib/" >> ~/.bashrc
fi


################
#   OPENMM     #
################

cd ${ROBOSAMPLE_HOME}/openmm/

git checkout ${OPENMM_GIT_BRANCH}
if [ ${OPENMM_GIT_PULL} -eq 1 ]
  then 
    git pull
fi

mkdir -p ${BUILD_DIRNAME}
mkdir -p ${INSTALL_DIRNAME}

cd ${BUILD_DIRNAME}

if [ ${OPENMM_REGENERATE} -eq 1 ]
  then
    rm -rf *;
    cmake \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_INSTALL_PREFIX=../${INSTALL_DIRNAME} ..
fi

make -j$(nproc)
make install


if [ ${OPENMM_EXPORT_PATH} -eq 1 ]
  then
    cd ../${INSTALL_DIRNAME}/
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$(pwd)/lib/" >> ~/.bashrc
fi



################
#   MOLMODEL   #
################


cd ${ROBOSAMPLE_HOME}/Molmodel/

git checkout ${MOLMODEL_GIT_BRANCH} 
if [ ${MOLMODEL_GIT_PULL} -eq 1 ]
  then 
    git pull;
fi

mkdir -p ${BUILD_DIRNAME}
mkdir -p ${INSTALL_DIRNAME}

cd ${BUILD_DIRNAME}

if [ ${MOLMODEL_REGENERATE} -eq 1 ];
  then;
    rm -rf *;
    cmake \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DSimTK_INSTALL_PREFIX=${ROBOSAMPLE_HOME}/Molmodel/${INSTALL_DIRNAME} \
      -DSimbody_DIR=${ROBOSAMPLE_HOME}/Molmodel/Simbody01/${INSTALL_DIRNAME}/lib/cmake/simbody/ \
      -DOpenMM_INCLUDE_DIR=${ROBOSAMPLE_HOME}/openmm/${INSTALL_DIRNAME}/include \
      -DOpenMM_LIBRARIES=${ROBOSAMPLE_HOME}/openmm/${INSTALL_DIRNAME}/lib/libOpenMM.so \
      -DOpenMM_LIBRARY=${ROBOSAMPLE_HOME}/openmm/${INSTALL_DIRNAME}/lib/libOpenMM.so ..
fi

make -j$(nproc)
make install



if [ ${MOLMODEL_EXPORT_PATH} -eq 1 ]
  then
    cd ../${INSTALL_DIRNAME}/
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$(pwd)/lib/" >> ~/.bashrc
    echo "OpenMMPlugin_PATH=${ROBOSAMPLE_HOME}/Molmodel/${INSTALL_DIRNAME}/lib/plugins" >> ~/.bashrc
fi



################
#  ROBOSAMPLE  #
################


cd ${ROBOSAMPLE_HOME}

git checkout ${ROBOSAMPLE_GIT_BRANCH}
if [ ${ROBOSAMPLE_GIT_PULL} -eq 1 ]
  then 
    git pull;
fi


mkdir -p ${BUILD_DIRNAME}
mkdir -p ${INSTALL_DIRNAME}

cd ${BUILD_DIRNAME}

if [ ${ROBOSAMPLE_REGENERATE} -eq 1 ]
  then
    rm -rf *;
    cmake \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_INSTALL_PREFIX=${ROBOSAMPLE_HOME}/${INSTALL_DIRNAME} \
      -DSimbody_DIR=${ROBOSAMPLE_HOME}/Molmodel/Simbody01/${INSTALL_DIRNAME} \
      -DMolmodel_DIR=${ROBOSAMPLE_HOME}/Molmodel/${INSTALL_DIRNAME} \
      -DOpenMM_INCLUDE_DIR=${ROBOSAMPLE_HOME}/openmm/${INSTALL_DIRNAME}/include \
      -DOpenMM_LIBRARIES=${ROBOSAMPLE_HOME}/openmm/${INSTALL_DIRNAME}/lib/libOpenMM.so \
      -DOpenMM_LIBRARY=${ROBOSAMPLE_HOME}/openmm/${INSTALL_DIRNAME}/lib/libOpenMM.so ..
fi

make -j$(nproc)


source ~/.bashrc
sudo ldconfig


# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# add tests
cp ../tests/test-memory.sh test-memory.sh

