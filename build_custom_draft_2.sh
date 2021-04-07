#!/bin/bash

# TODO: 
# update gitignore
# robosample section is unfinished !!!!!!!!!!



######################################################################
# CHOOSE THE FOLLOWING VARIABLES ACCORDING TO YOUR SYSTEM & DESIRES
######################################################################

# Robosample HOME DIR
export ROBOSAMPLE_HOME=$(pwd);

# Build type = "release" / "debug" (argv 1)
export BUILD_TYPE=$1

# 0/1 if you want to checkout to a specific branch/commit and pull updates 
# to the local repo

export SIMBODY_GIT_PULL=0
export SIMBODY_GIT_BRANCH="master"

export MOLMODEL_GIT_PULL=0
export MOLMODEL_GIT_BRANCH="master"

export OPENMM_GIT_PULL=0
export OPENMM_GIT_BRANCH="master"

export ROBOSAMPLE_GIT_PULL=0
export ROBOSAMPLE_GIT_BRANCH="master"


######################################################################


################
#   SIMBODY    #
################

cd ${ROBOSAMPLE_HOME}/Molmodel/Simbody01

if ${SIMBODY_GIT_PULL};
	then git checkout ${SIMBODY_GIT_BRANCH} && git pull;
fi;

mkdir -p build-${BUILD_TYPE}
mkdir -p install-${BUILD_TYPE}

cd build-${BUILD_TYPE}
rm -rf *

cmake \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_INSTALL_PREFIX=../install-${BUILD_TYPE} ..

make -j$((`nproc`*2))
make install



################
#   OPENMM     #
################

cd ${ROBOSAMPLE_HOME}/openmm/

if [ ${OPENMM_GIT_PULL} -eq 1 ];
then git checkout ${OPENMM_GIT_BRANCH} && git pull;
fi;


mkdir -p build-${BUILD_TYPE}
mkdir -p install-${BUILD_TYPE}

cd build-${BUILD_TYPE}
rm -rf *

cmake \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_INSTALL_PREFIX=../install-${BUILD_TYPE} ..

make -j$((`nproc`*2))
make install



################
#   MOLMODEL   #
################


cd ${ROBOSAMPLE_HOME}/Molmodel/

if [ ${MOLMODEL_GIT_PULL} -eq 1 ];
then git checkout ${MOLMODEL_GIT_BRANCH} && git pull;
fi;


mkdir -p build-${BUILD_TYPE}
mkdir -p install-${BUILD_TYPE}

cd build-${BUILD_TYPE}
rm -rf *

cmake \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DSimTK_INSTALL_PREFIX=${Robosample_HOME}/Molmodel/install-${BUILD_TYPE} \
	-DSimbody_DIR=${Robosample_HOME}/Molmodel/Simbody01/install-${BUILD_TYPE}/lib/cmake/simbody/ \
	-DOpenMM_INCLUDE_DIR=${Robosample_HOME}/openmm/install-${BUILD_TYPE}/include \
	-DOpenMM_LIBRARIES=${Robosample_HOME}/openmm/install-${BUILD_TYPE}/lib/libOpenMM.so \
	-DOpenMM_LIBRARY=${Robosample_HOME}/openmm/install-${BUILD_TYPE}/lib/libOpenMM.so ..

make -j$((`nproc`*2))
make install





################
#  ROBOSAMPLE  #
################


cd ${ROBOSAMPLE_HOME}

if [ ${ROBOSAMPLE_GIT_PULL} -eq 1 ];
then git checkout ${ROBOSAMPLE_GIT_BRANCH} && git pull;
fi;


mkdir -p build-${BUILD_TYPE}
mkdir -p install-${BUILD_TYPE}

cd build-${BUILD_TYPE}
rm -rf *

cmake \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	# -DCMAKE_INSTALL_PREFIX=${Robosample_HOME}/install-${BUILD_TYPE} \
	# -DSimbody_PATH=${Robosample_HOME}/Molmodel/Simbody01/install-${BUILD_TYPE} \
	# -DMolModel_PATH=${Robosample_HOME}/Molmodel/install-${BUILD_TYPE} \
	# -DOpenMM_PATH=${Robosample_HOME}/openmm/install-${BUILD_TYPE} ..


make -j$((`nproc`*2))

# make install




# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# add tests
cp ../tests/test-memory.sh test-memory.sh
