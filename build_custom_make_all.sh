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


######################################################################

if [ ! "${BUILD_TYPE}" = "Release" ] && [ ! "${BUILD_TYPE}" = "Debug" ] 
  then 
     echo "Specify build type: Debug/Release !"
     echo "Example:"
     echo "bash build_custom_1st_time.sh Release"
     exit 1;
fi

export BUILD_DIRNAME='build-'${BUILD_TYPE,,}
export INSTALL_DIRNAME='install-'${BUILD_TYPE,,}

echo "BUILD_DIRNAME is : ${BUILD_DIRNAME}"
echo "INSTALL_DIRNAME is : ${INSTALL_DIRNAME}"



################
#   SIMBODY    #
################

cd ${ROBOSAMPLE_HOME}/Molmodel/Simbody01
echo "Current directory is : $(pwd)"

git checkout ${SIMBODY_GIT_BRANCH};
if [ ${SIMBODY_GIT_PULL} -eq 1 ];
	then git pull;
fi

cd ${BUILD_DIRNAME}

make -j$(nproc)
make install


################
#   OPENMM     #
################

cd ${ROBOSAMPLE_HOME}/openmm/
echo "Current directory is : $(pwd)"

git checkout ${OPENMM_GIT_BRANCH}
if [ ${OPENMM_GIT_PULL} -eq 1 ]
  then 
    git pull
fi

cd ${BUILD_DIRNAME}

make -j$(nproc)
make install



################
#   MOLMODEL   #
################


cd ${ROBOSAMPLE_HOME}/Molmodel/
echo "Current directory is : $(pwd)"


git checkout ${MOLMODEL_GIT_BRANCH} 
if [ ${MOLMODEL_GIT_PULL} -eq 1 ]
  then 
    git pull;
fi

cd ${BUILD_DIRNAME}

make -j$(nproc)
make install




################
#  ROBOSAMPLE  #
################


cd ${ROBOSAMPLE_HOME}
echo "Current directory is : $(pwd)"

git checkout ${ROBOSAMPLE_GIT_BRANCH}
if [ ${ROBOSAMPLE_GIT_PULL} -eq 1 ]
  then 
    git pull;
fi


cd ${BUILD_DIRNAME}

make -j$(nproc)


