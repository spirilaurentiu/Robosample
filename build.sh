#!/bin/bash

#########################
# The command line help #
#########################
display_help() {
    echo "Usage: $0 -b [BUILD_CONFIG] [OPTIONS...]" >&2
    echo
	echo "   -h, --help                 Show help."
    echo "   -b, --build                debug/release (case insensitive)"
	echo "   -n, --nproc                How many processors to use when compiling (default is \$(nproc)*2)"
	echo "                              If not enough resources are available, compilation may fail."
	echo "                              When this is the case, use lower numbers (\$(nproc), 1, 2, ...)"
	echo "   -u, --ubsan                Use undefined behaviour sanitizer. If not specified, default to not using."
	echo "                              On by default with debug configuration (this argument is ignored)."
    echo
    exit 1
}

#####################
# Default arguments #
#####################
CPU=$((`nproc`*2))
ROBOSAMPLE_UBSAN="No"

################################
# Check command line arguments #
################################
while :
do
	case "$1" in
		-h | --help)
			display_help # show help
			exit 0
			;;
		-b | --build)
			BUILD_TYPE="$2"
			BUILD_TYPE=$(echo "$BUILD_TYPE" | awk '{print tolower($0)}')

			if [[ "$BUILD_TYPE" != "debug" && "$BUILD_TYPE" != "release" ]]; then
				echo "Error: Build configuration must be \"debug\" or \"release\" (case insensitive)."
				echo "Error: Build configuration was ${BUILD_TYPE}"
				echo
				display_help

				exit 1
			fi

			shift 2
			;;

		-n | --nproc)
			CPU="$2"

			re='^[0-9]+$'
			if ! [[ $CPU =~ $re ]] ; then
				echo "Error: nproc must be an integer."
				echo "Error: nproc was $CPU"
				echo
				display_help

				exit 1
			fi

			shift 2
			;;

		-n | --ubsan)
			ROBOSAMPLE_UBSAN="Yes"
			shift 1
			;;

		--) # End of all options
			shift
			break
			;;
		-*)
			echo "Error: Unknown option: $1" >&2
			display_help # show help
			exit 1 
			;;
		*)  # No more options
			break
			;;
	esac
done

if [ "$BUILD_TYPE" == "" ]; then
	echo "Error: No build (-b/--build) configuration specified."
	echo
	display_help

	exit 1
fi



######################################################################
# CHOOSE THE FOLLOWING VARIABLES ACCORDING TO YOUR SYSTEM & DESIRES
######################################################################

# # 0/1 if you want to checkout to a specific branch/commit and pull updates 
# # to the local repo

# export SIMBODY_GIT_PULL=0
# export SIMBODY_GIT_BRANCH="master"

# export MOLMODEL_GIT_PULL=0
# export MOLMODEL_GIT_BRANCH="master"

# export OPENMM_GIT_PULL=0
# export OPENMM_GIT_BRANCH="master"

# export ROBOSAMPLE_GIT_PULL=0
# export ROBOSAMPLE_GIT_BRANCH="master"


# Robosample directories
ROBOSAMPLE_HOME="$(pwd)"
BUILD_DIR="$(pwd)/build/${BUILD_TYPE}"
INSTALL_DIR="$(pwd)/install/${BUILD_TYPE}"

mkdir -p ${BUILD_DIR}
mkdir -p ${INSTALL_DIR}

###############################
# Write compilation warnings. #
###############################
if [ "$BUILD_TYPE" == "debug" ]; then
	out=/dev/stderr
else
	out=${BUILD_DIR}/out.txt
fi


######################################################################


################
#   SIMBODY    #
################

SIMBODY_SRC="${ROBOSAMPLE_HOME}/Molmodel/Simbody01/"
SIMBODY_BUILD_DIR="${BUILD_DIR}/Molmodel/Simbody01/"
SIMBODY_INSTALL_DIR="${INSTALL_DIR}/Molmodel/Simbody01/"

mkdir -p ${SIMBODY_BUILD_DIR}
mkdir -p ${SIMBODY_INSTALL_DIR}

# if ${SIMBODY_GIT_PULL}; then
# 	cd ${SIMBODY_SRC}
# 	git checkout ${SIMBODY_GIT_BRANCH} && git pull;
# fi;

cd ${SIMBODY_BUILD_DIR}
rm -rf *

cmake \
	${SIMBODY_SRC}
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_INSTALL_PREFIX=${SIMBODY_INSTALL_DIR} ..

make -j${CPU} 2> ${out}
make install



################
#   OPENMM     #
################

OPENMM_SRC="${ROBOSAMPLE_HOME}/openmm/"
OPENMM_BUILD_DIR="${BUILD_DIR}/openmm/"
OPENMM_INSTALL_DIR="${INSTALL_DIR}/openmm/"

mkdir -p ${OPENMM_BUILD_DIR}
mkdir -p ${OPENMM_INSTALL_DIR}

# if [ ${OPENMM_GIT_PULL} -eq 1 ]; then
# 	cd ${OPENMM_SRC}
# 	git checkout ${OPENMM_GIT_BRANCH} && git pull;
# fi;

cd ${OPENMM_BUILD_DIR}
rm -rf *

cmake \
	${OPENMM_SRC}
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_INSTALL_PREFIX=${OPENMM_INSTALL_DIR} ..

make -j${CPU} 2> ${out}
make install



################
#   MOLMODEL   #
################

MOLMODEL_SRC="${ROBOSAMPLE_HOME}/Molmodel/"
MOLMODEL_BUILD_DIR="${BUILD_DIR}/Molmodel/"
MOLMODEL_INSTALL_DIR="${INSTALL_DIR}/Molmodel/"

mkdir -p ${MOLMODEL_BUILD_DIR}
mkdir -p ${MOLMODEL_INSTALL_DIR}

# if [ ${MOLMODEL_GIT_PULL} -eq 1 ]; then
# 	cd ${MOLMODEL_SRC}
# 	git checkout ${MOLMODEL_GIT_BRANCH} && git pull;
# fi;

cd ${MOLMODEL_BUILD_DIR}
rm -rf *

cmake \
	${MOLMODEL_SRC}
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DSimTK_INSTALL_PREFIX=${MOLMODEL_INSTALL_DIR} \
	-DSimbody_DIR=${SIMBODY_INSTALL_DIR}/lib/cmake/simbody/ \
	-DOpenMM_INCLUDE_DIR=${OPENMM_INSTALL_DIR}/include \
	-DOpenMM_LIBRARIES=${OPENMM_INSTALL_DIR}/lib/libOpenMM.so \
	-DOpenMM_LIBRARY=${OPENMM_INSTALL_DIR}/lib/libOpenMM.so ..

make -j${CPU} 2> ${out}
make install



################
#  ROBOSAMPLE  #
################

ROBOSAMPLE_SRC="${ROBOSAMPLE_HOME}"
ROBOSAMPLE_BUILD_DIR="${BUILD_DIR}/robosample/"
ROBOSAMPLE_INSTALL_DIR="${INSTALL_DIR}/robosample/"

mkdir -p ${ROBOSAMPLE_BUILD_DIR}
mkdir -p ${ROBOSAMPLE_INSTALL_DIR}

# if [ ${ROBOSAMPLE_GIT_PULL} -eq 1 ]; then
# 	cd ${ROBOSAMPLE_SRC}
# 	git checkout ${ROBOSAMPLE_GIT_BRANCH} && git pull;
# fi;

cd ${ROBOSAMPLE_BUILD_DIR}
echo $(pwd)
rm -rf *

cmake \
	${ROBOSAMPLE_SRC}
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DROBO_UBSAN=${ROBOSAMPLE_UBSAN} \
	-DCMAKE_INSTALL_PREFIX=${ROBOSAMPLE_INSTALL_DIR} \
	-DSimbody_PATH=${SIMBODY_INSTALL_DIR} \
	-DMolModel_PATH=${MOLMODEL_INSTALL_DIR} \
	-DOpenMM_PATH=${OPENMM_INSTALL_DIR} ..


make -j${CPU} 2> ${out}
make install


##############################
#  Examples, tests and misc  #
##############################

# add test input files
cp -ri ${ROBOSAMPLE_HOME}/tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# add tests
cp ${ROBOSAMPLE_HOME}/tests/test-memory.sh test-memory.sh

# tell the user that we are done
echo
echo
echo
if [ "$BUILD_TYPE" == "debug" ]; then
	echo "Build Robosample and all submodules for debug configuration."
	echo "Build directory is ${ROBOSAMPLE_BUILD_DIR}."
else
	echo "Build Robosample and all submodules for release configuration."
	echo "Build directory is ${ROBOSAMPLE_BUILD_DIR}."
	echo "Any warnings or errors are written in ${out}."
fi
