#!/bin/bash

#########################
# The command line help #
#########################
display_help() {
    echo "Usage: $0 -b [BUILD_CONFIG] [C_COMPILER] [CPP_COMPILER] [OPTIONS...]" >&2
    echo
	echo "   -h, --help                 Show help."
    echo "   -b, --build                debug/release (case insensitive)"
	echo "   --cc                       C compiler."
	echo "   --cpp                      C++ compiler."
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

		--cc)
			C_COMPILER="$2"

			shift 2
			;;

		--cpp)
			CPP_COMPILER="$2"

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

if [ "${BUILD_TYPE}" == "" ]; then
	echo "Error: No build (-b/--build) configuration specified."
	echo
	display_help

	exit 1
fi

if [ "${C_COMPILER}" == "" ]; then
	echo "Error: No C compiler (--cc) specified."
	echo
	display_help

	exit 1
fi

if [ "${CPP_COMPILER}" == "" ]; then
	echo "Error: No C++ compiler (--cpp) specified."
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


######################################################################



# Robosample directories
ROBOSAMPLE_HOME="$(pwd)"
BUILD_DIR="${ROBOSAMPLE_HOME}/build/${BUILD_TYPE}"
INSTALL_DIR="${ROBOSAMPLE_HOME}/install/${BUILD_TYPE}"

mkdir -p ${BUILD_DIR}
mkdir -p ${INSTALL_DIR}

# create openmm paths and directories
OPENMM_SRC="${ROBOSAMPLE_HOME}/openmm/"
OPENMM_BUILD_DIR="${BUILD_DIR}/openmm/"
OPENMM_INSTALL_DIR="${INSTALL_DIR}/openmm/"

mkdir -p ${OPENMM_BUILD_DIR}
mkdir -p ${OPENMM_INSTALL_DIR}

# create simbody paths and directories
SIMBODY_SRC="${ROBOSAMPLE_HOME}/Molmodel/Simbody01/"
SIMBODY_BUILD_DIR="${BUILD_DIR}/Simbody01/"
SIMBODY_INSTALL_DIR="${INSTALL_DIR}/Simbody01/"

mkdir -p ${SIMBODY_BUILD_DIR}
mkdir -p ${SIMBODY_INSTALL_DIR}

# create molmodel paths and directories
# molmodel gets installed automatically in ${SIMBODY_INSTALL_DIR} 
MOLMODEL_SRC="${ROBOSAMPLE_HOME}/Molmodel/"
MOLMODEL_BUILD_DIR="${BUILD_DIR}/Molmodel/"

mkdir -p ${MOLMODEL_BUILD_DIR}

# create robosample paths and directories
# there is nothing to install here
ROBOSAMPLE_SRC="${ROBOSAMPLE_HOME}"
ROBOSAMPLE_BUILD_DIR="${BUILD_DIR}/robosample/"

mkdir -p ${ROBOSAMPLE_BUILD_DIR}



# On release builds, write warnings and errors to a file. Also, emit relocs for PGO.
if [ "${BUILD_TYPE}" == "debug" ]; then
	out=/dev/stderr
else
	# out file for compilation errors and warnings
	time_now=`date +"%b-%d-%Y-%H_%M_%S"`
	out=${BUILD_DIR}/out_${time_now}.txt
	touch ${out}

	# emit relocs for PGO
	LDFLAGS_TEMP=${LDFLAGS}
	export LDFLAGS="-Wl,-q"
fi


################
#   OPENMM     #
################

# if [ ${OPENMM_GIT_PULL} -eq 1 ]; then
# 	cd ${OPENMM_SRC}
# 	git checkout ${OPENMM_GIT_BRANCH} && git pull;
# fi;

cd ${OPENMM_INSTALL_DIR}
rm -rf *

cd ${OPENMM_BUILD_DIR}
rm -rf *

cmake \
	${OPENMM_SRC} \
	-DCMAKE_C_COMPILER=${C_COMPILER} \
	-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
	-DOPENMM_INSTALL_PREFIX=${OPENMM_INSTALL_DIR} \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_PREFIX_PATH=${OPENMM_INSTALL_DIR} \
	-DCMAKE_INSTALL_RPATH=${OPENMM_INSTALL_DIR}lib \
	-DCMAKE_INSTALL_PREFIX=${OPENMM_INSTALL_DIR}

make -j${CPU} 2>> ${out}
make install



################
#   SIMBODY    #
################

# if ${SIMBODY_GIT_PULL}; then
# 	cd ${SIMBODY_SRC}
# 	git checkout ${SIMBODY_GIT_BRANCH} && git pull;
# fi;

cd ${SIMBODY_INSTALL_DIR}
rm -rf *

cd ${SIMBODY_BUILD_DIR}
rm -rf *

cmake \
	${SIMBODY_SRC} \
	-DCMAKE_C_COMPILER=${C_COMPILER} \
	-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_INSTALL_LIBDIR="lib" \
	-DCMAKE_INSTALL_FULL_LIBDIR=${SIMBODY_INSTALL_DIR}lib \
	-DCMAKE_INSTALL_RPATH=${SIMBODY_INSTALL_DIR}lib \
	-DCMAKE_INSTALL_PREFIX=${SIMBODY_INSTALL_DIR}

make -j${CPU} 2>> ${out}
make install



# ################
# #   MOLMODEL   #
# ################

# if [ ${MOLMODEL_GIT_PULL} -eq 1 ]; then
# 	cd ${MOLMODEL_SRC}
# 	git checkout ${MOLMODEL_GIT_BRANCH} && git pull;
# fi;

cd ${MOLMODEL_BUILD_DIR}
rm -rf *

cmake \
	${MOLMODEL_SRC} \
	-DCMAKE_C_COMPILER=${C_COMPILER} \
	-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DROBO_OPENMM_SEARCH_PATH=${OPENMM_INSTALL_DIR} \
	-DSimTK_INSTALL_PREFIX=${SIMBODY_INSTALL_DIR} \
	-DCMAKE_INSTALL_RPATH="${OPENMM_INSTALL_DIR}lib;${SIMBODY_INSTALL_DIR}lib"
	
make -j${CPU} 2>> ${out}
make install



################
#  ROBOSAMPLE  #
################

# if [ ${ROBOSAMPLE_GIT_PULL} -eq 1 ]; then
# 	cd ${ROBOSAMPLE_SRC}
# 	git checkout ${ROBOSAMPLE_GIT_BRANCH} && git pull;
# fi;

cd ${ROBOSAMPLE_BUILD_DIR}
rm * -rf

cmake \
	${ROBOSAMPLE_SRC} \
	-DCMAKE_C_COMPILER=${C_COMPILER} \
	-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DROBO_UBSAN=${ROBOSAMPLE_UBSAN} \
	-DROBO_OPENMM_SEARCH_PATH=${OPENMM_INSTALL_DIR} \
	-DSimTK_INSTALL_DIR=${SIMBODY_INSTALL_DIR} \
	-DCMAKE_INSTALL_RPATH="${OPENMM_INSTALL_DIR}lib;${SIMBODY_INSTALL_DIR}lib"


make -j${CPU} 2>> ${out}



# On release builds, write warnings and errors to a file. Also, emit relocs for PGO.
if [ "${BUILD_TYPE}" == "release" ]; then
	LDFLAGS=${LDFLAGS_TEMP}
fi



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
