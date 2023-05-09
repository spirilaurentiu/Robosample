#!/bin/bash

#########################
# The command line help #
#########################
display_help() {
    echo "Usage: $0 -b [BUILD_CONFIG] --cc [C_COMPILER] --cpp [CPP_COMPILER] [OPTIONS...]" >&2
    echo
	echo "   -h, --help                 Show help."
    echo "   -b, --build                \"debug\"/\"release\" (case insensitive)"
	echo "   --cc                       C compiler."
	echo "   --cpp                      C++ compiler."
	echo "   --pgo                      Use PGO or not. Only applicable for release builds (ignored on debug builds)."
	echo "                              By default, release builds do not use PGO."
	echo "   -n, --njobs                How many processors to use when compiling (default is \$(nproc)*2)"
	echo "                              If not enough resources are available, compilation may fail."
	echo "                              When this is the case, use lower numbers (\$(nproc), 1, 2, ...)"
	echo "   -u, --ubsan                Use undefined behaviour sanitizer. If not specified, default to not using."
	echo "                              On by default with debug configuration (this argument is ignored)."
	echo "   --verbose                  Verbose compilation output. Either \"on\" or \"off\" (case insensitive)."
	echo "                              Off by default."
    echo
    exit 1
}

#####################
# Default arguments #
#####################
ROBO_PGO=false
ROBO_CPU=$((`nproc`*2))
ROBO_UBSAN="No"
ROBO_VERBOSE=false

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
			
			# lowercase all
			BUILD_TYPE=$(echo "$BUILD_TYPE" | awk '{print tolower($0)}')
			
			# capitalize first letter
			# this capitalization is required by all cmake configurations
			BUILD_TYPE="$(tr '[:lower:]' '[:upper:]' <<< ${BUILD_TYPE:0:1})${BUILD_TYPE:1}"

			if [[ "$BUILD_TYPE" != "Debug" && "$BUILD_TYPE" != "Release" ]]; then
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

		--pgo)
			ROBO_PGO=true
			shift 1
			;;

		-n | --njobs)
			ROBO_CPU="$2"

			re='^[0-9]+$'
			if ! [[ $ROBO_CPU =~ $re ]] ; then
				echo "Error: njobs must be an integer."
				echo "Error: njobs was $ROBO_CPU"
				echo
				display_help

				exit 1
			fi

			shift 2
			;;

		-n | --ubsan)
			ROBO_UBSAN="Yes"
			shift 1
			;;

		--verbose)
			ROBO_VERBOSE=true
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


# set build configuration dependant flags
if [ "${BUILD_TYPE}" == "Debug" ]; then
	# on debug builds, write everything to the screen
	out=/dev/stderr

	# do not use lto for debug builds
	ROBO_LTO=OFF

	# do not use pgo for debug builds
	ROBO_PGO_FLAGS_GENERATE=""
	ROBO_PGO_FLAGS_USE=""

	echo '' >> ~/.bashrc
	echo "export OPENMM_PLUGIN_DIR_DEBUG=${MOLMODEL_BUILD_DIR}" >> ~/.bashrc
	echo "export OPENMM_PLUGIN_PLATFORMS_DEBUG=${OPENMM_INSTALL_DIR}lib/plugins/" >> ~/.bashrc
else
	# get a lower-case string containing current date and time
	time_now=`date +"%b-%d-%Y-%H_%M_%S"`
	time_now=$(echo "$time_now" | awk '{print tolower($0)}')

	# on release builds, write warnings and errors to a file
	out=${BUILD_DIR}/out_${time_now}.txt
	touch ${out}

	# turn on lto
	ROBO_LTO=ON

	# profile guided optimization (pgo) compilation flags
	if ${ROBO_PGO}; then
		# -fprofile-generate enables -fprofile-arcs, -fprofile-values and -fvpt
		# -fprofile-use enables -fbranch-probabilities, -fvpt, -funroll-loops, -fpeel-loops and -ftracer
		ROBO_PGO_DIR="/tmp/robosample-pgo-${time_now}"
		ROBO_PGO_FLAGS_GENERATE="-fprofile-dir=${ROBO_PGO_DIR} -fprofile-generate=${ROBO_PGO_DIR}"
		ROBO_PGO_FLAGS_USE="-fprofile-dir=${ROBO_PGO_DIR} -fprofile-use=${ROBO_PGO_DIR} -fprofile-correction"

		# make the folder for pgo files
		mkdir -p ${ROBO_PGO_DIR}
		cd ${ROBO_PGO_DIR}
	else
		ROBO_PGO_FLAGS_GENERATE=""
		ROBO_PGO_FLAGS_USE=""
	fi

	echo '' >> ~/.bashrc
	echo "export OPENMM_PLUGIN_DIR_RELEASE=${MOLMODEL_BUILD_DIR}" >> ~/.bashrc
	echo "export OPENMM_PLUGIN_PLATFORMS_RELEASE=${OPENMM_INSTALL_DIR}lib/plugins/" >> ~/.bashrc
fi

source ~/.bashrc



# warning flags to be applied in all projects
ROBO_WARNINGS="-Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion -Wformat=2 -Wformat-truncation -Wundef -Wconversion -Wunused-parameter"
ROBO_WARNINGS=""

# sanitizers to use
ROBO_SANITIZER="-fsanitize=undefined -fsanitize=address -fsanitize-recover=address"

# instruction set
ROBO_BUILD_INSET="-march=native -mtune=native"

# GCC -flto-partition=one slows down the entire thing
# GCC -flto makes it go faster - done in cmake
# GCC -fwhole-program linker error on simbody
# GCC -fprefetch-loop-arrays slows down everything (a bit at least)
# now, lto is handled by cmake (no need to manually request)
# 
# https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gcc/Optimize-Options.html
# -fopt-info-vec-all prints everything (might be worth investigating)
# for dev only -fopt-info-vec-missed
ROBO_EXTRA_FLAGS="-pipe -fno-math-errno -Wl,-q,-emit-relocs"
ROBO_EXTRA_FLAGS=""

# flags to be applied in all projects
# we want to override some cmake defaults because some flags might be better suited
ROBO_DEBUG_FLAGS="-O0 -g ${ROBO_SANITIZER} ${ROBO_WARNINGS}"
ROBO_RELEASE_FLAGS="-O3 -DNDEBUG ${ROBO_BUILD_INSET} ${ROBO_EXTRA_FLAGS}"
ROBO_RELWITHDEBINFO_FLAGS="-O3 -DNDEBUG -g ${ROBO_BUILD_INSET} ${ROBO_EXTRA_FLAGS}"
ROBO_MINSIZEREL_FLAGS="-Os -DNDEBUG ${ROBO_BUILD_INSET} ${ROBO_EXTRA_FLAGS}"

# build tests and examples
ROBO_BUILD_EXAMPLES=OFF
ROBO_BUILD_TESTS=OFF
ROBO_BUILD_VISUALIZER=OFF



# if we use compiler pgo, then we must compile the code twice
if ${ROBO_PGO}; then
	ROUNDS=2
else
	ROUNDS=1
fi



# start building process
for ((i=0;i<ROUNDS;i++)); do

	# mark the projects for profiling if doing the first compilation round
	# compile using gathered data otherwise
	if ${ROBO_PGO}; then
		if [[ $i == 0 ]]; then
			ROBO_PGO_FLAGS=${ROBO_PGO_FLAGS_GENERATE}
		else
			ROBO_PGO_FLAGS=${ROBO_PGO_FLAGS_USE}
		fi
	fi

    # openmm
	cd ${OPENMM_INSTALL_DIR}
	rm -rf *

	cd ${OPENMM_BUILD_DIR}
	rm -rf *
	
	cmake \
		${OPENMM_SRC} \
		-DCMAKE_C_COMPILER=${C_COMPILER} \
		-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
		-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=${ROBO_LTO} \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_CXX_STANDARD_REQUIRED=Yes \
		-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
		-DCMAKE_CXX_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_C_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_PREFIX_PATH=${OPENMM_INSTALL_DIR} \
		-DCMAKE_INSTALL_RPATH=${OPENMM_INSTALL_DIR}lib \
		-DCMAKE_INSTALL_PREFIX=${OPENMM_INSTALL_DIR} \
		-DOPENMM_INSTALL_PREFIX=${OPENMM_INSTALL_DIR} \
		-DOPENMM_BUILD_EXAMPLES=${ROBO_BUILD_EXAMPLES} \
		-DBUILD_TESTING=${ROBO_BUILD_TESTS}

	# compile openmm
	if $ROBO_VERBOSE; then
		make VERBOSE=1 -j${ROBO_CPU} 2>> ${out}
	else
		make -j${ROBO_CPU} 2>> ${out}
	fi

	# install openmm
	# 'make install' does not work as we need to define some more variables in order to install locally
	cmake -DCMAKE_INSTALL_PREFIX=${OPENMM_INSTALL_DIR} -P cmake_install.cmake



	# simbody
	cd ${SIMBODY_INSTALL_DIR}
	rm -rf *

	cd ${SIMBODY_BUILD_DIR}
	rm -rf *

	cmake \
		${SIMBODY_SRC} \
		-DCMAKE_C_COMPILER=${C_COMPILER} \
		-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
		-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=${ROBO_LTO} \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_CXX_STANDARD_REQUIRED=Yes \
		-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
		-DCMAKE_CXX_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_C_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_INSTALL_LIBDIR="lib" \
		-DCMAKE_INSTALL_FULL_LIBDIR=${SIMBODY_INSTALL_DIR}lib \
		-DCMAKE_INSTALL_RPATH=${SIMBODY_INSTALL_DIR}lib \
		-DCMAKE_INSTALL_PREFIX=${SIMBODY_INSTALL_DIR} \
		-DBUILD_EXAMPLES=${ROBO_BUILD_EXAMPLES} \
		-DBUILD_TESTING=${ROBO_BUILD_TESTS} \
		-DBUILD_TESTS_AND_EXAMPLES_SHARED==${ROBO_BUILD_TESTS} \
		-DBUILD_TESTS_AND_EXAMPLES_STATIC=${ROBO_BUILD_TESTS} \
		-DBUILD_VISUALIZER=${ROBO_BUILD_VISUALIZER}

	# compile simbody
	if $ROBO_VERBOSE; then
		make VERBOSE=1 -j${ROBO_CPU} 2>> ${out}
	else
		make -j${ROBO_CPU} 2>> ${out}
	fi

	# install simbody
	make install



	# molmodel
	cd ${MOLMODEL_BUILD_DIR}
	rm -rf *

	cmake \
		${MOLMODEL_SRC} \
		-DCMAKE_C_COMPILER=${C_COMPILER} \
		-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
		-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=${ROBO_LTO} \
		-DCMAKE_POLICY_DEFAULT_CMP0069=NEW \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_CXX_STANDARD_REQUIRED=Yes \
		-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
		-DCMAKE_CXX_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_C_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DROBO_OPENMM_SEARCH_PATH=${OPENMM_INSTALL_DIR} \
		-DSimTK_INSTALL_PREFIX=${SIMBODY_INSTALL_DIR} \
		-DCMAKE_INSTALL_RPATH="${OPENMM_INSTALL_DIR}lib;${SIMBODY_INSTALL_DIR}lib" \
		-DBUILD_TESTING=${ROBO_BUILD_TESTS} \
		-DBUILD_TESTING_STATIC=${ROBO_BUILD_TESTS} \
		-DBUILD_TESTING_SHARED=${ROBO_BUILD_TESTS} \
		-DBUILD_EXAMPLES=${ROBO_BUILD_EXAMPLES}
		
	# compile molmodel
	if $ROBO_VERBOSE; then
		make VERBOSE=1 -j${ROBO_CPU} 2>> ${out}
	else
		make -j${ROBO_CPU} 2>> ${out}
	fi

	# install molmodel
	make install



	# robosample
	cd ${ROBOSAMPLE_BUILD_DIR}
	rm * -rf

	cmake \
		${ROBOSAMPLE_SRC} \
		-DCMAKE_C_COMPILER=${C_COMPILER} \
		-DCMAKE_CXX_COMPILER=${CPP_COMPILER} \
		-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=${ROBO_LTO} \
		-DCMAKE_POLICY_DEFAULT_CMP0069=NEW \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_CXX_STANDARD_REQUIRED=Yes \
		-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
		-DCMAKE_CXX_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_CXX_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_DEBUG="${ROBO_DEBUG_FLAGS}" \
		-DCMAKE_C_FLAGS_RELEASE="${ROBO_RELEASE_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_RELWITHDEBINFO="${ROBO_RELWITHDEBINFO_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DCMAKE_C_FLAGS_MINSIZEREL="${ROBO_MINSIZEREL_FLAGS} ${ROBO_PGO_FLAGS}" \
		-DROBO_UBSAN=${ROBO_UBSAN} \
		-DROBO_OPENMM_SEARCH_PATH=${OPENMM_INSTALL_DIR} \
		-DSimTK_INSTALL_DIR=${SIMBODY_INSTALL_DIR} \
		-DCMAKE_INSTALL_RPATH="${OPENMM_INSTALL_DIR}lib;${SIMBODY_INSTALL_DIR}lib"

	# compile robosample
	if $ROBO_VERBOSE; then
		make VERBOSE=1 -j${ROBO_CPU} 2>> ${out}
	else
		make -j${ROBO_CPU} 2>> ${out}
	fi



	# OpenMM plugin goes here
	echo "${SIMBODY_INSTALL_DIR}lib/plugins/" > "${ROBOSAMPLE_BUILD_DIR}src/openmmplugin"

	# add test input files
	cp -ri ${ROBOSAMPLE_HOME}/tests_inputs/* .
	mkdir -p temp
	mkdir -p temp/pdbs

	# add tests
	cp ${ROBOSAMPLE_HOME}/tests/test-memory.sh test-memory.sh



	# profile if using pgo and only after the first compilation round
	# this is to collect profile data and use it in the second round 
	if ${ROBO_PGO}; then
		if [[ $i == 0 ]]; then
			./src/GMOLMODEL_robo inp.aper.pgo
		fi
	fi
done

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
	echo "Any warnings or errors have been written to ${out}."
fi
