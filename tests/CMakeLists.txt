cmake_minimum_required(VERSION 3.0)

# Permit use of custom FindOpenMM module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

MESSAGE(" CMAKE_SOURCE_DIR -> ${CMAKE_SOURCE_DIR}")

set(OpenMM_ON OFF)
find_package(OpenMM)
if(OpenMM_FOUND)
	set(OpenMM_ON ON)
endif()

if ( OpenMM_ON )
	find_package(OpenMM REQUIRED)
	INCLUDE_DIRECTORIES(${OpenMM_INCLUDE_DIR})
	INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
	MESSAGE("Include OpenMM_INCLUDE_DIR: ${OpenMM_INCLUDE_DIR} ")
else( OpenMM_ON )
	MESSAGE( FATAL_ERROR "OpenMM_ON FALSE")
endif( OpenMM_ON )

# @TODO do we need 3.6 CONFIG?
find_package(Simbody REQUIRED) # you may have to set a path for config file location e.g.  modify the environment variable list specific for find_package..
INCLUDE_DIRECTORIES( ${Simbody_INCLUDE_DIR} )
LINK_DIRECTORIES( ${Simbody_LIB_DIR} )
MESSAGE("Simbody_INCLUDE_DIR: ${Simbody_INCLUDE_DIR}")
MESSAGE("Simbody_LIB_DIR: ${Simbody_LIB_DIR}")

#set( Integrators_PATH
#     "${CMAKE_CURRENT_LIST_DIR}/../simbody/molmodel_legacy/SimTKmath/Integrators/src"
#     CACHE FILEPATH "Path to Integrator.h from Simbody..." )
#INCLUDE_DIRECTORIES( ${Integrators_PATH} )
#MESSAGE("Integrators_PATH: ${Integrators_PATH}")

set( MOL_MODEL_INCLUDE
     "${Simbody_INCLUDE_DIR}/.."
     CACHE FILEPATH "header include folder path for Molmodel" )
INCLUDE_DIRECTORIES( ${MOL_MODEL_INCLUDE} )
MESSAGE("MOL_MODEL_INCLUDE: ${MOL_MODEL_INCLUDE}")

#find_package(Eigen3 REQUIRED)
#include_directories( ${EIGEN3_INCLUDE_DIRS} )
#MESSAGE("EIGEN3_INCLUDE_DIRS: ${EIGEN3_INCLUDE_DIRS}")


INCLUDE_DIRECTORIES(../include/gmolmodel)
LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/../build-debug/lib)
MESSAGE("Gmolmodel libs dir: ${CMAKE_CURRENT_LIST_DIR}/../build-debug/lib")

set(TESTS_SRCS
	TestSimbodyBasic.cpp
	TestMassSqrtAndDet.cpp
	TestMolmodelBasic.cpp
	TestOpenMMBasic.cpp
	TestParallelExecutor.cpp
	TestSampler.cpp
	TestAmberInput.cpp
	TestWriteDCD.cpp
	TestTopology.cpp
	TestSetupReader.cpp
	TestBallLongPendulum.cpp
	TestMatrix.cpp
	TestWorldReconstruction.cpp
	TestFixmanTorqueSim.cpp
	)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
	foreach(TEST_SRC ${TESTS_SRCS})
		STRING(REGEX REPLACE ".cpp" "" TEST_EXE ${TEST_SRC})
		MESSAGE(${TEST_EXE} ${TEST_SRC})
		add_executable(${TEST_EXE} ${TEST_SRC})
		TARGET_LINK_LIBRARIES( ${TEST_EXE}
			debug ${Simbody_LIBRARIES}
			SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
			# ${Boost_LIBRARIES}
			# ${BOOST_LIBS_BYHAND}
			GMOLMODEL_dynamic_d
			m )
	endforeach()
ELSE()
	foreach(TEST_SRC ${TESTS_SRCS})
		STRING(REGEX REPLACE ".cpp" "" TEST_EXE ${TEST_SRC})
		MESSAGE(${TEST_EXE} ${TEST_SRC})
		add_executable(${TEST_EXE} ${TEST_SRC})
		TARGET_LINK_LIBRARIES( ${TEST_EXE}
			optimized ${Simbody_LIBRARIES}
			SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
			# ${Boost_LIBRARIES}
			# ${BOOST_LIBS_BYHAND}
			GMOLMODEL_dynamic
			m )
	endforeach()
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

# MESSAGE("Simbody_LIBRARIES: ${Simbody_LIBRARIES}")
# MESSAGE("Boost_LIBRARIES: ${Boost_LIBRARIES}")
#MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
