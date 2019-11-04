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

find_package(Simbody 3.6 CONFIG REQUIRED) # you may have to set a path for config file location e.g.  modify the environment variable list specific for find_package..
INCLUDE_DIRECTORIES( ${Simbody_INCLUDE_DIR} )
LINK_DIRECTORIES( ${Simbody_LIB_DIR} )
MESSAGE("Simbody_INCLUDE_DIR: ${Simbody_INCLUDE_DIR}")
MESSAGE("Simbody_LIB_DIR: ${Simbody_LIB_DIR}")

set( Integrators_PATH
     "${CMAKE_CURRENT_LIST_DIR}/../simbody/molmodel_legacy/SimTKmath/Integrators/src"
     CACHE FILEPATH "Path to Integrator.h from Simbody..." )
INCLUDE_DIRECTORIES( ${Integrators_PATH} )
MESSAGE("Integrators_PATH: ${Integrators_PATH}")

set( MOL_MODEL_INCLUDE
     "${Simbody_INCLUDE_DIR}/.."
     CACHE FILEPATH "header include folder path for Molmodel" )
INCLUDE_DIRECTORIES( ${MOL_MODEL_INCLUDE} )
MESSAGE("MOL_MODEL_INCLUDE: ${MOL_MODEL_INCLUDE}")

#find_package(Eigen3 REQUIRED)
#include_directories( ${EIGEN3_INCLUDE_DIRS} )
#MESSAGE("EIGEN3_INCLUDE_DIRS: ${EIGEN3_INCLUDE_DIRS}")

#add_definitions(-DNUMPY=1)
#MESSAGE("ADD DEFINITION: DNUMPY=1")

#set(Python_ADDITIONAL_VERSIONS 2)
#FIND_PACKAGE(PythonInterp)
#FIND_PACKAGE(PythonLibs REQUIRED)
#include_directories( ${PYTHON_INCLUDE_DIRS} )
#link_directories( ${PYTHON_LIBRARIES} )
#MESSAGE("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
#MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

#set(mmtk_include
#    "${PYTHON_INCLUDE_DIRS}/MMTK"
#    CACHE FILEPATH "Adding header path to MMTK installed location")
#include_directories( ${mmtk_include} )
#MESSAGE("mmtk_include: ${mmtk_include}")

#execute_process(
#  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
#  OUTPUT_VARIABLE PYTHON_SITE_PATH
#  OUTPUT_STRIP_TRAILING_WHITESPACE )

#set(PYTHON_SITE-PACKAGES_PATH
#    ${PYTHON_SITE_PATH}
#    CACHE FILEPATH "side-packages python path")
#MESSAGE("PYTHON_SITE_PATH: ${PYTHON_SITE_PATH}")

#set(mmtk_libs
#    "${PYTHON_SITE-PACKAGES_PATH}/MMTK/linux2"
#    CACHE FILEPATH "Addinng path to MMTK libs folder")
#link_directories( ${mmtk_libs} )
#MESSAGE("mmtk_libs: ${mmtk_libs}")

#execute_process(
#  COMMAND ls "${mmtk_libs}"
#  OUTPUT_VARIABLE mmtk_LIBRARIES
#  OUTPUT_STRIP_TRAILING_WHITESPACE )

#set(numpy_include
#    "${PYTHON_SITE-PACKAGES_PATH}/numpy/core/include"
#    CACHE FILEPATH "Header path to numpy installed location")
#include_directories( ${numpy_include} )
#MESSAGE("numpy_include: ${numpy_include}")

#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

#find_package(Boost 1.55 REQUIRED) #I'm using 1.58
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
#LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
#MESSAGE("Boost_FOUND: ${Boost_FOUND}")
#MESSAGE("Boost_VERSION: ${Boost_VERSION}")
#MESSAGE("Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
#MESSAGE("Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")

#SET(BOOST_INCLUDE_DIR_BYHAND "/usr/include/boost/")
#INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR_BYHAND})
#MESSAGE("BOOST_INCLUDE_DIR_BYHAND: ${BOOST_INCLUDE_DIR_BYHAND}")

#LINK_DIRECTORIES("/usr/lib/x86_64-linux-gnu/")
#file(GLOB BOOST_LIBS_BYHAND "/usr/lib/x86_64-linux-gnu/*boost*")
#MESSAGE("BOOST_LIBS_BYHAND: ${BOOST_LIBS_BYHAND}")

INCLUDE_DIRECTORIES(../include/gmolmodel)
LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/../build-debug/lib)
MESSAGE("Gmolmodel libs dir: ${CMAKE_CURRENT_LIST_DIR}/../build-debug/lib")

set(TESTS_SRCS
        TestSimbodyBasic.cpp
	#TestMassSqrtAndDet.cpp
        TestMolmodelBasic.cpp
        TestOpenMMBasic.cpp
        TestParallelExecutor.cpp
        TestSampler.cpp
        TestAmberInput.cpp
        TestWriteDCD.cpp
        Robosample.cpp
        TestTopology.cpp
        TestSetupReader.cpp
        TestBallLongPendulum.cpp
        TestMatrix.cpp
        TestWorldReconstruction.cpp
	TestFixmanTorqueSim.cpp
	)

foreach(TEST_SRC ${TESTS_SRCS})
    STRING(REGEX REPLACE ".cpp" "" TEST_EXE ${TEST_SRC})
    MESSAGE(${TEST_EXE} ${TEST_SRC})
    add_executable(${TEST_EXE} ${TEST_SRC})
    TARGET_LINK_LIBRARIES( ${TEST_EXE}
                ${Simbody_LIBRARIES}
                SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
                ${Boost_LIBRARIES}
                ${BOOST_LIBS_BYHAND}
                #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense...
		GMOLMODEL_dynamic
                m )
endforeach()

MESSAGE("Simbody_LIBRARIES: ${Simbody_LIBRARIES}")
MESSAGE("Boost_LIBRARIES: ${Boost_LIBRARIES}")
#MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
