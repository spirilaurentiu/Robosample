cmake_minimum_required(VERSION 3.0)

# Permit use of custom FindOpenMM module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "Control volume of build output" )

set(OpenMM_ON OFF)
find_package(OpenMM)
if(OpenMM_FOUND)
    	set(OpenMM_ON ON)
        message("OpenMM found")
else()
        message("OpenMM not found")
endif()

if ( OpenMM_ON )
	find_package(OpenMM REQUIRED)
	INCLUDE_DIRECTORIES(${OpenMM_INCLUDE_DIR})
	INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
else( OpenMM_ON )

endif( OpenMM_ON )

find_package(Simbody 3.6 CONFIG REQUIRED) # you may have to set a path for config file location e.g.  modify the environment variable list specific for find_package..
INCLUDE_DIRECTORIES( ${Simbody_INCLUDE_DIR} )
LINK_DIRECTORIES( ${Simbody_LIB_DIR} )


#INCLUDE_DIRECTORIES("/home/nash_mit/git/refactoring/gmolmodel/simbody/SimTKmath/Integrators/src")
#this is a patch for the current structure of the project, it will be refactor after redesign of the gmodmodel integrator is done
set( Integrators_PATH
#     "${CMAKE_CURRENT_LIST_DIR}/../simbody/SimTKmath/Integrators/src"
     "${CMAKE_CURRENT_LIST_DIR}/../Molmodel/Simbody01/SimTKmath/Integrators/src"
     CACHE FILEPATH "Path to Integrator.h from Simbody..." )
INCLUDE_DIRECTORIES( ${Integrators_PATH} )


set( MOL_MODEL_INCLUDE
     "${Simbody_INCLUDE_DIR}/.."
     CACHE FILEPATH "header include folder path for Molmodel" )
INCLUDE_DIRECTORIES( ${MOL_MODEL_INCLUDE} )
#I should add link_directories for molmodel and add them to TARTGET_LINK!!!
#or maybe not.. at the moment I'm not convince that gmolmodel dependes by molmodel lib mostly because it's an extension that it's already compiled into gmolmode. but I need to look more into the code for that
#path for libs is the same as simbody's so, no need to add it!


#find_package(Eigen3 REQUIRED)
#find_package( PkgConfig )
#pkg_check_modules( EIGEN3 REQUIRED eigen3 )
#include_directories( ${EIGEN3_INCLUDE_DIRS} )


#this is done because /home/nash_mit/Anaconda/envs/mmtk_environment/include/python2.7/MMTK/arrayobject.h need to make a decison
#and we only have numpy/oldnumeric.h
#add_definitions(-DNUMPY=1)

#set(Python_ADDITIONAL_VERSIONS 2)
#FIND_PACKAGE(PythonInterp)
#FIND_PACKAGE(PythonLibs REQUIRED)
#include_directories( ${PYTHON_INCLUDE_DIRS} )
#link_directories( ${PYTHON_LIBRARIES} )

#set(mmtk_include
#    "${PYTHON_INCLUDE_DIRS}/MMTK"
#    CACHE FILEPATH "Adding header path to MMTK installed location")
#include_directories( ${mmtk_include} )

#execute_process(
#  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
#  OUTPUT_VARIABLE PYTHON_SITE_PATH
#  OUTPUT_STRIP_TRAILING_WHITESPACE )

#set(PYTHON_SITE-PACKAGES_PATH
#    ${PYTHON_SITE_PATH}
#    CACHE FILEPATH "side-packages python path")

#set(mmtk_libs "/home/nash_mit/git/refactoring/gmolmodel/mmtk_legacy/build/lib.linux-x86_64-2.6/MMTK/linux2")
#set(mmtk_libs
#    "${PYTHON_SITE-PACKAGES_PATH}/MMTK/linux2"
#    CACHE FILEPATH "Addinng path to MMTK libs folder")
#link_directories( ${mmtk_libs} )

#not finished yet... in the end it should export a list of libs without "lib" and without ".so"
#the resoult variable should be puted in TARGET_LINK...
#execute_process(
#  COMMAND ls "${mmtk_libs}"
#  OUTPUT_VARIABLE mmtk_LIBRARIES
#  OUTPUT_STRIP_TRAILING_WHITESPACE )
#message( ${mmtk_LIBRARIES} )

#set(numpy_include "/home/nash_mit/Anaconda/envs/mmtk_environment/lib/python2.7/site-packages/numpy/core/include" )
#set(numpy_include
#    "${PYTHON_SITE-PACKAGES_PATH}/numpy/core/include"
#    CACHE FILEPATH "Header path to numpy installed location")
#include_directories( ${numpy_include} )


#set (Boost_NO_SYSTEM_PATHS ON)
#set (Boost_USE_MULTITHREADED ON)
#set (Boost_USE_STATIC_LIBS ON)
#set (Boost_USE_STATIC_RUNTIME OFF)
#set (BOOST_ALL_DYN_LINK OFF)

find_package(Boost 1.55 REQUIRED) #I'm using 1.58
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
LINK_DIRECTORIES(${BOOST_LIBRARY_DIRS})


#it's for linking GMOLMODEL_static into GCHMC lib
LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/lib)

set(LIBGMOLMODEL_SRC
        TrivalentAtomTetra.cpp
        bSpecificAtom.cpp
        bBond.cpp
        Context.cpp # LAUR
        IState.cpp # LAUR
        MyState.cpp # LAUR
        Sampler.cpp # LAUR
        MonteCarloSampler.cpp # LAUR
        HMCSampler.cpp # LAUR
        LAHMCSampler.cpp # LAUR
        ConformationalSearch.cpp # LAUR
        GirolamiSampler.cpp # LAUR
	      #bAddParams.cpp
	      Topology.cpp
	      #bMoleculeReader.cpp
	      World.cpp
              FixmanTorque.cpp
        ParaMolecularDecorator.cpp
        readAmberInput.cpp
	      bgeneral.cpp
        SetupReader.cpp
        format/dcdplugin.c
        format/PDBObject.cpp
        format/TrajectoryObject.cpp
	)

include_directories(../include/gmolmodel)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
        add_library(GMOLMODEL_dynamic_d SHARED ${LIBGMOLMODEL_SRC})
        add_library(GMOLMODEL_static_d  STATIC ${LIBGMOLMODEL_SRC})
ELSE()
        add_library(GMOLMODEL_dynamic SHARED ${LIBGMOLMODEL_SRC})
        add_library(GMOLMODEL_static  STATIC ${LIBGMOLMODEL_SRC})
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

add_executable(GMOLMODEL_test ${LIBGMOLMODEL_SRC} testMain.cpp )

#simmain must be refactor since it contains mmtk data types and current version do not make use of mmtk anymore
#PYTHON_ADD_MODULE(
#		GCHMC
#		#${LIBGMOLMODEL_SRC}
#		simmain.cpp
#		bArgParser.cpp
#		)

#add_library( GCHMC SHARED ${LIBGMOLMODEL_SRC} simmain.cpp )

#since simmain.cpp contains mmtk data type this target it's not available until it's refactor
#TARGET_LINK_LIBRARIES( GCHMC
#		boost_python
#		GMOLMODEL_static
#                ${Simbody_LIBRARIES}
#                SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
#                ${BOOST_LIBRARIES}
#                ${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code..
#		${OpenMM_LIBRARIES}
#                m )
#                #${mmtk_libs} )
#

IF(CMAKE_BUILD_TYPE MATCHES Debug)
        message("Linking with debug libraries")

        TARGET_LINK_LIBRARIES( GMOLMODEL_test
                debug ${Simbody_LIBRARIES}
                SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
                ${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code..
		${OpenMM_LIBRARIES}
                m )
                #${mmtk_libs} )

        # link Simbody libraries
        TARGET_LINK_LIBRARIES(GMOLMODEL_dynamic_d
		debug ${Simbody_LIBRARIES}
		SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
		${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code...
		${OpenMM_LIBRARIES}
		m )
		#${mmtk_libs} )
                
        TARGET_LINK_LIBRARIES(GMOLMODEL_static_d
		debug ${Simbody_LIBRARIES}
		SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
		${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES}  # same as for dynamic lib... I don't think python shoud be part of this lib...
		${OpenMM_LIBRARIES}
		m )
                #${mmtk_libs} )
ELSE()
        message("Linking with release libraries")

        TARGET_LINK_LIBRARIES( GMOLMODEL_test
                optimized ${Simbody_LIBRARIES}
                SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
                ${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code..
		${OpenMM_LIBRARIES}
                m )
                #${mmtk_libs} )

        # link Simbody libraries
        TARGET_LINK_LIBRARIES(GMOLMODEL_dynamic
		optimized ${Simbody_LIBRARIES}
		SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
		${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code...
		${OpenMM_LIBRARIES}
		m )
		#${mmtk_libs} )

        TARGET_LINK_LIBRARIES(GMOLMODEL_static
		optimized ${Simbody_LIBRARIES}
		SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
		${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES}  # same as for dynamic lib... I don't think python shoud be part of this lib...
		${OpenMM_LIBRARIES}
		m )
                #${mmtk_libs} )
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

#set_target_properties(GMOLMODEL_dynamic_d PROPERTIES OUTPUT_NAME "gmolmodel")
#set_target_properties(GMOLMODEL_dynamic_d PROPERTIES VERSION 1.2 SOVERSION 1)

#set_target_properties(GMOLMODEL_static PROPERTIES OUTPUT_NAME "gmoldel")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)