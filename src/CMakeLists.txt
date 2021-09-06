cmake_minimum_required(VERSION 3.0)

# verbose ouput
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "Control volume of build output" )

# Permit use of custom FindOpenMM module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

find_package(OpenMM)
if(OpenMM_FOUND)
        message("OpenMM found")

        INCLUDE_DIRECTORIES(${OpenMM_INCLUDE_DIR})
        INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
else()
        message("OpenMM not found")
endif()

# find simbody
# @TODO do we need 3.6 CONFIG?
find_package(Simbody REQUIRED)
INCLUDE_DIRECTORIES( ${Simbody_INCLUDE_DIR} )
LINK_DIRECTORIES( ${Simbody_LIB_DIR} )

# find molmodel
set( MOL_MODEL_INCLUDE
     "${Simbody_INCLUDE_DIR}/.."
     CACHE FILEPATH "header include folder path for Molmodel" )
INCLUDE_DIRECTORIES( ${MOL_MODEL_INCLUDE} )


# #it's for linking GMOLMODEL_static into GCHMC lib
# LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/lib)

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
        Robosample.cpp
        )

include_directories(../include/gmolmodel)

# IF(CMAKE_BUILD_TYPE MATCHES Debug)
#          add_library(GMOLMODEL_dynamic_d SHARED ${LIBGMOLMODEL_SRC})
#          add_library(GMOLMODEL_static_d  STATIC ${LIBGMOLMODEL_SRC})
# ELSE()
#          add_library(GMOLMODEL_dynamic SHARED ${LIBGMOLMODEL_SRC})
#          add_library(GMOLMODEL_static  STATIC ${LIBGMOLMODEL_SRC})
# ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)


add_executable(GMOLMODEL_robo ${LIBGMOLMODEL_SRC} Robosample.cpp )



IF(CMAKE_BUILD_TYPE MATCHES Debug)
        message("Linking with debug libraries")

        TARGET_LINK_LIBRARIES( GMOLMODEL_robo
                debug ${Simbody_LIBRARIES}
                SimTKmolmodel_d
                ${OpenMM_LIBRARIES}
                m )

        # TARGET_LINK_LIBRARIES( GMOLMODEL_dynamic_d
	#  	debug ${Simbody_LIBRARIES}
	#  	SimTKmolmodel_d
	#  	${OpenMM_LIBRARIES}
	#  	m )
	
        # TARGET_LINK_LIBRARIES( GMOLMODEL_static_d
	#  	debug ${Simbody_LIBRARIES}
	#  	SimTKmolmodel_d
	#  	${OpenMM_LIBRARIES}
	#  	m )
        
ELSE()
        message("Linking with release libraries")

        TARGET_LINK_LIBRARIES( GMOLMODEL_robo
                optimized ${Simbody_LIBRARIES}
                SimTKmolmodel
                ${OpenMM_LIBRARIES}
                m )
                
        # TARGET_LINK_LIBRARIES( GMOLMODEL_dynamic
        #         optimized ${Simbody_LIBRARIES}
	#  	SimTKmolmodel
	#  	${OpenMM_LIBRARIES}
	#  	m )


        # TARGET_LINK_LIBRARIES( GMOLMODEL_static
	#  	optimized ${Simbody_LIBRARIES}
	#  	SimTKmolmodel
	#  	${OpenMM_LIBRARIES}
	#  	m )

ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

#set_target_properties(GMOLMODEL_dynamic_d PROPERTIES OUTPUT_NAME "gmolmodel")
#set_target_properties(GMOLMODEL_dynamic_d PROPERTIES VERSION 1.2 SOVERSION 1)

#set_target_properties(GMOLMODEL_static PROPERTIES OUTPUT_NAME "gmoldel")

# set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
