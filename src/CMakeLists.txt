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
#set( Integrators_PATH
#     "${CMAKE_CURRENT_LIST_DIR}/../simbody/SimTKmath/Integrators/src"
#     "${CMAKE_CURRENT_LIST_DIR}/../Molmodel/Simbody01/SimTKmath/Integrators/src"
#     CACHE FILEPATH "Path to Integrator.h from Simbody..." )
# INCLUDE_DIRECTORIES( ${Integrators_PATH} )


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
        Robosample.cpp
	)

include_directories(../include/gmolmodel)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
         add_library(GMOLMODEL_dynamic_d SHARED ${LIBGMOLMODEL_SRC})
         add_library(GMOLMODEL_static_d  STATIC ${LIBGMOLMODEL_SRC})
ELSE()
         add_library(GMOLMODEL_dynamic SHARED ${LIBGMOLMODEL_SRC})
         add_library(GMOLMODEL_static  STATIC ${LIBGMOLMODEL_SRC})
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)


add_executable(GMOLMODEL_robo ${LIBGMOLMODEL_SRC} Robosample.cpp )



IF(CMAKE_BUILD_TYPE MATCHES Debug)
         message("Linking with debug libraries")

         TARGET_LINK_LIBRARIES( GMOLMODEL_robo
                 debug ${Simbody_LIBRARIES}
                 SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
                 ${BOOST_LIBRARIES}
	  	 ${OpenMM_LIBRARIES}
                 m )

         # link Simbody libraries
         TARGET_LINK_LIBRARIES(GMOLMODEL_dynamic_d
	 	debug ${Simbody_LIBRARIES}
	 	SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
	 	${BOOST_LIBRARIES}
	 	${OpenMM_LIBRARIES}
	 	m )
	
                
        TARGET_LINK_LIBRARIES(GMOLMODEL_static_d
	 	debug ${Simbody_LIBRARIES}
	 	SimTKmolmodel_d #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
	 	${BOOST_LIBRARIES}
	 	${OpenMM_LIBRARIES}
	 	m )
        
ELSE()
        message("Linking with release libraries")

        # TARGET_LINK_LIBRARIES( GMOLMODEL_test
        #         optimized ${Simbody_LIBRARIES}
        #         SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
        #         ${BOOST_LIBRARIES}
        #         #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code..
	# 	${OpenMM_LIBRARIES}
        #         m )
        #         #${mmtk_libs} )

        TARGET_LINK_LIBRARIES( GMOLMODEL_robo
                optimized ${Simbody_LIBRARIES}
                SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
                ${BOOST_LIBRARIES}
                #${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code..
		${OpenMM_LIBRARIES}
                m )
                
        # link Simbody libraries
         TARGET_LINK_LIBRARIES(GMOLMODEL_dynamic
	 	optimized ${Simbody_LIBRARIES}
	 	SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
	 	${BOOST_LIBRARIES}
	 	${OpenMM_LIBRARIES}
	 	m )


         TARGET_LINK_LIBRARIES(GMOLMODEL_static
	 	optimized ${Simbody_LIBRARIES}
	 	SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
	 	${BOOST_LIBRARIES}
	 	${OpenMM_LIBRARIES}
	 	m )

ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

#set_target_properties(GMOLMODEL_dynamic_d PROPERTIES OUTPUT_NAME "gmolmodel")
#set_target_properties(GMOLMODEL_dynamic_d PROPERTIES VERSION 1.2 SOVERSION 1)

#set_target_properties(GMOLMODEL_static PROPERTIES OUTPUT_NAME "gmoldel")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
