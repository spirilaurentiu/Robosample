set(BUILD_OPTIMIZERS OFF)
set(BUILD_IMPULSE_SOLVER OFF)
set(BUILD_GEOMETRY OFF)
set(BUILD_XML OFF)
set(BUILD_RANDOM OFF)
set(BUILD_CPODES OFF)
set(BUILD_CONSTRAINTS ON)

set(SIMBODY_DIRS
    ${SIMBODY_DIRS}
    .
    # SimTKcommon
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/BigMatrix/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/BigMatrix/include/SimTKcommon/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/BigMatrix/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/include/SimTKcommon
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/include/SimTKcommon/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Mechanics/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Mechanics/include/SimTKcommon
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Mechanics/include/SimTKcommon/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Mechanics/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Scalar/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Scalar/include/SimTKcommon
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Scalar/include/SimTKcommon/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Scalar/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Simulation/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Simulation/include/SimTKcommon/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Simulation/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/SmallMatrix/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/SmallMatrix/include/SimTKcommon
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/SmallMatrix/include/SimTKcommon/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/SmallMatrix/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/src
    # SimTKmath
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/include/simmath
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/include/simmath/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/include
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/include/simmath
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/include/simmath/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/LinearAlgebra/src
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/src
    # Simbody
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src)

if(BUILD_OPTIMIZERS)
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Optimizers/src)
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Optimizers/src/IpOpt)
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Optimizers/src/c-cmaes)
endif()

if(BUILD_GEOMETRY)
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Polynomial/include)
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Polynomial/include/SimTKcommon/internal
  )
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Polynomial/src)

  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Geometry/include)
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Geometry/include/SimTKcommon/internal
  )
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Geometry/src)

  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Geometry/include)
  list(
    APPEND SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Geometry/include/simmath/internal)
  list(APPEND SIMBODY_DIRS ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Geometry/src)
endif()

if(BUILD_RANDOM)
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Random/include)
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Random/include/SimTKcommon/internal
  )
  list(APPEND SIMBODY_DIRS ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/Random/src)
endif()

if(BUILD_CPODES)
  list(APPEND SIMBODY_DIRS
       ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes)
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/include
  )
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/include/cpodes
  )
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/include/nvector
  )
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/include/sundials
  )
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/src/cpodes
  )
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/src/nvector
  )
  list(
    APPEND
    SIMBODY_DIRS
    ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodes/sundials/src/sundials
  )
endif()

# Get source files
foreach(subdir ${SIMBODY_DIRS})
  file(GLOB src_c_files ${subdir}/*.c)
  set(SIMBODY_SOURCE_C_FILES ${SIMBODY_SOURCE_C_FILES} ${src_c_files})

  file(GLOB src_cxx_files ${subdir}/*.cpp)
  set(SIMBODY_SOURCE_CXX_FILES ${SIMBODY_SOURCE_CXX_FILES} ${src_cxx_files})

  file(GLOB incl_files ${subdir}/*.h)
  set(SIMBODY_SOURCE_INCLUDE_FILES ${SIMBODY_SOURCE_INCLUDE_FILES}
                                   ${incl_files})
endforeach(subdir)

# # print all includes message(STATUS "SIMBODY_SOURCE_INCLUDE_FILES:")
# foreach(file ${SIMBODY_SOURCE_INCLUDE_FILES}) message(STATUS "  ${file}")
# endforeach()

# message(STATUS "SIMBODY_SOURCE_CXX_FILES:") foreach(file
# ${SIMBODY_SOURCE_CXX_FILES}) message(STATUS "  ${file}") endforeach()

if(NOT BUILD_OPTIMIZERS)
  set(ASSEMBLER_SOURCES
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Assembler.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/AssemblyCondition.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/AssemblyCondition_QValue.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Assembler.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/AssemblyCondition_Markers.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/AssemblyCondition_Markers.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/AssemblyCondition_OrientationSensors.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/AssemblyCondition_OrientationSensors.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/LocalEnergyMinimizer.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/LocalEnergyMinimizer.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/ObservedPointFitter.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ObservedPointFitter.cpp)

  list(REMOVE_ITEM SIMBODY_SOURCE_CXX_FILES ${ASSEMBLER_SOURCES})
  list(REMOVE_ITEM SIMBODY_SOURCE_INCLUDE_FILES ${ASSEMBLER_SOURCES})
endif()

if(NOT BUILD_IMPULSE_SOLVER)
  set(IMPULSE_SOLVER_SOURCES
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/ImpulseSolver.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/PGSImpulseSolver.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/PLUSImpulseSolver.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/SemiExplicitEulerTimeStepper.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ImpulseSolver.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/PGSImpulseSolver.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/PLUSImpulseSolver.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/SemiExplicitEulerTimeStepper.cpp
  )
  list(REMOVE_ITEM SIMBODY_SOURCE_CXX_FILES ${IMPULSE_SOLVER_SOURCES})
  list(REMOVE_ITEM SIMBODY_SOURCE_INCLUDE_FILES ${IMPULSE_SOLVER_SOURCES})
endif()

if(NOT BUILD_GEOMETRY)
  set(GEOMETRY_SOURCES
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/DecorationSubsystem.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/DecorationSubsystem.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/DecorationSubsystemRep.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/ElasticFoundationForce.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ElasticFoundationForceImpl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ElasticFoundationForce.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/ContactTrackerSubsystem.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/CompliantContactSubsystem.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/ContactSurface.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/HuntCrossleyContact.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/GeneralContactSubsystem.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/CompliantContactSubsystem.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/HuntCrossleyContact.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/GeneralContactSubsystem.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ContactTrackerSubsystem.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/CablePath.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/CableSpring.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/CableTrackerSubsystem.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/CablePath_Impl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/CableSpring.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/CableTrackerSubsystem.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/CableTrackerSubsystem_Impl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/CablePath.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/HuntCrossleyContact.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/HuntCrossleyForce.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/HuntCrossleyForce.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/HuntCrossleyContact.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/HuntCrossleyForceImpl.h)
  list(REMOVE_ITEM SIMBODY_SOURCE_CXX_FILES ${GEOMETRY_SOURCES})
  list(REMOVE_ITEM SIMBODY_SOURCE_INCLUDE_FILES ${GEOMETRY_SOURCES})
endif()

if(NOT BUILD_XML)
  set(XML_SOURCES
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/include/SimTKcommon/internal/Xml.h
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/src/tinyxml.h
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/src/Xml.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/src/tinyxmlparser.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/src/tinyxml.cpp)
  list(REMOVE_ITEM SIMBODY_SOURCE_CXX_FILES ${XML_SOURCES})
  list(REMOVE_ITEM SIMBODY_SOURCE_INCLUDE_FILES ${XML_SOURCES})
endif()

if(NOT BUILD_CPODES)
  set(CPODES_SOURCES
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/include/simmath/internal/SimTKcpodes.h
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/include/simmath/CPodesIntegrator.h
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodesIntegrator.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/Integrators/src/CPodesIntegratorRep.h
  )
  list(REMOVE_ITEM SIMBODY_SOURCE_CXX_FILES ${CPODES_SOURCES})
  list(REMOVE_ITEM SIMBODY_SOURCE_INCLUDE_FILES ${CPODES_SOURCES})
endif()

if(NOT BUILD_CONSTRAINTS)
  set(CONSTRAINT_SOURCES
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_Weld.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_LineOnLineContact.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_SphereOnSphereContact.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/ConditionalConstraint.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_BuiltIns.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_PointOnPlaneContact.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_SphereOnPlaneContact.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_PointInPlane.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_Rod.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal/Constraint_Ball.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ConstraintImpl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/OBSOLETE_LengthConstraints.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/OBSOLETE_LengthConstraints.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_SphereOnPlaneContactImpl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_LineOnLineContactImpl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_PointOnPlaneContact.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_LineOnLineContact.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_SphereOnPlaneContact.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_PointOnPlaneContactImpl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_SphereOnSphereContactImpl.h
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_SphereOnSphereContact.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/ConditionalConstraint.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_Rod.cpp
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/src/Constraint_RodImpl.h)
  list(REMOVE_ITEM SIMBODY_SOURCE_CXX_FILES ${CONSTRAINT_SOURCES})
  list(REMOVE_ITEM SIMBODY_SOURCE_INCLUDE_FILES ${CONSTRAINT_SOURCES})
endif()

set(SIMBODY_COMMON_DEFS
    SimTK_SimTKCOMMON_LIBRARY_NAME="SimTKcommon"
    SimTK_SimTKCOMMON_MAJOR_VERSION=3
    SimTK_SimTKCOMMON_MINOR_VERSION=8
    SimTK_SimTKCOMMON_PATCH_VERSION=0
    SimTK_SIMMATH_LIBRARY_NAME="SimTKmath"
    SimTK_SIMMATH_MAJOR_VERSION=3
    SimTK_SIMMATH_MINOR_VERSION=8
    SimTK_SIMMATH_PATCH_VERSION=0
    SimTK_SIMBODY_LIBRARY_NAME="SimTKsimbody"
    SimTK_SIMBODY_MAJOR_VERSION=3
    SimTK_SIMBODY_MINOR_VERSION=8
    SimTK_SIMBODY_PATCH_VERSION=0
    BUILD_OPTIMIZERS=${BUILD_OPTIMIZERS}
    BUILD_IMPULSE_SOLVER=${BUILD_IMPULSE_SOLVER}
    BUILD_GEOMETRY=${BUILD_GEOMETRY}
    BUILD_XML=${BUILD_XML}
    BUILD_RANDOM=${BUILD_RANDOM}
    BUILD_CPODES=${BUILD_CPODES}
    BUILD_CONSTRAINTS=1 # ${BUILD_CONSTRAINTS}
)

set_source_files_properties(
  ${SIMBODY_SOURCE_C_FILES} PROPERTIES COMPILE_DEFINITIONS
                                       "${SIMBODY_COMMON_DEFS}")
set_source_files_properties(
  ${SIMBODY_SOURCE_CXX_FILES} PROPERTIES COMPILE_DEFINITIONS
                                         "${SIMBODY_COMMON_DEFS}")
set_source_files_properties(
  ${SIMBODY_SOURCE_INCLUDE_FILES} PROPERTIES COMPILE_DEFINITIONS
                                             "${SIMBODY_COMMON_DEFS}")
