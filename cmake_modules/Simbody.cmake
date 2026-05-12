set(SIMBODY_DIRS)
set(SIMBODY_INCLUDE_DIRS)

set(SIMTK_COMMON_DIRS
    .
    Scalar
    SmallMatrix
    Mechanics
    BigMatrix
    Geometry
    Simulation
    Random
    Polynomial)
foreach(subdir ${SIMTK_COMMON_DIRS})
  set(SIMBODY_DIRS ${SIMBODY_DIRS}
                   ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir})
  set(SIMBODY_INCLUDE_DIRS
      ${SIMBODY_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include/SimTKcommon
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include/SimTKcommon/internal
  )
endforeach(subdir)

set(SIMTK_MATH_DIRS . LinearAlgebra Integrators Integrators/src/CPodes/sundials
                    Optimizers Geometry)
foreach(subdir ${SIMTK_MATH_DIRS})
  set(SIMBODY_DIRS ${SIMBODY_DIRS}
                   ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir})
  set(SIMBODY_INCLUDE_DIRS
      ${SIMBODY_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include/simmath
      ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include/simmath/internal
  )
endforeach(subdir)

set(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/Simbody)
if(BUILD_VISUALIZER)
  set(SIMBODY_DIRS ${SIMBODY_DIRS}
                   ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer)
endif(BUILD_VISUALIZER)

set(SIMBODY_INCLUDE_DIRS
    ${SIMBODY_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include/simbody
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include/simbody/internal)

if(BUILD_VISUALIZER)
  set(SIMBODY_INCLUDE_DIRS
      ${SIMBODY_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include/simbody
      ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include/simbody/internal)
endif(BUILD_VISUALIZER)

# find source and header files
set(SIMBODY_SOURCE_C_FILES)
set(SIMBODY_SOURCE_CXX_FILES)
set(SIMBODY_SOURCE_INCLUDE_FILES)

foreach(subdir ${SIMBODY_DIRS})
  file(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
  set(SIMBODY_SOURCE_C_FILES ${SIMBODY_SOURCE_C_FILES} ${src_c_files})

  file(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
  set(SIMBODY_SOURCE_CXX_FILES ${SIMBODY_SOURCE_CXX_FILES} ${src_cxx_files})

  # pimpl pattern is used and headers are stored in the src directory
  file(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
  set(SIMBODY_SOURCE_INCLUDE_FILES ${SIMBODY_SOURCE_INCLUDE_FILES}
                                   ${incl_files})
endforeach(subdir)

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
    BUILD_OPTIMIZERS=1 # ${BUILD_OPTIMIZERS}
    BUILD_IMPULSE_SOLVER=1 # ${BUILD_IMPULSE_SOLVER}
    BUILD_GEOMETRY=1 # ${BUILD_GEOMETRY}
    BUILD_XML=1 # ${BUILD_XML}
    BUILD_RANDOM=1 # ${BUILD_RANDOM}
    BUILD_CPODES=1 # ${BUILD_CPODES}
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
