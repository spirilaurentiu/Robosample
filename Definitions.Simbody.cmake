SET(SIMBODY_DIRS)
SET(SIMBODY_INCLUDE_DIRS)

SET(SIMTK_COMMON_DIRS . Scalar SmallMatrix Mechanics BigMatrix Geometry Simulation Random Polynomial)
foreach(subdir ${SIMTK_COMMON_DIRS})
    SET(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir})
    SET(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include/SimTKcommon
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include/SimTKcommon/internal)
endforeach(subdir)

SET(SIMTK_MATH_DIRS . LinearAlgebra Integrators Integrators/src/CPodes/sundials Optimizers Geometry)
foreach(subdir ${SIMTK_MATH_DIRS})
    SET(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir})
    SET(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include/simmath
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include/simmath/internal)
endforeach(subdir)

SET(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/Simbody)
IF(BUILD_VISUALIZER)
    SET(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer)
ENDIF(BUILD_VISUALIZER)

SET(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal)
IF(BUILD_VISUALIZER)
    SET(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include
        ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include/simbody
        ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/include/simbody/internal)
ENDIF(BUILD_VISUALIZER)

# find source and header files
SET(SIMBODY_SOURCE_C_FILES)
SET(SIMBODY_SOURCE_CXX_FILES)
SET(SIMBODY_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${SIMBODY_DIRS})
    FILE(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
    SET(SIMBODY_SOURCE_C_FILES ${SIMBODY_SOURCE_C_FILES} ${src_c_files})

    FILE(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    SET(SIMBODY_SOURCE_CXX_FILES ${SIMBODY_SOURCE_CXX_FILES} ${src_cxx_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(SIMBODY_SOURCE_INCLUDE_FILES ${SIMBODY_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

# set_source_files_properties(${SIMBODY_SOURCE_C_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

set(SIMBODY_COMPILE_DEFINITIONS
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
    SIMBODY_VISUALIZER_REL_INSTALL_DIR="${CMAKE_BINARY_DIR}"
    SIMBODY_PATH_FROM_LIBDIR_TO_VIZ_DIR="${CMAKE_BINARY_DIR}"
    SIMBODY_VISUALIZER_INSTALL_DIR="${CMAKE_BINARY_DIR}"
    SIMBODY_VISUALIZER_REL_INSTALL_DIR="${CMAKE_BINARY_DIR}"
)

foreach(compile_definition ${SIMBODY_COMPILE_DEFINITIONS})
    set_property(SOURCE ${SIMBODY_SOURCE_C_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    set_property(SOURCE ${SIMBODY_SOURCE_CXX_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    set_property(SOURCE ${SIMBODY_SOURCE_INCLUDE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
endforeach(compile_definition)
