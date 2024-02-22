SET(MOLMODEL_DIRS ${CMAKE_SOURCE_DIR}/Molmodel)
SET(MOLMODEL_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/Molmodel/include
    ${CMAKE_SOURCE_DIR}/Molmodel/include/molmodel
    ${CMAKE_SOURCE_DIR}/Molmodel/include/molmodel/internal)

SET(MOLMODEL_SOURCE_FILES)
SET(MOLMODEL_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${MOLMODEL_DIRS})
    FILE(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
    SET(MOLMODEL_C_SOURCE_FILES ${MOLMODEL_SOURCE_FILES} ${src_c_files})

    FILE(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    SET(MOLMODEL_CXX_SOURCE_FILES ${MOLMODEL_SOURCE_FILES} ${src_cxx_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(MOLMODEL_SOURCE_INCLUDE_FILES ${MOLMODEL_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

# set_source_files_properties(${MOLMODEL_C_SOURCE_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

set(MOLMODEL_COMPILE_DEFINITIONS
    MOLMODEL_COPYRIGHT_YEARS="2006-12"
    MOLMODEL_AUTHORS="Christopher.Bruns_Michael.Sherman"
    SimTK_MOLMODEL_LIBRARY_NAME="SimTKmolmodel"
    SimTK_MOLMODEL_MAJOR_VERSION=3
    SimTK_MOLMODEL_MINOR_VERSION=0
    SimTK_MOLMODEL_PATCH_VERSION=0
    OPENMM_PLATFORM_CPU=${USE_OPENMM_PLATFORM_CPU}
    OPENMM_PLATFORM_CUDA=${USE_OPENMM_PLATFORM_CUDA}
    OPENMM_PLATFORM_OPENCL=${USE_OPENMM_PLATFORM_OPENCL}
)

foreach(compile_definition ${MOLMODEL_COMPILE_DEFINITIONS})
    set_property(SOURCE ${MOLMODEL_C_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    set_property(SOURCE ${MOLMODEL_CXX_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    set_property(SOURCE ${MOLMODEL_SOURCE_INCLUDE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
endforeach(compile_definition)
