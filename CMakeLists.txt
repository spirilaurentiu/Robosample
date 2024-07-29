cmake_minimum_required(VERSION 3.17)

project(Robosample C CXX)

# robosample version and date
SET(ROBOSAMPLE_MAJOR_VERSION 2)
SET(ROBOSAMPLE_MINOR_VERSION 0)
SET(ROBOSAMPLE_PATCH_VERSION 0)
SET(ROBOSAMPLE_COPYRIGHT_YEARS "2021-01")

# underbar separated list of dotted authors, no spaces or commas
SET(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")

# OpenMM wants these set for some library
SET(COMMON_FLAGS "-msse4.1 -mavx")

# -O2 is marginally slower than -O3
SET(EXTRA_COMPILE_FLAGS_DEBUG   "${COMMON_FLAGS} -O1 -g -ggdb")
SET(EXTRA_COMPILE_FLAGS_RELEASE "${COMMON_FLAGS} -O3 -march=native -mtune=native")

# set the profile guided optimization (pgo) flags
SET(PGO_DIR ${CMAKE_BINARY_DIR}/profile-data)
SET(PGO_FILE ${CMAKE_BINARY_DIR}/profile-data/default.profraw)

IF(CMAKE_COMPILER_IS_GNUCXX)
    IF(CMAKE_BUILD_TYPE MATCHES "PGO_Train")
        SET(PGO_FLAGS "-fprofile-generate=${PGO_DIR} -fno-reorder-blocks-and-partition")
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Use")
        SET(PGO_FLAGS "-fprofile-use=${PGO_DIR} -fprofile-correction -fno-reorder-blocks-and-partition")
    ENDIF()
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    SET(PGO_CLANG_RAW ${CMAKE_BINARY_DIR}/profile-data/default.profraw)
    SET(PGO_CLANG_DATA ${CMAKE_BINARY_DIR}/profile-data/default.profdata)

    IF(CMAKE_BUILD_TYPE MATCHES "PGO_Train")
        # -fcoverage-mapping -ftest-coverage -fprofile-arcs
        SET(PGO_FLAGS "-fprofile-instr-generate=${PGO_CLANG_RAW}") # -fprofile-update=atomic
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Use")
        SET(PGO_FLAGS "-fprofile-instr-use=${PGO_CLANG_DATA}")
    ENDIF()
ENDIF()

set(BUILD_VISUALIZER ON)

# building tests requires us to rebuild the entire project as a shared object using different flags
# default builds require only static libraries (different flags)
# there are multiple build configurations
# - debug (no tests): build only as shared object for faster compilation (faster linking)
# - debug_tests: build only as shared object for faster compilation and compile tests
# - release (no tests): build only a standalone executable with ipo on from static objects
# - release_tests: build static and dynamic objects (for the executable and tests respectively)
# - pgo_train (no tests): same as release
# - pgo_use (no tests): same as release
IF(CMAKE_BUILD_TYPE MATCHES "Debug" OR CMAKE_BUILD_TYPE MATCHES "Debug_Tests" OR CMAKE_BUILD_TYPE MATCHES "Debug_Full")
    # set the sanitizer flags
    # https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html
    # https://clang.llvm.org/docs/AddressSanitizer.html
    #SET(SANITIZER_FLAGS "-fsanitize=undefined,leak,address,float-divide-by-zero -fsanitize-address-use-after-scope -fsanitize-recover=address")
    SET(SANITIZER_FLAGS "-fsanitize=leak,address,float-divide-by-zero -fsanitize-address-use-after-scope -fsanitize-recover=address")

    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_DEBUG} -fno-omit-frame-pointer")
    SET(SUFFIX "_d")

ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES "Release_Tests" OR CMAKE_BUILD_TYPE MATCHES "Release_Full")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(BOLT_FLAGS "-Wl,--emit-relocs")

ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Train")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(SUFFIX "_pgo_train")

ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Use")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(SUFFIX "_pgo_use")
    SET(BOLT_FLAGS "-Wl,--emit-relocs")
ENDIF()

separate_arguments(SANITIZER_FLAGS UNIX_COMMAND "${SANITIZER_FLAGS}")
separate_arguments(PGO_FLAGS UNIX_COMMAND "${PGO_FLAGS}")
separate_arguments(BUILD_FLAGS UNIX_COMMAND "${BUILD_FLAGS}")
separate_arguments(BOLT_FLAGS UNIX_COMMAND "${BOLT_FLAGS}")

# set names of the targets we will build
SET(ROBOSAMPLE_LIB_NAME "robosample${SUFFIX}")
SET(ROBOSAMPLE_EXE_NAME "robosample_exe${SUFFIX}")
SET(ROBOSAMPLE_NAME "robosample${SUFFIX}")

# we require ipo (interprocedural optimization) aka lto (link time optimization)
# this is fancy -flto basically
# on clang, it sets -flto=thin and this seems to work better than the default -flto(=full)
# even though the binary is larger, it runs faster
include(CheckIPOSupported)
check_ipo_supported()

IF(NOT CMAKE_BUILD_TYPE MATCHES "Debug" AND NOT CMAKE_BUILD_TYPE MATCHES "Debug_Tests")
    SET(USE_IPO TRUE)
ELSE()
    SET(USE_IPO FALSE)
ENDIF()

# SET(CMAKE_CXX_STANDARD 14)

# link prerequisites
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(DL_LIBRARY dl REQUIRED)
find_library(MATH_LIBRARY m REQUIRED)
find_library(RT_LIBRARY rt REQUIRED) # wants -lrt

SET(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

if(BUILD_VISUALIZER)
    set(OpenGL_GL_PREFERENCE LEGACY)
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)

    SET(VISUALIZER_LIBRARIES
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES}
    )
endif()

SET(LIBRARIES
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${RT_LIBRARY}
    Threads::Threads
    ${DL_LIBRARY}
    ${MATH_LIBRARY}
)

# # set the precompiled headers
# set(OPENMM_PCH pch/openmm_pch.h)
# set(SIMBODY_PCH pch/simbody_pch.h)
# set(MOLMODEL_PCH pch/molmodel_pch.h)
# set(ROBOSAMPLE_PCH pch/robosample_pch.h)

if (CMAKE_BUILD_TYPE MATCHES "Python_Debug" OR CMAKE_BUILD_TYPE MATCHES "Python_Release")
    # load pybind11
    # find_package(PythonLibs REQUIRED)
    # find_package(Python COMPONENTS Interpreter Development REQUIRED)
    add_subdirectory(pybind11)
endif()

include(Definitions.OpenMM.cmake)
include(Definitions.Simbody.cmake)
include(Definitions.Molmodel.cmake)

# build the robosample library
FILE(GLOB src_c_files ${CMAKE_SOURCE_DIR}/format/src/*.c)
SET(ROBOSAMPLE_C_SOURCE_FILES ${src_c_files})
set_source_files_properties(${ROBOSAMPLE_C_SOURCE_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

FILE(GLOB src_cxx_files ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/format/src/*.cpp)
SET(ROBOSAMPLE_CXX_SOURCE_FILES ${src_cxx_files})
LIST(REMOVE_ITEM ROBOSAMPLE_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)
LIST(REMOVE_ITEM ROBOSAMPLE_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp)

set(ROBOSAMPLE_SOURCE_FILES
    ${OPENMM_SOURCE_C_FILES} ${OPENMM_SOURCE_CXX_FILES} ${OPENMM_SOURCE_INCLUDE_FILES}
    ${SIMBODY_SOURCE_C_FILES} ${SIMBODY_SOURCE_CXX_FILES} ${SIMBODY_SOURCE_INCLUDE_FILES}
    ${MOLMODEL_C_SOURCE_FILES} ${MOLMODEL_CXX_SOURCE_FILES} ${MOLMODEL_SOURCE_INCLUDE_FILES}
    ${ROBOSAMPLE_C_SOURCE_FILES} ${ROBOSAMPLE_CXX_SOURCE_FILES}
)

set(ROBOSAMPLE_INCLUDE_DIRS
    ${OPENMM_INCLUDE_DIRS}
    ${SIMBODY_INCLUDE_DIRS}
    ${MOLMODEL_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/format/inc
    ${CMAKE_SOURCE_DIR}/format/inc/mdtraj
    ${CMAKE_SOURCE_DIR}/pcg-cpp/include
)

set(ROBOSAMPLE_COMPILE_DEFINITIONS
    OPENMM_PLATFORM_CPU=${USE_OPENMM_PLATFORM_CPU}
    OPENMM_PLATFORM_CUDA=${USE_OPENMM_PLATFORM_CUDA}
    OPENMM_PLATFORM_OPENCL=${USE_OPENMM_PLATFORM_OPENCL}
)

# target_precompile_headers(${SIMBODY_LIB_SHARED_NAME} PRIVATE ${SIMBODY_PCH})

# python bindings
if (CMAKE_BUILD_TYPE MATCHES "Python_Debug" OR CMAKE_BUILD_TYPE MATCHES "Python_Release")
    # we do not add sanitizer flags here because asan cannot be linked
    pybind11_add_module(${ROBOSAMPLE_LIB_NAME} ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp ${ROBOSAMPLE_SOURCE_FILES})

    target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE MODULE_NAME=${ROBOSAMPLE_LIB_NAME})
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    
    target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${BOLT_FLAGS} ${PGO_FLAGS})
    target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE ${LIBRARIES} ${OPENMM_LIBRARIES})

    # common for both
    target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
    set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

    foreach(compile_definition ${ROBOSAMPLE_COMPILE_DEFINITIONS})
        set_property(SOURCE ${ROBOSAMPLE_C_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
        set_property(SOURCE ${ROBOSAMPLE_CXX_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    endforeach(compile_definition)

    IF(NOT OPENMM_PLATFORM MATCHES "CPU")
        add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
    ENDIF(NOT OPENMM_PLATFORM MATCHES "CPU")

else()
    # robosample library
    add_library(${ROBOSAMPLE_LIB_NAME} SHARED ${ROBOSAMPLE_SOURCE_FILES})
    
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${BUILD_FLAGS} ${PGO_FLAGS})
    target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
    target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${LIBRARIES} ${OPENMM_LIBRARIES})

    # common for both
    target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
    set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

    foreach(compile_definition ${ROBOSAMPLE_COMPILE_DEFINITIONS})
        set_property(SOURCE ${ROBOSAMPLE_C_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
        set_property(SOURCE ${ROBOSAMPLE_CXX_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    endforeach(compile_definition)

    IF(NOT OPENMM_PLATFORM MATCHES "CPU")
        add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
    ENDIF(NOT OPENMM_PLATFORM MATCHES "CPU")

    # robosample executable
    add_executable(${ROBOSAMPLE_EXE_NAME} ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)
    target_include_directories(${ROBOSAMPLE_EXE_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
    target_compile_options(${ROBOSAMPLE_EXE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BUILD_FLAGS} ${PGO_FLAGS})
    set_target_properties(${ROBOSAMPLE_EXE_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${ROBOSAMPLE_EXE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
    target_link_libraries(${ROBOSAMPLE_EXE_NAME} PRIVATE ${ROBOSAMPLE_LIB_NAME})

    set_property(TARGET ${ROBOSAMPLE_EXE_NAME} PROPERTY OUTPUT_NAME ${ROBOSAMPLE_NAME})

    # build the simbody visualizer
    IF(BUILD_VISUALIZER)
        IF(NOT CMAKE_BUILD_TYPE MATCHES "Debug" AND NOT CMAKE_BUILD_TYPE MATCHES "Debug_Tests")
            SET(GUI_NAME "simbody-visualizer")
        ELSE()
            SET(GUI_NAME "simbody-visualizer_d")
        ENDIF()

        SET(VISUALIZER_DIR ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/simbody-visualizer)
        FILE(GLOB VISUALIZER_SOURCE_INCLUDE_FILES ${VISUALIZER_DIR}/*.h)
        FILE(GLOB VISUALIZER_SOURCE_CXX_FILES ${VISUALIZER_DIR}/*.cpp)

        add_executable(${GUI_NAME} ${VISUALIZER_SOURCE_INCLUDE_FILES} ${VISUALIZER_SOURCE_CXX_FILES})
        target_include_directories(${GUI_NAME} PRIVATE ${GLUT_INCLUDE_DIR} ${SIMBODY_INCLUDE_DIRS})

        # TODO do we need to link ${LIBRARIES}?
        target_link_libraries(${GUI_NAME} PRIVATE ${ROBOSAMPLE_LIB_NAME} ${LIBRARIES} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

        set_target_properties(${GUI_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
        target_link_options(${GUI_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    ENDIF(BUILD_VISUALIZER)

    # build tests
    IF (CMAKE_BUILD_TYPE MATCHES "Debug_Tests" OR CMAKE_BUILD_TYPE MATCHES "Release_Tests" OR
        CMAKE_BUILD_TYPE MATCHES "Debug_Full" OR CMAKE_BUILD_TYPE MATCHES "Release_Full")

        enable_testing()
        include(CTest)

        SET(TEST_NAMES)
        SET(TEST_SOURCES)

        # # find openmm tests
        # FILE(GLOB TESTS_SOURCES openmm/tests/*.cpp)
        # FOREACH(test_src ${TESTS_SOURCES})
        #     # get test name by removing file location and extension
        #     FILE(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/openmm/tests ${test_src})
        #     STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

        #     LIST(APPEND TEST_NAMES ${test_name})
        #     LIST(APPEND TEST_SOURCES ${test_src})
        # ENDFOREACH(test_src)

        # # find simbody tests
        # FILE(GLOB TESTS_SOURCES Simbody01/Simbody/tests/*.cpp)
        # FOREACH(test_src ${TESTS_SOURCES})
        #     # get test name by removing file location and extension
        #     FILE(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/tests ${test_src})
        #     STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

        #     LIST(APPEND TEST_NAMES ${test_name})
        #     LIST(APPEND TEST_SOURCES ${test_src})
        # ENDFOREACH(test_src)

        # find robosample tests
        FILE(GLOB TESTS_SOURCES tests/*.cpp)
        FOREACH(test_src ${TESTS_SOURCES})
            # get test name by removing file location and extension
            FILE(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/tests ${test_src})
            STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

            LIST(APPEND TEST_NAMES ${test_name})
            LIST(APPEND TEST_SOURCES ${test_src})
        ENDFOREACH(test_src)

        # add the tests
        FOREACH(test_name test_src IN ZIP_LISTS TEST_NAMES TEST_SOURCES)
            # the test executable is non-ipo
            add_executable(${test_name} ${test_src})
            target_compile_options(${test_name} PRIVATE ${BUILD_FLAGS})
            target_link_options(${test_name} PRIVATE ${SANITIZER_FLAGS})
            target_include_directories(${test_name} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
            target_link_libraries(${test_name} PRIVATE ${ROBOSAMPLE_LIB_NAME})

            add_test(NAME ${test_name} COMMAND ${test_name})
        ENDFOREACH()
    ENDIF ()
endif()

# add example inputs
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp/pdbs)
file(GLOB INPUTS ${CMAKE_SOURCE_DIR}/tests_inputs/*)
file(COPY ${INPUTS} DESTINATION ${CMAKE_BINARY_DIR})
