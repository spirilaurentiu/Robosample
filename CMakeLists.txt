cmake_minimum_required(VERSION 3.17)

# CUDA may be used as a platform for OpenMM
# Setting it here simplifies the logic even if CUDA is not used as the OpenMM platform
project(Robosample LANGUAGES C CXX)

# Define mutually exclusive hardware acceleration options
option(USE_CPU    "Enable CPU backend"    OFF)
option(USE_CUDA   "Enable CUDA backend"   OFF)
option(USE_OPENCL "Enable OpenCL backend" OFF)

# Tell CMake to enable CUDA
if(USE_CUDA)
    enable_language(CUDA)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "CUDA enabled but CMAKE_CUDA_ARCHITECTURES not set")
    endif()
endif()

# Robosample version and date
set(ROBOSAMPLE_MAJOR_VERSION 2)
set(ROBOSAMPLE_MINOR_VERSION 0)
set(ROBOSAMPLE_PATCH_VERSION 0)
set(ROBOSAMPLE_COPYRIGHT_YEARS "2021-01")
set(ROBOSAMPLE_LIB_NAME "robo_bindings")

# Underbar separated list of dotted authors, no spaces or commas
set(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")

# Generate compile_commands.json for tools like clangd and code analysis
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Ensure the Python extension uses relative RPATH ($ORIGIN)
# This prevents Conda's lib directory from shadowing system libraries
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN")

# Link prerequisites
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Threads REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# if(BUILD_VISUALIZER)
#     set(OpenGL_GL_PREFERENCE LEGACY)
#     find_package(OpenGL REQUIRED)
#     find_package(GLUT REQUIRED)

#     SET(VISUALIZER_LIBRARIES
#         ${OPENGL_LIBRARIES}
#         ${GLUT_LIBRARIES}
#     )
# endif()

include(${CMAKE_SOURCE_DIR}/cmake_modules/OpenMM.cmake)
include(${CMAKE_SOURCE_DIR}/cmake_modules/Simbody.cmake)
include(${CMAKE_SOURCE_DIR}/cmake_modules/Molmodel.cmake)

# Build the robosample library
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

# Compiler flags for Debug and Release builds
set(COMMON_FLAGS "-msse4.1 -mavx -Wno-c++17-extensions")
separate_arguments(COMMON_FLAGS UNIX_COMMAND "${COMMON_FLAGS}")

set(DEBUG_FLAGS "-O0 -g -ggdb -fno-omit-frame-pointer") #  -fsanitize=address
separate_arguments(DEBUG_FLAGS UNIX_COMMAND "${DEBUG_FLAGS}")

set(RELEASE_FLAGS "-O3 -march=native -mtune=native")
separate_arguments(RELEASE_FLAGS UNIX_COMMAND "${RELEASE_FLAGS}")

set(RELWITHDEBINFO_FLAGS "-O3 -g3 -ggdb3 -fno-lto -fvisibility=default -fno-omit-frame-pointer -march=native -mtune=native")
separate_arguments(RELWITHDEBINFO_FLAGS UNIX_COMMAND "${RELWITHDEBINFO_FLAGS}")

set(CMAKE_CXX_STANDARD 14)

# Define BOLT flags
set(BOLT_FLAGS "-Wl,--emit-relocs")
separate_arguments(BOLT_FLAGS UNIX_COMMAND "${BOLT_FLAGS}")

# Set sanitizer flags (only for Debug build type) FOR TESTS ONLY
set(SANITIZER_FLAGS "-fsanitize=address,leak,float-divide-by-zero")
separate_arguments(SANITIZER_FLAGS UNIX_COMMAND "${SANITIZER_FLAGS}")

# Profile-guided optimization (PGO) flags for different stages and compilers
set(PGO_DIR ${CMAKE_BINARY_DIR}/profile-data)
if(CMAKE_COMPILER_IS_GNUCXX)
    set(PGO_TRAIN_FLAGS "-fprofile-generate=${PGO_DIR} -fno-reorder-blocks-and-partition")
    set(PGO_USE_FLAGS "-fprofile-use=${PGO_DIR} -fprofile-correction -fno-reorder-blocks-and-partition")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(PGO_TRAIN_FLAGS "-fprofile-instr-generate=${PGO_DIR}/raw.profraw")
    set(PGO_USE_FLAGS "-fprofile-instr-use=${PGO_DIR}/data.profdata")
endif()

separate_arguments(PGO_TRAIN_FLAGS UNIX_COMMAND "${PGO_TRAIN_FLAGS}")
separate_arguments(PGO_USE_FLAGS UNIX_COMMAND "${PGO_USE_FLAGS} -g")

# Check if we can use IPO (Interprocedural Optimization) aka LTO (Link Time Optimization)
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED)
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "PGO_Train" OR CMAKE_BUILD_TYPE STREQUAL "PGO_Use")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

# Create the Python module
add_subdirectory(pybind11)
pybind11_add_module(${ROBOSAMPLE_LIB_NAME}
    ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp
    ${ROBOSAMPLE_SOURCE_FILES}
)

# Add preprocessor definitions
target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE MODULE_NAME=${ROBOSAMPLE_LIB_NAME})
foreach(backend CPU CUDA OPENCL)
    if(USE_${backend})
        target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE USE_${backend}=1)
    else()
        target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE USE_${backend}=0)
    endif()
endforeach()


# Apply common includes
target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})

# Set PIE
set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# # Apply base build flags
# target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${BUILD_FLAGS})

# Apply PGO and BOLT flags conditionally
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${DEBUG_FLAGS})
    add_compile_definitions(_GLIBCXX_DEBUG)
    # target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE -fsanitize=address)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${COMMON_FLAGS} ${RELEASE_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${RELWITHDEBINFO_FLAGS})
    target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${RELWITHDEBINFO_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "PGO_Train")
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${COMMON_FLAGS} ${RELEASE_FLAGS} ${PGO_TRAIN_FLAGS})
    target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${PGO_TRAIN_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "PGO_Use" OR CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${COMMON_FLAGS} ${RELEASE_FLAGS} ${PGO_USE_FLAGS})
    target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${PGO_USE_FLAGS} ${BOLT_FLAGS})
endif()

# # Apply sanitizers (Debug only; ASAN cannot link with Python runtime)
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS})
#     target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS})
#     target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ASAN_LIB})
# endif()

# Link with external libraries
target_link_libraries(${ROBOSAMPLE_LIB_NAME}
    PRIVATE
        BLAS::BLAS
        LAPACK::LAPACK
        Threads::Threads
)

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE CUDA::cuda_driver CUDA::cufft)
elseif(USE_OPENCL)
    find_package(OpenCL REQUIRED)
    target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE OpenCL::OpenCL)
endif()

# Optional: dependency on OpenMM build target
if(NOT USE_CPU)
    add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
endif()



## Install extension module into site-packages/robosample
#install(TARGETS ${ROBOSAMPLE_LIB_NAME} LIBRARY DESTINATION ${Python3_SITEARCH}/robosample) # ${CMAKE_INSTALL_LIBDIR}
#
## Install pure Python files
#install(DIRECTORY python/robosample DESTINATION ${Python3_SITEARCH})

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
# # # # # # # # # # # # # # # # # # # install(TARGETS ${ROBOSAMPLE_LIB_NAME} DESTINATION lib)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

#install(DIRECTORY ${CMAKE_SOURCE_DIR}/python/robosample
#        DESTINATION ${Python3_SITEARCH})









# # python bindings
# if (CMAKE_BUILD_TYPE MATCHES "Python_Debug" OR CMAKE_BUILD_TYPE MATCHES "Python_Release")
#     # we do not add sanitizer flags here because asan cannot be linked
#     pybind11_add_module(${ROBOSAMPLE_LIB_NAME} ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp ${ROBOSAMPLE_SOURCE_FILES})

#     target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE MODULE_NAME=${ROBOSAMPLE_LIB_NAME})
#     target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    
#     target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${BOLT_FLAGS} ${PGO_FLAGS})
#     target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE ${LIBRARIES} ${OPENMM_LIBRARIES})

#     # common for both
#     target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
#     set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

#     foreach(compile_definition ${ROBOSAMPLE_COMPILE_DEFINITIONS})
#         set_property(SOURCE ${ROBOSAMPLE_C_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
#         set_property(SOURCE ${ROBOSAMPLE_CXX_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
#     endforeach(compile_definition)

#     IF(NOT OPENMM_PLATFORM MATCHES "CPU")
#         add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
#     ENDIF(NOT OPENMM_PLATFORM MATCHES "CPU")

# else()
#     # robosample library
#     add_library(${ROBOSAMPLE_LIB_NAME} SHARED ${ROBOSAMPLE_SOURCE_FILES})
    
#     target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${BUILD_FLAGS} ${PGO_FLAGS})
#     target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
#     target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${LIBRARIES} ${OPENMM_LIBRARIES})

#     # common for both
#     target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
#     set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

#     foreach(compile_definition ${ROBOSAMPLE_COMPILE_DEFINITIONS})
#         set_property(SOURCE ${ROBOSAMPLE_C_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
#         set_property(SOURCE ${ROBOSAMPLE_CXX_SOURCE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
#     endforeach(compile_definition)

#     IF(NOT OPENMM_PLATFORM MATCHES "CPU")
#         add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
#     ENDIF(NOT OPENMM_PLATFORM MATCHES "CPU")

#     # robosample executable
#     add_executable(${ROBOSAMPLE_EXE_NAME} ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)
#     target_include_directories(${ROBOSAMPLE_EXE_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
#     target_compile_options(${ROBOSAMPLE_EXE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BUILD_FLAGS} ${PGO_FLAGS})
#     set_target_properties(${ROBOSAMPLE_EXE_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
#     target_link_options(${ROBOSAMPLE_EXE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
#     target_link_libraries(${ROBOSAMPLE_EXE_NAME} PRIVATE ${ROBOSAMPLE_LIB_NAME})

#     set_property(TARGET ${ROBOSAMPLE_EXE_NAME} PROPERTY OUTPUT_NAME ${ROBOSAMPLE_NAME})

#     # build the simbody visualizer
#     IF(BUILD_VISUALIZER)
#         IF(NOT CMAKE_BUILD_TYPE MATCHES "Debug" AND NOT CMAKE_BUILD_TYPE MATCHES "Debug_Tests")
#             SET(GUI_NAME "simbody-visualizer")
#         ELSE()
#             SET(GUI_NAME "simbody-visualizer_d")
#         ENDIF()

#         SET(VISUALIZER_DIR ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/simbody-visualizer)
#         FILE(GLOB VISUALIZER_SOURCE_INCLUDE_FILES ${VISUALIZER_DIR}/*.h)
#         FILE(GLOB VISUALIZER_SOURCE_CXX_FILES ${VISUALIZER_DIR}/*.cpp)

#         add_executable(${GUI_NAME} ${VISUALIZER_SOURCE_INCLUDE_FILES} ${VISUALIZER_SOURCE_CXX_FILES})
#         target_include_directories(${GUI_NAME} PRIVATE ${GLUT_INCLUDE_DIR} ${SIMBODY_INCLUDE_DIRS})

#         # TODO do we need to link ${LIBRARIES}?
#         target_link_libraries(${GUI_NAME} PRIVATE ${ROBOSAMPLE_LIB_NAME} ${LIBRARIES} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

#         set_target_properties(${GUI_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
#         target_link_options(${GUI_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
#     ENDIF(BUILD_VISUALIZER)

#     # build tests
#     IF (CMAKE_BUILD_TYPE MATCHES "Debug_Tests" OR CMAKE_BUILD_TYPE MATCHES "Release_Tests" OR
#         CMAKE_BUILD_TYPE MATCHES "Debug_Full" OR CMAKE_BUILD_TYPE MATCHES "Release_Full")

#         enable_testing()
#         include(CTest)

#         SET(TEST_NAMES)
#         SET(TEST_SOURCES)

#         # # find openmm tests
#         # FILE(GLOB TESTS_SOURCES openmm/tests/*.cpp)
#         # FOREACH(test_src ${TESTS_SOURCES})
#         #     # get test name by removing file location and extension
#         #     FILE(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/openmm/tests ${test_src})
#         #     STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

#         #     LIST(APPEND TEST_NAMES ${test_name})
#         #     LIST(APPEND TEST_SOURCES ${test_src})
#         # ENDFOREACH(test_src)

#         # # find simbody tests
#         # FILE(GLOB TESTS_SOURCES Simbody01/Simbody/tests/*.cpp)
#         # FOREACH(test_src ${TESTS_SOURCES})
#         #     # get test name by removing file location and extension
#         #     FILE(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/tests ${test_src})
#         #     STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

#         #     LIST(APPEND TEST_NAMES ${test_name})
#         #     LIST(APPEND TEST_SOURCES ${test_src})
#         # ENDFOREACH(test_src)

#         # find robosample tests
#         FILE(GLOB TESTS_SOURCES tests/*.cpp)
#         FOREACH(test_src ${TESTS_SOURCES})
#             # get test name by removing file location and extension
#             FILE(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/tests ${test_src})
#             STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

#             LIST(APPEND TEST_NAMES ${test_name})
#             LIST(APPEND TEST_SOURCES ${test_src})
#         ENDFOREACH(test_src)

#         # add the tests
#         FOREACH(test_name test_src IN ZIP_LISTS TEST_NAMES TEST_SOURCES)
#             # the test executable is non-ipo
#             add_executable(${test_name} ${test_src})
#             target_compile_options(${test_name} PRIVATE ${BUILD_FLAGS})
#             target_link_options(${test_name} PRIVATE ${SANITIZER_FLAGS})
#             target_include_directories(${test_name} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
#             target_link_libraries(${test_name} PRIVATE ${ROBOSAMPLE_LIB_NAME})

#             add_test(NAME ${test_name} COMMAND ${test_name})
#         ENDFOREACH()
#     ENDIF ()
# endif()

## add example inputs
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp/pdbs)
#file(GLOB INPUTS ${CMAKE_SOURCE_DIR}/tests_inputs/*)
#file(COPY ${INPUTS} DESTINATION ${CMAKE_BINARY_DIR})
#
## Copy Python scripts
#file(COPY ${CMAKE_SOURCE_DIR}/tools/simulate.py DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/tools/simulate_clusters.py DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/tools/simulate_decoys.py DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/tools/flexor.py DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/tools/openmm_hmc.py DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/tools/batstat.py DESTINATION ${CMAKE_BINARY_DIR})

# file(COPY ${CMAKE_SOURCE_DIR}/tools/roborun.py DESTINATION ${CMAKE_BINARY_DIR})
# file(COPY ${CMAKE_SOURCE_DIR}/tools/robobatstat.py DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/datasets/1APQ_correlations.npy DESTINATION ${CMAKE_BINARY_DIR})
