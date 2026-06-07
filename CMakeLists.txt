cmake_minimum_required(VERSION 3.17)

project(Robosample LANGUAGES C CXX)

# =============================================================================
# Options
# =============================================================================

# Mutually exclusive hardware acceleration backends
option(USE_CPU "Enable CPU backend" OFF)
option(USE_REFERENCE "Enable Reference backend" OFF)
option(USE_OPENCL "Enable OpenCL backend" OFF)
option(USE_CUDA "Enable CUDA backend" OFF)

if(USE_CUDA)
  enable_language(CUDA)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    message(FATAL_ERROR "CUDA enabled but CMAKE_CUDA_ARCHITECTURES not set")
  endif()
endif()

# =============================================================================
# Project metadata
# =============================================================================

set(ROBOSAMPLE_MAJOR_VERSION 2)
set(ROBOSAMPLE_MINOR_VERSION 0)
set(ROBOSAMPLE_PATCH_VERSION 0)
set(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")
set(ROBOSAMPLE_LIB_NAME "robo_bindings")

string(TIMESTAMP CURRENT_YEAR "%Y")
set(ROBOSAMPLE_COPYRIGHT_YEARS "2021-${CURRENT_YEAR}")

# =============================================================================
# Compile commands symlink
# =============================================================================

if(CMAKE_EXPORT_COMPILE_COMMANDS)
  add_custom_target(
    update_compile_commands_symlink ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/build/latest
    COMMENT "Updating `latest` symlink to ${CMAKE_BINARY_DIR}")
endif()

# =============================================================================
# Dependencies
# =============================================================================

# OpenBLAS should be the default BLAS implementation in the conda environment
# However, pytorch will pull in MKL
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(
  Python3
  COMPONENTS Interpreter Development
  REQUIRED)

include(${CMAKE_SOURCE_DIR}/cmake_modules/OpenMM.cmake)
include(${CMAKE_SOURCE_DIR}/cmake_modules/Simbody.cmake)
include(${CMAKE_SOURCE_DIR}/cmake_modules/Molmodel.cmake)

# =============================================================================
# Source files
# =============================================================================

file(GLOB ROBOSAMPLE_C_SOURCE_FILES ${CMAKE_SOURCE_DIR}/format/src/*.c)

file(GLOB ROBOSAMPLE_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp
     ${CMAKE_SOURCE_DIR}/format/src/*.cpp)

file(GLOB ROBOSAMPLE_SOURCE_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/include/*.h*
     ${CMAKE_SOURCE_DIR}/format/inc/*.h*
     ${CMAKE_SOURCE_DIR}/format/inc/mdtraj/*.h*)

# =============================================================================
# Build types
# =============================================================================

set(CMAKE_CONFIGURATION_TYPES
    Debug Release RelWithDebInfo PGO_Train PGO_Use
    CACHE STRING "" FORCE)

if(NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release
                                               RelWithDebInfo PGO_Train PGO_Use)
endif()

# =============================================================================
# Profile-Guided Optimization (PGO) flags
# =============================================================================

set(PGO_DATA_ROOT "${CMAKE_SOURCE_DIR}/profile-data")
file(MAKE_DIRECTORY ${PGO_DATA_ROOT})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(PGO_TRAIN_FLAGS
      -fprofile-generate=${PGO_DATA_ROOT}
      -fprofile-prefix-path=${CMAKE_BINARY_DIR} -fprofile-update=atomic
      -fno-reorder-blocks-and-partition)
  set(PGO_USE_FLAGS
      -fprofile-use=${PGO_DATA_ROOT} -fprofile-prefix-path=${CMAKE_BINARY_DIR}
      -fprofile-correction -fprofile-partial-training
      -fno-reorder-blocks-and-partition -fno-profile-reorder-functions)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(PGO_TRAIN_FLAGS -fprofile-instr-generate=${PGO_DATA_ROOT}/code.profraw)
  set(PGO_USE_FLAGS -fprofile-instr-use=${PGO_DATA_ROOT}/code.profdata)
endif()

# Disable pybind11's automatic strip post-build command (interferes with PGO)
set(CMAKE_STRIP
    ""
    CACHE FILEPATH "Strip disabled for PGO" FORCE)

# =============================================================================
# IPO / LTO check
# =============================================================================

include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT REASON)

if(NOT IPO_SUPPORTED)
  message(
    FATAL_ERROR
      "\n"
      "------------------------------------------------------------------\n"
      "FATAL: IPO/LTO is required but not supported.\n"
      "Reason: ${REASON}\n"
      "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\n"
      "Path:     ${CMAKE_CXX_COMPILER}\n"
      "Linker:   ${CMAKE_LINKER}\n"
      "------------------------------------------------------------------\n"
      "Common fixes:\n"
      "  1. Use a modern linker: lld or gold.\n"
      "  2. If using LLVM from Conda: -DCMAKE_LINKER=llvm-ld\n"
      "  3. Ensure llvm-ar and llvm-ranlib are in PATH.\n")
endif()

# =============================================================================
# Shared compile/link options - defined once, reused by all targets
# =============================================================================

# Common compile options (applied to robosample_objects, propagated to
# consumers)
set(ROBOSAMPLE_COMPILE_OPTIONS
    -march=native
    -mtune=native
    -Wno-deprecated-declarations # OpenMM double4 / double4_16a / double4_32a
    $<$<CONFIG:Debug>:
    -fno-omit-frame-pointer
    -O0
    -g>
    $<$<CONFIG:Release>:
    -O3>
    $<$<CONFIG:RelWithDebInfo>:
    -O3
    -g
    -fno-omit-frame-pointer
    --coverage>
    $<$<CONFIG:PGO_Train>:
    -O3
    -g
    -fno-omit-frame-pointer
    ${PGO_TRAIN_FLAGS}>
    $<$<CONFIG:PGO_Use>:
    -O3
    -g
    ${PGO_USE_FLAGS}>)

# Common link options (applied to final link targets: .so and test executables)
set(ROBOSAMPLE_LINK_OPTIONS
    $<$<CONFIG:RelWithDebInfo>:--coverage>
    $<$<CONFIG:PGO_Train>:-fprofile-generate=${PGO_DATA_ROOT}>
    $<$<CONFIG:PGO_Use>: -fprofile-use=${PGO_DATA_ROOT} -Wl,--emit-relocs>)

# =============================================================================
# Object library - compiled ONCE, linked into all consumers
#
# This is the key to fast incremental builds with many test targets: - N test
# executables x M source files = M compilations total (not NxM) - Each test only
# compiles its own Test*.cpp; everything else is reused
# =============================================================================

# All sources except the pybind11 entry point (not needed by test executables)
set(ROBOSAMPLE_CXX_NO_BINDINGS ${ROBOSAMPLE_CXX_SOURCE_FILES})
list(REMOVE_ITEM ROBOSAMPLE_CXX_NO_BINDINGS
     ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp)

add_library(
  robosample_objects OBJECT
  ${OPENMM_SOURCE_C_FILES}
  ${OPENMM_SOURCE_CXX_FILES}
  ${SIMBODY_SOURCE_C_FILES}
  ${SIMBODY_SOURCE_CXX_FILES}
  ${MOLMODEL_C_SOURCE_FILES}
  ${MOLMODEL_CXX_SOURCE_FILES}
  ${ROBOSAMPLE_C_SOURCE_FILES}
  ${ROBOSAMPLE_CXX_NO_BINDINGS}
  # Headers listed for IDE indexing; they are not compiled
  ${OPENMM_SOURCE_INCLUDE_FILES}
  ${SIMBODY_SOURCE_INCLUDE_FILES}
  ${MOLMODEL_SOURCE_INCLUDE_FILES}
  ${ROBOSAMPLE_SOURCE_INCLUDE_FILES})

target_include_directories(
  robosample_objects
  PUBLIC $ENV{CONDA_PREFIX}/targets/x86_64-linux/include
         ${OPENMM_INCLUDE_DIRS}
         ${SIMBODY_INCLUDE_DIRS}
         ${MOLMODEL_INCLUDE_DIRS}
         ${CMAKE_SOURCE_DIR}/include
         ${CMAKE_SOURCE_DIR}/format/inc
         ${CMAKE_SOURCE_DIR}/format/inc/mdtraj
         ${CMAKE_SOURCE_DIR}/pcg-cpp/include
         ${CMAKE_SOURCE_DIR}/pybind11/include)

# NOTE: MODULE_NAME is intentionally NOT set here. It is consumer-specific (the
# .so name and each test name differ), so each consuming target sets its own
# MODULE_NAME.
foreach(backend CPU REFERENCE OPENCL CUDA)
  if(USE_${backend})
    target_compile_definitions(robosample_objects PUBLIC USE_${backend}=1)
  else()
    target_compile_definitions(robosample_objects PUBLIC USE_${backend}=0)
  endif()
endforeach()

target_compile_options(robosample_objects PRIVATE ${ROBOSAMPLE_COMPILE_OPTIONS})

set_target_properties(
  robosample_objects
  PROPERTIES POSITION_INDEPENDENT_CODE ON # Required: objects are shared between
                                          # .so and executables
             CXX_STANDARD 17
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)

# =============================================================================
# Python module (pybind11)
# =============================================================================

add_subdirectory(pybind11)

pybind11_add_module(
  ${ROBOSAMPLE_LIB_NAME}
  ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp # the only file NOT in robosample_objects
  $<TARGET_OBJECTS:robosample_objects>)

target_include_directories(
  ${ROBOSAMPLE_LIB_NAME}
  PRIVATE $<TARGET_PROPERTY:robosample_objects,INCLUDE_DIRECTORIES>)

target_compile_definitions(
  ${ROBOSAMPLE_LIB_NAME}
  PRIVATE MODULE_NAME=${ROBOSAMPLE_LIB_NAME}
          $<TARGET_PROPERTY:robosample_objects,COMPILE_DEFINITIONS>)

# Remove the -s flag that pybind11_add_module injects (interferes with PGO/BOLT)
get_target_property(_link_opts ${ROBOSAMPLE_LIB_NAME} LINK_OPTIONS)
if(_link_opts)
  list(REMOVE_ITEM _link_opts "-s")
  set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES LINK_OPTIONS
                                                          "${_link_opts}")
endif()

target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_LINK_OPTIONS})

target_link_libraries(
  ${ROBOSAMPLE_LIB_NAME} PRIVATE BLAS::BLAS LAPACK::LAPACK OpenMP::OpenMP_CXX
                                 Threads::Threads)

if(USE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE CUDA::cuda_driver
                                                       CUDA::cufft CUDA::nvrtc)
elseif(USE_OPENCL)
  find_package(OpenCL REQUIRED)
  target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE OpenCL::OpenCL)
endif()

if(USE_CUDA OR USE_OPENCL)
  add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
endif()

set_target_properties(
  ${ROBOSAMPLE_LIB_NAME}
  PROPERTIES POSITION_INDEPENDENT_CODE ON
             CXX_STANDARD 17
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)

# IPO: enabled for Release / RelWithDebInfo / PGO_Use; NOT for PGO_Train
# (training binaries should not have cross-TU transformations that alter the
# call graph before profiling)
set_property(TARGET ${ROBOSAMPLE_LIB_NAME}
             PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
set_property(TARGET ${ROBOSAMPLE_LIB_NAME}
             PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
set_property(TARGET ${ROBOSAMPLE_LIB_NAME}
             PROPERTY INTERPROCEDURAL_OPTIMIZATION_PGO_Use TRUE)

# =============================================================================
# Python stubs
# =============================================================================

add_custom_command(
  TARGET ${ROBOSAMPLE_LIB_NAME}
  POST_BUILD
  COMMAND PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR} pybind11-stubgen
          ${ROBOSAMPLE_LIB_NAME} --output-dir ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating Python stubs for ${ROBOSAMPLE_LIB_NAME}")

# =============================================================================
# Install
# =============================================================================

install(TARGETS ${ROBOSAMPLE_LIB_NAME}
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${ROBOSAMPLE_LIB_NAME}.pyi"
        DESTINATION ${CMAKE_INSTALL_PREFIX})

# =============================================================================
# C++ Tests
#
# GoogleTest is resolved via find_package - install it externally: conda install
# -c conda-forge gtest
#
# To add a new test: drop tests/<folder>/Test<name>.cpp into the directory. It
# will be discovered automatically on the next CMake configure. No changes to
# this file are required.
# =============================================================================

include(CTest)

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(BUILD_TESTING ON)
else()
  set(BUILD_TESTING OFF)
endif()

if(BUILD_TESTING)
  enable_testing()

  # Helper: robosample_add_test(Test_<name> [EXTRA_SOURCES src1 src2 ...])
  #
  # Creates a test executable that reuses robosample_objects (.o files already
  # compiled). The only source compiled fresh per test is
  # tests/<folder>/Test<name>.cpp.
  #
  function(robosample_add_test FOLDER_NAME TEST_NAME)
    cmake_parse_arguments(ARG "" "" "EXTRA_SOURCES" ${ARGN})

    add_executable(
      ${TEST_NAME} tests/${FOLDER_NAME}/${TEST_NAME}.cpp ${ARG_EXTRA_SOURCES}
                   $<TARGET_OBJECTS:robosample_objects>)

    target_include_directories(
      ${TEST_NAME}
      PRIVATE $<TARGET_PROPERTY:robosample_objects,INCLUDE_DIRECTORIES>
              ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/tests)

    target_compile_definitions(
      ${TEST_NAME}
      PRIVATE MODULE_NAME=${TEST_NAME}
              $<TARGET_PROPERTY:robosample_objects,COMPILE_DEFINITIONS>)

    target_compile_options(${TEST_NAME} PRIVATE ${ROBOSAMPLE_COMPILE_OPTIONS})

    target_link_options(${TEST_NAME} PRIVATE ${ROBOSAMPLE_LINK_OPTIONS})

    target_link_libraries(
      ${TEST_NAME}
      PRIVATE BLAS::BLAS
              LAPACK::LAPACK
              OpenMP::OpenMP_CXX
              Threads::Threads
              GTest::gtest_main
              $<$<BOOL:${USE_CUDA}>:CUDA::cudart;CUDA::cublas;CUDA::cufft>
              $<$<BOOL:${USE_OPENCL}>:OpenCL::OpenCL>)

    # LTO disabled for tests: faster links, and profile data from tests
    # shouldn't pollute the production PGO corpus anyway
    set_target_properties(
      ${TEST_NAME}
      PROPERTIES POSITION_INDEPENDENT_CODE ON
                 CXX_STANDARD 17
                 CXX_STANDARD_REQUIRED ON
                 CXX_EXTENSIONS OFF
                 INTERPROCEDURAL_OPTIMIZATION FALSE)

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    gtest_discover_tests(${TEST_NAME})
  endfunction()

  # Locate GoogleTest (must be installed in the conda environment)
  find_package(GTest REQUIRED)
  include(GoogleTest)

  # Define the list of subdirectories to search
  set(TEST_FOLDERS SimTKMath SimTKCommon Simbody)

  foreach(folder ${TEST_FOLDERS})
    # Auto-discover all Test*.cpp files in the subdirectory
    file(GLOB _test_sources CONFIGURE_DEPENDS
         ${CMAKE_SOURCE_DIR}/tests/${folder}/Test*.cpp)

    foreach(_src IN LISTS _test_sources)
      cmake_path(GET _src STEM _test_name)
      robosample_add_test(${folder} ${_test_name})
    endforeach()
  endforeach()

endif(BUILD_TESTING)
