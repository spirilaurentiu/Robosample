# =============================================================================
# Robosample — top-level build
#
# WHAT THIS BUILDS One Python extension module (robo_bindings) and, in the
# RelWithDebInfo config only, a suite of C++ unit tests.
#
# Every C/C++ source — the bundled OpenMM / Simbody / Molmodel forks and
# Robosample's own code — is compiled ONCE into an OBJECT library
# (robosample_objects). Those pre-compiled objects are then linked into both the
# Python module and each test executable. Compiling once and reusing the objects
# is what keeps incremental builds fast: editing one file recompiles only that
# file, then relinks.
#
# BUILD TYPES (single-config generators such as Make / Ninja; pick with
# -DCMAKE_BUILD_TYPE=<type>) Debug          -O0 + full debug info, written as
# split DWARF and linked with mold for fast relinks. Day-to-day development
# build; step through any source, including the bundled forks. Release -O3 +
# LTO. Production build, no debug info. RelWithDebInfo -O3 + debug info + gcov
# coverage. The ONLY config that builds and runs the C++ tests. PGO_Train -O3,
# instrumented; writes profile data into ./profile-data. PGO_Use        -O3 +
# LTO; consumes profile data from a prior PGO_Train run, and emits relocations
# for an optional post-link BOLT pass.
#
# Example:  cmake -B build -DCMAKE_BUILD_TYPE=Debug -DUSE_CUDA=ON \
# -DCMAKE_CUDA_ARCHITECTURES=86
#
# HARDWARE BACKENDS Exactly one of USE_CPU / USE_REFERENCE / USE_OPENCL /
# USE_CUDA must be ON (enforced below). This selects the OpenMM platform that
# gets compiled in.
#
# REQUIREMENTS * A toolchain whose linker supports IPO/LTO (verified below). *
# ccache (required). * A conda environment providing OpenBLAS, OpenMP, Python3,
# and — for the tests — GoogleTest.
# =============================================================================

cmake_minimum_required(VERSION 3.17)

project(Robosample LANGUAGES C CXX)

# =============================================================================
# Options
# =============================================================================

# Hardware acceleration backends. These are mutually exclusive: choose the one
# matching the OpenMM platform you want. The check below enforces "exactly one".
option(USE_CPU "Enable CPU backend" OFF)
option(USE_REFERENCE "Enable Reference backend" OFF)
option(USE_OPENCL "Enable OpenCL backend" OFF)
option(USE_CUDA "Enable CUDA backend" OFF)

# Enforce mutual exclusivity. Without this guard a typo (two backends on, or
# none) configures cleanly and then misbehaves at runtime, which is hard to
# trace back to the build. Count how many are enabled and fail fast otherwise.
set(_enabled_backends 0)
foreach(backend CPU REFERENCE OPENCL CUDA)
  if(USE_${backend})
    math(EXPR _enabled_backends "${_enabled_backends} + 1")
  endif()
endforeach()
if(NOT _enabled_backends EQUAL 1)
  message(
    FATAL_ERROR
      "Exactly one hardware backend must be enabled, but ${_enabled_backends} "
      "are set. Pass exactly one of "
      "-DUSE_CPU=ON / -DUSE_REFERENCE=ON / -DUSE_OPENCL=ON / -DUSE_CUDA=ON.")
endif()

if(USE_CUDA)
  enable_language(CUDA)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    message(FATAL_ERROR "CUDA enabled but CMAKE_CUDA_ARCHITECTURES not set")
  endif()
endif()

# =============================================================================
# Project metadata
#
# NOTE: these values are informational only right now — they are not yet wired
# into project(... VERSION ...), configure_file(), or install metadata. Left in
# place for future use / tooling.
# =============================================================================

set(ROBOSAMPLE_MAJOR_VERSION 2)
set(ROBOSAMPLE_MINOR_VERSION 0)
set(ROBOSAMPLE_PATCH_VERSION 0)
set(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")
set(ROBOSAMPLE_LIB_NAME "robo_bindings")

string(TIMESTAMP CURRENT_YEAR "%Y")
set(ROBOSAMPLE_COPYRIGHT_YEARS "2021-${CURRENT_YEAR}")

# =============================================================================
# Build type selection
#
# This project uses custom build types (PGO_Train / PGO_Use) on top of the
# standard ones. How the build type is chosen depends on the generator: *
# Multi-config generators (Ninja Multi-Config, etc.) pick the config at BUILD
# time, so we populate CMAKE_CONFIGURATION_TYPES. * Single-config generators
# (Make / plain Ninja) pick it at CONFIGURE time via CMAKE_BUILD_TYPE, so we set
# a sensible default and offer the list as GUI dropdown choices.
# =============================================================================

set(ROBOSAMPLE_BUILD_TYPES Debug Release RelWithDebInfo PGO_Train PGO_Use)

get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(_is_multi_config)
  set(CMAKE_CONFIGURATION_TYPES
      "${ROBOSAMPLE_BUILD_TYPES}"
      CACHE STRING "" FORCE)
else()
  if(NOT CMAKE_BUILD_TYPE)
    # No build type given on the command line -> default to Release rather than
    # an empty/flagless build (an empty CMAKE_BUILD_TYPE matches none of the
    # generator expressions below, so you would get no -O, no -g, and no tests).
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                                               ${ROBOSAMPLE_BUILD_TYPES})
endif()

# =============================================================================
# Compile commands symlink
#
# When compile_commands.json export is on, keep build/latest pointing at the
# active binary dir so editors / clangd always find the newest database.
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

# ccache caches object files so unchanged translation units are not recompiled
# across reconfigures and branch switches. Required: the build assumes it.
find_program(CCACHE ccache REQUIRED)
message(STATUS "Using ccache: ${CCACHE}")
set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE})

# OpenBLAS should be the default BLAS implementation in the conda environment.
# However, pytorch will pull in MKL.
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(
  Python3
  COMPONENTS Interpreter Development
  REQUIRED)

# These modules define the OPENMM_* / SIMBODY_* / MOLMODEL_* source and include
# variables consumed by the object library below.
include(${CMAKE_SOURCE_DIR}/cmake_modules/OpenMM.cmake)
include(${CMAKE_SOURCE_DIR}/cmake_modules/Simbody.cmake)
include(${CMAKE_SOURCE_DIR}/cmake_modules/Molmodel.cmake)

# =============================================================================
# Source files
#
# Globbing is used for convenience. Note: a plain GLOB is only re-evaluated when
# CMake reconfigures, so adding/removing a source file here requires a manual
# reconfigure. (The test discovery further down uses CONFIGURE_DEPENDS, which
# re-globs automatically.)
# =============================================================================

file(GLOB ROBOSAMPLE_C_SOURCE_FILES ${CMAKE_SOURCE_DIR}/format/src/*.c)

file(GLOB ROBOSAMPLE_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp
     ${CMAKE_SOURCE_DIR}/format/src/*.cpp)

file(GLOB ROBOSAMPLE_SOURCE_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/include/*.h*
     ${CMAKE_SOURCE_DIR}/format/inc/*.h*
     ${CMAKE_SOURCE_DIR}/format/inc/mdtraj/*.h*)

# =============================================================================
# Profile-Guided Optimization (PGO) flags
#
# PGO is a two-pass workflow: 1. Build PGO_Train, run a representative workload
# -> profile data lands in ${PGO_DATA_ROOT}. 2. Build PGO_Use, which reads that
# data to guide optimization. The flag sets differ between GCC and Clang.
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

# Disable pybind11's automatic strip post-build command (interferes with PGO).
set(CMAKE_STRIP
    ""
    CACHE FILEPATH "Strip disabled for PGO" FORCE)

# =============================================================================
# IPO / LTO check
#
# Release / RelWithDebInfo / PGO_Use rely on link-time optimization, so a
# toolchain without working IPO is a hard error. Fail now with guidance rather
# than at link time with a cryptic plugin error.
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
      "  1. Use a linker with LTO-plugin support: mold, gold, or lld.\n"
      "  2. Ensure the matching archiver/ranlib are on PATH\n"
      "     (gcc-ar / gcc-ranlib for GCC, llvm-ar / llvm-ranlib for Clang).\n"
      "  3. For GCC, confirm the lto-wrapper / linker plugin is installed.\n")
endif()

# =============================================================================
# Shared compile/link options — defined once, reused by all targets
#
# Generator-expression note: each $<$<CONFIG:X>:...> block keeps its first flag
# on the same line as the colon. Starting the block with a newline would inject
# an empty leading argument into the flag list. Harmless in practice, but this
# form is cleaner and unambiguous.
# =============================================================================

# Common compile options (applied to robosample_objects, propagated to
# consumers). -march/-mtune=native : tune for the build host's CPU (all
# configs). Debug                : no optimization, full debug info as split
# DWARF (-gsplit-dwarf keeps the heavy .debug_* sections out of the link;
# -ggnu-pubnames feeds the gdb index). RelWithDebInfo       : optimized +
# debuggable + gcov instrumentation.
set(ROBOSAMPLE_COMPILE_OPTIONS
    -march=native
    -mtune=native
    -Wno-deprecated-declarations # OpenMM double4 / double4_16a / double4_32a
    $<$<CONFIG:Debug>:-fno-omit-frame-pointer
    -O0
    -g
    -gsplit-dwarf
    -ggnu-pubnames>
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O3
    -g
    -fno-omit-frame-pointer
    --coverage>
    $<$<CONFIG:PGO_Train>:-O3
    -g
    -fno-omit-frame-pointer
    ${PGO_TRAIN_FLAGS}>
    $<$<CONFIG:PGO_Use>:-O3
    -g
    ${PGO_USE_FLAGS}>)

# Common link options (applied to final link targets: .so and test executables).
# Debug          : link with mold (much faster on debug-info-heavy links) and
# build a gdb index for fast symbol lookup. RelWithDebInfo : gcov needs
# --coverage at link time too. PGO_*          : profile generate/use must also
# be passed to the linker; --emit-relocs preserves relocations for a later BOLT
# pass.
set(ROBOSAMPLE_LINK_OPTIONS
    $<$<CONFIG:Debug>:-fuse-ld=mold;-Wl,--gdb-index>
    $<$<CONFIG:RelWithDebInfo>:--coverage>
    $<$<CONFIG:PGO_Train>:-fprofile-generate=${PGO_DATA_ROOT}>
    $<$<CONFIG:PGO_Use>:-fprofile-use=${PGO_DATA_ROOT}
    -Wl,--emit-relocs>)

# =============================================================================
# Object library — compiled ONCE, linked into all consumers
#
# Why an OBJECT library: * N test executables share the same M source files. *
# As an object library, those M files compile a single time (M compiles), not
# once per executable (which would be N*M). * Each test then only compiles its
# own Test*.cpp; everything else is reused.
# =============================================================================

# All sources except the pybind11 entry point (not needed by test executables).
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
  # Headers are listed only so IDEs index them; they are not compiled.
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

# Expose every backend as a USE_<backend>=0/1 macro so the C++ can branch on all
# four unconditionally. (Exactly one is 1, enforced at the top.)
#
# MODULE_NAME is intentionally NOT set here: it is consumer-specific (the .so
# and each test use a different name), so each consuming target sets its own.
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

# Reuse the object library's include dirs and compile definitions so the entry
# point is compiled with the same configuration as everything else.
target_include_directories(
  ${ROBOSAMPLE_LIB_NAME}
  PRIVATE $<TARGET_PROPERTY:robosample_objects,INCLUDE_DIRECTORIES>)

target_compile_definitions(
  ${ROBOSAMPLE_LIB_NAME}
  PRIVATE MODULE_NAME=${ROBOSAMPLE_LIB_NAME}
          $<TARGET_PROPERTY:robosample_objects,COMPILE_DEFINITIONS>)

# Remove the -s (strip) flag that pybind11_add_module injects; stripping
# interferes with PGO and a later BOLT pass.
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

# IPO/LTO is enabled for Release / RelWithDebInfo / PGO_Use, but NOT PGO_Train
# (instrumented training binaries should not have cross-TU transformations that
# alter the call graph before profiling).
#
# IMPORTANT: the per-config property suffix is the configuration name in UPPER
# CASE — so PGO_Use -> INTERPROCEDURAL_OPTIMIZATION_PGO_USE. A mixed-case suffix
# silently never matches, leaving LTO off for that config.
set_property(TARGET ${ROBOSAMPLE_LIB_NAME}
             PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
set_property(TARGET ${ROBOSAMPLE_LIB_NAME}
             PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
set_property(TARGET ${ROBOSAMPLE_LIB_NAME}
             PROPERTY INTERPROCEDURAL_OPTIMIZATION_PGO_USE TRUE)

# =============================================================================
# Python stubs
#
# Generate a .pyi type stub next to the module after each build so editors get
# autocompletion for the compiled bindings.
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
# Tests build ONLY in the RelWithDebInfo configuration (see the guard below).
#
# GoogleTest is resolved via find_package — install it externally: conda install
# -c conda-forge gtest
#
# To add a new test: drop tests/<folder>/Test<name>.cpp into one of the folders
# listed in TEST_FOLDERS. It is discovered automatically on the next configure
# (the glob uses CONFIGURE_DEPENDS). No edits to this file are required.
# =============================================================================

include(CTest)

# CTest's include() above already defines a cached BUILD_TESTING option and
# calls enable_testing(). We deliberately override it to a normal (scope-local)
# variable so tests are gated strictly on the build type, regardless of the
# cache value.
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(BUILD_TESTING ON)
else()
  set(BUILD_TESTING OFF)
endif()

if(BUILD_TESTING)
  enable_testing()

  # robosample_add_test(<folder> <test_name> [EXTRA_SOURCES src1 src2 ...])
  #
  # Builds a test executable from tests/<folder>/<test_name>.cpp plus any
  # EXTRA_SOURCES, reusing the already-compiled robosample_objects. The only
  # source compiled fresh per test is the Test*.cpp itself.
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

    # LTO disabled for tests: faster links, and profile data from tests should
    # not pollute the production PGO corpus anyway.
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

  # Locate GoogleTest (must be installed in the conda environment).
  find_package(GTest REQUIRED)
  include(GoogleTest)

  # Subdirectories under tests/ to scan for Test*.cpp.
  set(TEST_FOLDERS SimTKMath SimTKCommon Simbody)

  foreach(folder ${TEST_FOLDERS})
    # Auto-discover all Test*.cpp files in the subdirectory. CONFIGURE_DEPENDS
    # re-runs the glob on build so new test files are picked up without a manual
    # reconfigure.
    file(GLOB _test_sources CONFIGURE_DEPENDS
         ${CMAKE_SOURCE_DIR}/tests/${folder}/Test*.cpp)

    foreach(_src IN LISTS _test_sources)
      cmake_path(GET _src STEM _test_name)
      robosample_add_test(${folder} ${_test_name})
    endforeach()
  endforeach()

endif(BUILD_TESTING)
