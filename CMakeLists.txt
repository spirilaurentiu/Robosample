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

set(ROBOSAMPLE_BUILD_TYPES Debug Release RelWithDebInfo PGO_Train PGO_Use Tests)

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

if(USE_CUDA)
  find_package(CUDAToolkit REQUIRED)
elseif(USE_OPENCL)
  find_package(OpenCL REQUIRED)
endif()

include(${CMAKE_SOURCE_DIR}/cmake_modules/OpenMM.cmake)

# =============================================================================
# Source files
#
# Globbing is used for convenience. Note: a plain GLOB is only re-evaluated when
# CMake reconfigures, so adding/removing a source file here requires a manual
# reconfigure. (The test discovery further down uses CONFIGURE_DEPENDS, which
# re-globs automatically.)
# =============================================================================

# file(GLOB ROBOSAMPLE_C_SOURCE_FILES ${CMAKE_SOURCE_DIR}/format/src/*.c)

file(GLOB ROBOSAMPLE_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)

file(GLOB ROBOSAMPLE_SOURCE_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/include/*.h*)

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
# Sanitizer flags for the Tests configuration.
#
# -fno-sanitize-recover=all turns every finding into a hard abort (so CI fails
# instead of logging and continuing); -fno-omit-frame-pointer +
# -fsanitize-address-use-after-scope give complete, accurate stack traces.
#
# These are applied ONLY to Robosample's own code (robosample_objects and the
# Robosample test TUs). OpenMM is third-party released code we consume; it is
# compiled into a SEPARATE, un-sanitized object library (openmm_objects) so that
# its known UBSan vptr quirks under hidden-visibility, the CUDA shadow-gap
# address-space OOM, and LeakSanitizer noise from the CUDA driver never reach
# us. See the object-library section below.
# =============================================================================

set(ROBOSAMPLE_SANITIZERS address,undefined)

set(ROBOSAMPLE_SANITIZER_COMPILE_FLAGS
    -fsanitize=${ROBOSAMPLE_SANITIZERS} -fno-sanitize-recover=all
    -fno-omit-frame-pointer -fsanitize-address-use-after-scope)

set(ROBOSAMPLE_SANITIZER_LINK_FLAGS -fsanitize=${ROBOSAMPLE_SANITIZERS})

# =============================================================================
# Compile / link option sets
#
# Split into a BASE set (everything that must be true of ALL compiled code,
# regardless of whether it is ours or OpenMM's) and the full ROBOSAMPLE_* set
# (BASE + sanitizers in the Tests config). OpenMM objects and the OpenMM test
# executables use the BASE sets; our code and our tests use the full sets.
#
# Debug                : no optimization, full debug info as split DWARF
# (-gsplit-dwarf keeps the heavy .debug_* sections out of the link;
# -ggnu-pubnames feeds the gdb index). RelWithDebInfo       : optimized +
# debuggable + gcov instrumentation. Tests                : light optimization
# keeps sanitizer reports fast and readable while preserving good line info.
# Sanitizer instrumentation (added only in the full set, and only for our code)
# is restricted to C/C++ TUs so it never reaches nvcc.
# =============================================================================

set(ROBOSAMPLE_BASE_COMPILE_OPTIONS
    -march=x86-64-v3
    # SIMD-only OpenMP: enables `#pragma omp simd` for the hot per-atom
    # transform loops (RobotEngine) WITHOUT pulling in the OpenMP threading
    # runtime (libgomp). No effect on TUs without the pragma. Robosample's code
    # only (not openmm_objects).
    -fopenmp-simd
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
    ${PGO_USE_FLAGS}>
    $<$<CONFIG:Tests>:-O1
    -g>)

# Robosample's own code additionally gets ASan/UBSan in the Tests config, gated
# to C/C++ TUs via COMPILE_LANGUAGE so the flags are never handed to nvcc.
set(ROBOSAMPLE_COMPILE_OPTIONS
    ${ROBOSAMPLE_BASE_COMPILE_OPTIONS}
    $<$<AND:$<CONFIG:Tests>,$<COMPILE_LANGUAGE:C,CXX>>:${ROBOSAMPLE_SANITIZER_COMPILE_FLAGS}>
)

# Common link options. Debug          : link with mold (much faster on
# debug-info-heavy links) and build a gdb index for fast symbol lookup.
# RelWithDebInfo : gcov needs --coverage at link time too. PGO_*          :
# profile generate/use must also be passed to the linker; --emit-relocs
# preserves relocations for a later BOLT pass.
set(ROBOSAMPLE_BASE_LINK_OPTIONS
    $<$<CONFIG:Debug>:-fuse-ld=mold;-Wl,--gdb-index>
    $<$<CONFIG:RelWithDebInfo>:--coverage>
    $<$<CONFIG:PGO_Train>:-fprofile-generate=${PGO_DATA_ROOT}>
    $<$<CONFIG:PGO_Use>:-fprofile-use=${PGO_DATA_ROOT}
    -Wl,--emit-relocs>)

set(ROBOSAMPLE_LINK_OPTIONS
    ${ROBOSAMPLE_BASE_LINK_OPTIONS}
    $<$<CONFIG:Tests>:${ROBOSAMPLE_SANITIZER_LINK_FLAGS}>)

# =============================================================================
# Object libraries - each compiled ONCE, linked into all consumers
#
# Two object libraries, deliberately separate so the sanitizer boundary falls
# exactly on the our-code / OpenMM line:
#
# robosample_objects : OUR code. Sanitized in the Tests config. Also the carrier
# of the PUBLIC include dirs and USE_<backend> macros, so the
# $<TARGET_PROPERTY:robosample_objects,...> queries below resolve for every
# consumer. openmm_objects     : OpenMM (third-party). NEVER sanitized. This is
# what keeps libasan out of the OpenMM test binaries entirely, which removes the
# CUDA shadow-gap OOM and the UBSan vptr report in CudaArray.cpp.
#
# Mixing sanitized and un-sanitized objects in one ASan binary is supported: the
# runtime intercepts allocation process-wide, so OpenMM allocating / our code
# freeing (and vice versa) is fine.
# =============================================================================

# --- Our code (sanitized) ----------------------------------------------------

# All our sources except the pybind11 entry point (not needed by test exes).
set(ROBOSAMPLE_CXX_NO_BINDINGS ${ROBOSAMPLE_CXX_SOURCE_FILES})
list(REMOVE_ITEM ROBOSAMPLE_CXX_NO_BINDINGS
     ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp)

add_library(robosample_objects OBJECT ${ROBOSAMPLE_CXX_NO_BINDINGS}
                                      ${ROBOSAMPLE_SOURCE_INCLUDE_FILES})

target_include_directories(
  robosample_objects
  PUBLIC $ENV{CONDA_PREFIX}/targets/x86_64-linux/include ${OPENMM_INCLUDE_DIRS}
         ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/pcg-cpp/include
         ${CMAKE_SOURCE_DIR}/pybind11/include ${CMAKE_SOURCE_DIR}/units/include)

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
  PROPERTIES POSITION_INDEPENDENT_CODE ON
             CXX_STANDARD 17
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)

# --- OpenMM (never sanitized) ------------------------------------------------

add_library(
  openmm_objects OBJECT ${OPENMM_SOURCE_C_FILES} ${OPENMM_SOURCE_CXX_FILES}
                        ${OPENMM_SOURCE_INCLUDE_FILES})

# Same include dirs and USE_<backend> macros as our code. Reusing the property
# queries keeps the two libraries in lockstep with no duplication.
target_include_directories(
  openmm_objects
  PRIVATE $<TARGET_PROPERTY:robosample_objects,INCLUDE_DIRECTORIES>)

target_compile_definitions(
  openmm_objects
  PRIVATE $<TARGET_PROPERTY:robosample_objects,COMPILE_DEFINITIONS>)

# BASE options only - no sanitizer flags, ever.
target_compile_options(openmm_objects
                       PRIVATE ${ROBOSAMPLE_BASE_COMPILE_OPTIONS})

set_target_properties(
  openmm_objects
  PROPERTIES POSITION_INDEPENDENT_CODE ON
             CXX_STANDARD 17
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)

# The generated kernel sources (part of OPENMM_SOURCE_CXX_FILES) are compiled
# here now, so the dependency on the generation targets lives on this library.
# Every consumer that references $<TARGET_OBJECTS:openmm_objects> therefore
# transitively waits for the kernels to be generated.
if(USE_CUDA OR USE_OPENCL)
  add_dependencies(openmm_objects ${OPENMM_DEPENDENCIES})
endif()

# =============================================================================
# Python module (pybind11)
#
# Skipped entirely in the Tests configuration: that config instruments our
# objects with sanitizers, and a sanitized extension module cannot be dlopen'd
# by a normal (un-instrumented) Python interpreter without an LD_PRELOAD'd
# runtime. The Tests config therefore builds ONLY the C++ test executables.
# =============================================================================

if(NOT CMAKE_BUILD_TYPE STREQUAL "Tests")

  add_subdirectory(pybind11)

  pybind11_add_module(
    ${ROBOSAMPLE_LIB_NAME}
    ${CMAKE_SOURCE_DIR}/src/PyBind11.cpp # the only file NOT in an object lib
    $<TARGET_OBJECTS:robosample_objects>
    $<TARGET_OBJECTS:openmm_objects>)

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
    target_link_libraries(${ROBOSAMPLE_LIB_NAME}
                          PRIVATE CUDA::cuda_driver CUDA::cufft CUDA::nvrtc)
  elseif(USE_OPENCL)
    target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE OpenCL::OpenCL)
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
  # (instrumented training binaries should not have cross-TU transformations
  # that alter the call graph before profiling).
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

endif() # NOT CMAKE_BUILD_TYPE STREQUAL "Tests"

# =============================================================================
# C++ Tests
# =============================================================================

include(CTest)

# CTest's include() above already defines a cached BUILD_TESTING option and
# calls enable_testing(). We deliberately override it to a normal (scope-local)
# variable so tests are gated strictly on the build type, regardless of the
# cache value.
if(CMAKE_BUILD_TYPE STREQUAL "Tests")
  set(BUILD_TESTING ON)
else()
  set(BUILD_TESTING OFF)
endif()

if(BUILD_TESTING)
  enable_testing()

  message(STATUS "C++ tests enabled (${CMAKE_BUILD_TYPE} build)")

  # Locate GoogleTest (must be installed in the conda environment).
  find_package(GTest REQUIRED)
  include(GoogleTest)

  # ===========================================================================
  # Robosample's own tests (GoogleTest)
  # ===========================================================================

  # robosample_add_test(<test_name> [EXTRA_SOURCES src1 src2 ...])
  #
  # Builds a test executable from tests/<test_name>.cpp plus any EXTRA_SOURCES,
  # reusing the already-compiled object libraries. The only source compiled
  # fresh per test is the Test*.cpp itself. Our code stays sanitized; OpenMM
  # (linked via openmm_objects) does not.
  function(robosample_add_test TEST_NAME)
    cmake_parse_arguments(ARG "" "" "EXTRA_SOURCES" ${ARGN})

    add_executable(
      ${TEST_NAME}
      tests/${TEST_NAME}.cpp ${ARG_EXTRA_SOURCES}
      $<TARGET_OBJECTS:robosample_objects> $<TARGET_OBJECTS:openmm_objects>)

    target_include_directories(
      ${TEST_NAME}
      PRIVATE $<TARGET_PROPERTY:robosample_objects,INCLUDE_DIRECTORIES>)

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
              $<$<BOOL:${USE_CUDA}>:CUDA::cuda_driver
              CUDA::cufft
              CUDA::nvrtc>
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

    # Our tests are ASan binaries AND run sampling on the GPU, so on a GPU
    # backend they need: * RESOURCE_LOCK "gpu" - one at a time, so ctest -j
    # never oversubscribes the single VRAM pool / set of CUDA contexts. *
    # ASAN_OPTIONS=protect_shadow_gap=0 - stop ASan reserving the shadow gap so
    # the CUDA driver can map its virtual address range (otherwise cuCtxCreate
    # reports CUDA_ERROR_OUT_OF_MEMORY even on an empty GPU). * detect_leaks=0 -
    # the CUDA driver's intentional allocations otherwise trip LSan at exit,
    # which -fno-sanitize-recover turns into a hard fail. (Drop this and use an
    # LSan suppressions file if you want to keep leak detection on our own
    # code.)
    set(_robo_test_props LABELS "robosample")
    if(USE_CUDA OR USE_OPENCL)
      list(
        APPEND
        _robo_test_props
        RESOURCE_LOCK
        "gpu"
        ENVIRONMENT
        "ASAN_OPTIONS=protect_shadow_gap=0:detect_leaks=0:new_delete_type_mismatch=0"
      )
    endif()

    gtest_discover_tests(${TEST_NAME} PROPERTIES ${_robo_test_props})
  endfunction()

  # Auto-discover all Test*.cpp files in the subdirectory. CONFIGURE_DEPENDS
  # re-runs the glob on build so new test files are picked up without a manual
  # reconfigure.
  file(GLOB _test_sources CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/Test*.cpp)

  foreach(_src IN LISTS _test_sources)
    cmake_path(GET _src STEM _test_name)
    robosample_add_test(${_test_name})
  endforeach()

  # ---------------------------------------------------------------------------
  # TestRoboticsOracle*-only wiring (docs/specs/robotics-oracle-differential.md
  # §9): the fixtures are <case>.npz + <case>.manifest.json on disk, loaded at
  # test runtime via the vendored cnpy (tests/third_party/cnpy) + a
  # CMake-defined fixture-dir macro. Test-scoped only -- no other test needs
  # cnpy/zlib, so this is added onto the TestRoboticsOracle* targets alone
  # rather than into robosample_add_test()/robosample_objects. Matches BOTH
  # TestRoboticsOracle (Scope A, tests/fixtures/robotics_oracle) and
  # TestRoboticsOracleMolecule (Scope B,
  # tests/fixtures/robotics_oracle_molecules) -- same cnpy/JSON machinery,
  # disjoint fixture directories (§4.2/§9).
  # ---------------------------------------------------------------------------
  foreach(_oracle_target TestRoboticsOracle TestRoboticsOracleMolecule)
    if(TARGET ${_oracle_target})
      find_package(ZLIB REQUIRED)
      target_sources(
        ${_oracle_target}
        PRIVATE ${CMAKE_SOURCE_DIR}/tests/third_party/cnpy/cnpy.cpp)
      target_include_directories(
        ${_oracle_target} PRIVATE ${CMAKE_SOURCE_DIR}/tests/third_party/cnpy
                                  ${CMAKE_SOURCE_DIR}/tests/third_party/sha256)
      target_link_libraries(${_oracle_target} PRIVATE ZLIB::ZLIB)
      target_compile_definitions(
        ${_oracle_target}
        PRIVATE
          ROBOTICS_ORACLE_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures/robotics_oracle"
          ROBOTICS_ORACLE_MOLECULES_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures/robotics_oracle_molecules"
      )
    endif()
  endforeach()

  # ===========================================================================
  # OpenMM's own test suite (platform-independent API + serialization + the one
  # selected backend's platform tests)
  #
  # These are NOT gtest programs: each *Test*.cpp has its own main() and uses
  # OpenMM's assertion macros; the GPU platforms take a precision arg
  # (single/mixed/double) on argv. So they are registered with add_test(), never
  # gtest_discover_tests, and never link gtest_main.
  #
  # They link openmm_objects ONLY (they don't need our code) and use the BASE
  # flag/link sets, so they carry no sanitizer instrumentation and no libasan -
  # which is exactly what removes the CUDA OOM and the UBSan vptr report.
  # ===========================================================================

  set(_omm ${CMAKE_SOURCE_DIR}/openmm)

  # Shared test headers (.h/.hpp) are scattered across these dirs and #included
  # by bare name. Per the project layout every OpenMM test gets the top-level
  # tests dir plus all four platform test dirs on its include path. Dirs not
  # actually #included from are inert, so listing all of them is harmless; we
  # still filter to existing dirs so a partial checkout doesn't hard-error.
  set(OPENMM_TEST_INCLUDE_DIRS)
  foreach(
    _d ${_omm}/tests ${_omm}/platforms/cpu/tests ${_omm}/platforms/cuda/tests
       ${_omm}/platforms/opencl/tests ${_omm}/platforms/reference/tests)
    if(IS_DIRECTORY ${_d})
      list(APPEND OPENMM_TEST_INCLUDE_DIRS ${_d})
    endif()
  endforeach()

  # robosample_add_openmm_test(<src> <target_prefix> <is_gpu>)
  function(robosample_add_openmm_test SRC PREFIX IS_GPU)
    cmake_path(GET SRC STEM _stem)
    set(_tgt "${PREFIX}_${_stem}") # prefix avoids clashes across dirs and with
                                   # our own tests/

    add_executable(${_tgt} ${SRC} $<TARGET_OBJECTS:openmm_objects>)

    target_include_directories(
      ${_tgt} PRIVATE $<TARGET_PROPERTY:robosample_objects,INCLUDE_DIRECTORIES>
                      ${OPENMM_TEST_INCLUDE_DIRS})

    target_compile_definitions(
      ${_tgt} PRIVATE $<TARGET_PROPERTY:robosample_objects,COMPILE_DEFINITIONS>)

    # BASE sets only: no sanitizer compile flags, no sanitizer runtime at link.
    target_compile_options(${_tgt} PRIVATE ${ROBOSAMPLE_BASE_COMPILE_OPTIONS})
    target_link_options(${_tgt} PRIVATE ${ROBOSAMPLE_BASE_LINK_OPTIONS})

    target_link_libraries(
      ${_tgt}
      PRIVATE BLAS::BLAS
              LAPACK::LAPACK
              OpenMP::OpenMP_CXX
              Threads::Threads
              $<$<BOOL:${USE_CUDA}>:CUDA::cuda_driver
              CUDA::cufft
              CUDA::nvrtc>
              $<$<BOOL:${USE_OPENCL}>:OpenCL::OpenCL>)
    # NO GTest::gtest_main - these bring their own main()

    set_target_properties(
      ${_tgt}
      PROPERTIES POSITION_INDEPENDENT_CODE ON
                 CXX_STANDARD 17
                 CXX_STANDARD_REQUIRED ON
                 CXX_EXTENSIONS OFF
                 INTERPROCEDURAL_OPTIMIZATION FALSE)

    if(IS_GPU)
      set(PRECISION mixed) # single double
      foreach(_prec ${PRECISION})
        add_test(NAME ${_tgt}_${_prec} COMMAND ${_tgt} ${_prec})
        set_tests_properties(${_tgt}_${_prec} PROPERTIES LABELS "openmm"
                                                         RESOURCE_LOCK "gpu")
      endforeach()
    endif()

  endfunction()

  # Source roots: always the platform-independent ones, plus exactly the one
  # backend USE_* selected. Encoded "dir|prefix|is_gpu". Only cuda/opencl are
  # GPU (they take the precision arg); cpu/reference/top-level/serialization do
  # not.
  set(_omm_test_roots "${_omm}/tests|openmm|OFF"
                      "${_omm}/serialization/tests|openmm_serialization|OFF")

  if(USE_REFERENCE)
    list(APPEND _omm_test_roots
         "${_omm}/platforms/reference/tests|openmm_reference|OFF")
  elseif(USE_CPU)
    list(APPEND _omm_test_roots "${_omm}/platforms/cpu/tests|openmm_cpu|OFF")
  elseif(USE_CUDA)
    list(APPEND _omm_test_roots "${_omm}/platforms/cuda/tests|openmm_cuda|ON")
  elseif(USE_OPENCL)
    list(APPEND _omm_test_roots
         "${_omm}/platforms/opencl/tests|openmm_opencl|ON")
  endif()

  foreach(_entry IN LISTS _omm_test_roots)
    string(REPLACE "|" ";" _p "${_entry}")
    list(GET _p 0 _dir)
    list(GET _p 1 _prefix)
    list(GET _p 2 _is_gpu)
    if(NOT IS_DIRECTORY ${_dir})
      message(STATUS "OpenMM tests: skipping missing ${_dir}")
      continue()
    endif()
    # OpenMM uses "*Test*.cpp" (Test anywhere), unlike our "Test*.cpp".
    file(GLOB _srcs CONFIGURE_DEPENDS ${_dir}/*Test*.cpp)
    foreach(_src IN LISTS _srcs)
      robosample_add_openmm_test(${_src} ${_prefix} ${_is_gpu})
    endforeach()
  endforeach()

endif(BUILD_TESTING)
