cmake_minimum_required(VERSION 3.17)

# CUDA may be used as a platform for OpenMM
# Setting it here simplifies the logic even if CUDA is not used as the OpenMM platform
project(Robosample LANGUAGES C CXX)

# Define mutually exclusive hardware acceleration options
option(USE_CPU          "Enable CPU backend"        OFF)
option(USE_REFERENCE    "Enable Reference backend"  OFF)
option(USE_OPENCL       "Enable OpenCL backend"     OFF)
option(USE_CUDA         "Enable CUDA backend"       OFF)

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

# Create a symlink from build/latest to the current binary directory
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    add_custom_target(update_compile_commands_symlink ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink 
            ${CMAKE_BINARY_DIR} 
            ${CMAKE_SOURCE_DIR}/build/latest
        COMMENT "Updating `latest` symlink to ${CMAKE_BINARY_DIR}"
    )
endif()

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

# Include headers (.h and .hpp)
file(GLOB include_files ${CMAKE_SOURCE_DIR}/include/*.h* ${CMAKE_SOURCE_DIR}/format/inc/*.h* ${CMAKE_SOURCE_DIR}/format/inc/mdtraj/*.h*)
set(ROBOSAMPLE_HEADER_FILES ${include_files})

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

# Allow custom build types
set(CMAKE_CONFIGURATION_TYPES Debug Release RelWithDebInfo PGO_Train PGO_Use CACHE STRING "" FORCE)

# For single-config generators (Ninja, Make)
if(NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      Debug Release RelWithDebInfo PGO_Train PGO_Use)
endif()

# Profile-guided optimization (PGO) flags for different stages and compilers
set(PGO_DATA_ROOT "${CMAKE_SOURCE_DIR}/profile-data")
file(MAKE_DIRECTORY ${PGO_DATA_ROOT})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(PGO_TRAIN_FLAGS
        -fprofile-generate=${PGO_DATA_ROOT}
        -fprofile-prefix-path=${CMAKE_BINARY_DIR}
        -fno-reorder-blocks-and-partition
    )
    set(PGO_USE_FLAGS
        -fprofile-use=${PGO_DATA_ROOT}
        -fprofile-prefix-path=${CMAKE_BINARY_DIR}
        -fprofile-correction
        -fno-reorder-blocks-and-partition
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(PGO_TRAIN_FLAGS
        -fprofile-instr-generate=${PGO_DATA_ROOT}/code.profraw
    )
    set(PGO_USE_FLAGS
        -fprofile-instr-use=${PGO_DATA_ROOT}/code.profdata
    )
endif()

# Disables pybind11's automatic strip post-build command
set(CMAKE_STRIP "" CACHE FILEPATH "Strip disabled for PGO" FORCE)

# Create the Python module
add_subdirectory(pybind11)
# add_library(${ROBOSAMPLE_LIB_NAME} SHARED ${ROBOSAMPLE_SOURCE_FILES})
# target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE pybind11::module)
pybind11_add_module(${ROBOSAMPLE_LIB_NAME} ${ROBOSAMPLE_SOURCE_FILES})

# Add preprocessor definitions
target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE MODULE_NAME=${ROBOSAMPLE_LIB_NAME})
foreach(backend CPU REFERENCE OPENCL CUDA)
    if(USE_${backend})
        target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE USE_${backend}=1)
    else()
        target_compile_definitions(${ROBOSAMPLE_LIB_NAME} PRIVATE USE_${backend}=0)
    endif()
endforeach()

# Must come AFTER pybind11_add_module (removes the -s it just added)
get_target_property(_link_opts ${ROBOSAMPLE_LIB_NAME} LINK_OPTIONS)
if(_link_opts)
    list(REMOVE_ITEM _link_opts "-s")
    set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES LINK_OPTIONS "${_link_opts}")
endif()

# Enable Interprocedural Optimization (IPO/LTO) only for the final optimized configurations.
# We deliberately DO NOT enable IPO for PGO_Train.
# Rationale:
# - The training build should collect representative execution profiles
#   without aggressive whole-program transformations.
# - IPO (LTO) performs cross-TU inlining and function importing, which
#   alters call graphs and control-flow structure.
# - If the training binary structure differs significantly from the
#   PGO_Use binary, profile quality and mapping accuracy may degrade.
# - Therefore:
#       PGO_Train  -> instrumentation only
#       PGO_Use    -> profile use + IPO
#       Release    -> IPO (no PGO)
#
# This mirrors LLVM's recommended PGO workflow:
#   1. Instrument
#   2. Run workload
#   3. Rebuild optimized (PGO + LTO)
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT REASON)

if(NOT IPO_SUPPORTED)
    message(FATAL_ERROR 
        "\n"
        "------------------------------------------------------------------\n"
        "FATAL: IPO/LTO (Interprocedural Optimization) is required but not supported.\n"
        "Reason: ${REASON}\n"
        "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\n"
        "Path:     ${CMAKE_CXX_COMPILER}\n"
        "Linker:   ${CMAKE_LINKER}\n"
        "------------------------------------------------------------------\n"
        "Common Fixes:\n"
        "1. Ensure you are using a modern linker (e.g., lld or gold).\n"
        "2. If using LLVM from Conda, set: -DCMAKE_LINKER=llvm-ld\n"
        "3. Check if the 'llvm-ar' and 'llvm-ranlib' wrappers are in your PATH.\n"
    )
endif()

set_property(TARGET ${ROBOSAMPLE_LIB_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
set_property(TARGET ${ROBOSAMPLE_LIB_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION_PGO_USE TRUE)

# Set properties for the library target
set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 14
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)

# Apply common includes
target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})

# Coverage
if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set_source_files_properties(${ROBOSAMPLE_HEADER_FILES} PROPERTIES COMPILE_FLAGS "--coverage -fprofile-arcs -ftest-coverage")
    set_source_files_properties(${ROBOSAMPLE_C_SOURCE_FILES} PROPERTIES COMPILE_FLAGS "--coverage -fprofile-arcs -ftest-coverage")
    set_source_files_properties(${ROBOSAMPLE_CXX_SOURCE_FILES} PROPERTIES COMPILE_FLAGS "--coverage -fprofile-arcs -ftest-coverage")
endif()

# Set configuration-specific compile options
target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE
    -march=native #x86-64-v3
    -mtune=native #generic
    -Wno-c++17-extensions

    # Debug
    $<$<CONFIG:Debug>:
        -O0
        -g
        -fno-omit-frame-pointer
    >

    # Release
    $<$<CONFIG:Release>:
        -O3
    >

    # Release with Debug Info
    $<$<CONFIG:RelWithDebInfo>:
        -O3
        -g3
        -fno-omit-frame-pointer
    >

    # PGO Train (GCC)
    $<$<CONFIG:PGO_Train>:
        -O3
        -g
        -fno-omit-frame-pointer
        ${PGO_TRAIN_FLAGS}
    >

    # PGO Use (GCC) - note that this binary will be optimized via LLVM-BOLT (perf)
    $<$<CONFIG:PGO_Use>:
        -O3
        -g
        -fno-omit-frame-pointer
        ${PGO_USE_FLAGS}
    >
)

# Set configuration-specific link options
target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE

    # Debug
    $<$<CONFIG:Debug>:
        $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:--coverage>
    >

    # Release
    $<$<CONFIG:Release>:
        # No extra flags
    >

    # Release with Debug Info
    $<$<CONFIG:RelWithDebInfo>:
        $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:--coverage>
    >

    # PGO Train
    $<$<CONFIG:PGO_TRAIN>:
        ${PGO_TRAIN_FLAGS}
    >

    # PGO Use - emit relocs for BOLT, keep debug info through LTO link
    $<$<CONFIG:PGO_USE>:
        ${PGO_USE_FLAGS}
        -g
        -Wl,--emit-relocs
    >
)

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
    target_link_libraries(${ROBOSAMPLE_LIB_NAME}
        PRIVATE
        CUDA::cuda_driver 
        CUDA::cufft
        CUDA::nvrtc
    )
elseif(USE_OPENCL)
    find_package(OpenCL REQUIRED)
    target_link_libraries(${ROBOSAMPLE_LIB_NAME}
        PRIVATE
        OpenCL::OpenCL
    )
endif()

# Optional: dependency on OpenMM build target
if(USE_CUDA OR USE_OPENCL)
    add_dependencies(${ROBOSAMPLE_LIB_NAME} ${OPENMM_DEPENDENCIES})
endif()

# Generate stubs for Python bindings using pybind11-stubgen
# pybind11-stubgen requires PYTHONPATH to find the built module, so we set it to build/*
add_custom_command(TARGET ${ROBOSAMPLE_LIB_NAME} POST_BUILD
    COMMAND PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR} pybind11-stubgen ${ROBOSAMPLE_LIB_NAME} --output-dir ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating python stubs for ${ROBOSAMPLE_LIB_NAME}"
)

# Install binaries and stubs
install(TARGETS ${ROBOSAMPLE_LIB_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${ROBOSAMPLE_LIB_NAME}.pyi" DESTINATION ${CMAKE_INSTALL_PREFIX})
