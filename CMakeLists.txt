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

# # Ensure the Python extension uses relative RPATH ($ORIGIN)
# # This prevents Conda's lib directory from shadowing system libraries
# set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
# set(CMAKE_INSTALL_RPATH "$ORIGIN")

# Link prerequisites
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Threads REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

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
pybind11_add_module(${ROBOSAMPLE_LIB_NAME} ${ROBOSAMPLE_SOURCE_FILES})

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

# Install Python bindings (.so)
install(TARGETS ${ROBOSAMPLE_LIB_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
