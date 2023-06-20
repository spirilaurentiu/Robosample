cmake_minimum_required(VERSION 3.10)

project(Robosample C CXX)

# robosample version and date
SET(ROBOSAMPLE_MAJOR_VERSION 2)
SET(ROBOSAMPLE_MINOR_VERSION 0)
SET(ROBOSAMPLE_PATCH_VERSION 0)
SET(ROBOSAMPLE_COPYRIGHT_YEARS "2021-01")

# underbar separated list of dotted authors, no spaces or commas
SET(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")

# set the sanitizer flags
# https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html
# https://clang.llvm.org/docs/AddressSanitizer.html
set(SANITIZER_FLAGS "-fsanitize=undefined,leak,address,float-divide-by-zero -fsanitize-address-use-after-scope -fsanitize-recover=address")

# OpenMM wants these set for some library
set(COMMON_FLAGS "-msse4.1 -mavx")

# -O2 is marginally slower than -O3
set(EXTRA_COMPILE_FLAGS_DEBUG   "${COMMON_FLAGS} -O0 -g")
set(EXTRA_COMPILE_FLAGS_RELEASE "${COMMON_FLAGS} -O3 -march=native -mtune=native")

# set the profile guided optimization (pgo) flags
set(PGO_DIR ${CMAKE_BINARY_DIR}/profile-data)
set(PGO_FILE ${CMAKE_BINARY_DIR}/profile-data/default.profraw)

if (CMAKE_COMPILER_IS_GNUCXX)
    SET(PGO_GENERATE "-fprofile-generate=${PGO_DIR} -fno-reorder-blocks-and-partition")
    SET(PGO_USE      "-fprofile-use=${PGO_DIR} -fprofile-correction -fno-reorder-blocks-and-partition")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(PGO_CLANG_RAW ${CMAKE_BINARY_DIR}/profile-data/default.profraw)
    set(PGO_CLANG_DATA ${CMAKE_BINARY_DIR}/profile-data/default.profdata)

    SET(PGO_GENERATE "-fprofile-instr-generate=${PGO_CLANG_RAW} -fcoverage-mapping") # -fprofile-arcs -ftest-coverage
    SET(PGO_USE      "-fprofile-instr-use=${PGO_CLANG_DATA}")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}   ${EXTRA_COMPILE_FLAGS_DEBUG} -fno-omit-frame-pointer ${SANITIZER_FLAGS}")
    set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_DEBUG} -fno-omit-frame-pointer ${SANITIZER_FLAGS}")
    set(PREFIX ".d")
    # set(INCREMENTAL_LINKING "-Wl,-i,-no-pie")
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}   ${EXTRA_COMPILE_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    set(BOLT_FLAGS "-Wl,--emit-relocs")
elseif(CMAKE_BUILD_TYPE MATCHES "PGO_Train")
    # note: using the sanitizer flags in debug builds make opencl fail on initialization
    set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}   ${EXTRA_COMPILE_FLAGS_DEBUG} ${PGO_GENERATE}")
    set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_DEBUG} ${PGO_GENERATE}")
    set(PREFIX ".pgo.train")
elseif(CMAKE_BUILD_TYPE MATCHES "PGO_Use")
    set(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS}   ${EXTRA_COMPILE_FLAGS_RELEASE} ${PGO_USE}")
    set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE} ${PGO_USE}")
    set(PREFIX ".pgo.use")
    set(BOLT_FLAGS "-Wl,--emit-relocs")
endif()

# set names of the targets we will build
set(OPENMM_LIB_NAME "openmm${PREFIX}")
set(SIMBODY_LIB_NAME "simbody${PREFIX}")
set(MOLMODEL_LIB_NAME "molmodel${PREFIX}")
set(ROBOSAMPLE_NAME "robosample${PREFIX}")

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(LIBRARY_TYPE SHARED)
else()
    set(LIBRARY_TYPE STATIC)
endif()

# we require ipo (interprocedural optimization) aka lto (link time optimization)
# this is fancy -flto basically
# on clang, it sets -flto=thin and this seems to work better than the default -flto(=full)
# even though the binary is larger, it runs faster
include(CheckIPOSupported)
check_ipo_supported()

if (CMAKE_BUILD_TYPE MATCHES "Release|PGO_Use")
    set(USE_IPO TRUE)
else()
    set(USE_IPO FALSE)
endif()

# set(CMAKE_CXX_STANDARD 14)


# link prerequisites
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(DL_LIBRARY dl REQUIRED)
find_library(MATH_LIBRARY m REQUIRED)
find_library(RT_LIBRARY rt REQUIRED) # wants -lrt

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

set(LIBRARIES
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${RT_LIBRARY}
    Threads::Threads
    ${DL_LIBRARY}
    ${MATH_LIBRARY}
)

# get openmm source directories
set(OPENMM_INCLUDE_DIRS)
SET(OPENMM_DIRS
    # ${CMAKE_SOURCE_DIR}/openmm/
    ${CMAKE_SOURCE_DIR}/openmm/openmmapi
    ${CMAKE_SOURCE_DIR}/openmm/olla
    ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit
    ${CMAKE_SOURCE_DIR}/openmm/libraries/jama
    ${CMAKE_SOURCE_DIR}/openmm/libraries/quern
    ${CMAKE_SOURCE_DIR}/openmm/libraries/lepton
    ${CMAKE_SOURCE_DIR}/openmm/libraries/sfmt
    ${CMAKE_SOURCE_DIR}/openmm/libraries/lbfgs
    ${CMAKE_SOURCE_DIR}/openmm/libraries/hilbert
    ${CMAKE_SOURCE_DIR}/openmm/libraries/csha1
    ${CMAKE_SOURCE_DIR}/openmm/libraries/irrxml
    ${CMAKE_SOURCE_DIR}/openmm/libraries/vecmath
    # ${CMAKE_SOURCE_DIR}/openmm/libraries/pthreads
    ${CMAKE_SOURCE_DIR}/openmm/platforms/reference
    ${CMAKE_SOURCE_DIR}/openmm/serialization
)

# all openmm platforms come with a function called registerPlatforms
# each platforms resides in its own shared object and this function initializes that platform when the .so is loaded
# because we build everything into one file, we can only build one platform at a time while the reference platform is excluded from this list
# thus, we can only roll with reference + cpu / cuda / opencl
# cpu will compile and run on any machine so it is a good starting point
# to change it, run `cmake -DOPENMM_PLATFORM=CUDA`
SET(OPENMM_PLATFORM "OPENCL" CACHE STRING "OpenMM platform: CPU, CUDA or OpenCL.")

IF(OPENMM_PLATFORM MATCHES "CPU")
    # build the default CPU platform
    message(STATUS "Building OpenMM CPU platform")
    SET(OPENMM_DIRS ${OPENMM_DIRS} ${CMAKE_SOURCE_DIR}/openmm/platforms/cpu)

    SET(OPENMM_PLATFORM_CPU 1)
    SET(OPENMM_PLATFORM_CUDA 0)
    SET(OPENMM_PLATFORM_OPENCL 0)

ELSEIF(OPENMM_PLATFORM MATCHES "CUDA" OR OPENMM_PLATFORM MATCHES "OPENCL")
    # building the common kernels is required for CUDA or OpenCL
    message(STATUS "Building OpenMM common platform")
    SET(OPENMM_PLATFORM_CPU 0)
    
    SET(COMMON_KERNEL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/src")
    SET(COMMON_KERNEL_SOURCE_CLASS CommonKernelSources)
    SET(COMMON_KERNELS_CPP ${COMMON_KERNEL_SOURCE_DIR}/${COMMON_KERNEL_SOURCE_CLASS}.cpp)
    SET(COMMON_KERNELS_H ${COMMON_KERNEL_SOURCE_DIR}/${COMMON_KERNEL_SOURCE_CLASS}.h)
    FILE(GLOB COMMON_KERNELS ${COMMON_KERNEL_SOURCE_DIR}/kernels/*.cc)
    ADD_CUSTOM_COMMAND(OUTPUT ${COMMON_KERNELS_CPP} ${COMMON_KERNELS_H}
        COMMAND ${CMAKE_COMMAND}
        ARGS -D KERNEL_SOURCE_DIR=${COMMON_KERNEL_SOURCE_DIR} -D KERNELS_CPP=${COMMON_KERNELS_CPP} -D KERNELS_H=${COMMON_KERNELS_H} -D KERNEL_SOURCE_CLASS=${COMMON_KERNEL_SOURCE_CLASS} -D KERNEL_FILE_EXTENSION=cc -P ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
        DEPENDS ${COMMON_KERNELS}
    )
    ADD_CUSTOM_TARGET(CommonKernels DEPENDS ${COMMON_KERNELS_CPP} ${COMMON_KERNELS_H})

    SET(OPENMM_DIRS ${OPENMM_DIRS} ${CMAKE_SOURCE_DIR}/openmm/platforms/common)
    SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/src)

    IF(OPENMM_PLATFORM MATCHES "CUDA")
        message(STATUS "Building OpenMM CUDA platform")
        SET(OPENMM_PLATFORM_CUDA 1)
        SET(OPENMM_PLATFORM_OPENCL 0)

        find_package(CUDAToolkit REQUIRED)
        SET(OPENMM_LIBRARIES CUDA::cuda_driver CUDA::cufft)

        # compile all cuda kernels into one single file
        SET(CUDA_KERNEL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src")
        SET(CUDA_KERNEL_SOURCE_CLASS CUDAKernelSources)
        SET(CUDA_KERNELS_CPP ${CUDA_KERNEL_SOURCE_DIR}/${CUDA_KERNEL_SOURCE_CLASS}.cpp)
        SET(CUDA_KERNELS_H ${CUDA_KERNEL_SOURCE_DIR}/${CUDA_KERNEL_SOURCE_CLASS}.h)
        FILE(GLOB CUDA_KERNELS ${CUDA_KERNEL_SOURCE_DIR}/kernels/*.cu)
        ADD_CUSTOM_COMMAND(OUTPUT ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H}
            COMMAND ${CMAKE_COMMAND}
            ARGS -D KERNEL_SOURCE_DIR=${CUDA_KERNEL_SOURCE_DIR} -D KERNELS_CPP=${CUDA_KERNELS_CPP} -D KERNELS_H=${CUDA_KERNELS_H} -D KERNEL_SOURCE_CLASS=${CUDA_KERNEL_SOURCE_CLASS} -D KERNEL_FILE_EXTENSION=cu -P ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
            DEPENDS ${CUDA_KERNELS}
        )
        ADD_CUSTOM_TARGET(CUDAKernels DEPENDS ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})

        SET(OPENMM_DIRS ${OPENMM_DIRS}
            ${CMAKE_SOURCE_DIR}/openmm/platforms/cuda
            # ${CMAKE_SOURCE_DIR}/openmm/plugins/cudacompiler
        )
        SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
            ${CUDAToolkit_INCLUDE_DIRS}
            ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src
            # ${CMAKE_CURRENT_SOURCE_DIR}/openmm/plugins/cudacompiler/src
        )
    ENDIF(OPENMM_PLATFORM MATCHES "CUDA")
        
    IF(OPENMM_PLATFORM MATCHES "OPENCL")
        message(STATUS "Building OpenMM OpenCL platform")
        SET(OPENMM_PLATFORM_CUDA 0)
        SET(OPENMM_PLATFORM_OPENCL 1)

        find_package(OpenCL REQUIRED)
        SET(OPENMM_LIBRARIES ${OpenCL_LIBRARY})

        # compile all opencl kernels into one single file
        SET(OPENCL_KERNEL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src")
        SET(OPENCL_KERNEL_SOURCE_CLASS OPENCLKernelSources)
        SET(OPENCL_KERNELS_CPP ${OPENCL_KERNEL_SOURCE_DIR}/${OPENCL_KERNEL_SOURCE_CLASS}.cpp)
        SET(OPENCL_KERNELS_H ${OPENCL_KERNEL_SOURCE_DIR}/${OPENCL_KERNEL_SOURCE_CLASS}.h)
        FILE(GLOB OPENCL_KERNELS ${OPENCL_KERNEL_SOURCE_DIR}/kernels/*.cl)
        ADD_CUSTOM_COMMAND(OUTPUT ${OPENCL_KERNELS_CPP} ${OPENCL_KERNELS_H}
            COMMAND ${CMAKE_COMMAND}
            ARGS -D KERNEL_SOURCE_DIR=${OPENCL_KERNEL_SOURCE_DIR} -D KERNELS_CPP=${OPENCL_KERNELS_CPP} -D KERNELS_H=${OPENCL_KERNELS_H} -D KERNEL_SOURCE_CLASS=${OPENCL_KERNEL_SOURCE_CLASS} -D KERNEL_FILE_EXTENSION=cl -P ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
            DEPENDS ${OPENCL_KERNELS}
        )
        ADD_CUSTOM_TARGET(OPENCLKernels DEPENDS ${OPENCL_KERNELS_CPP} ${OPENCL_KERNELS_H})

        SET(OPENMM_DIRS ${OPENMM_DIRS} ${CMAKE_SOURCE_DIR}/openmm/platforms/opencl)
        SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/src)
    ENDIF(OPENMM_PLATFORM MATCHES "OPENCL")

ELSE()
    message(ERROR "Unknown OPENMM_PLATFORM type " ${OPENMM_PLATFORM} ". Allowed types are CPU, CUDA or OPENCL")
ENDIF()

# get openmm include directories
FOREACH(subdir ${OPENMM_DIRS})
    SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
        ${subdir}
        ${subdir}/include
        ${subdir}/include/openmm
        ${subdir}/include/openmm/internal
        ${subdir}/include/openmm/common)
ENDFOREACH(subdir)

# find source and header files
set(OPENMM_SOURCE_C_FILES)
set(OPENMM_SOURCE_CXX_FILES)
set(OPENMM_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_DIRS})
    FILE(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
    SET(OPENMM_SOURCE_C_FILES ${OPENMM_SOURCE_C_FILES} ${src_c_files})

    FILE(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp ${subdir}/base/*.cpp ${subdir}/x86/*.cpp)
    SET(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${src_cxx_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(OPENMM_SOURCE_INCLUDE_FILES ${OPENMM_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

add_library(${OPENMM_LIB_NAME} ${LIBRARY_TYPE} ${OPENMM_SOURCE_C_FILES} ${OPENMM_SOURCE_CXX_FILES} ${OPENMM_SOURCE_INCLUDE_FILES})
target_include_directories(${OPENMM_LIB_NAME} PRIVATE ${OPENMM_INCLUDE_DIRS})
target_link_libraries(${OPENMM_LIB_NAME} PRIVATE ${OPENMM_LIBRARIES})
target_compile_definitions(${OPENMM_LIB_NAME}
    PRIVATE OPENMM_LIBRARY_NAME="OpenMM"
    PRIVATE OPENMM_MAJOR_VERSION=7
    PRIVATE OPENMM_MINOR_VERSION=5
    PRIVATE OPENMM_BUILD_VERSION=0
    PRIVATE HAVE_SSE2=1
    PRIVATE OPENMM_BUILD_STATIC_LIB=1
    PRIVATE IEEE_8087=1
    PRIVATE LEPTON_USE_JIT=1
    PRIVATE OPENMM_USE_STATIC_LIBRARIES=1
    PRIVATE LEPTON_USE_STATIC_LIBRARIES=1
    PRIVATE PTW32_STATIC_LIB=1
    PRIVATE OPENMM_BUILDING_STATIC_LIBRARY=1
    PRIVATE LEPTON_BUILDING_STATIC_LIBRARY=1
    PRIVATE OPENMM_COMMON_BUILDING_STATIC_LIBRARY=1
    # PRIVATE PTHREAD_BUILDING_STATIC_LIBRARY=1
)
set_target_properties(${OPENMM_LIB_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_precompile_headers(${OPENMM_LIB_NAME} PRIVATE pch/openmm_pch.h)
set_source_files_properties(${OPENMM_LIB_NAME} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)


# build simbody
set(SIMBODY_DIRS)
set(SIMBODY_INCLUDE_DIRS)

set(SIMTK_COMMON_DIRS . Scalar SmallMatrix Mechanics BigMatrix Geometry Simulation Random Polynomial)
foreach(subdir ${SIMTK_COMMON_DIRS})
    set(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir})
    set(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include/SimTKcommon
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKcommon/${subdir}/include/SimTKcommon/internal)
endforeach(subdir)

set(SIMTK_MATH_DIRS . LinearAlgebra Integrators Integrators/src/CPodes/sundials Optimizers Geometry)
foreach(subdir ${SIMTK_MATH_DIRS})
    set(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir})
    set(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include/simmath
        ${CMAKE_SOURCE_DIR}/Simbody01/SimTKmath/${subdir}/include/simmath/internal)
endforeach(subdir)

set(SIMBODY_DIRS ${SIMBODY_DIRS} ${CMAKE_SOURCE_DIR}/Simbody01/Simbody)
set(SIMBODY_INCLUDE_DIRS ${SIMBODY_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody
    ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/include/simbody/internal)

# find source and header files
set(SIMBODY_SOURCE_C_FILES)
set(SIMBODY_SOURCE_CXX_FILES)
set(SIMBODY_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${SIMBODY_DIRS})
    FILE(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
    SET(SIMBODY_SOURCE_C_FILES ${SIMBODY_SOURCE_C_FILES} ${src_c_files})

    FILE(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    SET(SIMBODY_SOURCE_CXX_FILES ${SIMBODY_SOURCE_CXX_FILES} ${src_cxx_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(SIMBODY_SOURCE_INCLUDE_FILES ${SIMBODY_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

add_library(${SIMBODY_LIB_NAME} ${LIBRARY_TYPE} ${SIMBODY_SOURCE_C_FILES} ${SIMBODY_SOURCE_CXX_FILES} ${SIMBODY_SOURCE_INCLUDE_FILES})
target_include_directories(${SIMBODY_LIB_NAME} PRIVATE ${SIMBODY_INCLUDE_DIRS})
target_compile_definitions(${SIMBODY_LIB_NAME}
    PRIVATE SimTK_SimTKCOMMON_LIBRARY_NAME="SimTKcommon"
    PRIVATE SimTK_SimTKCOMMON_MAJOR_VERSION=3
    PRIVATE SimTK_SimTKCOMMON_MINOR_VERSION=8
    PRIVATE SimTK_SimTKCOMMON_PATCH_VERSION=0
    PRIVATE SimTK_SIMMATH_LIBRARY_NAME="SimTKmath"
    PRIVATE SimTK_SIMMATH_MAJOR_VERSION=3
    PRIVATE SimTK_SIMMATH_MINOR_VERSION=8
    PRIVATE SimTK_SIMMATH_PATCH_VERSION=0
    PRIVATE SimTK_SIMBODY_LIBRARY_NAME="SimTKsimbody"
    PRIVATE SimTK_SIMBODY_MAJOR_VERSION=3
    PRIVATE SimTK_SIMBODY_MINOR_VERSION=8
    PRIVATE SimTK_SIMBODY_PATCH_VERSION=0
)
set_target_properties(${SIMBODY_LIB_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_precompile_headers(${SIMBODY_LIB_NAME} PRIVATE pch/simbody_pch.h)
set_source_files_properties(${SIMBODY_SOURCE_C_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

# build molmodel
set(MOLMODEL_DIRS ${CMAKE_SOURCE_DIR}/Molmodel)
set(MOLMODEL_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/Molmodel/include
    ${CMAKE_SOURCE_DIR}/Molmodel/include/molmodel
    ${CMAKE_SOURCE_DIR}/Molmodel/include/molmodel/internal)

set(MOLMODEL_SOURCE_FILES)
set(MOLMODEL_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${MOLMODEL_DIRS})
    FILE(GLOB src_files ${subdir}/src/*.cpp ${subdir}/src/*.c ${subdir}/src/*/*.cpp ${subdir}/src/*/*.c)
    SET(MOLMODEL_SOURCE_FILES ${MOLMODEL_SOURCE_FILES} ${src_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(MOLMODEL_SOURCE_INCLUDE_FILES ${MOLMODEL_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

add_library(${MOLMODEL_LIB_NAME} ${LIBRARY_TYPE} ${MOLMODEL_SOURCE_FILES} ${MOLMODEL_SOURCE_INCLUDE_FILES})
target_include_directories(${MOLMODEL_LIB_NAME} PRIVATE ${OPENMM_INCLUDE_DIRS} ${SIMBODY_INCLUDE_DIRS} ${MOLMODEL_INCLUDE_DIRS})
target_compile_definitions(${MOLMODEL_LIB_NAME}
    PRIVATE MOLMODEL_COPYRIGHT_YEARS="2006-12"
    PRIVATE MOLMODEL_AUTHORS="Christopher.Bruns_Michael.Sherman"
    PRIVATE SimTK_MOLMODEL_LIBRARY_NAME="SimTKmolmodel"
    PRIVATE SimTK_MOLMODEL_MAJOR_VERSION=3
    PRIVATE SimTK_MOLMODEL_MINOR_VERSION=0
    PRIVATE SimTK_MOLMODEL_PATCH_VERSION=0
    PRIVATE OPENMM_PLATFORM_CPU=${OPENMM_PLATFORM_CPU}
    PRIVATE OPENMM_PLATFORM_CUDA=${OPENMM_PLATFORM_CUDA}
    PRIVATE OPENMM_PLATFORM_OPENCL=${OPENMM_PLATFORM_OPENCL}
)
set_target_properties(${MOLMODEL_LIB_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_precompile_headers(${MOLMODEL_LIB_NAME} PRIVATE pch/molmodel_pch.h)
set_source_files_properties(Molmodel/src/vmdsock.c PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

# # build format
# FILE(GLOB FORMAT_SOURCE_FILES
#     ${CMAKE_SOURCE_DIR}/format/src/*.c
#     ${CMAKE_SOURCE_DIR}/format/src/*.cpp)
# SET(FORMAT_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/format/inc)
# add_library(format ${LIBRARY_TYPE} ${FORMAT_SOURCE_FILES})
# target_include_directories(format PRIVATE ${FORMAT_INCLUDE_DIRS})

# build robosample
FILE(GLOB ROBOSAMPLE_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)
set(ROBOSAMPLE_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

add_executable(${ROBOSAMPLE_NAME} ${ROBOSAMPLE_SOURCE_FILES})
target_include_directories(${ROBOSAMPLE_NAME} PRIVATE
    ${OPENMM_INCLUDE_DIRS}
    ${SIMBODY_INCLUDE_DIRS}
    ${MOLMODEL_INCLUDE_DIRS}
    # ${FORMAT_INCLUDE_DIRS}
    ${ROBOSAMPLE_INCLUDE_DIRS}
)

# the order of these libraries is extremely important
target_link_libraries(${ROBOSAMPLE_NAME} PRIVATE
    ${LIBRARIES}
    ${MOLMODEL_LIB_NAME}
    ${SIMBODY_LIB_NAME}
    ${OPENMM_LIB_NAME}
)
set_target_properties(${ROBOSAMPLE_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_precompile_headers(${ROBOSAMPLE_NAME} PRIVATE pch/robosample_pch.h)
target_link_options(${ROBOSAMPLE_NAME} PRIVATE ${BOLT_FLAGS})

# add example inputs
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp/pdbs)
file(GLOB INPUTS ${CMAKE_SOURCE_DIR}/tests_inputs/*)
file(COPY ${INPUTS} DESTINATION ${CMAKE_BINARY_DIR})
