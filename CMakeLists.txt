cmake_minimum_required(VERSION 3.10)

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
SET(EXTRA_COMPILE_FLAGS_DEBUG   "${COMMON_FLAGS} -O0 -g")
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
    SET(SANITIZER_FLAGS "-fsanitize=undefined,leak,address,float-divide-by-zero -fsanitize-address-use-after-scope -fsanitize-recover=address")

    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_DEBUG} -fno-omit-frame-pointer ${SANITIZER_FLAGS}")
    SET(PREFIX ".d")
    SET(BUILD_STATIC OFF)
    SET(BUILD_DYNAMIC ON)

ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES "Release_Tests" OR CMAKE_BUILD_TYPE MATCHES "Release_Full")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(BOLT_FLAGS "-Wl,--emit-relocs")

    SET(BUILD_STATIC ON)
    SET(BUILD_DYNAMIC OFF)

    IF (CMAKE_BUILD_TYPE MATCHES "Release_Tests")
        SET(BUILD_STATIC OFF)
        SET(BUILD_DYNAMIC ON)
    ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release_Full")
        SET(BUILD_STATIC ON)
        SET(BUILD_DYNAMIC ON)
    ENDIF()

ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Train")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(PREFIX ".pgo.train")
    SET(BUILD_STATIC ON)
    SET(BUILD_DYNAMIC OFF)

ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Use")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(PREFIX ".pgo.use")
    SET(BOLT_FLAGS "-Wl,--emit-relocs")
    SET(BUILD_STATIC ON)
    SET(BUILD_DYNAMIC OFF)
ENDIF()

separate_arguments(SANITIZER_FLAGS UNIX_COMMAND "${SANITIZER_FLAGS}")
separate_arguments(PGO_FLAGS UNIX_COMMAND "${PGO_FLAGS}")
separate_arguments(BUILD_FLAGS UNIX_COMMAND "${BUILD_FLAGS}")
separate_arguments(BOLT_FLAGS UNIX_COMMAND "${BOLT_FLAGS}")

# set names of the targets we will build
SET(OPENMM_LIB_SHARED_NAME "openmm.shared${PREFIX}")
SET(OPENMM_LIB_STATIC_NAME "openmm.static${PREFIX}")

SET(SIMBODY_LIB_SHARED_NAME "simbody.shared${PREFIX}")
SET(SIMBODY_LIB_STATIC_NAME "simbody.static${PREFIX}")

SET(MOLMODEL_LIB_SHARED_NAME "molmodel.shared${PREFIX}")
SET(MOLMODEL_LIB_STATIC_NAME "molmodel.static${PREFIX}")

SET(ROBOSAMPLE_LIB_SHARED_NAME "robosample.shared${PREFIX}")
SET(ROBOSAMPLE_LIB_STATIC_NAME "robosample.static${PREFIX}")

SET(ROBOSAMPLE_NAME "robosample${PREFIX}")

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

# set the precompiled headers
set(OPENMM_PCH pch/openmm_pch.h)
set(SIMBODY_PCH pch/simbody_pch.h)
set(MOLMODEL_PCH pch/molmodel_pch.h)
set(ROBOSAMPLE_PCH pch/robosample_pch.h)

# get openmm source directories
SET(OPENMM_INCLUDE_DIRS)
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
        COMMENT "Generating common kernel sources for OpenMM..."
    )

    # this command is executed when building, not when running CMakeLists.txt
    ADD_CUSTOM_TARGET(CommonKernels DEPENDS ${COMMON_KERNELS_CPP} ${COMMON_KERNELS_H})
    SET(OPENMM_DEPENDENCIES CommonKernels)

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
            COMMENT "Generating CUDA kernel sources for OpenMM..."
        )

        # this command is executed when building, not when running CMakeLists.txt
        ADD_CUSTOM_TARGET(CUDAKernels DEPENDS ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})
        SET(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} CUDAKernels)

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
        SET(OPENCL_KERNEL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/src")
        SET(OPENCL_KERNEL_SOURCE_CLASS OpenCLKernelSources)
        SET(OPENCL_KERNELS_CPP ${OPENCL_KERNEL_SOURCE_DIR}/${OPENCL_KERNEL_SOURCE_CLASS}.cpp)
        SET(OPENCL_KERNELS_H ${OPENCL_KERNEL_SOURCE_DIR}/${OPENCL_KERNEL_SOURCE_CLASS}.h)
        FILE(GLOB OPENCL_KERNELS ${OPENCL_KERNEL_SOURCE_DIR}/kernels/*.cl)
        ADD_CUSTOM_COMMAND(OUTPUT ${OPENCL_KERNELS_CPP} ${OPENCL_KERNELS_H}
            COMMAND ${CMAKE_COMMAND}
            ARGS -D KERNEL_SOURCE_DIR=${OPENCL_KERNEL_SOURCE_DIR} -D KERNELS_CPP=${OPENCL_KERNELS_CPP} -D KERNELS_H=${OPENCL_KERNELS_H} -D KERNEL_SOURCE_CLASS=${OPENCL_KERNEL_SOURCE_CLASS} -D KERNEL_FILE_EXTENSION=cl -P ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
            DEPENDS ${OPENCL_KERNELS}
            COMMENT "Generating OpenCL kernel sources for OpenMM..."
        )
        
        # this command is executed when building, not when running CMakeLists.txt
        ADD_CUSTOM_TARGET(OpenCLKernels DEPENDS ${OPENCL_KERNELS_CPP} ${OPENCL_KERNELS_H})
        SET(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} OpenCLKernels)

        SET(OPENMM_DIRS ${OPENMM_DIRS} ${CMAKE_SOURCE_DIR}/openmm/platforms/opencl)
        SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/src)
    ENDIF(OPENMM_PLATFORM MATCHES "OPENCL")

ELSE()
    message(ERROR "Unknown OPENMM_PLATFORM type " ${OPENMM_PLATFORM} ". Allowed types are CPU, CUDA or OPENCL")
ENDIF()

# set compile definitions for each library
set(OPENMM_COMPILE_DEFINITIONS
    OPENMM_LIBRARY_NAME="OpenMM"
    OPENMM_MAJOR_VERSION=7
    OPENMM_MINOR_VERSION=5
    OPENMM_BUILD_VERSION=0
    HAVE_SSE2=1
    OPENMM_BUILD_STATIC_LIB=1
    IEEE_8087=1
    LEPTON_USE_JIT=1
    OPENMM_USE_STATIC_LIBRARIES=1
    LEPTON_USE_STATIC_LIBRARIES=1
    PTW32_STATIC_LIB=1
    OPENMM_BUILDING_STATIC_LIBRARY=1
    LEPTON_BUILDING_STATIC_LIBRARY=1
    OPENMM_COMMON_BUILDING_STATIC_LIBRARY=1
    # PTHREAD_BUILDING_STATIC_LIBRARY=1
)

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

set(MOLMODEL_COMPILE_DEFINITIONS
    MOLMODEL_COPYRIGHT_YEARS="2006-12"
    MOLMODEL_AUTHORS="Christopher.Bruns_Michael.Sherman"
    SimTK_MOLMODEL_LIBRARY_NAME="SimTKmolmodel"
    SimTK_MOLMODEL_MAJOR_VERSION=3
    SimTK_MOLMODEL_MINOR_VERSION=0
    SimTK_MOLMODEL_PATCH_VERSION=0
    OPENMM_PLATFORM_CPU=${OPENMM_PLATFORM_CPU}
    OPENMM_PLATFORM_CUDA=${OPENMM_PLATFORM_CUDA}
    OPENMM_PLATFORM_OPENCL=${OPENMM_PLATFORM_OPENCL}
)

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
SET(OPENMM_SOURCE_C_FILES)
SET(OPENMM_SOURCE_CXX_FILES)
SET(OPENMM_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_DIRS})
    FILE(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
    SET(OPENMM_SOURCE_C_FILES ${OPENMM_SOURCE_C_FILES} ${src_c_files})

    FILE(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp ${subdir}/base/*.cpp ${subdir}/x86/*.cpp)
    SET(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${src_cxx_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(OPENMM_SOURCE_INCLUDE_FILES ${OPENMM_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

set_source_files_properties(${OPENMM_LIB_SHARED_NAME} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

IF (BUILD_DYNAMIC)
    add_library(${OPENMM_LIB_SHARED_NAME} SHARED ${OPENMM_SOURCE_C_FILES} ${OPENMM_SOURCE_CXX_FILES} ${OPENMM_SOURCE_INCLUDE_FILES})
    target_compile_options(${OPENMM_LIB_SHARED_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${OPENMM_LIB_SHARED_NAME} PRIVATE ${OPENMM_INCLUDE_DIRS})
    target_link_libraries(${OPENMM_LIB_SHARED_NAME} PRIVATE ${OPENMM_LIBRARIES})
    target_compile_definitions(${OPENMM_LIB_SHARED_NAME} PRIVATE ${OPENMM_COMPILE_DEFINITIONS})
    set_target_properties(${OPENMM_LIB_SHARED_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${OPENMM_LIB_SHARED_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    target_precompile_headers(${OPENMM_LIB_SHARED_NAME} PRIVATE ${OPENMM_PCH})
    set_property(TARGET ${OPENMM_LIB_SHARED_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF (BUILD_DYNAMIC)

IF (BUILD_STATIC)
    add_library(${OPENMM_LIB_STATIC_NAME} STATIC ${OPENMM_SOURCE_C_FILES} ${OPENMM_SOURCE_CXX_FILES} ${OPENMM_SOURCE_INCLUDE_FILES})
    target_compile_options(${OPENMM_LIB_STATIC_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${OPENMM_LIB_STATIC_NAME} PRIVATE ${OPENMM_INCLUDE_DIRS})
    target_link_libraries(${OPENMM_LIB_STATIC_NAME} PRIVATE ${OPENMM_LIBRARIES})
    target_compile_definitions(${OPENMM_LIB_STATIC_NAME} PRIVATE ${OPENMM_COMPILE_DEFINITIONS})
    set_target_properties(${OPENMM_LIB_STATIC_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${OPENMM_LIB_STATIC_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    target_precompile_headers(${OPENMM_LIB_STATIC_NAME} PRIVATE ${OPENMM_PCH})
ENDIF (BUILD_STATIC)

# this forces cmake to encode the kernels before compiling anything else
add_dependencies(${OPENMM_LIB_NAME} ${OPENMM_DEPENDENCIES})


# build simbody
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

set_source_files_properties(${SIMBODY_SOURCE_C_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

IF (BUILD_DYNAMIC)
    add_library(${SIMBODY_LIB_SHARED_NAME} SHARED ${SIMBODY_SOURCE_C_FILES} ${SIMBODY_SOURCE_CXX_FILES} ${SIMBODY_SOURCE_INCLUDE_FILES})
    target_compile_options(${SIMBODY_LIB_SHARED_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${SIMBODY_LIB_SHARED_NAME} PRIVATE ${SIMBODY_INCLUDE_DIRS})
    target_compile_definitions(${SIMBODY_LIB_SHARED_NAME} PRIVATE ${SIMBODY_COMPILE_DEFINITIONS})
    set_target_properties(${SIMBODY_LIB_SHARED_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${SIMBODY_LIB_SHARED_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    target_precompile_headers(${SIMBODY_LIB_SHARED_NAME} PRIVATE ${SIMBODY_PCH})
    set_property(TARGET ${SIMBODY_LIB_SHARED_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF (BUILD_DYNAMIC)

IF (BUILD_STATIC)
    add_library(${SIMBODY_LIB_STATIC_NAME} STATIC ${SIMBODY_SOURCE_C_FILES} ${SIMBODY_SOURCE_CXX_FILES} ${SIMBODY_SOURCE_INCLUDE_FILES})
    target_compile_options(${SIMBODY_LIB_STATIC_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${SIMBODY_LIB_STATIC_NAME} PRIVATE ${SIMBODY_INCLUDE_DIRS})
    target_compile_definitions(${SIMBODY_LIB_STATIC_NAME} PRIVATE ${SIMBODY_COMPILE_DEFINITIONS})
    set_target_properties(${SIMBODY_LIB_STATIC_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${SIMBODY_LIB_STATIC_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    target_precompile_headers(${SIMBODY_LIB_STATIC_NAME} PRIVATE ${SIMBODY_PCH})
ENDIF (BUILD_STATIC)

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

    IF (BUILD_DYNAMIC)
        target_link_libraries(${GUI_NAME} ${SIMBODY_LIB_SHARED_NAME} ${LIBRARIES} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})
    ELSE()
        target_link_libraries(${GUI_NAME} ${SIMBODY_LIB_STATIC_NAME} ${LIBRARIES} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})
    ENDIF()
    
    set_target_properties(${GUI_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${GUI_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
ENDIF(BUILD_VISUALIZER)

# build molmodel
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

set_source_files_properties(${MOLMODEL_C_SOURCE_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

IF (BUILD_DYNAMIC)
    add_library(${MOLMODEL_LIB_SHARED_NAME} SHARED ${MOLMODEL_C_SOURCE_FILES} ${MOLMODEL_CXX_SOURCE_FILES} ${MOLMODEL_SOURCE_INCLUDE_FILES})
    target_compile_options(${MOLMODEL_LIB_SHARED_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${MOLMODEL_LIB_SHARED_NAME} PRIVATE ${OPENMM_INCLUDE_DIRS} ${SIMBODY_INCLUDE_DIRS} ${MOLMODEL_INCLUDE_DIRS})
    target_compile_definitions(${MOLMODEL_LIB_SHARED_NAME} PRIVATE ${MOLMODEL_COMPILE_DEFINITIONS})
    set_target_properties(${MOLMODEL_LIB_SHARED_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${MOLMODEL_LIB_SHARED_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    target_precompile_headers(${MOLMODEL_LIB_SHARED_NAME} PRIVATE ${MOLMODEL_PCH})
    set_property(TARGET ${MOLMODEL_LIB_SHARED_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF (BUILD_DYNAMIC)

IF (BUILD_STATIC)
    add_library(${MOLMODEL_LIB_STATIC_NAME} STATIC ${MOLMODEL_C_SOURCE_FILES} ${MOLMODEL_CXX_SOURCE_FILES} ${MOLMODEL_SOURCE_INCLUDE_FILES})
    target_compile_options(${MOLMODEL_LIB_STATIC_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${MOLMODEL_LIB_STATIC_NAME} PRIVATE ${OPENMM_INCLUDE_DIRS} ${SIMBODY_INCLUDE_DIRS} ${MOLMODEL_INCLUDE_DIRS})
    target_compile_definitions(${MOLMODEL_LIB_STATIC_NAME} PRIVATE ${MOLMODEL_COMPILE_DEFINITIONS})
    set_target_properties(${MOLMODEL_LIB_STATIC_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_link_options(${MOLMODEL_LIB_STATIC_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS})
    target_precompile_headers(${MOLMODEL_LIB_STATIC_NAME} PRIVATE ${MOLMODEL_PCH})
ENDIF (BUILD_STATIC)


# build the robosample library
FILE(GLOB src_c_files ${CMAKE_SOURCE_DIR}/format/src/*.c)
SET(ROBOSAMPLE_C_SOURCE_FILES ${src_c_files})
set_source_files_properties(${ROBOSAMPLE_C_SOURCE_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

FILE(GLOB src_cxx_files ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/format/src/*.cpp)
SET(ROBOSAMPLE_CXX_SOURCE_FILES ${src_cxx_files})
LIST(REMOVE_ITEM ROBOSAMPLE_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)

SET(ROBOSAMPLE_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/format/inc
    ${CMAKE_SOURCE_DIR}/format/inc/mdtraj
    ${OPENMM_INCLUDE_DIRS}
    ${SIMBODY_INCLUDE_DIRS}
    ${MOLMODEL_INCLUDE_DIRS}
)

# the order of these libraries is extremely important
SET(ROBOSAMPLE_LIBRARIES_SHARED
    ${LIBRARIES}
    ${MOLMODEL_LIB_SHARED_NAME}
    ${SIMBODY_LIB_SHARED_NAME}
    ${OPENMM_LIB_SHARED_NAME}
)

SET(ROBOSAMPLE_LIBRARIES_STATIC
    ${LIBRARIES}
    ${MOLMODEL_LIB_STATIC_NAME}
    ${SIMBODY_LIB_STATIC_NAME}
    ${OPENMM_LIB_STATIC_NAME}
)

IF (BUILD_DYNAMIC)
    add_library(${ROBOSAMPLE_LIB_SHARED_NAME} SHARED ${ROBOSAMPLE_C_SOURCE_FILES} ${ROBOSAMPLE_CXX_SOURCE_FILES})
    target_compile_options(${ROBOSAMPLE_LIB_SHARED_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${ROBOSAMPLE_LIB_SHARED_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
    target_link_libraries(${ROBOSAMPLE_LIB_SHARED_NAME} PRIVATE ${ROBOSAMPLE_LIBRARIES_SHARED})
    set_target_properties(${ROBOSAMPLE_LIB_SHARED_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_precompile_headers(${ROBOSAMPLE_LIB_SHARED_NAME} PRIVATE ${ROBOSAMPLE_PCH})
    target_link_options(${ROBOSAMPLE_LIB_SHARED_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
    set_property(TARGET ${ROBOSAMPLE_LIB_SHARED_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF (BUILD_DYNAMIC)

IF (BUILD_STATIC)
    add_library(${ROBOSAMPLE_LIB_STATIC_NAME} STATIC ${ROBOSAMPLE_C_SOURCE_FILES} ${ROBOSAMPLE_CXX_SOURCE_FILES})
    target_compile_options(${ROBOSAMPLE_LIB_STATIC_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS})
    target_include_directories(${ROBOSAMPLE_LIB_STATIC_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
    target_link_libraries(${ROBOSAMPLE_LIB_STATIC_NAME} PRIVATE ${ROBOSAMPLE_LIBRARIES_STATIC})
    set_target_properties(${ROBOSAMPLE_LIB_STATIC_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
    target_precompile_headers(${ROBOSAMPLE_LIB_STATIC_NAME} PRIVATE ${ROBOSAMPLE_PCH})
    target_link_options(${ROBOSAMPLE_LIB_STATIC_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
ENDIF (BUILD_STATIC)

# add the actual robosample executable
add_executable(${ROBOSAMPLE_NAME} ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)
target_compile_options(${ROBOSAMPLE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})
target_include_directories(${ROBOSAMPLE_NAME} PRIVATE ${ROBOSAMPLE_INCLUDE_DIRS})
set_target_properties(${ROBOSAMPLE_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_link_options(${ROBOSAMPLE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS})

IF(BUILD_DYNAMIC)
    target_link_libraries(${ROBOSAMPLE_NAME} PRIVATE ${ROBOSAMPLE_LIB_SHARED_NAME} ${MOLMODEL_LIB_SHARED_NAME} ${SIMBODY_LIB_SHARED_NAME})
ELSE()
    target_link_libraries(${ROBOSAMPLE_NAME} PRIVATE ${ROBOSAMPLE_LIB_STATIC_NAME} ${MOLMODEL_LIB_STATIC_NAME} ${SIMBODY_LIB_STATIC_NAME})
ENDIF(BUILD_DYNAMIC)


# add tests
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
        target_link_libraries(${test_name} PRIVATE ${ROBOSAMPLE_LIB_SHARED_NAME} ${MOLMODEL_LIB_SHARED_NAME} ${SIMBODY_LIB_SHARED_NAME} ${LIBRARIES})

        add_test(NAME ${test_name} COMMAND ${test_name})
    ENDFOREACH()
ENDIF ()


# add example inputs
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp/pdbs)
file(GLOB INPUTS ${CMAKE_SOURCE_DIR}/tests_inputs/*)
file(COPY ${INPUTS} DESTINATION ${CMAKE_BINARY_DIR})
