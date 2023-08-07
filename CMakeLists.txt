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

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    # set the sanitizer flags
    # https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html
    # https://clang.llvm.org/docs/AddressSanitizer.html
    SET(SANITIZER_FLAGS "-fsanitize=undefined,leak,address,float-divide-by-zero -fsanitize-address-use-after-scope -fsanitize-recover=address")

    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_DEBUG} -fno-omit-frame-pointer ${SANITIZER_FLAGS}")
    SET(PREFIX ".d")
ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(BOLT_FLAGS "-Wl,--emit-relocs")
ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Train")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(PREFIX ".pgo.train")
ELSEIF(CMAKE_BUILD_TYPE MATCHES "PGO_Use")
    SET(BUILD_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILE_FLAGS_RELEASE}")
    SET(PREFIX ".pgo.use")
    SET(BOLT_FLAGS "-Wl,--emit-relocs")
ENDIF()

separate_arguments(SANITIZER_FLAGS UNIX_COMMAND "${SANITIZER_FLAGS}")
separate_arguments(PGO_FLAGS UNIX_COMMAND "${PGO_FLAGS}")
separate_arguments(BUILD_FLAGS UNIX_COMMAND "${BUILD_FLAGS}")
separate_arguments(BOLT_FLAGS UNIX_COMMAND "${BOLT_FLAGS}")

# set names of the targets we will build
SET(OPENMM_LIB_NAME "openmm${PREFIX}")
SET(SIMBODY_LIB_NAME "simbody${PREFIX}")
SET(MOLMODEL_LIB_NAME "molmodel${PREFIX}")
SET(FORMATS_LIB_NAME "formats${PREFIX}")
SET(ROBOSAMPLE_LIB_NAME "robosample.lib${PREFIX}")
SET(ROBOSAMPLE_NAME "robosample${PREFIX}")

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    SET(LIBRARY_TYPE SHARED)
ELSE()
    SET(LIBRARY_TYPE STATIC)
ENDIF()

# we require ipo (interprocedural optimization) aka lto (link time optimization)
# this is fancy -flto basically
# on clang, it sets -flto=thin and this seems to work better than the default -flto(=full)
# even though the binary is larger, it runs faster
include(CheckIPOSupported)
check_ipo_supported()

IF(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    SET(USE_IPO TRUE)
ELSE()
    SET(USE_IPO FALSE)
ENDIF()

# SET(CMAKE_CXX_STANDARD 14)

set(BUILD_VISUALIZER ON)

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

add_library(${OPENMM_LIB_NAME} ${LIBRARY_TYPE} ${OPENMM_SOURCE_C_FILES} ${OPENMM_SOURCE_CXX_FILES} ${OPENMM_SOURCE_INCLUDE_FILES})
target_compile_options(${OPENMM_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
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
target_link_options(${OPENMM_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
target_precompile_headers(${OPENMM_LIB_NAME} PRIVATE ${OPENMM_PCH})
set_source_files_properties(${OPENMM_LIB_NAME} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

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

add_library(${SIMBODY_LIB_NAME} ${LIBRARY_TYPE} ${SIMBODY_SOURCE_C_FILES} ${SIMBODY_SOURCE_CXX_FILES} ${SIMBODY_SOURCE_INCLUDE_FILES})
target_compile_options(${SIMBODY_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
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
    PRIVATE SIMBODY_VISUALIZER_REL_INSTALL_DIR="${CMAKE_BINARY_DIR}"
    PRIVATE SIMBODY_PATH_FROM_LIBDIR_TO_VIZ_DIR="${CMAKE_BINARY_DIR}"
    PRIVATE SIMBODY_VISUALIZER_INSTALL_DIR="${CMAKE_BINARY_DIR}"
    PRIVATE SIMBODY_VISUALIZER_REL_INSTALL_DIR="${CMAKE_BINARY_DIR}"
)
set_target_properties(${SIMBODY_LIB_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_link_options(${SIMBODY_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
target_precompile_headers(${SIMBODY_LIB_NAME} PRIVATE ${SIMBODY_PCH})
set_source_files_properties(${SIMBODY_SOURCE_C_FILES} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

IF(BUILD_VISUALIZER)
    IF(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
        SET(GUI_NAME "simbody-visualizer")
    ELSE()
        SET(GUI_NAME "simbody-visualizer_d")
    ENDIF()
    
    SET(VISUALIZER_DIR ${CMAKE_SOURCE_DIR}/Simbody01/Simbody/Visualizer/simbody-visualizer)
    FILE(GLOB VISUALIZER_SOURCE_INCLUDE_FILES ${VISUALIZER_DIR}/*.h)
    FILE(GLOB VISUALIZER_SOURCE_CXX_FILES ${VISUALIZER_DIR}/*.cpp)

    add_executable(${GUI_NAME} ${VISUALIZER_SOURCE_INCLUDE_FILES} ${VISUALIZER_SOURCE_CXX_FILES})
    target_include_directories(${GUI_NAME} PRIVATE ${GLUT_INCLUDE_DIR} ${SIMBODY_INCLUDE_DIRS})
    target_link_libraries(${GUI_NAME} ${SIMBODY_LIB_NAME} ${LIBRARIES} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})
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
    FILE(GLOB src_files ${subdir}/src/*.cpp ${subdir}/src/*.c ${subdir}/src/*/*.cpp ${subdir}/src/*/*.c)
    SET(MOLMODEL_SOURCE_FILES ${MOLMODEL_SOURCE_FILES} ${src_files})

    # pimpl pattern is used and headers are stored in the src directory
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(MOLMODEL_SOURCE_INCLUDE_FILES ${MOLMODEL_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

add_library(${MOLMODEL_LIB_NAME} ${LIBRARY_TYPE} ${MOLMODEL_SOURCE_FILES} ${MOLMODEL_SOURCE_INCLUDE_FILES})
target_compile_options(${MOLMODEL_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
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
target_link_options(${MOLMODEL_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
target_precompile_headers(${MOLMODEL_LIB_NAME} PRIVATE ${MOLMODEL_PCH})
set_source_files_properties(Molmodel/src/vmdsock.c PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

# # build format - this does not have ipo or pch yet
# FILE(GLOB FORMAT_SOURCE_FILES
#     ${CMAKE_SOURCE_DIR}/format/src/*.c
#     ${CMAKE_SOURCE_DIR}/format/src/*.cpp)
# SET(FORMAT_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/format/inc)
# add_library(${FORMATS_LIB_NAME} ${LIBRARY_TYPE} ${FORMAT_SOURCE_FILES})
# target_compile_options(${FORMATS_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}
# target_include_directories(${FORMATS_LIB_NAME} PRIVATE ${FORMAT_INCLUDE_DIRS})
# target_link_options(${FORMATS_LIB_NAME} PRIVATE ${BOLT_FLAGS} ${SANITIZER_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}

# build robosample
FILE(GLOB ROBOSAMPLE_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM ROBOSAMPLE_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)
SET(ROBOSAMPLE_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

add_library(${ROBOSAMPLE_LIB_NAME} ${LIBRARY_TYPE} ${ROBOSAMPLE_SOURCE_FILES})
target_compile_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${BUILD_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}

target_include_directories(${ROBOSAMPLE_LIB_NAME} PRIVATE
    ${OPENMM_INCLUDE_DIRS}
    ${SIMBODY_INCLUDE_DIRS}
    ${MOLMODEL_INCLUDE_DIRS}
    # ${FORMAT_INCLUDE_DIRS}
    ${ROBOSAMPLE_INCLUDE_DIRS}
)

# the order of these libraries is extremely important
target_link_libraries(${ROBOSAMPLE_LIB_NAME} PRIVATE
    ${LIBRARIES}
    # ${FORMATS_LIB_NAME}
    ${MOLMODEL_LIB_NAME}
    ${SIMBODY_LIB_NAME}
    ${OPENMM_LIB_NAME}
)
set_target_properties(${ROBOSAMPLE_LIB_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_precompile_headers(${ROBOSAMPLE_LIB_NAME} PRIVATE ${ROBOSAMPLE_PCH})
target_link_options(${ROBOSAMPLE_LIB_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}

# add the actual robosample executable
add_executable(${ROBOSAMPLE_NAME} ${CMAKE_SOURCE_DIR}/src/Robosample.cpp)
target_compile_options(${ROBOSAMPLE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}

target_include_directories(${ROBOSAMPLE_NAME} PRIVATE
    ${OPENMM_INCLUDE_DIRS}
    ${SIMBODY_INCLUDE_DIRS}
    ${MOLMODEL_INCLUDE_DIRS}
    # ${FORMAT_INCLUDE_DIRS}
    ${ROBOSAMPLE_INCLUDE_DIRS}
)

# the order of these libraries is extremely important
target_link_libraries(${ROBOSAMPLE_NAME} PRIVATE
    ${ROBOSAMPLE_LIB_NAME}
    # ${FORMATS_LIB_NAME}
    ${MOLMODEL_LIB_NAME}
    ${SIMBODY_LIB_NAME}
)
set_target_properties(${ROBOSAMPLE_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${USE_IPO})
target_link_options(${ROBOSAMPLE_NAME} PRIVATE ${SANITIZER_FLAGS} ${BOLT_FLAGS} ${PGO_FLAGS}) # ${PGO_FLAGS}

# # add tests
# enable_testing()
# include(CTest)
# file(GLOB TESTS_SOURCES tests/*.cpp)
# foreach(test_src ${TESTS_SOURCES})
#     # get test name by removing file location and extension
#     file(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR}/tests ${test_src})
#     STRING(REGEX REPLACE ".cpp" "" test_name ${test_name})

#     # needs refactoring
#     add_executable(${test_name} ${test_src})
#     target_compile_options(${test_name} PRIVATE ${BUILD_FLAGS})
#     target_link_options(${test_name} PRIVATE ${SANITIZER_FLAGS})

#     target_include_directories(${test_name} PRIVATE
#         ${OPENMM_INCLUDE_DIRS}
#         ${SIMBODY_INCLUDE_DIRS}
#         ${MOLMODEL_INCLUDE_DIRS}
#         # ${FORMAT_INCLUDE_DIRS}
#         ${ROBOSAMPLE_INCLUDE_DIRS}
#     )
#     target_link_libraries(${test_name} PRIVATE
#         ${ROBOSAMPLE_LIB_NAME}
#         # ${FORMATS_LIB_NAME}
#         ${MOLMODEL_LIB_NAME}
#         ${SIMBODY_LIB_NAME}
#     )

#     add_test(NAME ${test_name} COMMAND ${test_name})

# endforeach(test_src )


# add example inputs
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp/pdbs)
file(GLOB INPUTS ${CMAKE_SOURCE_DIR}/tests_inputs/*)
file(COPY ${INPUTS} DESTINATION ${CMAKE_BINARY_DIR})
