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
# SET(OPENMM_PLATFORM "CPU" CACHE STRING "OpenMM platform: CPU, CUDA or OpenCL.")
SET(OPENMM_PLATFORM "CUDA")
# SET(OPENMM_PLATFORM "OPENCL")

IF(OPENMM_PLATFORM MATCHES "CPU")
    # build the default CPU platform
    message(STATUS "Building OpenMM CPU platform")
    SET(OPENMM_DIRS ${OPENMM_DIRS} ${CMAKE_SOURCE_DIR}/openmm/platforms/cpu)

    SET(USE_OPENMM_PLATFORM_CPU 1)
    SET(USE_OPENMM_PLATFORM_CUDA 0)
    SET(USE_OPENMM_PLATFORM_OPENCL 0)

ELSEIF(OPENMM_PLATFORM MATCHES "CUDA" OR OPENMM_PLATFORM MATCHES "OPENCL")
    # building the common kernels is required for CUDA or OpenCL
    message(STATUS "Building OpenMM common platform")
    SET(USE_OPENMM_PLATFORM_CPU 0)
    
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
        SET(USE_OPENMM_PLATFORM_CUDA 1)
        SET(USE_OPENMM_PLATFORM_OPENCL 0)

        find_package(CUDAToolkit REQUIRED)
        SET(OPENMM_LIBRARIES CUDA::cuda_driver CUDA::cufft)

        # compile all cuda kernels into one single file
        SET(CUDA_KERNEL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src")
        SET(CUDA_KERNEL_SOURCE_CLASS CudaKernelSources)
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
        ADD_CUSTOM_TARGET(CudaKernels DEPENDS ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})
        SET(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} CudaKernels)

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
        SET(USE_OPENMM_PLATFORM_CUDA 0)
        SET(USE_OPENMM_PLATFORM_OPENCL 1)

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

# set_source_files_properties(${OPENMM_LIB_SHARED_NAME} PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

# set compile definitions for each library
set(OPENMM_COMPILE_DEFINITIONS
    OPENMM_LIBRARY_NAME="OpenMM"
    OPENMM_MAJOR_VERSION=7
    OPENMM_MINOR_VERSION=5
    OPENMM_BUILD_VERSION=0
    HAVE_SSE2=1
    OPENMM_BUILD_STATIC_LIB=0
    IEEE_8087=1
    LEPTON_USE_JIT=1
    OPENMM_USE_STATIC_LIBRARIES=0
    LEPTON_USE_STATIC_LIBRARIES=0
    PTW32_STATIC_LIB=0
    OPENMM_BUILDING_STATIC_LIBRARY=0
    LEPTON_BUILDING_STATIC_LIBRARY=0
    OPENMM_COMMON_BUILDING_STATIC_LIBRARY=0
    # PLATFORM_CPU_NOSNAODIWEAD=1
    OPENMM_PLATFORM_CPU=${USE_OPENMM_PLATFORM_CPU}
    OPENMM_PLATFORM_CUDA=${USE_OPENMM_PLATFORM_CUDA}
    OPENMM_PLATFORM_OPENCL=${USE_OPENMM_PLATFORM_OPENCL}
)

foreach(compile_definition ${OPENMM_COMPILE_DEFINITIONS})
    set_property(SOURCE ${OPENMM_SOURCE_CXX_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
    set_property(SOURCE ${OPENMM_SOURCE_INCLUDE_FILES} APPEND_STRING PROPERTY COMPILE_DEFINITIONS ${compile_definition})
endforeach(compile_definition)
