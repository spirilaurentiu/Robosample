SET(OPENMM_GENERATED_CXX_FILES)

# Get openmm source directories
SET(OPENMM_INCLUDE_DIRS)
SET(OPENMM_SOURCE_SUBDIRS
    ${CMAKE_SOURCE_DIR}/openmm/
    ${CMAKE_SOURCE_DIR}/openmm/openmmapi
    ${CMAKE_SOURCE_DIR}/openmm/olla

    ${CMAKE_SOURCE_DIR}/openmm/libraries/jama
    ${CMAKE_SOURCE_DIR}/openmm/libraries/quern
    ${CMAKE_SOURCE_DIR}/openmm/libraries/lepton
    ${CMAKE_SOURCE_DIR}/openmm/libraries/sfmt
    ${CMAKE_SOURCE_DIR}/openmm/libraries/lbfgs
    ${CMAKE_SOURCE_DIR}/openmm/libraries/hilbert
    ${CMAKE_SOURCE_DIR}/openmm/libraries/csha1
    ${CMAKE_SOURCE_DIR}/openmm/libraries/pocketfft
    ${CMAKE_SOURCE_DIR}/openmm/libraries/vkfft
    ${CMAKE_SOURCE_DIR}/openmm/libraries/irrxml
    ${CMAKE_SOURCE_DIR}/openmm/libraries/vecmath

    # Needed by custom kernels and forces
    ${CMAKE_SOURCE_DIR}/openmm/platforms/reference
    
    ${CMAKE_SOURCE_DIR}/openmm/serialization
)

IF(USE_CPU)
    SET(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
        ${CMAKE_SOURCE_DIR}/openmm/platforms/cpu
    )
    SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cpu/include
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cpu/src
    )

ELSEIF(USE_CUDA OR USE_OPENCL)
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

    # This command is executed when building, not when running CMakeLists.txt
    ADD_CUSTOM_TARGET(CommonKernels DEPENDS ${COMMON_KERNELS_CPP} ${COMMON_KERNELS_H})
    SET(OPENMM_GENERATED_CXX_FILES ${OPENMM_GENERATED_CXX_FILES} ${COMMON_KERNELS_CPP})
    SET(OPENMM_DEPENDENCIES CommonKernels)

    SET(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
        ${CMAKE_SOURCE_DIR}/openmm/platforms/common
    )
    SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/include
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/src
    )

    IF(USE_CUDA)
        # Compile all CUDA kernels into one single file
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

        ADD_CUSTOM_TARGET(CudaKernels DEPENDS ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})
        SET(OPENMM_GENERATED_CXX_FILES ${OPENMM_GENERATED_CXX_FILES} ${CUDA_KERNELS_CPP})
        SET(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} CudaKernels)

        SET(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
            ${CMAKE_SOURCE_DIR}/openmm/platforms/cuda
        )
        SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
            ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/include
            ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src
        )
    ELSEIF(USE_OPENCL)
        # Compile all OpenCL kernels into one single file
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
        
        ADD_CUSTOM_TARGET(OpenCLKernels DEPENDS ${OPENCL_KERNELS_CPP} ${OPENCL_KERNELS_H})
        SET(OPENMM_GENERATED_CXX_FILES ${OPENMM_GENERATED_CXX_FILES} ${OPENCL_KERNELS_CPP})
        SET(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} OpenCLKernels)

        SET(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
            ${CMAKE_SOURCE_DIR}/openmm/platforms/opencl
        )
        SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
            ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/include
            ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/src
        )
    ENDIF()
    
ENDIF()

# Set generated files as generated to avoid warnings about missing headers
set_source_files_properties(${OPENMM_GENERATED_CXX_FILES} PROPERTIES GENERATED TRUE)

# Find which include directories actually exist
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    SET(potential_paths
        ${subdir}
        ${subdir}/include
        ${subdir}/include/openmm
        ${subdir}/include/openmm/internal
        ${subdir}/include/openmm/common
    )

    FOREACH(path ${potential_paths})
        IF(IS_DIRECTORY ${path})
            LIST(APPEND OPENMM_INCLUDE_DIRS ${path})
        ENDIF()
    ENDFOREACH()
ENDFOREACH()

SET(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit)

# Find source and header files
SET(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${OPENMM_GENERATED_CXX_FILES})

FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
    SET(OPENMM_SOURCE_C_FILES ${OPENMM_SOURCE_C_FILES} ${src_c_files})

    FILE(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp ${subdir}/base/*.cpp ${subdir}/x86/*.cpp)
    SET(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${src_cxx_files})

    FILE(GLOB incl_files ${subdir}/*.h ${subdir}/src/*.h ${subdir}/src/*/*.h)
    SET(OPENMM_SOURCE_INCLUDE_FILES ${OPENMM_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

FILE(GLOB src_files ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit/asmjit/*/*.cpp)
SET(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${src_files})

FILE(GLOB incl_files ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit/*.h)
SET(OPENMM_SOURCE_INCLUDE_FILES ${OPENMM_SOURCE_INCLUDE_FILES} ${incl_files})

# Set compile definitions for each library
set(OPENMM_COMMON_DEFS
    OPENMM_MAJOR_VERSION=8
    OPENMM_MINOR_VERSION=5
    OPENMM_BUILD_VERSION=0
    HAVE_SSE2=1
    HAVE_EMMINTRIN_H=1
    IEEE_8087=1
    LEPTON_USE_JIT=1
    OPENMM_BUILDING_STATIC_LIBRARY=0
    LEPTON_BUILDING_STATIC_LIBRARY=0
    OPENMM_COMMON_BUILDING_STATIC_LIBRARY=0
)

set_source_files_properties(${OPENMM_SOURCE_INCLUDE_FILES}
    PROPERTIES COMPILE_DEFINITIONS "${OPENMM_COMMON_DEFS}"
)
set_source_files_properties(${OPENMM_SOURCE_CXX_FILES}
    PROPERTIES COMPILE_DEFINITIONS "${OPENMM_COMMON_DEFS}"
)
set_source_files_properties(${OPENMM_SOURCE_INCLUDE_FILES}
    PROPERTIES COMPILE_DEFINITIONS "${OPENMM_COMMON_DEFS}"
)
