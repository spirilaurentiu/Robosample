set(OPENMM_GENERATED_CXX_FILES)

# Get openmm source directories
set(OPENMM_INCLUDE_DIRS)
set(OPENMM_SOURCE_SUBDIRS
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
    ${CMAKE_SOURCE_DIR}/openmm/platforms/reference
    ${CMAKE_SOURCE_DIR}/openmm/serialization)

if(USE_CPU)
  set(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
                            ${CMAKE_SOURCE_DIR}/openmm/platforms/cpu)
  set(OPENMM_INCLUDE_DIRS
      ${OPENMM_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cpu/include
      ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cpu/src)

elseif(USE_CUDA OR USE_OPENCL)
  set(COMMON_KERNEL_SOURCE_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/src")
  set(COMMON_KERNEL_SOURCE_CLASS CommonKernelSources)
  set(COMMON_KERNELS_CPP
      ${COMMON_KERNEL_SOURCE_DIR}/${COMMON_KERNEL_SOURCE_CLASS}.cpp)
  set(COMMON_KERNELS_H
      ${COMMON_KERNEL_SOURCE_DIR}/${COMMON_KERNEL_SOURCE_CLASS}.h)
  file(GLOB COMMON_KERNELS ${COMMON_KERNEL_SOURCE_DIR}/kernels/*.cc)
  set(ENCODE_KERNELS_DEPS
      ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
      ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/strip_comments.py)
  add_custom_command(
    OUTPUT ${COMMON_KERNELS_CPP} ${COMMON_KERNELS_H}
    COMMAND
      ${CMAKE_COMMAND} ARGS -D KERNEL_SOURCE_DIR=${COMMON_KERNEL_SOURCE_DIR} -D
      KERNELS_CPP=${COMMON_KERNELS_CPP} -D KERNELS_H=${COMMON_KERNELS_H} -D
      KERNEL_SOURCE_CLASS=${COMMON_KERNEL_SOURCE_CLASS} -D
      KERNEL_FILE_EXTENSION=cc -P
      ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
    DEPENDS ${COMMON_KERNELS} ${ENCODE_KERNELS_DEPS}
    COMMENT "Generating common kernel sources for OpenMM...")

  # This command is executed when building, not when running CMakeLists.txt
  add_custom_target(CommonKernels DEPENDS ${COMMON_KERNELS_CPP}
                                          ${COMMON_KERNELS_H})
  set(OPENMM_GENERATED_CXX_FILES ${OPENMM_GENERATED_CXX_FILES}
                                 ${COMMON_KERNELS_CPP})
  set(OPENMM_DEPENDENCIES CommonKernels)

  set(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
                            ${CMAKE_SOURCE_DIR}/openmm/platforms/common)
  set(OPENMM_INCLUDE_DIRS
      ${OPENMM_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/include
      ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/common/src)

  if(USE_CUDA)
    # Compile all CUDA kernels into one single file
    set(CUDA_KERNEL_SOURCE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src")
    set(CUDA_KERNEL_SOURCE_CLASS CudaKernelSources)
    set(CUDA_KERNELS_CPP
        ${CUDA_KERNEL_SOURCE_DIR}/${CUDA_KERNEL_SOURCE_CLASS}.cpp)
    set(CUDA_KERNELS_H ${CUDA_KERNEL_SOURCE_DIR}/${CUDA_KERNEL_SOURCE_CLASS}.h)
    file(GLOB CUDA_KERNELS ${CUDA_KERNEL_SOURCE_DIR}/kernels/*.cu)
    add_custom_command(
      OUTPUT ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H}
      COMMAND
        ${CMAKE_COMMAND} ARGS -D KERNEL_SOURCE_DIR=${CUDA_KERNEL_SOURCE_DIR} -D
        KERNELS_CPP=${CUDA_KERNELS_CPP} -D KERNELS_H=${CUDA_KERNELS_H} -D
        KERNEL_SOURCE_CLASS=${CUDA_KERNEL_SOURCE_CLASS} -D
        KERNEL_FILE_EXTENSION=cu -P
        ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
      DEPENDS ${CUDA_KERNELS} ${ENCODE_KERNELS_DEPS}
      COMMENT "Generating CUDA kernel sources for OpenMM...")

    add_custom_target(CudaKernels DEPENDS ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})
    set(OPENMM_GENERATED_CXX_FILES ${OPENMM_GENERATED_CXX_FILES}
                                   ${CUDA_KERNELS_CPP})
    set(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} CudaKernels)

    set(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
                              ${CMAKE_SOURCE_DIR}/openmm/platforms/cuda)
    set(OPENMM_INCLUDE_DIRS
        ${OPENMM_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/include
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/cuda/src)
  elseif(USE_OPENCL)
    # Compile all OpenCL kernels into one single file
    set(OPENCL_KERNEL_SOURCE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/src")
    set(OPENCL_KERNEL_SOURCE_CLASS OpenCLKernelSources)
    set(OPENCL_KERNELS_CPP
        ${OPENCL_KERNEL_SOURCE_DIR}/${OPENCL_KERNEL_SOURCE_CLASS}.cpp)
    set(OPENCL_KERNELS_H
        ${OPENCL_KERNEL_SOURCE_DIR}/${OPENCL_KERNEL_SOURCE_CLASS}.h)
    file(GLOB OPENCL_KERNELS ${OPENCL_KERNEL_SOURCE_DIR}/kernels/*.cl)
    add_custom_command(
      OUTPUT ${OPENCL_KERNELS_CPP} ${OPENCL_KERNELS_H}
      COMMAND
        ${CMAKE_COMMAND} ARGS -D KERNEL_SOURCE_DIR=${OPENCL_KERNEL_SOURCE_DIR}
        -D KERNELS_CPP=${OPENCL_KERNELS_CPP} -D KERNELS_H=${OPENCL_KERNELS_H} -D
        KERNEL_SOURCE_CLASS=${OPENCL_KERNEL_SOURCE_CLASS} -D
        KERNEL_FILE_EXTENSION=cl -P
        ${CMAKE_SOURCE_DIR}/openmm/cmake_modules/EncodeKernelFiles.cmake
      DEPENDS ${OPENCL_KERNELS} ${ENCODE_KERNELS_DEPS}
      COMMENT "Generating OpenCL kernel sources for OpenMM...")

    add_custom_target(OpenCLKernels DEPENDS ${OPENCL_KERNELS_CPP}
                                            ${OPENCL_KERNELS_H})
    set(OPENMM_GENERATED_CXX_FILES ${OPENMM_GENERATED_CXX_FILES}
                                   ${OPENCL_KERNELS_CPP})
    set(OPENMM_DEPENDENCIES ${OPENMM_DEPENDENCIES} OpenCLKernels)

    set(OPENMM_SOURCE_SUBDIRS ${OPENMM_SOURCE_SUBDIRS}
                              ${CMAKE_SOURCE_DIR}/openmm/platforms/opencl)
    set(OPENMM_INCLUDE_DIRS
        ${OPENMM_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/include
        ${CMAKE_CURRENT_SOURCE_DIR}/openmm/platforms/opencl/src)
  endif()

endif()

# Set generated files as generated to avoid warnings about missing headers
set_source_files_properties(${OPENMM_GENERATED_CXX_FILES} PROPERTIES GENERATED
                                                                     TRUE)

# Find which include directories actually exist
foreach(subdir ${OPENMM_SOURCE_SUBDIRS})
  set(potential_paths
      ${subdir} ${subdir}/include ${subdir}/include/openmm
      ${subdir}/include/openmm/internal ${subdir}/include/openmm/common)

  foreach(path ${potential_paths})
    if(IS_DIRECTORY ${path})
      list(APPEND OPENMM_INCLUDE_DIRS ${path})
    endif()
  endforeach()
endforeach()

set(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIRS}
                        ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit)

# Find source and header files
set(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES}
                            ${OPENMM_GENERATED_CXX_FILES})

foreach(subdir ${OPENMM_SOURCE_SUBDIRS})
  file(GLOB src_c_files ${subdir}/src/*.c ${subdir}/src/*/*.c)
  set(OPENMM_SOURCE_C_FILES ${OPENMM_SOURCE_C_FILES} ${src_c_files})

  file(GLOB src_cxx_files ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp
       ${subdir}/base/*.cpp ${subdir}/x86/*.cpp)
  set(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${src_cxx_files})

  file(GLOB incl_files ${subdir}/*.h ${subdir}/src/*.h ${subdir}/src/*/*.h)
  set(OPENMM_SOURCE_INCLUDE_FILES ${OPENMM_SOURCE_INCLUDE_FILES} ${incl_files})
endforeach(subdir)

file(GLOB src_files ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit/asmjit/*/*.cpp)
set(OPENMM_SOURCE_CXX_FILES ${OPENMM_SOURCE_CXX_FILES} ${src_files})

file(GLOB incl_files ${CMAKE_SOURCE_DIR}/openmm/libraries/asmjit/*.h)
set(OPENMM_SOURCE_INCLUDE_FILES ${OPENMM_SOURCE_INCLUDE_FILES} ${incl_files})

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
    OPENMM_COMMON_BUILDING_STATIC_LIBRARY=0)

set_source_files_properties(
  ${OPENMM_SOURCE_INCLUDE_FILES} PROPERTIES COMPILE_DEFINITIONS
                                            "${OPENMM_COMMON_DEFS}")
set_source_files_properties(
  ${OPENMM_SOURCE_CXX_FILES} PROPERTIES COMPILE_DEFINITIONS
                                        "${OPENMM_COMMON_DEFS}")
set_source_files_properties(
  ${OPENMM_SOURCE_INCLUDE_FILES} PROPERTIES COMPILE_DEFINITIONS
                                            "${OPENMM_COMMON_DEFS}")
