cmake_minimum_required(VERSION 3.0)

set( CMAKE_C_COMPILER "gcc" )
set( CMAKE_CXX_COMPILER "g++" )

project(Robosample)
set(ROBOSAMPLE_MAJOR_VERSION 2)
set(ROBOSAMPLE_MINOR_VERSION 0)
set(ROBOSAMPLE_PATCH_VERSION 0)
set(ROBOSAMPLE_COPYRIGHT_YEARS "2021-01")

# underbar separated list of dotted authors, no spaces or commas
set(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")

# set up compilation flags
# TODO how to set for .c? we have one such file: src/format/dcdplugin.c
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC -g" )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set( ROBO_WARNINGS "-Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion -Wformat=2 -Wformat-truncation -Wundef -Wconversion -Wunused-parameter" )
    # set( ROBO_SANITIZER "-fsanitize=undefined -fsanitize=address -fsanitize-recover=address")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROBO_WARNINGS} -O0 ${ROBO_SANITIZER}" )
ELSE()
    # set up profile guided optimization
    set(ROBO_PGO)

    IF(ROBO_PGO MATCHES Generate)
        set(ROBO_PGO "-fprofile-dir=${ROBO_PGO_PATH} -fprofile-generate")
    ENDIF(ROBO_PGO MATCHES Generate)

    IF(ROBO_PGO MATCHES Use)
        set(ROBO_PGO "-fprofile-dir=${ROBO_PGO_PATH} -fprofile-use -fprofile-correction")
    ENDIF(ROBO_PGO MATCHES Use)

    set( ROBO_WARNINGS "-Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion -Wformat=2 -Wformat-truncation -Wundef -Wconversion -Wunused-parameter" )

    # https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gcc/Optimize-Options.html
    # -fopt-info-vec-all prints everything (might be worth investigating)
    # for dev only -fopt-info-vec-missed
    # -flto
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native ${ROBO_PGO}" )

    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROBO_PGO}" )
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ROBO_PGO}" )
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${ROBO_PGO}" )

ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
add_subdirectory(tests) # Laur




#########
# molmodel vs gmolmodel
# what does src/CMakeLists.txt do? why build both static and dynamic?
