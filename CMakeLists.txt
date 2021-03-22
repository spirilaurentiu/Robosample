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

# set up profile guided optimization
set(PGO)

IF(ROBO_PGO MATCHES Generate)
    set(PGO "-fprofile-dir=pgo -fprofile-generate=pgo")
ENDIF(ROBO_PGO MATCHES Generate)

IF(ROBO_PGO MATCHES Profile)
    set(PGO "-fprofile-dir=pgo -fprofile-use=pgo -fprofile-correction")
ENDIF(ROBO_PGO MATCHES Profile)


# set up compilation flags
set( CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -Wshadow -Wdouble-promotion -Wformat=2 -Wformat-truncation -Wundef -Wconversion =Wunused-parameter -fPIC -g" )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
ELSE()
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -flto -march=native ${PGO}" )
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
add_subdirectory(tests) # Laur

