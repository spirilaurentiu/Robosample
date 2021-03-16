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

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set( CMAKE_CXX_FLAGS "-std=c++17 -Wall -fPIC -g" )
ELSE()
    set( CMAKE_CXX_FLAGS "-std=c++17 -Wall -fPIC -g -O2" )
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
add_subdirectory(tests) # Laur

