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

set( CMAKE_CXX_FLAGS "-std=c++11 -Wall -fPIC -g" )

add_subdirectory(src)
add_subdirectory(tests) # Laur

# These are used in Doxyfile.in and SimbodyConfig.cmake.in.
set(ROBOSAMPLE_INSTALL_DOXYGENDIR   "${CMAKE_INSTALL_DOCDIR}/api")
set(ROBOSAMPLE_DOXYGEN_TAGFILE_NAME "RobosampleDoxygenTagfile")
add_subdirectory(doc)

