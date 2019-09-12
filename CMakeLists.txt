cmake_minimum_required(VERSION 3.0)

set( CMAKE_C_COMPILER "gcc" )
set( CMAKE_CXX_COMPILER "g++" )

project( gmolmodel VERSION 0.1.1 LANGUAGES C CXX )

set( CMAKE_CXX_FLAGS "-std=c++11 -Wall -fPIC -g" )


add_subdirectory(src)
add_subdirectory(tests) # Laur
