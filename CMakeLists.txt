# use 3.9.4 for cmake lto flags
cmake_minimum_required(VERSION 3.9.4)

# define project languages here
project(Robosample C CXX)

# robosample version and date
SET(ROBOSAMPLE_MAJOR_VERSION 2)
SET(ROBOSAMPLE_MINOR_VERSION 0)
SET(ROBOSAMPLE_PATCH_VERSION 0)
SET(ROBOSAMPLE_COPYRIGHT_YEARS "2021-01")

# underbar separated list of dotted authors, no spaces or commas
SET(ROBOSAMPLE_AUTHORS "Laurentiu.Spiridon")

# C++
IF( ROBO_UBSAN MATCHES Yes )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined" )
ENDIF( ROBO_UBSAN MATCHES Yes )

# C
IF( ROBO_UBSAN MATCHES Yes )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined" )
ENDIF( ROBO_UBSAN MATCHES Yes )



# add sources
add_subdirectory(src)
# add_subdirectory(tests) # Laur