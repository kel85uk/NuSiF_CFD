cmake_minimum_required( VERSION 2.6 ) 

project ( NuSiF )

#set( CMAKE_CXX_COMPILER "icpc")
set( CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -Wall -Wconversion -Wshadow" )

# Default build type
if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Debug CACHE STRING "Build Types: Debug Release"  FORCE )
endif()
set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release )

add_subdirectory( src   )
add_subdirectory( tests )
