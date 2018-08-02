  cmake_minimum_required(VERSION 2.8)
project(ngmlr)

include(CheckCXXCompilerFlag)

option(STATIC "Build static binary" OFF)

set( NGM_VERSION_MAJOR 0 )
set( NGM_VERSION_MINOR 2 )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode!")
    set( NGM_VERSION_BUILD 8-debug )
else()
  set( NGM_VERSION_BUILD 8 )
ENDIF()

find_package( ZLIB REQUIRED )
find_package ( Threads REQUIRED )

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    IF(NOT APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    ELSE()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
    ENDIF()
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/ngmlr-${NGM_VERSION_MAJOR}.${NGM_VERSION_MINOR}.${NGM_VERSION_BUILD}/)
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

# Set a default build type for single-configuration
# CMake generators if no build type is set.
# RELEASE, RELWITHDEBINFO, DEBUG
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Using 'RELWITHDEBINFO'")
  SET(CMAKE_BUILD_TYPE RELWITHDEBINFO)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

add_subdirectory(src)