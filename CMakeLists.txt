cmake_minimum_required(VERSION 2.8)
project(NextGenMap-LR)

set( NGM_VERSION_MAJOR 0 )
set( NGM_VERSION_MINOR 1 )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode!")
    set( NGM_VERSION_BUILD 4-debug )
else()
  set( NGM_VERSION_BUILD 4 )
ENDIF()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/ngmlr-${NGM_VERSION_MAJOR}.${NGM_VERSION_MINOR}.${NGM_VERSION_BUILD}/)
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Using 'release'")
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/lib/zlib-1.2.7")
include_directories("${PROJECT_SOURCE_DIR}/lib/zlib-1.2.7")

add_subdirectory(lib/zlib-1.2.7)
add_subdirectory(lib/bamtools-2.3.0)

add_subdirectory(src)