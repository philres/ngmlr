cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

set( NGM_VERSION_MAJOR 0 )
set( NGM_VERSION_MINOR 4 )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode!")
    set( NGM_VERSION_BUILD 12-debug )
else()
	set( NGM_VERSION_BUILD 12 )
ENDIF()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/ngm-${NGM_VERSION_MAJOR}.${NGM_VERSION_MINOR}.${NGM_VERSION_BUILD}/)
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

#add_definitions(-D INSTANCE_COUNTING)
#VERBOSE OUTPUT 
#add_definitions(-D VERBOSE)
#DEBUG OUTPUT FOR CS (FILE: cs-results.txt)
#add_definitions(-D _DEBUGCS)
#DEBUG OUTPUT FOR CS (CONSOLE)
#add_definitions(-D _DEBUGCSVERBOSE)
#DEBUG OUTPUT FOR SCORE COMPTATION
#add_definitions(-D _DEBUGSW)
#add_definitions(-D _DEBUGOUT)
#add_definitions(-D _DEBUGCSKMERS)
#add_definitions(-D _DEBUGRP)
#Print bed file with all CMRs
#add_definitions(-D _DEBUGCMRS)

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type specified. Using 'release'")
	SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

#SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -L${EXECUTABLE_OUTPUT_PATH}/opencl/lib -Wl,-rpath,'\$ORIGIN/amd-app-sdk/'")
#SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0 -Wall -L${EXECUTABLE_OUTPUT_PATH}/opencl/lib -Wl,-rpath,'\$ORIGIN/amd-app-sdk/'")

SET (CMAKE_CXX_FLAGS_RELEASE "-L${EXECUTABLE_OUTPUT_PATH}/opencl/lib")
SET (CMAKE_CXX_FLAGS_DEBUG "-L${EXECUTABLE_OUTPUT_PATH}/opencl/lib")

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_subdirectory(opencl-sdk)

set(OPENCL_LIBRARIES " -lOpenCL")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/lib/zlib-1.2.7")
include_directories("${PROJECT_SOURCE_DIR}/lib/zlib-1.2.7")

add_subdirectory(lib/zlib-1.2.7)
add_subdirectory(lib/bamtools-2.3.0)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/opencl-sdk/AMD-APP-SDK-v2.8.1.0-RC-lnx64/include")

add_subdirectory(lib/mason/opencl)
add_subdirectory(src)

add_subdirectory(utils)
