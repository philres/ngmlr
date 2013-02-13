cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

set( NGM_VERSION_MAJOR 0 )
set( NGM_VERSION_MINOR 4 )
set( NGM_VERSION_BUILD 4 )

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

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -g -L${EXECUTABLE_OUTPUT_PATH}/opencl/lib -Wl,-rpath,'\$ORIGIN/amd-app-sdk/'")
SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -L${EXECUTABLE_OUTPUT_PATH}/opencl/lib -Wl,-rpath,'\$ORIGIN/amd-app-sdk/'")
SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0 -Wall -L${EXECUTABLE_OUTPUT_PATH}/opencl/lib -Wl,-rpath,'\$ORIGIN/amd-app-sdk/'")

#include(FindOpenMP)
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#FIND_PACKAGE(OpenCL REQUIRED)
#if(OPENCL_FOUND)
#	MESSAGE(STATUS "OpenCL include dir: ${OPENCL_INCLUDE_DIRS}")
#endif()
add_subdirectory(opencl-sdk)
set(OPENCL_LIBRARIES " -lOpenCL")

FIND_PACKAGE(ZLIB REQUIRED)
#FIND_PACKAGE(Bamtools QUIET COMPONENTS)

#ADD_DEFINITIONS(-g -O3)

    #SET (CMAKE_EXE_LINKER_FLAGS
    #"${CMAKE_EXE_LINKER_FLAGS} -L${PROJECT_SOURCE_DIR}/bin/openecl/lib -Wl,-rpath,'$ORIGIN/amd-app-sdk/'" )

    #SET (CMAKE_SHARED_LINKER_FLAGS
    #"${CMAKE_SHARED_LINKER_FLAGS} -L${PROJECT_SOURCE_DIR}/bin/opencl/lib -Wl,-rpath,'$ORIGIN/amd-app-sdk'" )

#SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN'" )
#SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN:$ORIGIN/amd-app-sdk'" )
#SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN/mason:$ORIGIN/amd-app/lib/x86_64'" )

add_subdirectory(lib/bamtools)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/opencl-sdk/AMD-APP-SDK-v2.8-RC-lnx64/include")
#LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/bin/amd-app-sdk")
#include_directories(${OPENCL_INCLUDE_DIRS})

add_subdirectory(lib/mason/opencl)
add_subdirectory(src)

add_subdirectory(utils)
