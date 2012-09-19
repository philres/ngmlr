cmake_minimum_required(VERSION 2.8)
project(MASon)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/src/opencl/tools")

INCLUDE(CPack)
INCLUDE(FindCUDA)

ADD_DEFINITIONS(-g -O3 -Wall)

find_package(OpenCL COMPONENTS)
if(OPENCL_FOUND)
	add_subdirectory(src/opencl)
else()
	MESSAGE(STATUS "Skipping OpenCL library.")
endif()


find_package(CUDA COMPONENTS)
if(CUDA_FOUND)
	add_subdirectory(src/cuda)
else()
	MESSAGE(STATUS "Skipping CUDA library.")
endif()

#add_subdirectory(src/sse)

add_subdirectory(src/example)

install (FILES "${PROJECT_SOURCE_DIR}/include/IAlignment.h" "${PROJECT_SOURCE_DIR}/include/ILog.h" "${PROJECT_SOURCE_DIR}/include/IConfig.h" DESTINATION include)