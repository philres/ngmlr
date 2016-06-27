cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

set( NGM_VERSION_MAJOR 0 )
set( NGM_VERSION_MINOR 1 )
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode!")
    set( NGM_VERSION_BUILD 13-debug )
else()
  set( NGM_VERSION_BUILD 1 )
ENDIF()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/ngmlr-${NGM_VERSION_MAJOR}.${NGM_VERSION_MINOR}.${NGM_VERSION_BUILD}/)
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

IF(NOT APPLE)
	SET (CMAKE_CXX_FLAGS_RELEASE "-L${EXECUTABLE_OUTPUT_PATH}/opencl/lib")
	SET (CMAKE_CXX_FLAGS_DEBUG "-L${EXECUTABLE_OUTPUT_PATH}/opencl/lib")
endif()



IF(NOT APPLE)
	add_subdirectory(opencl-sdk)
endif()

IF(APPLE) 
  set(GDB "ggdb")
  set(OPENCL_LIBRARIES "-framework OpenCL")
ELSE()
  set(GDB "gdb")
  set(OPENCL_LIBRARIES "-lOpenCL")
ENDIF()

include_directories("${CMAKE_CURRENT_BINARY_DIR}/lib/zlib-1.2.7")
include_directories("${PROJECT_SOURCE_DIR}/lib/zlib-1.2.7")

add_subdirectory(lib/zlib-1.2.7)
add_subdirectory(lib/bamtools-2.3.0)

IF(NOT APPLE)
	include_directories("${CMAKE_CURRENT_BINARY_DIR}/opencl-sdk/AMD-APP-SDK-v2.8.1.0-RC-lnx64/include")
endif()

add_subdirectory(lib/mason/opencl)
add_subdirectory(src)

add_subdirectory(utils)


#Writing wrapper scripts for NextGenMape: ngm, ngm-debug, ngm-log

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ngm" "#!/bin/bash

pushd . > /dev/null
SCRIPT_PATH=\"\${BASH_SOURCE[0]}\";
  while([ -h \"\${SCRIPT_PATH}\" ]) do 
    cd \"`dirname \"\${SCRIPT_PATH}\"`\"
    SCRIPT_PATH=\"$(readlink \"`basename \"\${SCRIPT_PATH}\"`\")\"; 
  done
cd \"`dirname \"\${SCRIPT_PATH}\"`\" > /dev/null
SCRIPT_PATH=\"`pwd`\";
popd  > /dev/null
#echo \"srcipt=[\${SCRIPT_PATH}]\"
#echo \"pwd   =[`pwd`]\"


OPENCL_VENDOR_PATH=\"$SCRIPT_PATH/opencl/vendor\" LD_LIBRARY_PATH=\"$SCRIPT_PATH/opencl/lib\" exec $SCRIPT_PATH/ngm-core $@
")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ngm-log" "#!/bin/bash

pushd . > /dev/null
SCRIPT_PATH=\"\${BASH_SOURCE[0]}\";
  while([ -h \"\${SCRIPT_PATH}\" ]) do 
    cd \"`dirname \"\${SCRIPT_PATH}\"`\"
    SCRIPT_PATH=\"$(readlink \"`basename \"\${SCRIPT_PATH}\"`\")\"; 
  done
cd \"`dirname \"\${SCRIPT_PATH}\"`\" > /dev/null
SCRIPT_PATH=\"`pwd`\";
popd  > /dev/null
#echo \"srcipt=[\${SCRIPT_PATH}]\"
#echo \"pwd   =[`pwd`]\"


OPENCL_VENDOR_PATH=\"$SCRIPT_PATH/opencl/vendor\" LD_LIBRARY_PATH=\"$SCRIPT_PATH/opencl/lib\" exec $SCRIPT_PATH/ngm-core-debug $@
")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ngm-debug" "#!/bin/bash

pushd . > /dev/null
SCRIPT_PATH=\"\${BASH_SOURCE[0]}\";
  while([ -h \"\${SCRIPT_PATH}\" ]) do 
    cd \"`dirname \"\${SCRIPT_PATH}\"`\"
    SCRIPT_PATH=\"$(readlink \"`basename \"\${SCRIPT_PATH}\"`\")\"; 
  done
cd \"`dirname \"\${SCRIPT_PATH}\"`\" > /dev/null
SCRIPT_PATH=\"`pwd`\";
popd  > /dev/null
#echo \"srcipt=[\${SCRIPT_PATH}]\"
#echo \"pwd   =[`pwd`]\"


OPENCL_VENDOR_PATH=\"$SCRIPT_PATH/opencl/vendor\" LD_LIBRARY_PATH=\"$SCRIPT_PATH/opencl/lib\" exec ${GDB} --args $SCRIPT_PATH/ngm-core-debug $@
")

file(COPY "${CMAKE_CURRENT_BINARY_DIR}/ngm-debug" DESTINATION "${EXECUTABLE_OUTPUT_PATH}" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ NO_SOURCE_PERMISSIONS)
file(COPY "${CMAKE_CURRENT_BINARY_DIR}/ngm" DESTINATION "${EXECUTABLE_OUTPUT_PATH}" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ NO_SOURCE_PERMISSIONS)
file(COPY "${CMAKE_CURRENT_BINARY_DIR}/ngm-log" DESTINATION "${EXECUTABLE_OUTPUT_PATH}" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ NO_SOURCE_PERMISSIONS)


