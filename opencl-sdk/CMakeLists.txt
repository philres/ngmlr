cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")



if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8.1.0-RC-lnx64" AND EXISTS "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	MESSAGE(STATUS "AMD-APP-SDK 2.8.1 found.")
else()
	MESSAGE(STATUS "Downloading AMD-APP-SDK 2.8.1. This might take a while (~100mb)")
	
	
	#AMD-APP-SDK 2.8.1
	file(DOWNLOAD "http://developer.amd.com/wordpress/media/download.php?f=QU1ELUFQUC1TREstdjIuOC4xLjAtbG54NjQudGd6=" "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8.1.0-lnx64.tgz" 
	    STATUS status 
	    EXPECTED_MD5 269b323e4b5f64735b3f8309000929c4)
	#AMD-APP-SDK 2.8
	#file(DOWNLOAD "http://developer.amd.com/wordpress/media/download.php?f=QU1ELUFQUC1TREstdjIuOC1sbng2NC50Z3o=" "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-lnx64.tgz" EXPECTED_MD5 cf8ae9eecc5e3573486e1bc73e0b316c)
    list(GET status 0 status_code) 
    list(GET status 1 status_string) 
    
    if(NOT status_code EQUAL 0)
        message( FATAL_ERROR "Failed to download libStatGen. Please check whether ${url} is available. If not, please try again later. Error: ${error_message}" ) 
        #message(FATAL_ERROR "error: downloading 'http://zlib.net/zlib-1.2.7.tar.gz' failed status_code: ${status_code}  status_string: ${status_string}") 
    endif() 

	MESSAGE(STATUS "Extracting AMD-APP-SDK 2.8.1")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xzf "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8.1.0-lnx64.tgz" "${CMAKE_CURRENT_BINARY_DIR}/"
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
	
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xzf "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8.1.0-RC-lnx64.tgz"
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
	
	file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8.1.0-RC-lnx64/lib/x86_64/" "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	execute_process(COMMAND chmod go+rx "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	
	#MESSAGE(STATUS "Cleaning up")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8.1.0-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/icd-registration.tgz")
	#file(REMOVE_RECURSE  "default-install_lnx_64.pl" "Install-AMD-APP.sh")	
endif()



file(WRITE "${EXECUTABLE_OUTPUT_PATH}/opencl/vendor/amdocl64.icd" "libamdocl64.so")
file(WRITE "${EXECUTABLE_OUTPUT_PATH}/opencl/vendor/nvidiaold.icd" "libcuda.so")
file(WRITE "${EXECUTABLE_OUTPUT_PATH}/opencl/vendor/nvidia.icd" "libnvidia-opencl.so.1")

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


OPENCL_VENDOR_PATH=\"$SCRIPT_PATH/opencl/vendor\" LD_LIBRARY_PATH=\"$SCRIPT_PATH/opencl/lib\" exec gdb --args $SCRIPT_PATH/ngm-core $@
")

file(COPY "${CMAKE_CURRENT_BINARY_DIR}/ngm-debug" DESTINATION "${EXECUTABLE_OUTPUT_PATH}" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ NO_SOURCE_PERMISSIONS)
file(COPY "${CMAKE_CURRENT_BINARY_DIR}/ngm" DESTINATION "${EXECUTABLE_OUTPUT_PATH}" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ NO_SOURCE_PERMISSIONS)

