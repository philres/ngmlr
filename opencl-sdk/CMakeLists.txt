cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")



if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64")
	MESSAGE(STATUS "AMD-APP-SDK 2.8 found.")
else()
	MESSAGE(STATUS "Downloading AMD-APP-SDK 2.8. This might take a while (~100mb)")
	file(DOWNLOAD "http://developer.amd.com/wordpress/media/download.php?f=QU1ELUFQUC1TREstdjIuOC1sbng2NC50Z3o=" "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-lnx64.tgz" EXPECTED_MD5 cf8ae9eecc5e3573486e1bc73e0b316c)

	MESSAGE(STATUS "Extracting AMD-APP-SDK 2.8")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xzf "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-lnx64.tgz" "${CMAKE_CURRENT_BINARY_DIR}/"
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
	
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xzf "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64.tgz"
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
	
	file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64/lib/x86_64/" "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	
	#MESSAGE(STATUS "Cleaning up")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/icd-registration.tgz")
	#file(REMOVE_RECURSE  "default-install_lnx_64.pl" "Install-AMD-APP.sh")	
endif()


file(WRITE "${EXECUTABLE_OUTPUT_PATH}/opencl/vendor/amdocl64.icd" "libamdocl64.so")
file(WRITE "${EXECUTABLE_OUTPUT_PATH}/opencl/vendor/nvidia.icd" "libcuda.so")

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

file(COPY "${CMAKE_CURRENT_BINARY_DIR}/ngm" DESTINATION "${EXECUTABLE_OUTPUT_PATH}" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ NO_SOURCE_PERMISSIONS)