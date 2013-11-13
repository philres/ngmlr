cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")



if(EXISTS "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	MESSAGE(STATUS "AMD-APP-SDK 2.9 found.")
else()
	
	
	file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	file(COPY "x86_64/" DESTINATION "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
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

