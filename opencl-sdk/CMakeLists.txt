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

