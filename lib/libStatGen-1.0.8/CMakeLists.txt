cmake_minimum_required(VERSION 2.8)
project(NextGenMap)



if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/libStatGen-1.0.8/libStatGen.a")
	MESSAGE(STATUS "libStatGen-1.0.8 found.")
else()
	MESSAGE(STATUS "Downloading libStatGen-1.0.8.")
	file(DOWNLOAD "https://github.com/statgen/libStatGen/archive/v1.0.8.tar.gz" "${CMAKE_CURRENT_BINARY_DIR}/v1.0.8.tar.gz")

	MESSAGE(STATUS "Extracting libStatGen-1.0.8")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xzf "${CMAKE_CURRENT_BINARY_DIR}/v1.0.8.tar.gz" "${CMAKE_CURRENT_BINARY_DIR}/"
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
		
	execute_process(
		COMMAND make
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libStatGen-1.0.8/
	)
	
	#file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	#file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	#file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64/lib/x86_64/" "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	
	#MESSAGE(STATUS "Cleaning up")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/icd-registration.tgz")
	#file(REMOVE_RECURSE  "default-install_lnx_64.pl" "Install-AMD-APP.sh")	
endif()


