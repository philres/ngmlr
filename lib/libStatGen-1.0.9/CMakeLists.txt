cmake_minimum_required(VERSION 2.8)
project(NextGenMap)



#if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/libStatGen-1.0.8/libStatGen.a")
#	MESSAGE(STATUS "libStatGen-1.0.8 found.")
#else()
#	MESSAGE(STATUS "Downloading libStatGen-1.0.8.")
#	set(url  "http://github.com/statgen/libStatGen/archive/v1.0.8.tar.gz")
#	file(DOWNLOAD ${url} "${CMAKE_CURRENT_BINARY_DIR}/v1.0.8.tar.gz" 
#	    STATUS status  
#	    EXPECTED_MD5 4499f46291ff5cc25ace20e8c2fb741c)

#    list(GET status 0 status_code) 
#    list(GET status 1 status_string) 
    
#    if(NOT status_code EQUAL 0)
#        message( FATAL_ERROR "Failed to download libStatGen. Please check whether ${url} is available. If not, please try again later. Error: ${error_message}" ) 
#        #message(FATAL_ERROR "error: downloading 'http://zlib.net/zlib-1.2.7.tar.gz' failed status_code: ${status_code}  status_string: ${status_string}") 
#    endif() 


	MESSAGE(STATUS "Extracting libStatGen-1.0.9")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xzf "${PROJECT_SOURCE_DIR}/LibStatGen.1.0.9.tgz" "${CMAKE_CURRENT_BINARY_DIR}/"
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
		

	execute_process(
		COMMAND make		
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libStatGen_1.0.9/
	)
	
	#file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	#file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}opencl/lib)
	#file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64/lib/x86_64/" "${EXECUTABLE_OUTPUT_PATH}opencl/lib")
	
	#MESSAGE(STATUS "Cleaning up")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/AMD-APP-SDK-v2.8-RC-lnx64.tgz")
	#file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/icd-registration.tgz")
	#file(REMOVE_RECURSE  "default-install_lnx_64.pl" "Install-AMD-APP.sh")	
#endif()


