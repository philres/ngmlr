cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
project(NextGenMap)

include_directories (../lib/bamtools/src/api)
include_directories (../lib/mason/opencl)
include_directories (../include)
include_directories (parser)
include_directories (writer)
include_directories (core)
include_directories (misc)
include_directories (log)
include_directories (config)
include_directories (.)

configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )

add_executable(ngm-core 	AlignmentDispatcher.cpp
					parser/BamParser.cpp
					writer/BAMWriter.cpp					
					core/Buffer.cpp										
					config/Config.cpp					
					CS.cpp					
					misc/Debug.cpp					
					writer/FFormatWriter.cpp					
					log/Logging.cpp
					MappedRead.cpp
					NGM_main.cpp
					NGM.cpp
					core/NGMStats.cpp
					core/NGMTask.cpp
					Output.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					parser/SamParser.cpp
					writer/SAMWriter.cpp
					SequenceProvider.cpp
					SW.cpp
					core/unix.cpp
					core/unix_threads.cpp
					core/windows_threads.cpp
					core/windows.cpp)


ADD_CUSTOM_TARGET(recompile_always ALL touch "${CMAKE_SOURCE_DIR}/src/NGM_main.cpp")

	#STATIC
    #set(CMAKE_EXE_LINKER_FLAGS -static) 
    
    #set(CMAKE_EXE_LINKER_FLAGS " -static-libgcc -static-libstdc++")
    #set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 
    #set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic 
    #set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS) 
    #set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC 
    #set(CMAKE_SHARED_LIBRARY_CXX_FLAGS) 
    #set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic 
    #set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
     
    # Maybe this works as well, haven't tried yet. 
	#set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE) 
	#target_link_libraries(ngm-core /scratch/philipp/src/zlib-1.2.7/zlib/lib/libz.a)
	

#target_link_libraries(ngm-core z)					
#target_link_libraries(ngm-core dl)
target_link_libraries(ngm-core pthread)
target_link_libraries(ngm-core MASonOpenCl)
#target_link_libraries(ngm-core /software/ngm/src/bamtools/lib/libbamtools.a)	

#if( ${BAMTOOLS_FOUND} )		
	#ADD_LIBRARY(bamtools STATIC)
	#SET_TARGET_PROPERTIES(bamtools PROPERTIES IMPORTED_LOCATION ${BAMTOOLS_LIBRARY})
	#include_directories(${BAMTOOLS_INCLUDE_DIRS})
	TARGET_LINK_LIBRARIES(ngm-core BamTools-static)
	TARGET_LINK_LIBRARIES(ngm-core zlibstatic)
	#TARGET_LINK_LIBRARIES(ngm-core zlib)
	
	add_definitions(-D _BAM) 
#endif()

#find_library(HOARD hoard)
#if(HOARD) 
#	message(" --Found libhoard: ${HOARD}")	
#	target_link_libraries(ngm-core ${HOARD})
#else()
#	message(" --Libhoard not found!")	
#endif()		

