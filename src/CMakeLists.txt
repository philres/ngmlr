cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

include_directories (../lib/seqan-library-1.4.1/include)
#include_directories (../lib/bamtools-2.3.0/src/api)
include_directories (../lib/bamtools-2.3.0/src)
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

add_executable(ngm-core 	
					parser/BamParser.cpp
					parser/VcfParser.cpp
					writer/BAMWriter.cpp																	
					config/Config.cpp					
					CS.cpp					
					CSstatic.cpp
					misc/Debug.cpp													
					log/Logging.cpp
					MappedRead.cpp
					NGM_main.cpp
					NGM.cpp
					core/NGMStats.cpp
					core/NGMTask.cpp
					AlignmentBuffer.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					parser/SamParser.cpp
					writer/SAMWriter.cpp
					writer/ScoreWriter.cpp
					seqan/EndToEndAffine.cpp
					SequenceProvider.cpp	
					OutputReadBuffer.cpp	
					ScoreBuffer.cpp
					core/unix.cpp
					core/unix_threads.cpp
					core/windows_threads.cpp
					core/windows.cpp)										
					
ADD_CUSTOM_TARGET(recompile_always ALL touch "${CMAKE_SOURCE_DIR}/src/NGM_main.cpp")

SET_TARGET_PROPERTIES(ngm-core PROPERTIES COMPILE_FLAGS "-g -O3 -DNDEBUG")

target_link_libraries(ngm-core pthread)
target_link_libraries(ngm-core MASonOpenCl)	
TARGET_LINK_LIBRARIES(ngm-core BamTools-static)
TARGET_LINK_LIBRARIES(ngm-core zlibstatic)
	
add_definitions(-D _BAM) 


add_executable(ngm-core-debug	
					parser/BamParser.cpp
					writer/BAMWriter.cpp																	
					config/Config.cpp					
					CS.cpp					
					CSstatic.cpp
					misc/Debug.cpp													
					log/Logging.cpp
					MappedRead.cpp
					NGM_main.cpp
					NGM.cpp
					core/NGMStats.cpp
					core/NGMTask.cpp
					AlignmentBuffer.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					parser/SamParser.cpp
					writer/SAMWriter.cpp
					writer/ScoreWriter.cpp
					seqan/EndToEndAffine.cpp
					SequenceProvider.cpp	
					OutputReadBuffer.cpp				
					ScoreBuffer.cpp
					core/unix.cpp
					core/unix_threads.cpp
					core/windows_threads.cpp
					core/windows.cpp)



set_target_properties(ngm-core-debug PROPERTIES COMPILE_DEFINITIONS "DEBUGLOG")

SET_TARGET_PROPERTIES(ngm-core-debug PROPERTIES COMPILE_FLAGS "-g3 -O0")

target_link_libraries(ngm-core-debug pthread)
target_link_libraries(ngm-core-debug MASonOpenCl-debug)	
TARGET_LINK_LIBRARIES(ngm-core-debug BamTools-static)
TARGET_LINK_LIBRARIES(ngm-core-debug zlibstatic)

