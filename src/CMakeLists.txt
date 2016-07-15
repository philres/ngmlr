cmake_minimum_required(VERSION 2.8)
project(NextGenMap)

include_directories (../lib/seqan-library-1.4.1/include)
#include_directories (../lib/bamtools-2.3.0/src/api)
include_directories (../lib/bamtools-2.3.0/src)
include_directories (../lib/mason/cpu)
include_directories (../lib/mason/opencl)
include_directories (../lib/Complete-Striped-Smith-Waterman-Library/src)
include_directories (../include)
include_directories (../lib)
include_directories (parser)
include_directories (writer)
include_directories (core)
include_directories (misc)
include_directories (log)
include_directories (config)
include_directories (.)

configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )
 	
add_executable(ngm-core
					../lib/Complete-Striped-Smith-Waterman-Library/src/ssw.c
					parser/BamParser.cpp
					writer/BAMWriter.cpp
					BitpalAligner.cpp
					parser/VcfParser.cpp																	
					config/Config.cpp					
					CS.cpp					
					CSstatic.cpp
					misc/Debug.cpp
					misc/LinearRegression.cpp													
					log/Logging.cpp
					MappedRead.cpp
					NGM_main.cpp
					NGM.cpp
					UpdateCheck.cpp
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
					StrippedSW.cpp
					SWCPU.cpp
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
TARGET_LINK_LIBRARIES(ngm-core minimapstatic)
	
add_definitions(-D _BAM) 


add_executable(ngm-core-debug
					../lib/Complete-Striped-Smith-Waterman-Library/src/ssw.c
					parser/BamParser.cpp
					writer/BAMWriter.cpp
					BitpalAligner.cpp
					parser/VcfParser.cpp																		
					config/Config.cpp					
					CS.cpp					
					CSstatic.cpp
					misc/Debug.cpp
					misc/LinearRegression.cpp													
					log/Logging.cpp
					MappedRead.cpp
					NGM_main.cpp
					NGM.cpp
					UpdateCheck.cpp
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
					StrippedSW.cpp
					SWCPU.cpp
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
TARGET_LINK_LIBRARIES(ngm-core-debug minimapstatic)

