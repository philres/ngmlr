cmake_minimum_required(VERSION 2.8)
project(NextGenMap-LR)

include_directories (../lib/seqan-library-1.4.1/include)
include_directories (../lib/bamtools-2.3.0/src)
include_directories (../lib/Complete-Striped-Smith-Waterman-Library/src)
include_directories (../lib)

configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )

#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")
 	
add_executable(ngmlr
                    AlignmentMatrix.cpp
                    AlignmentMatrixFast.cpp
                    ArgParser.cpp
					../lib/Complete-Striped-Smith-Waterman-Library/src/ssw.c
					BamParser.cpp
					BAMWriter.cpp
					BitpalAligner.cpp
					VcfParser.cpp																					
					CS.cpp					
					CSstatic.cpp
					ConvexAlign.cpp
					ConvexAlignFast.cpp
					LinearRegression.cpp													
					Logging.cpp
					MappedRead.cpp
					main.cpp
					NGM.cpp
					UpdateCheck.cpp
					NGMTask.cpp
					AlignmentBuffer.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					SamParser.cpp
					SAMWriter.cpp
					EndToEndAffine.cpp
					SequenceProvider.cpp	
					OutputReadBuffer.cpp	
					ScoreBuffer.cpp
					StrippedSW.cpp
					unix.cpp
					unix_threads.cpp)										
					
ADD_CUSTOM_TARGET(recompile_always ALL touch "${CMAKE_SOURCE_DIR}/src/main.cpp")

SET_TARGET_PROPERTIES(ngmlr PROPERTIES COMPILE_FLAGS " -g -O3 -DNDEBUG")

TARGET_LINK_LIBRARIES(ngmlr pthread)	
TARGET_LINK_LIBRARIES(ngmlr BamTools-static)
TARGET_LINK_LIBRARIES(ngmlr zlibstatic)
	
add_definitions(-D _BAM) 

add_executable(ngmlr-debug
                    AlignmentMatrix.cpp
                    AlignmentMatrixFast.cpp
                    ArgParser.cpp
					../lib/Complete-Striped-Smith-Waterman-Library/src/ssw.c
					BamParser.cpp
					BAMWriter.cpp
					BitpalAligner.cpp
					VcfParser.cpp																							
					CS.cpp					
					CSstatic.cpp
					ConvexAlign.cpp
					ConvexAlignFast.cpp
					LinearRegression
					Logging.cpp
					MappedRead.cpp
					main.cpp
					NGM.cpp
					UpdateCheck.cpp
					NGMTask.cpp
					AlignmentBuffer.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					SamParser.cpp
					SAMWriter.cpp
					EndToEndAffine.cpp
					SequenceProvider.cpp	
					OutputReadBuffer.cpp				
					ScoreBuffer.cpp
					StrippedSW.cpp
					unix.cpp
					unix_threads.cpp)

SET_TARGET_PROPERTIES(ngmlr-debug PROPERTIES COMPILE_DEFINITIONS "DEBUGLOG")
SET_TARGET_PROPERTIES(ngmlr-debug PROPERTIES COMPILE_FLAGS " -g3 -O0")

TARGET_LINK_LIBRARIES(ngmlr-debug pthread)	
TARGET_LINK_LIBRARIES(ngmlr-debug BamTools-static)
TARGET_LINK_LIBRARIES(ngmlr-debug zlibstatic)

