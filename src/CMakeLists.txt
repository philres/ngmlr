cmake_minimum_required(VERSION 2.8)
project(ngmlr)

include_directories (../lib/seqan-library-2.2.0/include)
include_directories (../lib/Complete-Striped-Smith-Waterman-Library/src)
include_directories (../lib)
include_directories( ${ZLIB_INCLUDE_DIRS} )

configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )

IF(STATIC)
    IF(APPLE)
        message(FATAL_ERROR "MacOSX does not support static linking. Please use -DSTATIC=OFF")
    ENDIF()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")    
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
ENDIF()

add_executable(ngmlr
                    AlignmentMatrix.cpp
                    AlignmentMatrixFast.cpp
                    ArgParser.cpp																					
					CS.cpp					
					CSstatic.cpp
					ConvexAlign.cpp
					ConvexAlignFast.cpp
					LinearRegression.cpp													
					Logging.cpp
					MappedRead.cpp
					main.cpp
					NGM.cpp
					NGMTask.cpp
					AlignmentBuffer.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					SamParser.cpp
					SAMWriter.cpp
					SequenceProvider.cpp	
					OutputReadBuffer.cpp	
					ScoreBuffer.cpp
					StrippedSW.cpp
					../lib/Complete-Striped-Smith-Waterman-Library/src/ssw.c
					unix.cpp
					unix_threads.cpp)										
					
ADD_CUSTOM_TARGET(recompile_always ALL touch "${CMAKE_SOURCE_DIR}/src/main.cpp")

#SET_TARGET_PROPERTIES(ngmlr PROPERTIES COMPILE_FLAGS " -g -O3 -DNDEBUG")
 
TARGET_LINK_LIBRARIES(ngmlr ${ZLIB_LIBRARIES})
TARGET_LINK_LIBRARIES(ngmlr ${CMAKE_THREAD_LIBS_INIT})

install (TARGETS ngmlr DESTINATION bin)