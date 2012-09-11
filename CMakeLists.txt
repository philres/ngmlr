cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
project(NGM)

#add_definitions(-D INSTANCE_COUNTING) 
#add_definitions(-D VERBOSE)
#add_definitions(-D _DEBUGCS)

include(FindOpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${PROJECT_SOURCE_DIR}/lib/)

ADD_DEFINITIONS(-g -O3 -Wall)

SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN'" )
#SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN/mason:$ORIGIN/amd-app/lib/x86_64'" )

#SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN/mason:$ORIGIN/amd-app/lib/x86_64'" )

SET (CMAKE_CXX_FLAGS " -g -Wall")
SET (CMAKE_CXX_FLAGS_DEBUG "-D _DEBUG -g3 -Wall")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/)

add_executable(NGM 	AlignmentDispatcher.cpp
					BamParser.cpp
					BAMWriter.cpp					
					Buffer.cpp										
					Config.cpp					
					CS.cpp					
					Debug.cpp					
					FFormatWriter.cpp					
					Logging.cpp
					MappedRead.cpp
					NGM_main.cpp
					NGM.cpp
					NGMStats.cpp
					NGMTask.cpp
					Output.cpp
					PrefixTable.cpp														
					ReadProvider.cpp
					SamParser.cpp
					SAMWriter.cpp
					SequenceProvider.cpp
					SW.cpp
					unix.cpp
					unix_threads.cpp
					windows_threads.cpp
					windows.cpp)

#FIND_PACKAGE(ZLIB REQUIRED)

FIND_PACKAGE(Bamtools QUIET COMPONENTS)
if( ${BAMTOOLS_FOUND} )		
	ADD_LIBRARY(bamtools STATIC IMPORTED)
	SET_TARGET_PROPERTIES(bamtools PROPERTIES IMPORTED_LOCATION ${BAMTOOLS_LIBRARY})
	include_directories(${BAMTOOLS_INCLUDE_DIRS})
	TARGET_LINK_LIBRARIES(NGM bamtools)
	
	add_definitions(-D _BAM) 
endif()

#find_library(HOARD hoard)
#if(HOARD) 
#	message(" --Found libhoard: ${HOARD}")	
#	target_link_libraries(NGM ${HOARD})
#else()
#	message(" --Libhoard not found!")	
#endif()

#SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mtune=amdfam10")
#SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=amdfam10")


	#STATIC
    #set(CMAKE_EXE_LINKER_FLAGS -static) 
    #set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 
    #set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic 
    #set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS) 
    #set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC 
    #set(CMAKE_SHARED_LIBRARY_CXX_FLAGS) 
    #set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic 
    #set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS) 
    # Maybe this works as well, haven't tried yet. 
	# set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE) 
	#target_link_libraries(NGM /scratch/philipp/src/zlib-1.2.7/zlib/lib/libz.a)
	#target_link_libraries(NGM /software/ngm/src/bamtools/lib/libbamtools.a)

target_link_libraries(NGM z)					
target_link_libraries(NGM dl)
target_link_libraries(NGM pthread)			

