cmake_minimum_required(VERSION 2.8)

#include_directories (${CMAKE_BINARY_DIR}/lib/libStatGen_1.0.9/include)
include_directories (../include)
include_directories (../src)
include_directories (../src/parser)
include_directories (../src/writer)
include_directories (../src/misc)
include_directories (../src/log)
include_directories (../src/core)
include_directories (../src/config)
include_directories (writer)
include_directories (paired)


include_directories (.)

add_executable(ngm-utils main.cpp
								kmers/kmer-distribution.cpp
								paired/interleave-pairs.cpp
								writer/Writer.cpp								
								../src/parser/SamParser.cpp
								../src/log/Logging.cpp
								../src/core/unix.cpp
								../src/core/unix_threads.cpp
								../src/MappedRead.cpp
								../src/SequenceProvider.cpp
								../src/config/Config.cpp
								../src/PrefixTable.cpp
								../src/CSstatic.cpp
 					)		

SET_TARGET_PROPERTIES(ngm-utils PROPERTIES COMPILE_FLAGS "-g -O3 -DNDEBUG")

TARGET_LINK_LIBRARIES(ngm-utils zlibstatic)
