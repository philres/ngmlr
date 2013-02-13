cmake_minimum_required(VERSION 2.8)

include_directories (../include)
include_directories (../src)
include_directories (../src/parser)
include_directories (../src/writer)
include_directories (../src/misc)
include_directories (../src/log)
include_directories (../src/core)
include_directories (writer)
include_directories (paired)

include_directories (.)

add_executable(ngm-utils main.cpp
								paired/interleave-pairs.cpp
								writer/Writer.cpp
								../src/parser/SamParser.cpp
								../src/log/Logging.cpp
								../src/core/unix.cpp
								../src/core/unix_threads.cpp
								../src/MappedRead.cpp
								
 					)		

TARGET_LINK_LIBRARIES(ngm-utils zlibstatic)
#TARGET_LINK_LIBRARIES(ngm-utils zlib)
#target_link_libraries(ngm-utils z)