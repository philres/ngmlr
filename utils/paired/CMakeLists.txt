cmake_minimum_required(VERSION 2.8)

include_directories (../../include)
include_directories (../../src)
include_directories (../../src/parser)
include_directories (../../src/writer)
include_directories (../../src/misc)
include_directories (../../src/log)
include_directories (../../src/core)
#include_directories (../../src/config)

include_directories (.)

add_executable(ngm-join-paired main.cpp
								Writer.cpp
								../../src/parser/SamParser.cpp
								../../src/log/Logging.cpp
								../../src/core/unix.cpp
								../../src/core/unix_threads.cpp
								../../src/MappedRead.cpp
 					)		

target_link_libraries(ngm-join-paired z)