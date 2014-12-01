cmake_minimum_required(VERSION 2.8)

include_directories (../include)
include_directories (../src)
include_directories (../lib/bamtools-2.3.0/src)
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
								reformat_fasta/reformat_fasta.cpp
								cout_reads/cout-reads.cpp
								kmers/kmer-distribution.cpp
								paired/interleave-pairs.cpp
								writer/Writer.cpp								
								../src/parser/SamParser.cpp
								../src/parser/BamParser.cpp
								../src/log/Logging.cpp
								../src/core/unix.cpp
								../src/core/unix_threads.cpp
								../src/MappedRead.cpp
								../src/SequenceProvider.cpp
								../src/config/Config.cpp
								../src/PrefixTable.cpp
								../src/CSstatic.cpp
 					)		
 
 add_executable(ngm-utils-debug main.cpp
								reformat_fasta/reformat_fasta.cpp
								cout_reads/cout-reads.cpp
								kmers/kmer-distribution.cpp
								paired/interleave-pairs.cpp
								writer/Writer.cpp								
								../src/parser/SamParser.cpp
								../src/parser/BamParser.cpp
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

SET_TARGET_PROPERTIES(ngm-utils-debug PROPERTIES COMPILE_FLAGS "-g3 -O0")

add_definitions(-D _BAM)

TARGET_LINK_LIBRARIES(ngm-utils pthread)
TARGET_LINK_LIBRARIES(ngm-utils-debug pthread)

TARGET_LINK_LIBRARIES(ngm-utils BamTools-static)
TARGET_LINK_LIBRARIES(ngm-utils-debug BamTools-static)

TARGET_LINK_LIBRARIES(ngm-utils zlibstatic)
TARGET_LINK_LIBRARIES(ngm-utils-debug zlibstatic)
