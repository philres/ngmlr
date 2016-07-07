cmake_minimum_required(VERSION 2.4.4)

project(minimap)

set(VERSION "0.2-r124-dirty")

#SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -Wall -O2 -Wc++-compat -Wno-unused-function")

add_library(minimapstatic STATIC bseq.c
index.c
kthread.c
map.c
misc.c
sdust.c
sketch.c)


