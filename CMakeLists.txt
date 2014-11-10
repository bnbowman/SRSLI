cmake_minimum_required(VERSION 2.8.0)
project(SRSLI)

list(APPEND CMAKE_MODULE_PATH .)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/local/lib/seqan-1.4.1/util/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(SeqAn 1.41 REQUIRED)

message( STATUS "SEQAN_INCLUDE_DIRS:     " ${SEQAN_INCLUDE_DIRS} )
include_directories(${CMAKE_BINARY_DIR}/src)

add_subdirectory(src)