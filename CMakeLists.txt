cmake_minimum_required(VERSION 2.8.0)
project(SRSLI)

# Set this as a Clang-compiled C++ project
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

# Set this as a C++11 project
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Find and include the SeqAn library of templates
list(APPEND CMAKE_MODULE_PATH ./cmake)
find_package(SeqAn 1.4.1 REQUIRED)

add_subdirectory(src/C++)
