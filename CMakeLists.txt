cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(PolyMaker)
enable_language(CXX)

add_executable(poly main.cpp)

find_package(Eigen3 3.0.0)

if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "-- Could not find any EIGEN3 library")
endif()

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall -lm")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

