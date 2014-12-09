cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(PolyMaker)
enable_language(CXX)

find_package(Eigen3)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_PATH})
    add_definitions(-DENABLE_EIGEN3)
else()
    message(FATAL_ERROR "-- Could not find any EIGEN3 library")
endif()

add_executable(poly main.c)
