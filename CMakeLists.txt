# CMakeLists for the tagged_dis_he4 generator repository

# Set up
cmake_minimum_required (VERSION 3.0)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (tagged_dis_he4)

# Includes
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Generator executable
add_executable(tagged_dis_gen tagged_dis_gen.cpp)
target_link_libraries(tagged_dis_gen ${ROOT_LIBRARIES})
