cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project (SuperTrack VERSION 1.0.0 LANGUAGES CXX CUDA)


set(CMAKE_CXX_STANDARD "17")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 ")
#set(CMAKE_CUDA_STANDARD "17")
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++") 

add_library(cuda_test STATIC src/testCUDA.cu include/testCUDA.cuh)
#set_target_properties(cuda_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)


# Finding the ROOT package
find_package(ROOT 6.25 CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
message("ROOT Version: ${ROOT_VERSION}")
#message("${ROOT_LIBARIES}")

#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh )
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cu)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

#Add an executable program and link to ROOT libraries
add_executable(SuperTrack SuperTrack.cu)
target_link_libraries(SuperTrack PUBLIC ${ROOT_LIBRARIES})
#set_property(TARGET SuperimposeTracks PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(SuperTrack PRIVATE cuda_test)


