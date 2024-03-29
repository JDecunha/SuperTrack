#CMakeLists file for SuperTrack

cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project (SuperTrack VERSION 1.0.0 LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD "17")
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G") #have host and device symbols in debug code

###
### Find and include ROOT libraries
###

find_package(ROOT 6.25 CONFIG REQUIRED)
include("${ROOT_USE_FILE}") #include the ROOT CMakelists which defines $ROOT_LIBRARIES

###
### Add the .ini reader header only library
###

#interface because header only
add_library(inih INTERFACE)
target_include_directories(inih INTERFACE ${CMAKE_SOURCE_DIR}/extern/inih)

###
### Add googletest
###

add_subdirectory(extern/googletest)

###
### Configure src
###

add_executable(SuperTrack SuperTrack.cu)

#locate our project source files and make a library
file(GLOB_RECURSE sources ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/*.cu)
add_library(src STATIC ${sources})

#public below means the main SuperTrack.cu inherits the include directories of src
#hard coding the include directories is sloppy, but simple
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include/Data_output)
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include/GPU_helper_types)
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include/Simulation_method)
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include/Simulation_method/VoxelConstrainedSphereMethod)
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include/Thread_allocation)

#Link src to the required libraries
target_link_libraries(src PUBLIC ${ROOT_LIBRARIES} inih -lcurand)
set_target_properties(src PROPERTIES CUDA_SEPARABLE_COMPILATION ON)


#Link my executable to src, which itself is linked to all of the required libraries
target_link_libraries(SuperTrack PRIVATE src)

###
### Add tests
###

add_subdirectory(tst)






