#CMakeLists file for SuperTrack

cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project (SuperTrack VERSION 1.0.0 LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD "17")

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
### Configure executable and in-project libraries
###

add_executable(SuperTrack SuperTrack.cu)

#locate our project source files and make a library
	file(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/*.cu)
	add_library(src STATIC ${sources})

	#public below means the main SuperTrack.cu inherits the include directories of src
	target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include)
	set_target_properties(src PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

#Now link my libraries
	target_link_libraries(SuperTrack PRIVATE src ${ROOT_LIBRARIES} inih -lcurand)


