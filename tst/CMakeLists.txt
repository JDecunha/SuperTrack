#Set up src_tst
file(GLOB_RECURSE src_tst ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cu)

#Set up executable
add_executable(SuperTrack_tst ${src_tst})
target_link_libraries(SuperTrack_tst PUBLIC gtest_main src)