#not recommended but, easy enough
file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.hh *.cc *.cu)

add_executable(SuperTrack_tst ${TEST_SOURCES})

target_link_libraries(SuperTrack_tst PUBLIC src gtest)