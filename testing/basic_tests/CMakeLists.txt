include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runBasicTests check_complements.cpp)

target_link_libraries(runBasicTests gtest_main)
target_link_libraries(runBasicTests Common)