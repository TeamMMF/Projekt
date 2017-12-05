include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories("${LIBRARY_SOURCES_DIRECTORY}/Dynamic")

add_executable(runTests LCSTests.cpp)

target_link_libraries(runTests gtest_main)
target_link_libraries(runTests Dynamic)