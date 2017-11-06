include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories("${LIBRARY_SOURCES_DIRECTORY}/Common")

add_executable(runBasicTests check_complements.cpp)

target_link_libraries(runBasicTests gtest_main)
target_link_libraries(runBasicTests Common)