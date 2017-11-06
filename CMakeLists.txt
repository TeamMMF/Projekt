# minimum CMake version required
cmake_minimum_required (VERSION 3.2)

# The project's name
project (SequenceOverlaping)

# The project's version number.
set (SequenceOverlaping_VERSION_MAJOR 1)
set (SequenceOverlaping_VERSION_MINOR 0)
 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/ProjectConfig.h.in"
  "${PROJECT_BINARY_DIR}/ProjectConfig.h"
  )

# ----------------- VARIABLES USED----------------------------#
# CMAKE_SOURCE_DIR - the directory from which CMake was started
# PROJECT_NAME - the name of the project ste by the project() command 
# PROJECT_SOURCE_DIR - the full path to the project source directory root
# Set the output program where your program will be created
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# The location of our source files
set(SOURCES "${PROJECT_SOURCE_DIR}/src")

# The directory all libraries should be stored in
set(LIBRARY_SOURCES_DIRECTORY "${PROJECT_SOURCE_DIR}/libraries")

set( CMAKE_CXX_FLAGS  "-std=c++11" )

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${LIBRARY_SOURCES_DIRECTORY}/Common")
add_subdirectory(testing)
add_subdirectory(${LIBRARY_SOURCES_DIRECTORY}/Common)

# including bioparser
add_subdirectory(libraries/bioparser EXCLUDE_FROM_ALL)

# add the executable
add_executable(sample ${SOURCES}/sample.cpp)
target_link_libraries(sample Common)

add_executable(vectors ${SOURCES}/vectorSample.cpp)

add_executable(fastaSample ${SOURCES}/fastaSample.cpp )
target_link_libraries(fastaSample bioparser)