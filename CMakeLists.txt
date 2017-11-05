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

 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h

include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# The location of our source files
set(SOURCES "${PROJECT_SOURCE_DIR}/src")

# The directory all libraries should be stored in
set(LIBRARY_SOURCES_DIRECTORY "${PROJECT_SOURCE_DIR}/libraries")

SET( CMAKE_CXX_FLAGS  "-Iinclude/ -std=c++11" )

# my simple test library
include_directories("${LIBRARY_SOURCES_DIRECTORY}/Common")
add_subdirectory(${LIBRARY_SOURCES_DIRECTORY}/Common)
 
# including bioparser
add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

# add the executable
add_executable(sample ${SOURCES}/sample.cpp)
target_link_libraries(sample Common)

add_executable(vectors ${SOURCES}/vectorSample.cpp)

# 
add_executable(fastaSample ${SOURCES}/fastaSample.cpp)
target_link_libraries(fastaSample bioparser)