# minimum CMake version required
cmake_minimum_required (VERSION 2.8)

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

# sets the output program where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)

# sets the top level directory of the build tree, separating CMake files from our code
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# The location of our source files
set(SOURCES "${PROJECT_SOURCE_DIR}/src")

# The directory all libraries should be stored in
set(LIBRARIES "${PROJECT_SOURCE_DIR}/libraries")

include_directories("${LIBRARIES}/Common")
add_subdirectory(${LIBRARIES}/Common)

# sets the location of binary library files
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
 
# add the executable
add_executable(sample ${SOURCES}/sample.cpp)
target_link_libraries(sample Common)