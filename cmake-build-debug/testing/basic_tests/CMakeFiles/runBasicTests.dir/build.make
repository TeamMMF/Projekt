# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2017.2.3\bin\cmake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2017.2.3\bin\cmake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\cmake-build-debug"

# Include any dependencies generated for this target.
include testing/basic_tests/CMakeFiles/runBasicTests.dir/depend.make

# Include the progress variables for this target.
include testing/basic_tests/CMakeFiles/runBasicTests.dir/progress.make

# Include the compile flags for this target's objects.
include testing/basic_tests/CMakeFiles/runBasicTests.dir/flags.make

testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj: testing/basic_tests/CMakeFiles/runBasicTests.dir/flags.make
testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj: testing/basic_tests/CMakeFiles/runBasicTests.dir/includes_CXX.rsp
testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj: ../testing/basic_tests/check_complements.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj"
	cd /d C:\Faks\55956~1.SEM\PROJEK~1\Projekt\CMAKE-~1\testing\BASIC_~1 && C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\runBasicTests.dir\check_complements.cpp.obj -c "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\testing\basic_tests\check_complements.cpp"

testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runBasicTests.dir/check_complements.cpp.i"
	cd /d C:\Faks\55956~1.SEM\PROJEK~1\Projekt\CMAKE-~1\testing\BASIC_~1 && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\testing\basic_tests\check_complements.cpp" > CMakeFiles\runBasicTests.dir\check_complements.cpp.i

testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runBasicTests.dir/check_complements.cpp.s"
	cd /d C:\Faks\55956~1.SEM\PROJEK~1\Projekt\CMAKE-~1\testing\BASIC_~1 && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\testing\basic_tests\check_complements.cpp" -o CMakeFiles\runBasicTests.dir\check_complements.cpp.s

testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.requires:

.PHONY : testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.requires

testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.provides: testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.requires
	$(MAKE) -f testing\basic_tests\CMakeFiles\runBasicTests.dir\build.make testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.provides.build
.PHONY : testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.provides

testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.provides.build: testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj


# Object files for target runBasicTests
runBasicTests_OBJECTS = \
"CMakeFiles/runBasicTests.dir/check_complements.cpp.obj"

# External object files for target runBasicTests
runBasicTests_EXTERNAL_OBJECTS =

bin/runBasicTests.exe: testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj
bin/runBasicTests.exe: testing/basic_tests/CMakeFiles/runBasicTests.dir/build.make
bin/runBasicTests.exe: lib/libgtest_maind.a
bin/runBasicTests.exe: lib/libCommond.a
bin/runBasicTests.exe: lib/libgtestd.a
bin/runBasicTests.exe: testing/basic_tests/CMakeFiles/runBasicTests.dir/linklibs.rsp
bin/runBasicTests.exe: testing/basic_tests/CMakeFiles/runBasicTests.dir/objects1.rsp
bin/runBasicTests.exe: testing/basic_tests/CMakeFiles/runBasicTests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ..\..\bin\runBasicTests.exe"
	cd /d C:\Faks\55956~1.SEM\PROJEK~1\Projekt\CMAKE-~1\testing\BASIC_~1 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\runBasicTests.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
testing/basic_tests/CMakeFiles/runBasicTests.dir/build: bin/runBasicTests.exe

.PHONY : testing/basic_tests/CMakeFiles/runBasicTests.dir/build

testing/basic_tests/CMakeFiles/runBasicTests.dir/requires: testing/basic_tests/CMakeFiles/runBasicTests.dir/check_complements.cpp.obj.requires

.PHONY : testing/basic_tests/CMakeFiles/runBasicTests.dir/requires

testing/basic_tests/CMakeFiles/runBasicTests.dir/clean:
	cd /d C:\Faks\55956~1.SEM\PROJEK~1\Projekt\CMAKE-~1\testing\BASIC_~1 && $(CMAKE_COMMAND) -P CMakeFiles\runBasicTests.dir\cmake_clean.cmake
.PHONY : testing/basic_tests/CMakeFiles/runBasicTests.dir/clean

testing/basic_tests/CMakeFiles/runBasicTests.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt" "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\testing\basic_tests" "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\cmake-build-debug" "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\cmake-build-debug\testing\basic_tests" "C:\Faks\5. semestar\Projekt iz programske potpore\Projekt\cmake-build-debug\testing\basic_tests\CMakeFiles\runBasicTests.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : testing/basic_tests/CMakeFiles/runBasicTests.dir/depend
