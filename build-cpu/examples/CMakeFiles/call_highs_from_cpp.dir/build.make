# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ivet/code/HiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivet/code/HiGHS/build-cpu

# Include any dependencies generated for this target.
include examples/CMakeFiles/call_highs_from_cpp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/CMakeFiles/call_highs_from_cpp.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/call_highs_from_cpp.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/call_highs_from_cpp.dir/flags.make

examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o: examples/CMakeFiles/call_highs_from_cpp.dir/flags.make
examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o: /home/ivet/code/HiGHS/examples/call_highs_from_cpp.cpp
examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o: examples/CMakeFiles/call_highs_from_cpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/ivet/code/HiGHS/build-cpu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o"
	cd /home/ivet/code/HiGHS/build-cpu/examples && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o -MF CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o.d -o CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o -c /home/ivet/code/HiGHS/examples/call_highs_from_cpp.cpp

examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.i"
	cd /home/ivet/code/HiGHS/build-cpu/examples && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivet/code/HiGHS/examples/call_highs_from_cpp.cpp > CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.i

examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.s"
	cd /home/ivet/code/HiGHS/build-cpu/examples && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivet/code/HiGHS/examples/call_highs_from_cpp.cpp -o CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.s

# Object files for target call_highs_from_cpp
call_highs_from_cpp_OBJECTS = \
"CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o"

# External object files for target call_highs_from_cpp
call_highs_from_cpp_EXTERNAL_OBJECTS =

bin/call_highs_from_cpp: examples/CMakeFiles/call_highs_from_cpp.dir/call_highs_from_cpp.cpp.o
bin/call_highs_from_cpp: examples/CMakeFiles/call_highs_from_cpp.dir/build.make
bin/call_highs_from_cpp: lib/libhighs.so.1.8.1
bin/call_highs_from_cpp: /usr/lib/x86_64-linux-gnu/libz.so
bin/call_highs_from_cpp: examples/CMakeFiles/call_highs_from_cpp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/ivet/code/HiGHS/build-cpu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/call_highs_from_cpp"
	cd /home/ivet/code/HiGHS/build-cpu/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/call_highs_from_cpp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/call_highs_from_cpp.dir/build: bin/call_highs_from_cpp
.PHONY : examples/CMakeFiles/call_highs_from_cpp.dir/build

examples/CMakeFiles/call_highs_from_cpp.dir/clean:
	cd /home/ivet/code/HiGHS/build-cpu/examples && $(CMAKE_COMMAND) -P CMakeFiles/call_highs_from_cpp.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/call_highs_from_cpp.dir/clean

examples/CMakeFiles/call_highs_from_cpp.dir/depend:
	cd /home/ivet/code/HiGHS/build-cpu && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivet/code/HiGHS /home/ivet/code/HiGHS/examples /home/ivet/code/HiGHS/build-cpu /home/ivet/code/HiGHS/build-cpu/examples /home/ivet/code/HiGHS/build-cpu/examples/CMakeFiles/call_highs_from_cpp.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : examples/CMakeFiles/call_highs_from_cpp.dir/depend

