# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jens/Uni/data_topology/Aleph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jens/Uni/data_topology/Aleph

# Include any dependencies generated for this target.
include examples/CMakeFiles/1D.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/1D.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/1D.dir/flags.make

examples/CMakeFiles/1D.dir/1D.cc.o: examples/CMakeFiles/1D.dir/flags.make
examples/CMakeFiles/1D.dir/1D.cc.o: examples/1D.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/1D.dir/1D.cc.o"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/1D.dir/1D.cc.o -c /home/jens/Uni/data_topology/Aleph/examples/1D.cc

examples/CMakeFiles/1D.dir/1D.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/1D.dir/1D.cc.i"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Uni/data_topology/Aleph/examples/1D.cc > CMakeFiles/1D.dir/1D.cc.i

examples/CMakeFiles/1D.dir/1D.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/1D.dir/1D.cc.s"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Uni/data_topology/Aleph/examples/1D.cc -o CMakeFiles/1D.dir/1D.cc.s

examples/CMakeFiles/1D.dir/1D.cc.o.requires:

.PHONY : examples/CMakeFiles/1D.dir/1D.cc.o.requires

examples/CMakeFiles/1D.dir/1D.cc.o.provides: examples/CMakeFiles/1D.dir/1D.cc.o.requires
	$(MAKE) -f examples/CMakeFiles/1D.dir/build.make examples/CMakeFiles/1D.dir/1D.cc.o.provides.build
.PHONY : examples/CMakeFiles/1D.dir/1D.cc.o.provides

examples/CMakeFiles/1D.dir/1D.cc.o.provides.build: examples/CMakeFiles/1D.dir/1D.cc.o


# Object files for target 1D
1D_OBJECTS = \
"CMakeFiles/1D.dir/1D.cc.o"

# External object files for target 1D
1D_EXTERNAL_OBJECTS =

examples/1D: examples/CMakeFiles/1D.dir/1D.cc.o
examples/1D: examples/CMakeFiles/1D.dir/build.make
examples/1D: examples/CMakeFiles/1D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 1D"
	cd /home/jens/Uni/data_topology/Aleph/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/1D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/1D.dir/build: examples/1D

.PHONY : examples/CMakeFiles/1D.dir/build

examples/CMakeFiles/1D.dir/requires: examples/CMakeFiles/1D.dir/1D.cc.o.requires

.PHONY : examples/CMakeFiles/1D.dir/requires

examples/CMakeFiles/1D.dir/clean:
	cd /home/jens/Uni/data_topology/Aleph/examples && $(CMAKE_COMMAND) -P CMakeFiles/1D.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/1D.dir/clean

examples/CMakeFiles/1D.dir/depend:
	cd /home/jens/Uni/data_topology/Aleph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/examples /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/examples /home/jens/Uni/data_topology/Aleph/examples/CMakeFiles/1D.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/1D.dir/depend

