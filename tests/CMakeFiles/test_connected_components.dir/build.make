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
include tests/CMakeFiles/test_connected_components.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/test_connected_components.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/test_connected_components.dir/flags.make

tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o: tests/CMakeFiles/test_connected_components.dir/flags.make
tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o: tests/test_connected_components.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o"
	cd /home/jens/Uni/data_topology/Aleph/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_connected_components.dir/test_connected_components.cc.o -c /home/jens/Uni/data_topology/Aleph/tests/test_connected_components.cc

tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_connected_components.dir/test_connected_components.cc.i"
	cd /home/jens/Uni/data_topology/Aleph/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Uni/data_topology/Aleph/tests/test_connected_components.cc > CMakeFiles/test_connected_components.dir/test_connected_components.cc.i

tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_connected_components.dir/test_connected_components.cc.s"
	cd /home/jens/Uni/data_topology/Aleph/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Uni/data_topology/Aleph/tests/test_connected_components.cc -o CMakeFiles/test_connected_components.dir/test_connected_components.cc.s

tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.requires:

.PHONY : tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.requires

tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.provides: tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.requires
	$(MAKE) -f tests/CMakeFiles/test_connected_components.dir/build.make tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.provides.build
.PHONY : tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.provides

tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.provides.build: tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o


# Object files for target test_connected_components
test_connected_components_OBJECTS = \
"CMakeFiles/test_connected_components.dir/test_connected_components.cc.o"

# External object files for target test_connected_components
test_connected_components_EXTERNAL_OBJECTS =

tests/test_connected_components: tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o
tests/test_connected_components: tests/CMakeFiles/test_connected_components.dir/build.make
tests/test_connected_components: tests/CMakeFiles/test_connected_components.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_connected_components"
	cd /home/jens/Uni/data_topology/Aleph/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_connected_components.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/test_connected_components.dir/build: tests/test_connected_components

.PHONY : tests/CMakeFiles/test_connected_components.dir/build

tests/CMakeFiles/test_connected_components.dir/requires: tests/CMakeFiles/test_connected_components.dir/test_connected_components.cc.o.requires

.PHONY : tests/CMakeFiles/test_connected_components.dir/requires

tests/CMakeFiles/test_connected_components.dir/clean:
	cd /home/jens/Uni/data_topology/Aleph/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_connected_components.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/test_connected_components.dir/clean

tests/CMakeFiles/test_connected_components.dir/depend:
	cd /home/jens/Uni/data_topology/Aleph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/tests /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/tests /home/jens/Uni/data_topology/Aleph/tests/CMakeFiles/test_connected_components.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/test_connected_components.dir/depend
