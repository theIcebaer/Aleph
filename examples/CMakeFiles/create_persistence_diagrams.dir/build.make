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
include examples/CMakeFiles/create_persistence_diagrams.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/create_persistence_diagrams.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/create_persistence_diagrams.dir/flags.make

examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o: examples/CMakeFiles/create_persistence_diagrams.dir/flags.make
examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o: examples/create_persistence_diagrams.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o -c /home/jens/Uni/data_topology/Aleph/examples/create_persistence_diagrams.cc

examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.i"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Uni/data_topology/Aleph/examples/create_persistence_diagrams.cc > CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.i

examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.s"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Uni/data_topology/Aleph/examples/create_persistence_diagrams.cc -o CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.s

examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.requires:

.PHONY : examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.requires

examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.provides: examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.requires
	$(MAKE) -f examples/CMakeFiles/create_persistence_diagrams.dir/build.make examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.provides.build
.PHONY : examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.provides

examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.provides.build: examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o


# Object files for target create_persistence_diagrams
create_persistence_diagrams_OBJECTS = \
"CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o"

# External object files for target create_persistence_diagrams
create_persistence_diagrams_EXTERNAL_OBJECTS =

examples/create_persistence_diagrams: examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o
examples/create_persistence_diagrams: examples/CMakeFiles/create_persistence_diagrams.dir/build.make
examples/create_persistence_diagrams: examples/CMakeFiles/create_persistence_diagrams.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable create_persistence_diagrams"
	cd /home/jens/Uni/data_topology/Aleph/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/create_persistence_diagrams.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/create_persistence_diagrams.dir/build: examples/create_persistence_diagrams

.PHONY : examples/CMakeFiles/create_persistence_diagrams.dir/build

examples/CMakeFiles/create_persistence_diagrams.dir/requires: examples/CMakeFiles/create_persistence_diagrams.dir/create_persistence_diagrams.cc.o.requires

.PHONY : examples/CMakeFiles/create_persistence_diagrams.dir/requires

examples/CMakeFiles/create_persistence_diagrams.dir/clean:
	cd /home/jens/Uni/data_topology/Aleph/examples && $(CMAKE_COMMAND) -P CMakeFiles/create_persistence_diagrams.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/create_persistence_diagrams.dir/clean

examples/CMakeFiles/create_persistence_diagrams.dir/depend:
	cd /home/jens/Uni/data_topology/Aleph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/examples /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/examples /home/jens/Uni/data_topology/Aleph/examples/CMakeFiles/create_persistence_diagrams.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/create_persistence_diagrams.dir/depend

