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
include examples/CMakeFiles/vietoris_rips.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/vietoris_rips.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/vietoris_rips.dir/flags.make

examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o: examples/CMakeFiles/vietoris_rips.dir/flags.make
examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o: examples/vietoris_rips.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o -c /home/jens/Uni/data_topology/Aleph/examples/vietoris_rips.cc

examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.i"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Uni/data_topology/Aleph/examples/vietoris_rips.cc > CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.i

examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.s"
	cd /home/jens/Uni/data_topology/Aleph/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Uni/data_topology/Aleph/examples/vietoris_rips.cc -o CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.s

examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.requires:

.PHONY : examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.requires

examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.provides: examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.requires
	$(MAKE) -f examples/CMakeFiles/vietoris_rips.dir/build.make examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.provides.build
.PHONY : examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.provides

examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.provides.build: examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o


# Object files for target vietoris_rips
vietoris_rips_OBJECTS = \
"CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o"

# External object files for target vietoris_rips
vietoris_rips_EXTERNAL_OBJECTS =

examples/vietoris_rips: examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o
examples/vietoris_rips: examples/CMakeFiles/vietoris_rips.dir/build.make
examples/vietoris_rips: examples/CMakeFiles/vietoris_rips.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vietoris_rips"
	cd /home/jens/Uni/data_topology/Aleph/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vietoris_rips.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/vietoris_rips.dir/build: examples/vietoris_rips

.PHONY : examples/CMakeFiles/vietoris_rips.dir/build

examples/CMakeFiles/vietoris_rips.dir/requires: examples/CMakeFiles/vietoris_rips.dir/vietoris_rips.cc.o.requires

.PHONY : examples/CMakeFiles/vietoris_rips.dir/requires

examples/CMakeFiles/vietoris_rips.dir/clean:
	cd /home/jens/Uni/data_topology/Aleph/examples && $(CMAKE_COMMAND) -P CMakeFiles/vietoris_rips.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/vietoris_rips.dir/clean

examples/CMakeFiles/vietoris_rips.dir/depend:
	cd /home/jens/Uni/data_topology/Aleph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/examples /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/examples /home/jens/Uni/data_topology/Aleph/examples/CMakeFiles/vietoris_rips.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/vietoris_rips.dir/depend

