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
include src/tools/CMakeFiles/wicked_triangulations.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/wicked_triangulations.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/wicked_triangulations.dir/flags.make

src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o: src/tools/CMakeFiles/wicked_triangulations.dir/flags.make
src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o: src/tools/wicked_triangulations.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o"
	cd /home/jens/Uni/data_topology/Aleph/src/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o -c /home/jens/Uni/data_topology/Aleph/src/tools/wicked_triangulations.cc

src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.i"
	cd /home/jens/Uni/data_topology/Aleph/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Uni/data_topology/Aleph/src/tools/wicked_triangulations.cc > CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.i

src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.s"
	cd /home/jens/Uni/data_topology/Aleph/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Uni/data_topology/Aleph/src/tools/wicked_triangulations.cc -o CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.s

src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.requires:

.PHONY : src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.requires

src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.provides: src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.requires
	$(MAKE) -f src/tools/CMakeFiles/wicked_triangulations.dir/build.make src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.provides.build
.PHONY : src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.provides

src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.provides.build: src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o


# Object files for target wicked_triangulations
wicked_triangulations_OBJECTS = \
"CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o"

# External object files for target wicked_triangulations
wicked_triangulations_EXTERNAL_OBJECTS =

src/tools/wicked_triangulations: src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o
src/tools/wicked_triangulations: src/tools/CMakeFiles/wicked_triangulations.dir/build.make
src/tools/wicked_triangulations: src/tools/CMakeFiles/wicked_triangulations.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jens/Uni/data_topology/Aleph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable wicked_triangulations"
	cd /home/jens/Uni/data_topology/Aleph/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wicked_triangulations.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/wicked_triangulations.dir/build: src/tools/wicked_triangulations

.PHONY : src/tools/CMakeFiles/wicked_triangulations.dir/build

src/tools/CMakeFiles/wicked_triangulations.dir/requires: src/tools/CMakeFiles/wicked_triangulations.dir/wicked_triangulations.cc.o.requires

.PHONY : src/tools/CMakeFiles/wicked_triangulations.dir/requires

src/tools/CMakeFiles/wicked_triangulations.dir/clean:
	cd /home/jens/Uni/data_topology/Aleph/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/wicked_triangulations.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/wicked_triangulations.dir/clean

src/tools/CMakeFiles/wicked_triangulations.dir/depend:
	cd /home/jens/Uni/data_topology/Aleph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/src/tools /home/jens/Uni/data_topology/Aleph /home/jens/Uni/data_topology/Aleph/src/tools /home/jens/Uni/data_topology/Aleph/src/tools/CMakeFiles/wicked_triangulations.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/wicked_triangulations.dir/depend

