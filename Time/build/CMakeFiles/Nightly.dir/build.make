# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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
CMAKE_COMMAND = /afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin/cmake

# The command to remove a file.
RM = /afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build

# Utility rule file for Nightly.

CMakeFiles/Nightly:
	/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin/ctest -D Nightly

Nightly: CMakeFiles/Nightly
Nightly: CMakeFiles/Nightly.dir/build.make
.PHONY : Nightly

# Rule to build all files generated by this target.
CMakeFiles/Nightly.dir/build: Nightly
.PHONY : CMakeFiles/Nightly.dir/build

CMakeFiles/Nightly.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Nightly.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Nightly.dir/clean

CMakeFiles/Nightly.dir/depend:
	cd /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build/CMakeFiles/Nightly.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Nightly.dir/depend

