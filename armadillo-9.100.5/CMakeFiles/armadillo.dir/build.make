# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5

# Include any dependencies generated for this target.
include CMakeFiles/armadillo.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/armadillo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/armadillo.dir/flags.make

CMakeFiles/armadillo.dir/src/wrapper.cpp.o: CMakeFiles/armadillo.dir/flags.make
CMakeFiles/armadillo.dir/src/wrapper.cpp.o: src/wrapper.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/armadillo.dir/src/wrapper.cpp.o"
	/ssoft/spack/paien/v2/opt/spack/linux-rhel7-x86_E5v2_IntelIB/gcc-4.8.5/gcc-6.4.0-5lz2sine3ujgoc7aswaexuinafuepktx/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/armadillo.dir/src/wrapper.cpp.o -c /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5/src/wrapper.cpp

CMakeFiles/armadillo.dir/src/wrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armadillo.dir/src/wrapper.cpp.i"
	/ssoft/spack/paien/v2/opt/spack/linux-rhel7-x86_E5v2_IntelIB/gcc-4.8.5/gcc-6.4.0-5lz2sine3ujgoc7aswaexuinafuepktx/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5/src/wrapper.cpp > CMakeFiles/armadillo.dir/src/wrapper.cpp.i

CMakeFiles/armadillo.dir/src/wrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armadillo.dir/src/wrapper.cpp.s"
	/ssoft/spack/paien/v2/opt/spack/linux-rhel7-x86_E5v2_IntelIB/gcc-4.8.5/gcc-6.4.0-5lz2sine3ujgoc7aswaexuinafuepktx/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5/src/wrapper.cpp -o CMakeFiles/armadillo.dir/src/wrapper.cpp.s

CMakeFiles/armadillo.dir/src/wrapper.cpp.o.requires:
.PHONY : CMakeFiles/armadillo.dir/src/wrapper.cpp.o.requires

CMakeFiles/armadillo.dir/src/wrapper.cpp.o.provides: CMakeFiles/armadillo.dir/src/wrapper.cpp.o.requires
	$(MAKE) -f CMakeFiles/armadillo.dir/build.make CMakeFiles/armadillo.dir/src/wrapper.cpp.o.provides.build
.PHONY : CMakeFiles/armadillo.dir/src/wrapper.cpp.o.provides

CMakeFiles/armadillo.dir/src/wrapper.cpp.o.provides.build: CMakeFiles/armadillo.dir/src/wrapper.cpp.o

# Object files for target armadillo
armadillo_OBJECTS = \
"CMakeFiles/armadillo.dir/src/wrapper.cpp.o"

# External object files for target armadillo
armadillo_EXTERNAL_OBJECTS =

libarmadillo.so.9.100.5: CMakeFiles/armadillo.dir/src/wrapper.cpp.o
libarmadillo.so.9.100.5: CMakeFiles/armadillo.dir/build.make
libarmadillo.so.9.100.5: /ssoft/spack/paien/v2/opt/spack/linux-rhel7-x86_E5v2_IntelIB/gcc-6.4.0/openblas-0.2.20-gjpsmh7vkmgwtqrodsfzsdl324vpv6y7/lib/libopenblas.so
libarmadillo.so.9.100.5: CMakeFiles/armadillo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library libarmadillo.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/armadillo.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library libarmadillo.so.9.100.5 libarmadillo.so.9 libarmadillo.so

libarmadillo.so.9: libarmadillo.so.9.100.5

libarmadillo.so: libarmadillo.so.9.100.5

# Rule to build all files generated by this target.
CMakeFiles/armadillo.dir/build: libarmadillo.so
.PHONY : CMakeFiles/armadillo.dir/build

CMakeFiles/armadillo.dir/requires: CMakeFiles/armadillo.dir/src/wrapper.cpp.o.requires
.PHONY : CMakeFiles/armadillo.dir/requires

CMakeFiles/armadillo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/armadillo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/armadillo.dir/clean

CMakeFiles/armadillo.dir/depend:
	cd /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5 /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5 /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5 /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5 /work/lcvmm/crea/cgFrames_SC/armadillo-9.100.5/CMakeFiles/armadillo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/armadillo.dir/depend

