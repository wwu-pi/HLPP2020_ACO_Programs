# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /home/schredder/Programs/clion-2020.1.1/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/schredder/Programs/clion-2020.1.1/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/schredder/Research/Musket/src-gen/BPP/CUDA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/schredder/Research/Musket/src-gen/BPP/CUDA

# Include any dependencies generated for this target.
include CMakeFiles/BPP_0.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/BPP_0.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/BPP_0.dir/flags.make

CMakeFiles/BPP_0.dir/src/BPP_0.cu.o: CMakeFiles/BPP_0.dir/flags.make
CMakeFiles/BPP_0.dir/src/BPP_0.cu.o: src/BPP_0.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/schredder/Research/Musket/src-gen/BPP/CUDA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/BPP_0.dir/src/BPP_0.cu.o"
	/usr/local/cuda-10.2/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/schredder/Research/Musket/src-gen/BPP/CUDA/src/BPP_0.cu -o CMakeFiles/BPP_0.dir/src/BPP_0.cu.o

CMakeFiles/BPP_0.dir/src/BPP_0.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/BPP_0.dir/src/BPP_0.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/BPP_0.dir/src/BPP_0.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/BPP_0.dir/src/BPP_0.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target BPP_0
BPP_0_OBJECTS = \
"CMakeFiles/BPP_0.dir/src/BPP_0.cu.o"

# External object files for target BPP_0
BPP_0_EXTERNAL_OBJECTS =

bin/BPP_0: CMakeFiles/BPP_0.dir/src/BPP_0.cu.o
bin/BPP_0: CMakeFiles/BPP_0.dir/build.make
bin/BPP_0: CMakeFiles/BPP_0.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/schredder/Research/Musket/src-gen/BPP/CUDA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA executable bin/BPP_0"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BPP_0.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/BPP_0.dir/build: bin/BPP_0

.PHONY : CMakeFiles/BPP_0.dir/build

CMakeFiles/BPP_0.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/BPP_0.dir/cmake_clean.cmake
.PHONY : CMakeFiles/BPP_0.dir/clean

CMakeFiles/BPP_0.dir/depend:
	cd /home/schredder/Research/Musket/src-gen/BPP/CUDA && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/schredder/Research/Musket/src-gen/BPP/CUDA /home/schredder/Research/Musket/src-gen/BPP/CUDA /home/schredder/Research/Musket/src-gen/BPP/CUDA /home/schredder/Research/Musket/src-gen/BPP/CUDA /home/schredder/Research/Musket/src-gen/BPP/CUDA/CMakeFiles/BPP_0.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/BPP_0.dir/depend

