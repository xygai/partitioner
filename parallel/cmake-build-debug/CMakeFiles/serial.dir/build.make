# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/xgai/Final Project/serial"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/xgai/Final Project/serial/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/serial.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/serial.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/serial.dir/flags.make

CMakeFiles/serial.dir/main.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/serial.dir/main.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/main.c.o   -c "/Users/xgai/Final Project/serial/main.c"

CMakeFiles/serial.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/main.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/main.c" > CMakeFiles/serial.dir/main.c.i

CMakeFiles/serial.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/main.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/main.c" -o CMakeFiles/serial.dir/main.c.s

CMakeFiles/serial.dir/graph.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/graph.c.o: ../graph.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/serial.dir/graph.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/graph.c.o   -c "/Users/xgai/Final Project/serial/graph.c"

CMakeFiles/serial.dir/graph.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/graph.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/graph.c" > CMakeFiles/serial.dir/graph.c.i

CMakeFiles/serial.dir/graph.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/graph.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/graph.c" -o CMakeFiles/serial.dir/graph.c.s

CMakeFiles/serial.dir/refine.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/refine.c.o: ../refine.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/serial.dir/refine.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/refine.c.o   -c "/Users/xgai/Final Project/serial/refine.c"

CMakeFiles/serial.dir/refine.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/refine.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/refine.c" > CMakeFiles/serial.dir/refine.c.i

CMakeFiles/serial.dir/refine.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/refine.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/refine.c" -o CMakeFiles/serial.dir/refine.c.s

CMakeFiles/serial.dir/initpart.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/initpart.c.o: ../initpart.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/serial.dir/initpart.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/initpart.c.o   -c "/Users/xgai/Final Project/serial/initpart.c"

CMakeFiles/serial.dir/initpart.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/initpart.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/initpart.c" > CMakeFiles/serial.dir/initpart.c.i

CMakeFiles/serial.dir/initpart.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/initpart.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/initpart.c" -o CMakeFiles/serial.dir/initpart.c.s

CMakeFiles/serial.dir/coarsen.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/coarsen.c.o: ../coarsen.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/serial.dir/coarsen.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/coarsen.c.o   -c "/Users/xgai/Final Project/serial/coarsen.c"

CMakeFiles/serial.dir/coarsen.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/coarsen.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/coarsen.c" > CMakeFiles/serial.dir/coarsen.c.i

CMakeFiles/serial.dir/coarsen.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/coarsen.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/coarsen.c" -o CMakeFiles/serial.dir/coarsen.c.s

CMakeFiles/serial.dir/multilevel.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/multilevel.c.o: ../multilevel.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/serial.dir/multilevel.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/multilevel.c.o   -c "/Users/xgai/Final Project/serial/multilevel.c"

CMakeFiles/serial.dir/multilevel.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/multilevel.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/multilevel.c" > CMakeFiles/serial.dir/multilevel.c.i

CMakeFiles/serial.dir/multilevel.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/multilevel.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/multilevel.c" -o CMakeFiles/serial.dir/multilevel.c.s

CMakeFiles/serial.dir/utils.c.o: CMakeFiles/serial.dir/flags.make
CMakeFiles/serial.dir/utils.c.o: ../utils.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/serial.dir/utils.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/serial.dir/utils.c.o   -c "/Users/xgai/Final Project/serial/utils.c"

CMakeFiles/serial.dir/utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/serial.dir/utils.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/xgai/Final Project/serial/utils.c" > CMakeFiles/serial.dir/utils.c.i

CMakeFiles/serial.dir/utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/serial.dir/utils.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/xgai/Final Project/serial/utils.c" -o CMakeFiles/serial.dir/utils.c.s

# Object files for target serial
serial_OBJECTS = \
"CMakeFiles/serial.dir/main.c.o" \
"CMakeFiles/serial.dir/graph.c.o" \
"CMakeFiles/serial.dir/refine.c.o" \
"CMakeFiles/serial.dir/initpart.c.o" \
"CMakeFiles/serial.dir/coarsen.c.o" \
"CMakeFiles/serial.dir/multilevel.c.o" \
"CMakeFiles/serial.dir/utils.c.o"

# External object files for target serial
serial_EXTERNAL_OBJECTS =

serial: CMakeFiles/serial.dir/main.c.o
serial: CMakeFiles/serial.dir/graph.c.o
serial: CMakeFiles/serial.dir/refine.c.o
serial: CMakeFiles/serial.dir/initpart.c.o
serial: CMakeFiles/serial.dir/coarsen.c.o
serial: CMakeFiles/serial.dir/multilevel.c.o
serial: CMakeFiles/serial.dir/utils.c.o
serial: CMakeFiles/serial.dir/build.make
serial: CMakeFiles/serial.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_8) "Linking C executable serial"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/serial.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/serial.dir/build: serial

.PHONY : CMakeFiles/serial.dir/build

CMakeFiles/serial.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/serial.dir/cmake_clean.cmake
.PHONY : CMakeFiles/serial.dir/clean

CMakeFiles/serial.dir/depend:
	cd "/Users/xgai/Final Project/serial/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/xgai/Final Project/serial" "/Users/xgai/Final Project/serial" "/Users/xgai/Final Project/serial/cmake-build-debug" "/Users/xgai/Final Project/serial/cmake-build-debug" "/Users/xgai/Final Project/serial/cmake-build-debug/CMakeFiles/serial.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/serial.dir/depend

