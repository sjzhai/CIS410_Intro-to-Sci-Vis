# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/shengjiezhai/Desktop/project8

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/shengjiezhai/Desktop/project8

# Include any dependencies generated for this target.
include CMakeFiles/proj8.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/proj8.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/proj8.dir/flags.make

CMakeFiles/proj8.dir/proj8.cxx.o: CMakeFiles/proj8.dir/flags.make
CMakeFiles/proj8.dir/proj8.cxx.o: proj8.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/shengjiezhai/Desktop/project8/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/proj8.dir/proj8.cxx.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/proj8.dir/proj8.cxx.o -c /Users/shengjiezhai/Desktop/project8/proj8.cxx

CMakeFiles/proj8.dir/proj8.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proj8.dir/proj8.cxx.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/shengjiezhai/Desktop/project8/proj8.cxx > CMakeFiles/proj8.dir/proj8.cxx.i

CMakeFiles/proj8.dir/proj8.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proj8.dir/proj8.cxx.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/shengjiezhai/Desktop/project8/proj8.cxx -o CMakeFiles/proj8.dir/proj8.cxx.s

CMakeFiles/proj8.dir/proj8.cxx.o.requires:

.PHONY : CMakeFiles/proj8.dir/proj8.cxx.o.requires

CMakeFiles/proj8.dir/proj8.cxx.o.provides: CMakeFiles/proj8.dir/proj8.cxx.o.requires
	$(MAKE) -f CMakeFiles/proj8.dir/build.make CMakeFiles/proj8.dir/proj8.cxx.o.provides.build
.PHONY : CMakeFiles/proj8.dir/proj8.cxx.o.provides

CMakeFiles/proj8.dir/proj8.cxx.o.provides.build: CMakeFiles/proj8.dir/proj8.cxx.o


# Object files for target proj8
proj8_OBJECTS = \
"CMakeFiles/proj8.dir/proj8.cxx.o"

# External object files for target proj8
proj8_EXTERNAL_OBJECTS =

proj8.app/Contents/MacOS/proj8: CMakeFiles/proj8.dir/proj8.cxx.o
proj8.app/Contents/MacOS/proj8: CMakeFiles/proj8.dir/build.make
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkDomainsChemistryOpenGL2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersFlowPaths-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersGeneric-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersHyperTree-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersParallelImaging-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersPoints-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersProgrammable-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersSMP-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersSelection-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersTexture-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersTopology-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersVerdict-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkGeovisCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOAMR-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOEnSight-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOExodus-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOExportOpenGL2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOImport-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOInfovis-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOLSDyna-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOMINC-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOMovie-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOPLY-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOParallel-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOParallelXML-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOSQL-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOTecplotTable-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOVideo-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingMorphological-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingStatistics-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingStencil-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkInteractionImage-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingContextOpenGL2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingImage-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingLOD-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingVolumeOpenGL2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkViewsContext2D-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkViewsInfovis-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkDomainsChemistry-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkverdict-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkproj4-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersAMR-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOExport-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingGL2PSOpenGL2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkgl2ps-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtklibharu-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtklibxml2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkoggtheora-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersParallel-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkexoIIc-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOGeometry-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIONetCDF-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtknetcdfcpp-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkNetCDF-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkhdf5_hl-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkhdf5-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkjsoncpp-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkParallelCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOLegacy-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtksqlite-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingOpenGL2-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkglew-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingMath-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkChartsCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingContext2D-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersImaging-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkInfovisLayout-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkInfovisCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkViewsCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkInteractionWidgets-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersHybrid-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingGeneral-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingSources-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersModeling-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingHybrid-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOImage-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkDICOMParser-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkmetaio-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkpng-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtktiff-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkjpeg-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /usr/lib/libm.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkInteractionStyle-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersExtraction-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersStatistics-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingFourier-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkalglib-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingAnnotation-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingColor-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingVolume-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkImagingCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOXML-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOXMLParser-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkIOCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtklz4-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkexpat-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingLabel-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingFreeType-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkRenderingCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonColor-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersGeometry-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersSources-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersGeneral-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonComputationalGeometry-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkFiltersCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonExecutionModel-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonDataModel-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonMisc-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonSystem-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtksys-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonTransforms-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonMath-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkCommonCore-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkfreetype-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: /Users/shengjiezhai/Develop/VTK/VTKBuild/lib/libvtkzlib-8.1.1.dylib
proj8.app/Contents/MacOS/proj8: CMakeFiles/proj8.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/shengjiezhai/Desktop/project8/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable proj8.app/Contents/MacOS/proj8"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/proj8.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/proj8.dir/build: proj8.app/Contents/MacOS/proj8

.PHONY : CMakeFiles/proj8.dir/build

CMakeFiles/proj8.dir/requires: CMakeFiles/proj8.dir/proj8.cxx.o.requires

.PHONY : CMakeFiles/proj8.dir/requires

CMakeFiles/proj8.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/proj8.dir/cmake_clean.cmake
.PHONY : CMakeFiles/proj8.dir/clean

CMakeFiles/proj8.dir/depend:
	cd /Users/shengjiezhai/Desktop/project8 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/shengjiezhai/Desktop/project8 /Users/shengjiezhai/Desktop/project8 /Users/shengjiezhai/Desktop/project8 /Users/shengjiezhai/Desktop/project8 /Users/shengjiezhai/Desktop/project8/CMakeFiles/proj8.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/proj8.dir/depend

