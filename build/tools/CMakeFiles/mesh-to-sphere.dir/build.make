# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/benuix/Setup/SphericalMesh

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/benuix/Setup/SphericalMesh/build

# Include any dependencies generated for this target.
include tools/CMakeFiles/mesh-to-sphere.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/mesh-to-sphere.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/mesh-to-sphere.dir/flags.make

tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o: tools/CMakeFiles/mesh-to-sphere.dir/flags.make
tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o: ../tools/mesh-to-sphere.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/benuix/Setup/SphericalMesh/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o"
	cd /home/benuix/Setup/SphericalMesh/build/tools && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o -c /home/benuix/Setup/SphericalMesh/tools/mesh-to-sphere.cc

tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.i"
	cd /home/benuix/Setup/SphericalMesh/build/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/benuix/Setup/SphericalMesh/tools/mesh-to-sphere.cc > CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.i

tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.s"
	cd /home/benuix/Setup/SphericalMesh/build/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/benuix/Setup/SphericalMesh/tools/mesh-to-sphere.cc -o CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.s

tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.requires:

.PHONY : tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.requires

tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.provides: tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.requires
	$(MAKE) -f tools/CMakeFiles/mesh-to-sphere.dir/build.make tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.provides.build
.PHONY : tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.provides

tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.provides.build: tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o


# Object files for target mesh-to-sphere
mesh__to__sphere_OBJECTS = \
"CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o"

# External object files for target mesh-to-sphere
mesh__to__sphere_EXTERNAL_OBJECTS =

bin/mesh-to-sphere: tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o
bin/mesh-to-sphere: tools/CMakeFiles/mesh-to-sphere.dir/build.make
bin/mesh-to-sphere: lib/libsphericalMesh.a
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersGeneric-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersParallelImaging-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingMath-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOExodus-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOMovie-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOLSDyna-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingImage-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkViewsInfovis-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkChartsCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOEnSight-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingStencil-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingLIC-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersHyperTree-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOParallelXML-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkInteractionImage-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingFreeTypeOpenGL-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersProgrammable-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOPLY-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersSelection-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOMINC-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOAMR-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOImport-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOInfovis-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtklibxml2-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersVerdict-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingLOD-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersFlowPaths-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOExport-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingStatistics-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersTexture-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkDomainsChemistry-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOSQL-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingVolumeOpenGL-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkGeovisCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOVideo-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingMorphological-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersSMP-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOParallel-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIONetCDF-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkViewsContext2D-6.2.so.1
bin/mesh-to-sphere: lib/libfastMarchingGeodesic.a
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkoggtheora-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersImaging-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersAMR-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkverdict-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingGL2PS-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingContextOpenGL-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkgl2ps-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingLabel-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtksqlite-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingOpenGL-6.2.so.1
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libGLU.so
bin/mesh-to-sphere: /usr/lib/libXNVCtrl.a
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libSM.so
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libICE.so
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libX11.so
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libXext.so
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libXt.so
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkInfovisLayout-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkInfovisCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkproj4-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersParallel-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkParallelCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOLegacy-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOXML-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOGeometry-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkjsoncpp-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOXMLParser-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkexpat-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkexoIIc-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkNetCDF_cxx-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkNetCDF-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkhdf5_hl-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkhdf5-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingContext2D-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkViewsCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkInteractionWidgets-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersHybrid-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersModeling-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingGeneral-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingSources-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingHybrid-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOImage-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkDICOMParser-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkIOCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkmetaio-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkpng-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtktiff-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkjpeg-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkInteractionStyle-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingAnnotation-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingColor-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingFreeType-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkftgl-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkfreetype-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkzlib-6.2.so.1
bin/mesh-to-sphere: /usr/lib/x86_64-linux-gnu/libGL.so
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingVolume-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkRenderingCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonColor-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersExtraction-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersStatistics-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingFourier-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkImagingCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkalglib-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersGeometry-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersSources-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersGeneral-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkFiltersCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonExecutionModel-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonComputationalGeometry-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonDataModel-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonMisc-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonTransforms-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonMath-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonSystem-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtkCommonCore-6.2.so.1
bin/mesh-to-sphere: /home/benuix/Setup/VTK-6.2.0/build/lib/libvtksys-6.2.so.1
bin/mesh-to-sphere: tools/CMakeFiles/mesh-to-sphere.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/benuix/Setup/SphericalMesh/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/mesh-to-sphere"
	cd /home/benuix/Setup/SphericalMesh/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mesh-to-sphere.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/mesh-to-sphere.dir/build: bin/mesh-to-sphere

.PHONY : tools/CMakeFiles/mesh-to-sphere.dir/build

tools/CMakeFiles/mesh-to-sphere.dir/requires: tools/CMakeFiles/mesh-to-sphere.dir/mesh-to-sphere.cc.o.requires

.PHONY : tools/CMakeFiles/mesh-to-sphere.dir/requires

tools/CMakeFiles/mesh-to-sphere.dir/clean:
	cd /home/benuix/Setup/SphericalMesh/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/mesh-to-sphere.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/mesh-to-sphere.dir/clean

tools/CMakeFiles/mesh-to-sphere.dir/depend:
	cd /home/benuix/Setup/SphericalMesh/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/benuix/Setup/SphericalMesh /home/benuix/Setup/SphericalMesh/tools /home/benuix/Setup/SphericalMesh/build /home/benuix/Setup/SphericalMesh/build/tools /home/benuix/Setup/SphericalMesh/build/tools/CMakeFiles/mesh-to-sphere.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/mesh-to-sphere.dir/depend

