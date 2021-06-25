cmake_minimum_required(VERSION 3.12)
include(CMakeDependentOption)

file(READ VERSION FULL_VERSION_STRING)
string(STRIP "${FULL_VERSION_STRING}" FULL_VERSION_STRING)
string(REGEX MATCH "^[0-9]+(\\.[0-9]+)?(\\.[0-9]+)?(\\.[0-9]+)?" numeric_version "${FULL_VERSION_STRING}")

project(arbor VERSION ${numeric_version})
enable_language(CXX)

include(GNUInstallDirs)

# Turn on this option to force the compilers to produce color output when output is
# redirected from the terminal (e.g. when using ninja or a pager).

option(ARBDEV_COLOR "Always produce ANSI-colored output (GNU/Clang only)." OFF)

#----------------------------------------------------------
# Configure-time build options for Arbor:
#----------------------------------------------------------

# Specify target archiecture.

set(ARB_ARCH "native" CACHE STRING "Target architecture for arbor libraries")

# Perform explicit vectorization?

option(ARB_VECTORIZE "use explicit SIMD code in generated mechanisms" OFF)

# Use externally built modcc?

set(ARB_MODCC "" CACHE STRING "path to external modcc NMODL compiler")

# Use libunwind to generate stack traces on errors?

option(ARB_UNWIND "Use libunwind for stack trace printing if available" OFF)

# Specify GPU build type

set(ARB_GPU "none" CACHE STRING "GPU backend and compiler configuration")
set_property(CACHE PROPERTY STRINGS "none" "cuda" "cuda-clang" "hip")

# Use bundled 3rd party libraries

option(ARB_USE_BUNDLED_LIBS "Use bundled 3rd party libraries" OFF)

# Optional additional CXX Flags used for all code that will run on the target
# CPU architecture. Recorded in installed target, for downstream dependencies
# to use.
# Useful, for example, when a user wants to compile with target-specific
# optimization flags.
set(ARB_CXX_FLAGS_TARGET "" CACHE STRING "Optional additional flags for compilation")

#----------------------------------------------------------
# Debug support
#----------------------------------------------------------

# Print builtin catalogue configuration while building
option(ARB_CAT_VERBOSE "Print catalogue build information" OFF)
mark_as_advanced(ARB_CAT_VERBOSE)

#----------------------------------------------------------
# Configure-time features for Arbor:
#----------------------------------------------------------

option(ARB_WITH_MPI "build with MPI support" OFF)

option(ARB_WITH_PROFILING "use built-in profiling" OFF)

option(ARB_WITH_ASSERTIONS "enable arb_assert() assertions in code" OFF)

#----------------------------------------------------------
# NeuroML support library:
#----------------------------------------------------------

option(ARB_WITH_NEUROML "build NeuroML support library" OFF)

#----------------------------------------------------------
# Python front end for Arbor:
#----------------------------------------------------------

option(ARB_WITH_PYTHON "enable Python front end" OFF)

#----------------------------------------------------------
# Global CMake configuration
#----------------------------------------------------------

# Include own CMake modules in search path, load common modules.

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(GitSubmodule) # required for check_git_submodule
include(ErrorTarget)  # reguired for add_error_target
include(FindThreadsCudaFix) # bug work around

# Set release as the default build type (CMake default is debug.)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "debug" "release")
endif()

# When we find threads, prefer not to use the -pthread option
# in order to avoid a CMake 3.12 issue.
set(THREADS_PREFER_PTHREAD_FLAG OFF)

# Add CUDA as a language if GPU support requested.
# (This has to be set early so as to enable CUDA tests in generator
# expressions.)

if(ARB_GPU STREQUAL "cuda")
    set(ARB_WITH_NVCC TRUE)

    # CMake 18 and later set the default CUDA architecture for
    # each target according to CMAKE_CUDA_ARCHITECTURES.
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 60 70 80)
    endif()

    # This fixes nvcc picking up a wrong host compiler for linking, causing issues
    # with outdated libraries, eg libstdc++ and std::filesystem. Must happend before
    # all calls to enable_language(CUDA)
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

    enable_language(CUDA)

    # Despite native CUDA support, the CUDA package is still required to export
    # the cuda library dependencies from the installed target.
    find_package(CUDA 10 REQUIRED)

elseif(ARB_GPU STREQUAL "cuda-clang")
    set(ARB_WITH_CUDA_CLANG TRUE)

    # The CUDA package is needed for clang compilation for the same reasons as above.
    # enable_langaue(CUDA) has a bug with clang
    find_package(CUDA 10 REQUIRED)
elseif(ARB_GPU STREQUAL "hip")
    set(ARB_WITH_HIP_CLANG TRUE)
endif()

if(ARB_WITH_NVCC OR ARB_WITH_CUDA_CLANG OR ARB_WITH_HIP_CLANG)
    set(ARB_WITH_GPU TRUE)
endif()

# Build paths.

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Generate a .json file with full compilation command for each file.

set(CMAKE_EXPORT_COMPILE_COMMANDS "YES")

# Detect and deprecate xlC.

include("CheckCompilerXLC")

# Compiler options common to library, examples, tests, etc.

include("CompilerOptions")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${CXXOPT_WALL}>")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#----------------------------------------------------------
# Set up flags and dependencies:
#----------------------------------------------------------

# Note: any target dependency of arbor needs to be explicitly added
# to the 'export set', even the private ones, and this must be done
# in the same CMakeLists.txt in which the target is defined.

# Data and internal scripts go here
set(ARB_INSTALL_DATADIR ${CMAKE_INSTALL_FULL_DATAROOTDIR}/arbor)

# Interface library `arbor-config-defs` collects configure-time defines
# for arbor, arborenv, arborio, of the form ARB_HAVE_XXX. These
# defines should _not_ be used in any installed public headers.

add_library(arbor-config-defs INTERFACE)
install(TARGETS arbor-config-defs EXPORT arbor-targets)

# Interface library `arbor-private-deps` collects dependencies, options etc.
# for the arbor library.
add_library(arbor-private-deps INTERFACE)
target_link_libraries(arbor-private-deps INTERFACE arbor-config-defs ext-random123 ${CMAKE_DL_LIBS})
install(TARGETS arbor-private-deps EXPORT arbor-targets)

# Interface library `arborenv-private-deps` collects dependencies, options etc.
# for the arborenv library.

add_library(arborenv-private-deps INTERFACE)
target_link_libraries(arborenv-private-deps INTERFACE arbor-config-defs)
install(TARGETS arborenv-private-deps EXPORT arbor-targets)

# Interface library `arborio-private-deps` collects dependencies, options etc.
# for the arborio library.

add_library(arborio-private-deps INTERFACE)
target_link_libraries(arborio-private-deps INTERFACE arbor-config-defs)
install(TARGETS arborio-private-deps EXPORT arbor-targets)

# Interface library `arbor-public-deps` collects requirements for the
# users of the arbor library (e.g. mpi) that will become part
# of arbor's PUBLIC interface.

add_library(arbor-public-deps INTERFACE)
install(TARGETS arbor-public-deps EXPORT arbor-targets)

# Interface library `arborio-public-deps` collects requirements for the
# users of the arborio library (e.g. xml libs) that will become part
# of arborio's PUBLIC interface.

add_library(arborio-public-deps INTERFACE)
install(TARGETS arborio-public-deps EXPORT arborio-targets)

# Add scripts and supporting CMake for setting up external catalogues

configure_file(scripts/build-catalogue.in ${CMAKE_CURRENT_BINARY_DIR}/build-catalogue @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/build-catalogue DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
install(FILES mechanisms/BuildModules.cmake DESTINATION ${ARB_INSTALL_DATADIR})
install(FILES mechanisms/generate_catalogue DESTINATION ${ARB_INSTALL_DATADIR} PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
# External libraries in `ext` sub-directory: json, tinyopt and randon123.
# Creates interface libraries `ext-json`, `ext-tinyopt` and `ext-random123`

cmake_dependent_option(ARB_USE_BUNDLED_JSON "Use bundled Niels Lohmann's json library." ON "ARB_USE_BUNDLED_LIBS" OFF)
if(NOT ARB_USE_BUNDLED_JSON)
  find_package(nlohmann_json)
  set(json_library_name nlohmann_json::nlohmann_json)
else()
  unset(nlohmann_json_DIR CACHE)
endif()
add_subdirectory(ext)

cmake_dependent_option(ARB_USE_BUNDLED_FMT "Use bundled FMT lib." ON "ARB_USE_BUNDLED_LIBS" OFF)

# Keep track of packages we need to add to the generated CMake config
# file for arbor.

set(arbor_export_dependencies)

# Keep track of which 'components' of arbor are included (this is
# currently just 'MPI' support and 'neuroml' for NeuroML support in
# libarborio.)

set(arbor_supported_components)

# Target microarchitecture for building arbor libraries, tests and examples
#---------------------------------------------------------------------------

# Set the full set of target flags in ARB_CXX_FLAGS_TARGET_FULL, which
# will include target-specific -march flags if ARB_ARCH is not "none".
if(ARB_ARCH STREQUAL "none")
    set(ARB_CXX_FLAGS_TARGET_FULL ${ARB_CXX_FLAGS_TARGET})
else()
    set_arch_target(ARB_CXXOPT_ARCH ${ARB_ARCH})
    set(ARB_CXX_FLAGS_TARGET_FULL ${ARB_CXX_FLAGS_TARGET} ${ARB_CXXOPT_ARCH})
endif()
separate_arguments(ARB_CXX_FLAGS_TARGET_FULL)

target_compile_options(arbor-private-deps INTERFACE ${ARB_CXX_FLAGS_TARGET_FULL})
target_compile_options(arborenv-private-deps INTERFACE ${ARB_CXX_FLAGS_TARGET_FULL})

# Profiling and test features
#-----------------------------

if(ARB_WITH_PROFILING)
    target_compile_definitions(arbor-config-defs INTERFACE ARB_HAVE_PROFILING)
endif()
if(ARB_WITH_ASSERTIONS)
    target_compile_definitions(arbor-config-defs INTERFACE ARB_HAVE_ASSERTIONS)
endif()

# Python bindings
#----------------------------------------------------------

# The minimum version of Python supported by Arbor.
set(arb_py_version 3.6.0)

if(DEFINED PYTHON_EXECUTABLE)
    set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()

if(ARB_WITH_PYTHON)
    cmake_dependent_option(ARB_USE_BUNDLED_PYBIND11 "Use bundled pybind11" ON "ARB_WITH_PYTHON;ARB_USE_BUNDLED_LIBS" OFF)

    if(DEFINED ENV{CIBUILDWHEEL})
        find_package(Python3 ${arb_py_version} COMPONENTS Interpreter Development.Module REQUIRED)
    else()
        find_package(Python3 ${arb_py_version} COMPONENTS Interpreter Development REQUIRED)
    endif()

    # Required to link the dynamic libraries for python modules.
    # Effectively adds '-fpic' flag to CXX_FLAGS.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
else()
    # If not building the Python module, the interpreter is still required
    # to build some targets, e.g. when building the documentation.
    find_package(Python3 ${arb_py_version} COMPONENTS Interpreter)
endif()

if(${Python3_FOUND})
    set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
    message(STATUS "PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}")
endif()


# Threading model
#-----------------

find_package(Threads REQUIRED)
find_threads_cuda_fix()
target_link_libraries(arbor-private-deps INTERFACE Threads::Threads)

list(APPEND arbor_export_dependencies "Threads")

# MPI support
#-------------------

if(ARB_WITH_MPI)
    find_package(MPI REQUIRED CXX)
    target_compile_definitions(arbor-config-defs INTERFACE ARB_HAVE_MPI)

    # target_compile_definitions(MPI::MPI_CXX INTERFACE MPICH_SKIP_MPICXX=1 OMPI_SKIP_MPICXX=1)
    # target_link_libraries(arbor-public-deps INTERFACE MPI::MPI_CXX)

    # CMake 3.9 does not allow us to add definitions to an import target.
    # so wrap MPI::MPI_CXX in an interface library 'mpi-wrap' instead.

    add_library(mpi-wrap INTERFACE)
    target_link_libraries(mpi-wrap INTERFACE MPI::MPI_CXX)
    target_compile_definitions(mpi-wrap INTERFACE MPICH_SKIP_MPICXX=1 OMPI_SKIP_MPICXX=1)

    target_link_libraries(arbor-public-deps INTERFACE mpi-wrap)
    install(TARGETS mpi-wrap EXPORT arbor-targets)

    list(APPEND arbor_export_dependencies "MPI\;COMPONENTS\;CXX")
    list(APPEND arbor_supported_components "MPI")
endif()

# CUDA support
#--------------

if(ARB_WITH_GPU)
    target_compile_definitions(arbor-config-defs INTERFACE ARB_HAVE_GPU)

    if(ARB_WITH_NVCC OR ARB_WITH_CUDA_CLANG)
        target_include_directories(arborenv-private-deps INTERFACE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

        add_compile_options(
                "$<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe=--diag_suppress=integer_sign_change>"
                "$<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe=--diag_suppress=unsigned_compare_with_zero>")
    endif()


    if(ARB_WITH_NVCC)
        target_compile_options(arbor-private-deps INTERFACE
            $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_60,code=sm_60>)
        target_compile_options(arbor-private-deps INTERFACE
            $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_70,code=sm_70>)
        if (${CUDA_VERSION_MAJOR} GREATER 10)
            target_compile_options(arbor-private-deps INTERFACE
                $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_80,code=sm_80>)
        endif()

        target_compile_definitions(arbor-private-deps INTERFACE ARB_CUDA)
        target_compile_definitions(arborenv-private-deps INTERFACE ARB_CUDA)

    elseif(ARB_WITH_CUDA_CLANG)
        set(clang_options_ -DARB_CUDA -xcuda --cuda-gpu-arch=sm_60 --cuda-gpu-arch=sm_70 --cuda-gpu-arch=sm_80 --cuda-path=${CUDA_TOOLKIT_ROOT_DIR})
        target_compile_options(arbor-private-deps INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${clang_options_}>)
        target_compile_options(arborenv-private-deps INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${clang_options_}>)

    elseif(ARB_WITH_HIP_CLANG)
        set(clang_options_ -DARB_HIP -xhip --amdgpu-target=gfx906 --amdgpu-target=gfx900)
        target_compile_options(arbor-private-deps INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${clang_options_}>)
        target_compile_options(arborenv-private-deps INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${clang_options_}>)
    endif()
endif()

# Use libunwind if requested for pretty printing stack traces
#-------------------------------------------------------------

if (ARB_UNWIND)
    find_package(Unwind REQUIRED)
    if(Unwind_FOUND)
        target_link_libraries(arbor-private-deps INTERFACE Unwind::unwind)
        target_compile_definitions(arbor-private-deps INTERFACE WITH_UNWIND)

        list(APPEND arbor_export_dependencies "Unwind")
    endif()
endif()

# Build and use modcc unless explicit path given
#------------------------------------------------

if(ARB_MODCC)
    find_program(modcc NAMES ${ARB_MODCC} NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
    if(NOT modcc)
        message(FATAL_ERROR "Unable to find modcc executable.")
    endif()
    set(ARB_WITH_EXTERNAL_MODCC TRUE)
else()
    set(modcc $<TARGET_FILE:modcc>)
    set(ARB_WITH_EXTERNAL_MODCC FALSE)
endif()

set(ARB_MODCC_FLAGS)
if(ARB_VECTORIZE)
    list(APPEND ARB_MODCC_FLAGS "--simd")
endif()
if(ARB_WITH_PROFILING)
    list(APPEND ARB_MODCC_FLAGS "--profile")
endif()

#----------------------------------------------------------
# Set up install paths, permissions.
#----------------------------------------------------------

# Set up install paths according to GNU conventions.
#
# GNUInstallDirs picks (e.g.) `lib64` for the library install path on some
# systems where this is definitely not correct (e.g. Arch Linux). If there
# are cases where `lib` is inappropriate, we will have to incorporate special
# case behaviour here.

if(NOT CMAKE_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR lib)
endif()
include(GNUInstallDirs)

# Implicitly created directories require permissions to be set explicitly
# via this CMake variable.
#
# Note that this has no effect until CMake version 3.11.

set(CMAKE_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
    WORLD_READ
    WORLD_EXECUTE)

# CMake versions 3.11 and 3.12 ignore this variable for directories
# implicitly created by install(DIRECTORY ...), which for us corresponds
# to our doc and include directories. Work-around by trying to install
# a non-existant file to these locations.

foreach(directory "${CMAKE_INSTALL_DOCDIR}" "${CMAKE_INSTALL_INCLUDEDIR}")
    install(FILES _no_such_file_ OPTIONAL DESTINATION "${directory}")
endforeach()

#----------------------------------------------------------
# Configure targets in sub-directories.
#----------------------------------------------------------

# arbor-public-headers:
add_subdirectory(arbor/include)

# arbor-sup:
add_subdirectory(sup)

# modcc, libmodcc:
add_subdirectory(modcc)

# arbor, arbor-private-headers:
add_subdirectory(arbor)

# arborenv, arborenv-public-headers:
add_subdirectory(arborenv)

# arborio, arborio-public-headers:
add_subdirectory(arborio)

# unit, unit-mpi, unit-local, unit-modcc
add_subdirectory(test)

# self contained examples:
add_subdirectory(example)

# html:
add_subdirectory(doc)

# python interface:
if(ARB_WITH_PYTHON)
    add_subdirectory(python)
endif()

#----------------------------------------------------------
# Generate CMake config/version files for install.
#----------------------------------------------------------

# Note: each dependency for the arbor library target, private or otherwise,
# needs to add itself to the arbor-exports EXPORT target in the subdirectory
# in which they are defined, or none of this will work.

set(cmake_config_dir "${CMAKE_INSTALL_LIBDIR}/cmake/arbor")
install(EXPORT arbor-targets NAMESPACE arbor:: DESTINATION "${cmake_config_dir}")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/arbor-config-version.cmake"
    COMPATIBILITY SameMajorVersion)

# Template file will use contents of arbor_export_dependencies to include the
# required `find_dependency` statements, and arbor_supported_components will
# be used to check feature support.
#
# To avoid CMake users of the installed arbor library conditionally requiring
# that they add CUDA to their project language, explicitly munge the import
# language and library dependencies on the installed target if ARB_WITH_GPU
# is set, via the variables arbor_override_import_lang and arbor_add_import_libs.
# arbor_build_config records our build type in a way compatible with the
# generated export cmake files.

set(arbor_build_config NOCONFIG)
if(CMAKE_BUILD_TYPE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" arbor_build_config)
endif()

set(arbor_override_import_lang)
set(arbor_add_import_libs)
set(arborenv_add_import_libs)
set(arborio_add_import_libs)

if(ARB_WITH_GPU)
    set(arbor_override_import_lang CXX)
    set(arbor_add_import_libs ${CUDA_LIBRARIES})
    set(arborenv_add_import_libs ${CUDA_LIBRARIES})
endif()

# (We remove old generated one so that the generation happens every time we run cmake.)
file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/arbor-config.cmake")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/arbor-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/arbor-config.cmake"
    @ONLY)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/arbor-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/arbor-config-version.cmake"
        cmake/FindUnwind.cmake
    DESTINATION "${cmake_config_dir}")

