# Works with 3.14 and tested through 3.18
cmake_minimum_required(VERSION 3.14...3.18)

# Project name and a few useful settings. Other commands can pick up the results
project(
  MeshKernel
  VERSION 0.10.0
  DESCRIPTION "Library for creating and editing meshes."
  LANGUAGES CXX C)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O2 -std=c++17")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -std=c++17")
  else()
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} /EHsc /MP /std:c++17")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /EHsc /MP /std:c++17")
  endif()

  # Create position independent binaries
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  # Disable compiler specific extensions
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Note this needs to be done in the main CMakeLists since it calls
  # enable_testing, which must be in the main CMakeLists.
  include(CTest)

  # Docs only available if this is the main app
  add_subdirectory(docs)

  # Run packaging scripts
  add_subdirectory(package)
endif()

# Add third party directory
add_subdirectory(extern)

# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)

# Boost
set(Boost_USE_STATIC_LIBS ON)

find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIR})

# NetCDF
if(WIN32)
  find_path(
    NETCDF_BIN_PATH
    NAMES "${CMAKE_SHARED_LIBRARY_PREFIX}netcdf${CMAKE_SHARED_LIBRARY_SUFFIX}"
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/netcdf/netCDF 4.6.1/bin")
  message(STATUS "netcdf bin dir is ${NETCDF_BIN_PATH}")
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# The static library
add_subdirectory(src/MeshKernel)

# The dynamic library
add_subdirectory(src/MeshKernelApi)

# Testing only available if this is the main app.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(tests)
endif()
