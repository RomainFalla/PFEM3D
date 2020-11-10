include(CheckCXXCompilerFlag)

cmake_minimum_required(VERSION 3.13)

project(PFEM VERSION 1.0
             DESCRIPTION "A Particle Finite Element Method solver for various problems."
             LANGUAGES CXX)

# build type is "" by default in Linux
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)

message(STATUS "Build Type is " ${CMAKE_BUILD_TYPE})

# Use this only if you do not want to make portable libraries !
if(NOT USE_MARCH)
	set(USE_MARCH 0 CACHE STRING "" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
	if(USE_MARCH)
		check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
		if(COMPILER_SUPPORTS_MARCH_NATIVE)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
			message(STATUS "Using -march=native")
		endif()
	endif()
endif()

# put all generated files in a single folder (required in windows for exe linked to dlls)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/ CACHE PATH "")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/ CACHE PATH "")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# enable C++17
set(CMAKE_CXX_STANDARD 17) # newer way to set C++17 (requires cmake>3.1)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_WARNINGS)
    add_definitions(-D_USE_MATH_DEFINES) # for M_PI
endif()

#add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON) # for Valgrind

if(APPLE)
    # on macOS, do not give priority to frameworks/apps
    set(CMAKE_FIND_APPBUNDLE LAST)
    set(CMAKE_FIND_FRAMEWORK LAST)
endif()

# find gmsh-sdk
# gmsh.h
find_path(GMSH_INCLUDE_DIRS NAMES "gmsh.h")
if(NOT GMSH_INCLUDE_DIRS)
    message(FATAL_ERROR "gmsh.h not found!")
else()
    message(STATUS "Found GMSH: " ${GMSH_INCLUDE_DIRS})
endif()

# libgmsh.so
find_library(GMSH_LIBRARIES gmsh)
if(NOT GMSH_LIBRARIES)
    message(FATAL_ERROR "GMSH library not found!")
else()
    message(STATUS "Found GMSH: " ${GMSH_LIBRARIES})
endif()

find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "OpenMP not found!")
else()
    message(STATUS "Found OpenMP: " ${OpenMP_CXX_FLAGS})
    message(STATUS "Found OpenMP: " ${OpenMP_CXX_LIBRARIES})
endif()

add_subdirectory(srcs)
enable_testing()
