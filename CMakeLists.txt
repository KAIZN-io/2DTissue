# CMake, our favorite Build System Generator
# PLEASE SCROLL TO THE BOTTOM OF THIS FILE FOR READING THE INSTRUCTIONS ON HOW TO ADD YOUR C++ SCRIPT TO THE PROJECT

# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.22)

# Set the project name and version
project(2DTissue VERSION 1.0)

# Make the project foulder path available to the source code
add_definitions("-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"")
add_definitions("-DMeshCartographyLib_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}/MeshCartographyLib\"")

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)

# Set macOS RPATH
set(CMAKE_MACOSX_RPATH 1)

# Set library output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Find required packages
find_package(Boost COMPONENTS filesystem REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)

# Include directories for local header files
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GTEST_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/simulation")

# Collect all source files in the 'src/simulation' directory
file(GLOB_RECURSE SIMULATION_SOURCES src/simulation/*.cpp)

# Add subdirectories
add_subdirectory(MeshCartographyLib)

# Add library with simulation sources
add_library(simulation_lib STATIC ${SIMULATION_SOURCES})
target_link_libraries(simulation_lib PRIVATE MeshCartographyLib Boost::boost Boost::filesystem)

# Add C++ source files
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE
        MeshCartographyLib
        simulation_lib
)

# Install the target
install(TARGETS main
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Collect all test source files
file(GLOB_RECURSE TEST_SOURCES tests/simulation/*.cpp)

# Create a test executable for all available tests
add_executable(all_tests ${TEST_SOURCES})

# Link your libraries to the test executable
target_link_libraries(all_tests PRIVATE simulation_lib MeshCartographyLib ${GTEST_BOTH_LIBRARIES} pthread)
