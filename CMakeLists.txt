# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.22)

# Set the project name and version
project(
  2DTissue
  VERSION 1.0
  DESCRIPTION "Simulating 3D motion in 2D"
  LANGUAGES CXX
)

# Make the project folder path available to the source code
add_definitions("-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"")
add_definitions("-DMeshCartographyLib_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}/MeshCartographyLib\"")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set library output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Find and link the necessary libraries
find_package(GTest REQUIRED)
find_package(PkgConfig REQUIRED)

# Add paths for Kafka and SQLite
pkg_check_modules(LIBRDKAFKA REQUIRED IMPORTED_TARGET librdkafka)
find_package(SQLite3 REQUIRED)

# Include directories
include_directories(${GTEST_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/simulation")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/utils")
include_directories(${SQLite3_INCLUDE_DIRS})

# Collect all source files in the 'src/simulation' directory
file(GLOB_RECURSE SIMULATION_SOURCES src/simulation/*.cpp)

# Collect the source files in the 'src/utils' directory
file(GLOB_RECURSE UTILS_SOURCES src/utils/*.cpp)

# Add subdirectories
add_subdirectory(MeshCartographyLib)

# Add library with simulation sources
add_library(simulation_lib STATIC ${SIMULATION_SOURCES} ${UTILS_SOURCES})
target_link_libraries(simulation_lib PRIVATE MeshCartographyLib PkgConfig::LIBRDKAFKA SQLite::SQLite3)

# Add C++ source files for main executable
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE
        MeshCartographyLib
        simulation_lib
        PkgConfig::LIBRDKAFKA
        SQLite::SQLite3
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
target_link_libraries(all_tests PRIVATE simulation_lib MeshCartographyLib ${GTEST_BOTH_LIBRARIES} pthread PkgConfig::LIBRDKAFKA SQLite::SQLite3)
