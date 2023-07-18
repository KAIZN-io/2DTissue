# CMake, our favorite Build System Generator
# PLEASE SCROLL TO THE BOTTOM OF THIS FILE FOR READING THE INSTRUCTIONS ON HOW TO ADD YOUR C++ SCRIPT TO THE PROJECT

# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.14)

# Set the project name and version
project(2DTissue VERSION 1.0)

# Make the project foulder path available to the source code
add_definitions("-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"")

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)

# Set macOS RPATH
set(CMAKE_MACOSX_RPATH 1)

# Set library output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Find required packages
find_package(assimp CONFIG REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost COMPONENTS filesystem REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Set the Assimp path
if(assimp_FOUND)
  message(STATUS "Found Assimp library")
else()
  message(FATAL_ERROR "Assimp library not found. Please install Assimp library.")
endif()

# Check if Boost is found
if(NOT Boost_FOUND)
  message(FATAL_ERROR "This project requires the Boost library, and will not be compiled.")
endif()

if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(WARNING "OpenMP not found. Some parts of the project might run slower without OpenMP support.")
endif()

# Include directories for local header files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include(CGAL_Eigen3_support)

# ##########################################################
# Creating entries for all C++ files with "main" routine
# Use target_link_libraries() to express direct dependencies!
# PUBLIC, if you use the libraries in the headers
# ##########################################################

# Add C++ source files
create_single_source_cgal_program("src/simulation/main.cpp")


add_library(utilities_lib STATIC
    src/simulation/2DTissue.cpp
    src/simulation/Simulator.cpp
    src/simulation/IO.cpp
    src/simulation/Cell.cpp
    src/simulation/GeometryProcessing.cpp
    src/simulation/utilities/2D_mapping_fixed_border.cpp
    src/simulation/utilities/analytics.cpp
    src/simulation/utilities/angles_to_unit_vectors.cpp
    src/simulation/utilities/check_boundary.cpp
    src/simulation/utilities/check_validity.cpp
    src/simulation/utilities/dye_particle.cpp
    src/simulation/utilities/error_checking.cpp
    src/simulation/utilities/matrix_algebra.cpp
    src/simulation/utilities/sim_structs.cpp
    src/simulation/utilities/update.cpp
)
target_include_directories(utilities_lib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(utilities_lib PRIVATE CGAL::Eigen3_support Boost::boost Boost::filesystem assimp::assimp)

# Link required libraries to the targets
target_link_libraries(main PRIVATE CGAL::Eigen3_support utilities_lib)

# Install the target
install(TARGETS main
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Collect all test source files
file(GLOB TEST_SOURCES tests/simulation/**/*.cpp)
file(GLOB TEST_SOURCES_CLASS tests/simulation/test_Simulator.cpp tests/simulation/test_GeometryProcessing.cpp)

# Create a test executable for all available tests
add_executable(all_tests ${TEST_SOURCES} ${TEST_SOURCES_CLASS})

# Link your libraries to the test executable
target_link_libraries(all_tests PRIVATE CGAL::Eigen3_support utilities_lib ${GTEST_BOTH_LIBRARIES} pthread)

