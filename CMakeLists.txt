# CMake, our favorite Build System Generator
# PLEASE SCROLL TO THE BOTTOM OF THIS FILE FOR READING THE INSTRUCTIONS ON HOW TO ADD YOUR C++ SCRIPT TO THE PROJECT

# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.22)

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
# find_package(Sundials REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
#include_directories(
#    "${CMAKE_CURRENT_SOURCE_DIR}/libsbml/src"
#    "${CMAKE_CURRENT_SOURCE_DIR}/libsbml/build/src"
#)

set(ROADRUNNER_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/roadrunner/install-release")
set(LLVM_INSTALL_PREFIX "../../llvm-13.x/install-release")
set(CMAKE_PREFIX_PATH
        "${ROADRUNNER_INSTALL_PREFIX}/lib/cmake/Poco" # needed because this is where poco installs config files to
        "${ROADRUNNER_INSTALL_PREFIX}/lib/cmake" # all other normal packages install config files to here
        )
list(APPEND CMAKE_MODULE_PATH "${ROADRUNNER_INSTALL_PREFIX}/lib/cmake")

include(ImportRoadRunnerAndDependencies)
include_directories("${ROADRUNNER_INSTALL_PREFIX}/include")

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
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/simulation")
include(CGAL_Eigen3_support)

# Collect all source files in the 'src/simulation' directory
file(GLOB SIMULATION_SOURCES src/simulation/*.cpp)

# Add library with simulation sources
add_library(simulation_lib STATIC ${SIMULATION_SOURCES})
target_link_libraries(simulation_lib PRIVATE CGAL::Eigen3_support Boost::boost Boost::filesystem assimp::assimp)

# Add C++ source files
create_single_source_cgal_program("src/main.cpp")
target_link_libraries(main PRIVATE
        CGAL::Eigen3_support
        simulation_lib
        "${CMAKE_CURRENT_SOURCE_DIR}/libsbml/build/src/libsbml.dylib"
        roadrunner-static::roadrunner-static
)

# Install the target
install(TARGETS main
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Collect all test source files
file(GLOB TEST_SOURCES tests/simulation/*.cpp)

# Create a test executable for all available tests
add_executable(all_tests ${TEST_SOURCES})

# Link your libraries to the test executable
target_link_libraries(all_tests PRIVATE CGAL::Eigen3_support simulation_lib ${GTEST_BOTH_LIBRARIES} pthread)
