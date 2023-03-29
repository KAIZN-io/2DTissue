# CMake, our favorite Build System Generator
# PLEASE SCROLL TO THE BOTTOM OF THIS FILE FOR READING THE INSTRUCTIONS ON HOW TO ADD YOUR C++ SCRIPT TO THE PROJECT

# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.14)

# Set the project name and version
project(PhysicsEngineBuilds VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set macOS RPATH
set(CMAKE_MACOSX_RPATH 1)

# Set library output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Find required packages
find_package(CGAL 5.5.1 REQUIRED)
find_package(Boost 1.80.0 REQUIRED)
find_package(Eigen3 3.4.0 REQUIRED)
find_package(JlCxx REQUIRED)

# Check if Boost is found
if(NOT Boost_FOUND)
  message(FATAL_ERROR "This project requires the Boost library, and will not be compiled.")
endif()

# Get JlCxx location and set install RPATH
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")

# Include directories for local header files
include_directories("${PROJECT_DIR}/include")
include(CGAL_Eigen3_support)


# ##########################################################
# Creating entries for all C++ files with "main" routine
# Use target_link_libraries() to express direct dependencies!
# PUBLIC, if you use the libraries in the headers
# ##########################################################

# Add C++ source files
create_single_source_cgal_program("src/GeomProcessing/geodesic_distance.cpp")
create_single_source_cgal_program("src/GeomProcessing/create_uv_surface.cpp")
create_single_source_cgal_program("src/GeomProcessing/basic.cpp")
create_single_source_cgal_program("src/GeomProcessing/particle_simulation.cpp")

# Link required libraries to the targets
target_link_libraries(geodesic_distance PRIVATE CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
target_link_libraries(create_uv_surface PRIVATE CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
target_link_libraries(basic PRIVATE CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
target_link_libraries(particle_simulation PRIVATE CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)

# Install the target
install(TARGETS geodesic_distance create_uv_surface basic particle_simulation
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
