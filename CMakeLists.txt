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

# Add the installation prefix of "JlCxx" to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH "${HOME}/.julia/artifacts/434257acf70ce7407e3946f2c8452e814d7cca64")

# Find CGAL and required components
find_package(CGAL 5.5.1 REQUIRED)

# Include here your dependend libraries (e.g. CGAL, Eigen, Boost, etc.)
find_package(Boost 1.80.0 REQUIRED)

if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled."
    )
  return()
endif()

# Find Eigen
find_package(Eigen3 3.4.0 REQUIRED)

include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message(
    STATUS
      "This project requires the Eigen library (3.3 or greater), and will not be compiled."
  )
  return()
endif()

# Find JlCxx
find_package(JlCxx REQUIRED)

# Get JlCxx location and set install RPATH
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")

# Include local header files
include_directories("${PROJECT_DIR}/include")


# ##########################################################
# Creating entries for all C++ files with "main" routine
# Use target_link_libraries() to express direct dependencies!
# PUBLIC, if you use the libraries in the headers
# ##########################################################

# Add C++ source files
create_single_source_cgal_program("src/GeomProcessing/geodesic_distance.cpp")
target_link_libraries(geodesic_distance PUBLIC CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
create_single_source_cgal_program("src/GeomProcessing/create_uv_surface.cpp")
target_link_libraries(create_uv_surface PUBLIC CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)

# Install the target
install(TARGETS geodesic_distance create_uv_surface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)
