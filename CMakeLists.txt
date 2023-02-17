project(CppHello)

cmake_minimum_required(VERSION 3.14)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_STANDARD 14)  # Sets the C++ standard to C++17
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


# `using Pkg; Pkg.instantiate()`
# Add the installation prefix of "JlCxx" to CMAKE_PREFIX_PATH or set
# "JlCxx_DIR" to a directory containing one of the above files.  If "JlCxx"
# provides a separate development package or SDK, be sure it has been
# installed.
# ! Next, you can build your own code against this by setting the JlCxx_DIR CMake variable to the build directory 
# (libcxxwrap-julia-build) used above, or add it to the CMAKE_PREFIX_PATH CMake variable.

list(APPEND CMAKE_PREFIX_PATH "${HOME}/.julia/artifacts/434257acf70ce7407e3946f2c8452e814d7cca64")

find_package(CGAL 5.5.1 REQUIRED)

# Include here your dependend libraries (e.g. CGAL, Eigen, Boost, etc.)
find_package(Boost 1.80.0 REQUIRED)  # find_package(Boost 1.56 REQUIRED COMPONENTS date_time filesystem iostreams)

if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled."
    )
  return()
endif()

find_package(Eigen3 3.4.0 REQUIRED)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message(
    STATUS
      "This project requires the Eigen library (3.3 or greater), and will not be compiled."
  )
  return()
endif()

find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")

# include for local directory
include_directories(BEFORE include)

# Creating entries for all C++ files with "main" routine
# Use target_link_libraries() to express direct dependencies!
# PUBLIC, if you use the libraries in the headers
# ##########################################################

create_single_source_cgal_program("src/geodesic_distance.cpp")
target_link_libraries(geodesic_distance PUBLIC CGAL::Eigen3_support JlCxx::cxxwrap_julia)
create_single_source_cgal_program("src/create_uv_surface.cpp")
target_link_libraries(create_uv_surface PUBLIC CGAL::Eigen3_support JlCxx::cxxwrap_julia)

# Compile
install(TARGETS geodesic_distance create_uv_surface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)
