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

# Include here your dependend libraries (e.g. CGAL, Eigen, Boost, etc.)
find_package(Boost 1.80.0 REQUIRED)  # find_package(Boost 1.56 REQUIRED COMPONENTS date_time filesystem iostreams)
find_package(CGAL 5.5.1 REQUIRED)
find_package(Eigen3 3.4.0 REQUIRED)
find_package(JlCxx REQUIRED)

get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")

# Add our library
add_library(hello SHARED src/hello.cpp)
add_library(func SHARED src/call_func.cpp)
add_executable(geodesic_distance src/geodesic_distance.cpp)
add_library(geo_distance SHARED src/geodesic_distance.cpp)

# Link to JlCxx
# PUBLIC, if you use the libraries in the header
target_link_libraries(hello JlCxx::cxxwrap_julia)
target_link_libraries(func JlCxx::cxxwrap_julia)
target_sources(geodesic_distance PUBLIC Boost CGAL Eigen3 JlCxx::cxxwrap_julia)

# Install
install(TARGETS hello func geo_distance
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)
