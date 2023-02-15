project(CppHello)

cmake_minimum_required(VERSION 3.14)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


# `using Pkg; Pkg.instantiate()`
# Add the installation prefix of "JlCxx" to CMAKE_PREFIX_PATH or set
# "JlCxx_DIR" to a directory containing one of the above files.  If "JlCxx"
# provides a separate development package or SDK, be sure it has been
# installed.
# ! Next, you can build your own code against this by setting the JlCxx_DIR CMake variable to the build directory 
# (libcxxwrap-julia-build) used above, or add it to the CMAKE_PREFIX_PATH CMake variable.
list(APPEND CMAKE_PREFIX_PATH "/Users/jan-piotraschke/.julia/artifacts/434257acf70ce7407e3946f2c8452e814d7cca64")

# Find JlCxx
find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")

# Add our library
add_library(libhello SHARED src/libhello.cpp)

# Link to JlCxx
target_link_libraries(libhello JlCxx::cxxwrap_julia)
# target_link_libraries(libhello JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)

# Install
install(TARGETS libhello
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)
