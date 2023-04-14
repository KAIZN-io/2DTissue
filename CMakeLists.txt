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
find_package(OpenMP)

# Set the Assimp path
set(ASSIMP_ROOT "/opt/homebrew/Cellar/assimp/5.2.5")
set(ASSIMP_INCLUDE_DIRS "${ASSIMP_ROOT}/include")
set(ASSIMP_LIBRARY_DIRS "${ASSIMP_ROOT}/lib")

# Include the Assimp headers and libraries
include_directories(${ASSIMP_INCLUDE_DIRS})
link_directories(${ASSIMP_LIBRARY_DIRS})


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


add_library(particle_forces_lib STATIC
    src/GeomProcessing/particle_forces.cpp
)
target_include_directories(particle_forces_lib PRIVATE src)
target_link_libraries(particle_forces_lib PRIVATE CGAL::Eigen3_support Boost::boost)

add_library(matrix_algebra_lib STATIC
    src/GeomProcessing/matrix_algebra.cpp
)
target_include_directories(matrix_algebra_lib PRIVATE src)
target_link_libraries(matrix_algebra_lib PRIVATE CGAL::Eigen3_support Boost::boost)

add_library(uv_surface_lib STATIC
    src/GeomProcessing/dye_particle.cpp
    src/GeomProcessing/flight_of_the_particle.cpp
    src/GeomProcessing/geo_distance.cpp
    src/GeomProcessing/mesh_analysis.cpp
    src/GeomProcessing/mesh_loader.cpp
    src/GeomProcessing/particle_vector.cpp
    src/GeomProcessing/julia_handler.cpp
    src/GeomProcessing/uv_surface.cpp
    src/GeomProcessing/analytics.cpp
    src/GeomProcessing/csv_loader.cpp
    src/GeomProcessing/uv_operations.cpp
)
target_include_directories(uv_surface_lib PRIVATE src)
target_link_libraries(uv_surface_lib PRIVATE CGAL::Eigen3_support Boost::boost JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)


# Link required libraries to the targets
target_link_libraries(geodesic_distance PRIVATE uv_surface_lib CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
target_link_libraries(create_uv_surface PRIVATE CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
target_link_libraries(basic PRIVATE CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl)
target_link_libraries(particle_simulation PRIVATE uv_surface_lib particle_forces_lib matrix_algebra_lib CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl assimp)

# Install the target
install(TARGETS geodesic_distance create_uv_surface basic particle_simulation
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
