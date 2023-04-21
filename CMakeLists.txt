# CMake, our favorite Build System Generator
# PLEASE SCROLL TO THE BOTTOM OF THIS FILE FOR READING THE INSTRUCTIONS ON HOW TO ADD YOUR C++ SCRIPT TO THE PROJECT

# Set the minimum version of CMake required
cmake_minimum_required(VERSION 3.14)

# Set the Assimp path
find_program(ASSIMP_EXECUTABLE NAMES assimp)
if(ASSIMP_EXECUTABLE)
  get_filename_component(ASSIMP_ROOT_DIR "${ASSIMP_EXECUTABLE}" DIRECTORY)
  get_filename_component(ASSIMP_ROOT_DIR "${ASSIMP_ROOT_DIR}/.." ABSOLUTE)
  message(STATUS "Found Assimp executable: ${ASSIMP_EXECUTABLE}")
  message(STATUS "Assimp root directory: ${ASSIMP_ROOT_DIR}")
else()
  message(FATAL_ERROR "Assimp executable not found. Please install Assimp library.")
endif()

set(ASSIMP_ROOT "${ASSIMP_ROOT_DIR}")
set(ASSIMP_INCLUDE_DIRS "${ASSIMP_ROOT}/include")
set(ASSIMP_LIBRARY_DIRS "${ASSIMP_ROOT}/lib")

# Set the project name and version
project(2DTissue VERSION 1.0)

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
find_package(OpenCV REQUIRED)

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
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include(CGAL_Eigen3_support)

# Include the Assimp headers and libraries
include_directories(${ASSIMP_INCLUDE_DIRS})
link_directories(${ASSIMP_LIBRARY_DIRS})

include_directories(${OpenCV_INCLUDE_DIRS})

# ##########################################################
# Creating entries for all C++ files with "main" routine
# Use target_link_libraries() to express direct dependencies!
# PUBLIC, if you use the libraries in the headers
# ##########################################################

# Add C++ source files
create_single_source_cgal_program("src/simulation/main.cpp")

add_library(io_lib STATIC
    src/simulation/io/csv.cpp
    src/simulation/io/mesh_loader.cpp
    src/simulation/utilities/sim_structs.cpp
)
target_include_directories(io_lib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(io_lib PRIVATE CGAL::Eigen3_support Boost::boost)

add_library(particle_simulation_lib STATIC
    src/simulation/particle_simulation/motion.cpp
    src/simulation/particle_simulation/forces.cpp
    src/simulation/particle_simulation/particle_vector.cpp
    src/simulation/particle_simulation/process_invalid_particle.cpp
)
target_include_directories(particle_simulation_lib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(particle_simulation_lib PRIVATE CGAL::Eigen3_support Boost::boost)

add_library(utilities_lib STATIC
    src/simulation/utilities/analytics.cpp
    src/simulation/utilities/nearest_map.cpp
    src/simulation/utilities/distance.cpp
    src/simulation/utilities/init_particle.cpp
    src/simulation/utilities/julia_handler.cpp
    src/simulation/utilities/matrix_algebra.cpp
    src/simulation/utilities/mesh_analysis.cpp
    src/simulation/utilities/sim_structs.cpp
    src/simulation/utilities/splay_state.cpp
    src/simulation/utilities/update.cpp
    src/simulation/utilities/boundary_check.cpp
    src/simulation/utilities/uv_surface.cpp
    src/simulation/utilities/validity_check.cpp
    src/simulation/utilities/dye_particle.cpp
)
target_include_directories(utilities_lib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(utilities_lib PRIVATE CGAL::Eigen3_support Boost::boost JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl ${OpenCV_LIBS})

# Link required libraries to the targets
target_link_libraries(main PRIVATE io_lib particle_simulation_lib utilities_lib CGAL::Eigen3_support JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl assimp)

# Install the target
install(TARGETS main
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

