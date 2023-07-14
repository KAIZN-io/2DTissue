// author: @Jan-Piotraschke
// date: 2023-06-13
// license: Apache License 2.0
// version: 0.1.0

// Eigen
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_USE_THREADS
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

// Standard libraries
#include <atomic>
#include <cmath>
#include <cstddef>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <stdexcept>
#include <thread>
#include <vector>
#include <limits>

#include <io/csv.h>
#include <io/mesh_loader.h>

#include <particle_simulation/motion.h>

#include <utilities/analytics.h>
#include <utilities/dye_particle.h>
#include <utilities/2D_mapping_fixed_border.h>
// #include <utilities/2D_mapping_free_border.h>
#include <utilities/2D_surface.h>
#include <utilities/error_checking.h>

#include <particle_simulation/simulation.h>


void perform_particle_simulation(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
    Eigen::VectorXd& n,
    Eigen::VectorXi& particles_color,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd distance_matrix_v,
    Eigen::VectorXd& v_order,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double step_size,
    int current_step,
    int num_part,
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    double plotstep
){
    // Get the original mesh from the dictionary
    auto mesh_struct = vertices_2DTissue_map[0];
    Eigen::MatrixXd halfedges_uv = mesh_struct.mesh;
    std::vector<int64_t> h_v_mapping = mesh_struct.h_v_mapping;
    Eigen::MatrixXd vertices_UV = mesh_struct.vertices_UV;
    Eigen::MatrixXd vertices_3D = mesh_struct.vertices_3D;
    std::string mesh_file_path = mesh_struct.mesh_file_path;
    Eigen::MatrixXi faces_uv;
    loadMeshFaces(mesh_file_path, faces_uv);

    // 1. Simulate the flight of the particle on the UV mesh
    auto dist_length = simulate_flight(r_UV, r_dot, n, vertices_3D_active, distance_matrix_v, v0, k, σ, μ, r_adh, k_adh, step_size);

    // Map the new UV coordinates back to the UV mesh
    auto mesh_UV_name = get_mesh_name(mesh_file_path);

    // ! TODO: try to find out why the mesh parametrization can result in different UV mapping logics
    // ? is it because of the seam edge cut line?
    if (mesh_UV_name == "sphere_uv"){
        opposite_seam_edges_square_border(r_UV);
    }
    else {
        diagonal_seam_edges_square_border(r_UV_old, r_UV, n);
    }

    /*
    Error checkings
    */
    error_lost_particles(r_UV, num_part);  // 1. Check if we lost particles
    error_invalid_values(r_UV);  // 2. Check if there are invalid values like NaN or Inf in the output

    // Dye the particles based on their distance
    particles_color = dye_particles(dist_length, σ);

    // The new UV coordinates are the old ones for the next step
    r_UV_old = r_UV;

    // Calculate the order parameter
    calculate_order_parameter(v_order, r_UV, r_dot, current_step);
}
