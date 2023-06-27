// author: @Jan-Piotraschke
// date: 2023-06-13
// license: Apache License 2.0
// version: 0.1.0

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

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
#include <filesystem>
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
#include <particle_simulation/particle_vector.h>
#include <particle_simulation/simulation.h>

#include <utilities/process_points.h>
#include <utilities/analytics.h>
#include <utilities/dye_particle.h>
#include <utilities/distance.h>
#include <utilities/init_particle.h>
#include <utilities/matrix_algebra.h>
#include <utilities/2D_3D_mapping.h>
#include <utilities/2D_mapping.h>
#include <utilities/sim_structs.h>
#include <utilities/splay_state.h>
#include <utilities/update.h>
#include <utilities/boundary_check.h>
#include <utilities/2D_surface.h>
#include <utilities/validity_check.h>
#include <utilities/process_invalid_particle.h>

// CGAL type aliases
using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;


void error_unvalid_vertices(
    std::vector<VertexData> vertex_data
){
    if (!are_all_valid(vertex_data)) {
        throw std::runtime_error("There are still particles outside the mesh");
    }
}


void error_invalid_values(
    Eigen::MatrixXd r_new
){
    if (checkForInvalidValues(r_new)) {
        std::exit(1);  // stop script execution
    }
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd> perform_particle_simulation(
    Eigen::MatrixXd& r_UV,
    Eigen::MatrixXd& n,
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
    Eigen::MatrixXi faces_uv = loadMeshFaces(mesh_file_path);

    // 1. Simulate the flight of the particle on the UV mesh
    auto [r_UV_new, r_dot, dist_length] = simulate_flight(r_UV, n, vertices_3D_active, distance_matrix_v, v0, k, σ, μ, r_adh, k_adh, step_size);

    // Map the new UV coordinates back to the UV mesh
    opposite_seam_edges(r_UV_new);

    /*
    Error checkings
    */
    // 1. Check if we lost particles
    if (find_inside_uv_vertices_id(r_UV_new).size() != num_part) {
        throw std::runtime_error("We lost particles after getting the original UV mesh coord");
    }

    // 2. Check if there are invalid values like NaN or Inf in the output
    error_invalid_values(r_UV_new);

    // Dye the particles based on their distance
    Eigen::VectorXd particles_color = dye_particles(dist_length, σ);

    // Calculate the order parameter
    calculate_order_parameter(v_order, r_UV, r_dot, current_step);

    return std::make_tuple(r_UV_new, r_dot, dist_length, n, particles_color);
}
