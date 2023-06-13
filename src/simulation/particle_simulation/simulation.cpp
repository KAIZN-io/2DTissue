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

#include <utilities/analytics.h>
#include <utilities/dye_particle.h>
#include <utilities/distance.h>
#include <utilities/init_particle.h>
#include <utilities/matrix_algebra.h>
#include <utilities/2D_3D_mapping.h>
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
        std::cout << "Invalid vertices:\n";
        for (const VertexData& vd : vertex_data) {
            if (!vd.valid) {
                std::cout << "Old ID: " << vd.old_id << ", Next ID: " << vd.next_id << ", Valid: " << vd.valid << ", UV Mesh ID: " << vd.uv_mesh_id << '\n';
            }
        }
        throw std::runtime_error("There are still particles outside the mesh");
    }
}


void validate_vertices(
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    std::vector<VertexData>& vertex_data,
    int num_part,
    Eigen::MatrixXd distance_matrix_v,
    Eigen::MatrixXd n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
){
    if (!are_all_valid(vertex_data)) {
        process_if_not_valid(vertices_2DTissue_map, vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);
    }

    // Throw an error if there are still invalid vertices
    error_unvalid_vertices(vertex_data);
}


std::vector<int64_t> get_next_ids(const std::vector<VertexData>& vertex_data){
    std::vector<int64_t> vertices_next_id(vertex_data.size());
    for (size_t i = 0; i < vertex_data.size(); ++i) {
        vertices_next_id[i] = vertex_data[i].next_id;
    }

    return vertices_next_id;
}


std::vector<int64_t> find_next_position(
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    std::vector<VertexData>& vertex_data,
    int num_part,
    Eigen::MatrixXd distance_matrix_v,
    Eigen::MatrixXd n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
){
    validate_vertices(vertices_2DTissue_map, vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);

    return get_next_ids(vertex_data);
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd> perform_particle_simulation(
    Eigen::MatrixXd& r,
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
    double dt,
    int tt,
    int num_part,
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    double plotstep
){
    // Simulate the flight of the particle
    auto [r_new, r_dot, dist_length] = simulate_flight(r, n, vertices_3D_active, distance_matrix_v, v0, k, σ, μ, r_adh, k_adh, dt);

    // Find the particles which landed inside the mesh
    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_new);
    // Find the particles which landed outside the mesh
    std::vector<int> outside_uv_ids = set_difference(num_part, inside_uv_ids);

    // Get the original mesh from the dictionary
    auto [halfedges_uv, h_v_mapping] = get_mesh_data(vertices_2DTissue_map, 0);

    // Find the new suggested 3D vertex ids, even if they are outside the 2D mesh
    Eigen::VectorXd vertice_3D_id = get_vertice_id(r_new, halfedges_uv, h_v_mapping);

    // Update our struct for this time step for the particles which landed Inside the mesh
    std::vector<VertexData> vertex_data = update_vertex_data(vertices_3D_active, vertice_3D_id, inside_uv_ids, 0);

    // Check and process invalid particles, which landed Outside the mesh
    auto vertices_next_id = find_next_position(vertices_2DTissue_map, vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);

    // Update the data for the previous particles which landed Outside
    for (int i : outside_uv_ids) {
        std::vector<int64_t> single_vertex_next_id = {vertices_next_id[i]};
        std::vector<int64_t> halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(single_vertex_next_id, h_v_mapping);

        Eigen::MatrixXd r_new_temp_single_row = get_r_from_halfedge_id(halfedge_id, halfedges_uv);
        r_new.row(i) = r_new_temp_single_row.row(0);
    }

    // Check if we lost particles
    if (find_inside_uv_vertices_id(r_new).size() != num_part) {
        throw std::runtime_error("We lost particles after getting the original mesh halfedges coord");
    }

    // Check if there are invalid values like NaN or Inf in the output
    if (checkForInvalidValues(r_new)) {
        std::cout << "Invalid values found in r: " << std::endl;
        std::cout << r_new << std::endl;
        std::exit(1);  // stop script execution
    }

    // Dye the particles based on their distance
    Eigen::VectorXd particles_color = dye_particles(dist_length, σ);

    // Calculate the order parameter
    calculate_order_parameter(v_order, r, r_dot, tt);

    return std::make_tuple(r_new, r_dot, dist_length, n, particles_color);
}
