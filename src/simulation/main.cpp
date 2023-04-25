// author: @Jan-Piotraschke
// date: 2023-04-14
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
    if (!are_all_valid(vertex_data)) {
        process_if_not_valid(vertices_2DTissue_map, vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);
    }

    if (!are_all_valid(vertex_data)) {
        std::cout << "Invalid vertices:\n";
        for (const VertexData& vd : vertex_data) {
            if (!vd.valid) {
                std::cout << "Old ID: " << vd.old_id << ", Next ID: " << vd.next_id << ", Valid: " << vd.valid << ", UV Mesh ID: " << vd.uv_mesh_id << '\n';
            }
        }
        throw std::runtime_error("There are still particles outside the mesh");
    }

    std::vector<int64_t> vertices_next_id(vertex_data.size());
    for (size_t i = 0; i < vertex_data.size(); ++i) {
        vertices_next_id[i] = vertex_data[i].next_id;
    }

    return vertices_next_id;
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd> perform_particle_simulation(
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
    double plotstep = 0.1
){
    // Simulate the flight of the particle
    auto [r_new, r_dot, dist_length] = simulate_flight(r, n, vertices_3D_active, distance_matrix_v, v0, k, σ, μ, r_adh, k_adh, dt);

    // Find the particles which landed inside the mesh
    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_new);
    // Find the particles which landed outside the mesh
    std::vector<int> outside_uv_ids = set_difference(num_part, inside_uv_ids);

    // Specify the file path of the 3D model you want to load
    Eigen::MatrixXd vertices_3D = loadMeshVertices("/Users/jan-piotraschke/git_repos/2DTissue/meshes/ellipsoid_x4.off");

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

    // Calculate the particle vectors
    auto [ntest, nr_dot] = calculate_particle_vectors(r_dot, n, dt);
    ntest = n;
    // Calculate the output vector v_order
    calculate_order_parameter(v_order, r, r_dot, tt);

    return std::make_tuple(r_new, r_dot, dist_length, ntest, nr_dot, particles_color);
}


int main()
{
    auto v0 = 0.1;
    auto k = 10;
    auto k_next = 10;
    auto v0_next = 0.1;
    auto σ = 0.4166666666666667;
    auto μ = 1;
    auto r_adh = 1;
    auto k_adh = 0.75;
    auto dt = 0.001;
    int num_part = 200;
    int num_frames = 10;

    static std::unordered_map<int, Mesh_UV_Struct> vertices_2DTissue_map;

    auto result = create_uv_surface_intern("Ellipsoid", 0);
    std::vector<int64_t> h_v_mapping_vector = std::get<0>(result);  // halfedge-vertice mapping
    std::string mesh_file_path = std::get<1>(result);

    Eigen::MatrixXd halfedge_uv = loadMeshVertices(mesh_file_path);
    Eigen::MatrixXi faces_uv = loadMeshFaces(mesh_file_path);
    Eigen::MatrixXd r(num_part, 3);
    Eigen::MatrixXd n(num_part, 1);  // degree based vector

    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_uv, h_v_mapping_vector};

    init_particle_position(faces_uv, halfedge_uv, num_part, r, n);

    Eigen::VectorXd vertices_3D_active_eigen = get_vertice_id(r, halfedge_uv, h_v_mapping_vector);
    std::vector<int> vertices_3D_active(vertices_3D_active_eigen.data(), vertices_3D_active_eigen.data() + vertices_3D_active_eigen.size());

    const Eigen::MatrixXd distance_matrix = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/meshes/data/ellipsoid_x4_distance_matrix_static.csv");

    // Prefill the vertices_2DTissue_map
    auto [splay_state_coord, splay_state_halfedges] = get_splay_state_vertices(faces_uv, halfedge_uv, 3);
    auto splay_state_vertices = get_vertice_id(splay_state_coord, halfedge_uv, h_v_mapping_vector);

    for (int i = 0; i < splay_state_vertices.size(); ++i) {
        int splay_state_v = splay_state_vertices[i];

        auto result_virtual = create_uv_surface_intern("Ellipsoid", splay_state_v);
        std::vector<int64_t> h_v_mapping_vector_virtual = std::get<0>(result_virtual);
        std::string mesh_file_path_virtual = std::get<1>(result_virtual);
        Eigen::MatrixXd halfedge_uv_virtual = loadMeshVertices(mesh_file_path);

        // Store the virtual meshes
        vertices_2DTissue_map[splay_state_v] = Mesh_UV_Struct{splay_state_v, halfedge_uv_virtual, h_v_mapping_vector_virtual};
    }

    // Start the simulation
    std::clock_t start = std::clock();

    Eigen::VectorXd v_order(num_frames);

    for (int tt = 1; tt <= num_frames; ++tt) {
        auto [r_new, r_dot, dist_length, ntest, nr_dot, particles_color] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v_order, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, num_part, vertices_2DTissue_map);
        r = r_new;
        n = ntest;

        auto new_vertices_3D_active_eigen = get_vertice_id(r, halfedge_uv, h_v_mapping_vector);
        std::vector<int> new_vertices_3D_active(new_vertices_3D_active_eigen.data(), new_vertices_3D_active_eigen.data() + new_vertices_3D_active_eigen.size());
        vertices_3D_active = new_vertices_3D_active;

        // std::string file_name = "r_data_" + std::to_string(tt) + ".csv";
        // save_matrix_to_csv(r, file_name);
        // std::string file_name_n = "n_data_" + std::to_string(tt) + ".csv";
        // Eigen::MatrixXi n_int = n.cast<int>();
        // save_matrix_to_csv(n_int, file_name_n);
    }
    std::cout << "Order parameter: " << v_order << std::endl;

    std::clock_t end = std::clock();
    double duration = (end - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    return 0;
}
