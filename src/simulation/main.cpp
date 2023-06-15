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
#include <utilities/barycentric_coord.h>


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
    int num_frames = 1;

    static std::unordered_map<int, Mesh_UV_Struct> vertices_2DTissue_map;

    // Initialize the simulation
    auto [h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path] = create_uv_surface_intern("Ellipsoid", 0);

    // Check if the distance matrix of the static 3D mesh already exists
    if (!std::filesystem::exists("/Users/jan-piotraschke/git_repos/2DTissue/meshes/data/ellipsoid_x4_distance_matrix_static.csv")) {

        // Calculate the distance matrix of the static 3D mesh
        get_all_distances();
    }

    Eigen::MatrixXd halfedge_uv = loadMeshVertices(mesh_file_path);
    Eigen::MatrixXi faces_uv = loadMeshFaces(mesh_file_path);

    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_uv, h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path};

    // Load the distance matrix
    const Eigen::MatrixXd distance_matrix = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/meshes/data/ellipsoid_x4_distance_matrix_static.csv");

    /*
    Prefill the vertices_2DTissue_map with the virtual meshes
    */
    // Get the vertices that are selected for the splay state in 3D
    auto splay_state_vertices_id = get_3D_splay_vertices(distance_matrix, 30);

    // auto [splay_state_UV_coord, splay_state_halfedges] = get_splay_state_vertices(faces_uv, halfedge_uv, 3);
    // auto [splay_state_3D_coord, splay_state_vertices_id] = get_r3d(splay_state_UV_coord, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping_vector);

    for (int i = 0; i < splay_state_vertices_id.size(); ++i) {
        int splay_state_v = splay_state_vertices_id[i];

        auto [h_v_mapping_vector_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual] = create_uv_surface_intern("Ellipsoid", splay_state_v);
        Eigen::MatrixXd halfedge_uv_virtual = loadMeshVertices(mesh_file_path_virtual);

        // Store the virtual meshes
        vertices_2DTissue_map[splay_state_v] = Mesh_UV_Struct{splay_state_v, halfedge_uv_virtual, h_v_mapping_vector_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual};
    }

    std::string mesh_3D_file_path = "/Users/jan-piotraschke/git_repos/2DTissue/meshes/ellipsoid_x4.off";

    /*
    Run the simulation
    */
    for (int num_part = 1; num_part <= 1; num_part += 100) {

        // Repeat the loop 5 times for each num_part
        for (int repeat = 0; repeat < 1; ++repeat) {

            // Initialize the particles in 2D
            Eigen::MatrixXd r(num_part, 3);
            Eigen::MatrixXd n(num_part, 1);
            init_particle_position(faces_uv, halfedge_uv, num_part, r, n);

            // Map the 2D coordinates to their 3D vertices counterparts
            auto [start_3D_points, vertices_3D_active] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping_vector);
            std::cout << "start_3D_points: " << start_3D_points << std::endl;
            std::cout << "r: " << r << std::endl;

            // Find the nearest point for each start_3D_points in vertices_3D

            for (int i = 0;  i < start_3D_points.rows(); ++i)
            {
                Eigen::MatrixXi faces_3D_static = loadMeshFaces(mesh_3D_file_path);
                std::vector<std::pair<double, int>> distances(faces_3D_static.rows());

                for (int j = 0; j < faces_3D_static.rows(); ++j) {

                    // find first occurence of vertice_a in h_v_mapping_vector and get the row index
                    int halfedge_a = std::find(h_v_mapping_vector.begin(), h_v_mapping_vector.end(), faces_3D_static(j, 0)) - h_v_mapping_vector.begin();
                    int halfedge_b = std::find(h_v_mapping_vector.begin(), h_v_mapping_vector.end(), faces_3D_static(j, 1)) - h_v_mapping_vector.begin();
                    int halfedge_c = std::find(h_v_mapping_vector.begin(), h_v_mapping_vector.end(), faces_3D_static(j, 2)) - h_v_mapping_vector.begin();

                    Eigen::Vector3d uv_a = vertices_3D.row(halfedge_a);
                    Eigen::Vector3d uv_b = vertices_3D.row(halfedge_b);
                    Eigen::Vector3d uv_c = vertices_3D.row(halfedge_c);

                    distances[j] = {pointTriangleDistance(start_3D_points.row(i).head(2), uv_a, uv_b, uv_c), j};
                }

                std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

                int vertice_a = faces_3D_static(min_distance.second, 0);
                int vertice_b = faces_3D_static(min_distance.second, 1);
                int vertice_c = faces_3D_static(min_distance.second, 2);

                int nearest_halfedge_a = std::find(h_v_mapping_vector.begin(), h_v_mapping_vector.end(), vertice_a) - h_v_mapping_vector.begin();
                int nearest_halfedge_b = std::find(h_v_mapping_vector.begin(), h_v_mapping_vector.end(), vertice_b) - h_v_mapping_vector.begin();
                int nearest_halfedge_c = std::find(h_v_mapping_vector.begin(), h_v_mapping_vector.end(), vertice_c) - h_v_mapping_vector.begin();

                // Get the 3D coordinates of the 3 halfedges
                Eigen::Vector3d a = vertices_3D.row(nearest_halfedge_a);
                Eigen::Vector3d b = vertices_3D.row(nearest_halfedge_b);
                Eigen::Vector3d c = vertices_3D.row(nearest_halfedge_c);

            }




            // auto r_UV_check = get_r2d(start_3D_points, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping_vector);
            // std::cout << "r_UV_check: " << r_UV_check << std::endl;

            // // Start the simulation
            // std::clock_t start = std::clock();

            // Eigen::VectorXd v_order(num_frames);

            // for (int tt = 1; tt <= num_frames; ++tt) {

            //     // Simulate the particles on the 2D surface
            //     auto [r_new, r_dot, dist_length, n_new, particles_color] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v_order, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, num_part, vertices_2DTissue_map);
            //     r = r_new;
            //     n = n_new;

            //     // Get the 3D vertices coordinates from the 2D particle position coordinates
            //     auto [new_3D_points, new_vertices_3D_active] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping_vector);
            //     vertices_3D_active = new_vertices_3D_active;

            //     // // Save the data
            //     // std::string file_name = "r_data_" + std::to_string(tt) + ".csv";
            //     // save_matrix_to_csv(r, file_name, num_part);
            //     // std::string file_name_color = "color_data_" + std::to_string(tt) + ".csv";
            //     // save_matrix_to_csv(particles_color, file_name_color, num_part);
            //     // std::string file_name_3D = "r_3D_data_" + std::to_string(tt) + ".csv";
            //     // save_matrix_to_csv(new_3D_points, file_name_3D, num_part);
            // }

            // std::clock_t end = std::clock();
            // double duration = (end - start) / (double) CLOCKS_PER_SEC;
            // std::cout << "Time taken: " << duration << " seconds" << std::endl;
        }
    }

    return 0;
}
