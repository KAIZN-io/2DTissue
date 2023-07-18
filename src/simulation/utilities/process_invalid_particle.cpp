// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include <IO.h>
#include <utilities/2D_3D_mapping.h>
#include <utilities/nearest_map.h>
#include <utilities/sim_structs.h>
#include <utilities/splay_state.h>
#include <utilities/update.h>
#include <utilities/check_boundary.h>
#include <utilities/2D_surface.h>
#include <utilities/check_validity.h>

#include <particle_simulation/motion.h>
#include <utilities/process_invalid_particle.h>


void process_invalid_particle(
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    int old_id,
    std::vector<int> old_ids,
    std::vector<VertexData>& vertex_struct,
    const VertexData& particle,
    int num_part,
    const Eigen::MatrixXd& distance_matrix,
    Eigen::MatrixXd& n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double current_step
) {
    // Get the nearest vertice map
    auto [halfedges_uv, h_v_mapping, vertices_UV, vertices_3D, mesh_file_path] = find_nearest_vertice_map(old_id, distance_matrix, vertices_2DTissue_map);

    // Find the new row indices of the used vertices
    auto row_indices = find_vertice_rows_index(h_v_mapping, old_ids);

    // Get the coordinates of the vertices based on the row indices
    auto r_active = get_coordinates(row_indices, vertices_UV);

    // Simulate the flight of the particle
    auto [r_UV_virtual, r_dot, dist_length] = simulate_flight(r_active, n, old_ids, distance_matrix, v0, k, σ, μ, r_adh, k_adh, dt);

    // Get the new vertice id
    Eigen::MatrixXi faces_uv;
    loadMeshFaces(mesh_file_path, faces_uv);

    // Map them to the 3D coordinates
    auto [r_3D_virtual, vertices_3D_active] = get_r3d(r_UV_virtual, halfedges_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

    update_if_valid(vertex_struct, r_UV_virtual, r_3D_virtual, old_id);
}


void process_if_not_valid(
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    std::vector<int> old_vertices_3D,
    std::vector<VertexData>& vertex_struct,
    int num_part,
    Eigen::MatrixXd& distance_matrix_v,
    Eigen::MatrixXd& n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double step_size,
    double current_step
) {
    std::vector<int> invalid_ids;

    for (int i = 0; i < vertex_struct.size(); ++i) {
        if (!vertex_struct[i].valid) {
            invalid_ids.push_back(i);
        }
    }

    for (int invalid_id : invalid_ids) {
        process_invalid_particle(vertices_2DTissue_map, invalid_id, old_vertices_3D, vertex_struct, vertex_struct[invalid_id], num_part, distance_matrix_v, n, v0, k, v0_next, k_next, σ, μ, r_adh, k_adh, step_size, current_step);

        if (are_all_valid(vertex_struct)) {
            break;
        }
    }

    std::vector<int> still_invalid_ids;

    for (int i = 0; i < vertex_struct.size(); ++i) {
        if (!vertex_struct[i].valid) {
            still_invalid_ids.push_back(i);
        }
    }
    if (still_invalid_ids.size() > 0) {

        // Create new 2D surfaces for the still invalid ids
        for (int i = 0; i < still_invalid_ids.size(); ++i) {

            int invalid_particle = still_invalid_ids[i];

            Eigen::VectorXd particle_distance = distance_matrix_v.row(invalid_particle);  // 0-based indexing, so the "fifth" row is at index 4

            Eigen::VectorXd::Index maxIndex;
            double max_distance = particle_distance.maxCoeff(&maxIndex);

            // transfrom the maxIndex to int
            int maxIndex_int = static_cast<int>(maxIndex);

            std::cout << "Creating new 2D surface for particle " << invalid_particle << " with the distance " << max_distance << " in the row " << maxIndex_int << std::endl;

            // because it is on the Seam Edge line of its own mesh !!
            auto [h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path] = create_uv_surface("Ellipsoid", maxIndex_int);
            Eigen::MatrixXd halfedge_uv;
            loadMeshVertices(mesh_file_path, halfedge_uv);

            // Store the new meshes
            vertices_2DTissue_map[maxIndex_int] = Mesh_UV_Struct{maxIndex_int, halfedge_uv, h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path};
            process_invalid_particle(vertices_2DTissue_map, maxIndex_int, old_vertices_3D, vertex_struct, vertex_struct[maxIndex_int], num_part, distance_matrix_v, n, v0, k, v0_next, k_next, σ, μ, r_adh, k_adh, step_size, current_step);

            if (are_all_valid(vertex_struct)) {
                break;
            }
        }
    
    }
}
