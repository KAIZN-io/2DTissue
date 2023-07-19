// author: @Jan-Piotraschke
// date: 2023-07-19
// license: Apache License 2.0
// version: 0.1.0

#include <iostream>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <cstdint>

// #include <utilities/2D_3D_mapping.h>
// #include <utilities/nearest_map.h>
// #include <utilities/sim_structs.h>
// #include <utilities/update.h>
// #include <utilities/check_boundary.h>
// #include <utilities/2D_surface.h>
// #include <utilities/check_validity.h>

// #include <particle_simulation/motion.h>

#include <IO.h>
#include <VirtualMesh.h>
#include <Validation.h>

using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;


std::tuple<Eigen::MatrixXd, std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> find_nearest_vertice_map(
    int target_vertex,
    const Eigen::MatrixXd distance_matrix,
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map
) {
    // Get all the availabe 2D maps
    std::vector<int> vertices_2DTissue_map_keys;
    for (auto const& [key, val] : vertices_2DTissue_map) {
        vertices_2DTissue_map_keys.push_back(key);
    }

    // double min_distance = std::numeric_limits<double>::max();
    // int nearest_vertex = -1;

    // for (int vertex : vertices_2DTissue_map_keys) {
    //     double distance = distance_matrix(target_vertex, vertex);
    //     if (distance < min_distance) {
    //         min_distance = distance;
    //         nearest_vertex = vertex;
    //     }
    // }

    double max_distance = std::numeric_limits<double>::lowest();
    int furthest_vertex = -1;

    for (int vertex : vertices_2DTissue_map_keys) {
        double distance = distance_matrix(target_vertex, vertex);
        if (distance > max_distance) {
            max_distance = distance;
            furthest_vertex = vertex;
        }
    }

    Eigen::MatrixXd halfedges_uv;
    std::vector<int64_t> h_v_mapping;
    Eigen::MatrixXd vertices_UV;
    Eigen::MatrixXd vertices_3D;
    std::string mesh_file_path;

    auto it = vertices_2DTissue_map.find(furthest_vertex);
    if (it != vertices_2DTissue_map.end()) {
        // Load the mesh
        halfedges_uv = it->second.mesh;
        h_v_mapping = it->second.h_v_mapping;
        vertices_UV = it->second.vertices_UV;
        vertices_3D = it->second.vertices_3D;
        mesh_file_path = it->second.mesh_file_path;
    }

    return std::tuple(halfedges_uv, h_v_mapping, vertices_UV, vertices_3D, mesh_file_path);
}


std::set<std::pair<int, int>> find_border_edges(const Matrix3Xi& mesh) {
    std::set<std::pair<int, int>> boundary_edges;

    for (int i = 0; i < mesh.rows(); ++i) {
        std::array<std::pair<int, int>, 3> edges = {{
            {mesh(i, 0), mesh(i, 1)},
            {mesh(i, 1), mesh(i, 2)},
            {mesh(i, 2), mesh(i, 0)}
        }};

        for (const auto& edge : edges) {
            auto rev_edge = std::make_pair(edge.second, edge.first);
            auto it = boundary_edges.find(rev_edge);
            if (it != boundary_edges.end()) {
                boundary_edges.erase(it);
            } else {
                boundary_edges.insert(edge);
            }
        }
    }

    return boundary_edges;
}

std::vector<int> order_of_edges_vertices(const Eigen::MatrixXi& border_edges_array) {
    std::vector<int> border_edges_vec(border_edges_array.rows());
    int start = border_edges_array(0, 1);

    for (int i = 0; i < border_edges_array.rows(); ++i) {
        int next = -1;
        for (int j = 0; j < border_edges_array.rows(); ++j) {
            if (border_edges_array(j, 0) == start) {
                next = j;
                break;
            }
        }
        start = border_edges_array(next, 1);
        border_edges_vec[i] = border_edges_array(next, 0);
    }

    return border_edges_vec;
}

std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 2>, std::vector<int>> get_splay_state_vertices(const Matrix3Xi& mesh_loaded_uv,
                                                      const Eigen::MatrixXd& halfedges_uv,
                                                      int modula_mode) {
    std::set<std::pair<int, int>> border_edges = find_border_edges(mesh_loaded_uv);

    Eigen::MatrixXi border_edges_array(border_edges.size(), 2);
    int idx = 0;
    for (const auto& edge : border_edges) {
        border_edges_array(idx, 0) = edge.first;
        border_edges_array(idx, 1) = edge.second;
        ++idx;
    }

    std::vector<int> border_edges_order = order_of_edges_vertices(border_edges_array);

    std::vector<int> border_vertices;
    for (size_t i = 0; i < border_edges_order.size(); ++i) {
        if ((i + 1) % modula_mode == 0) {
            border_vertices.push_back(border_edges_order[i]);
        }
    }

    Eigen::MatrixXd splay_state_vertices_test(border_vertices.size(), 3);
    for (size_t i = 0; i < border_vertices.size(); ++i) {
        splay_state_vertices_test.row(i) = halfedges_uv.row(border_vertices[i]);
    }

    Eigen::Matrix<double, Eigen::Dynamic, 2> splay_state_vertices;
    // copy the first two columns of splay_state_vertices_test for each row into splay_state_vertices
    splay_state_vertices = splay_state_vertices_test.block(0, 0, splay_state_vertices_test.rows(), 2);

    return {splay_state_vertices, border_vertices};
}


std::vector<int> get_3D_splay_vertices(
    Eigen::MatrixXd distance_matrix,
    int number_vertices
){
    std::vector<int> selected_vertices;

    // Start at a random vertex
    int start_vertex = rand() % distance_matrix.rows();
    selected_vertices.push_back(start_vertex);

    while (selected_vertices.size() < number_vertices) {
        double max_distance = -1;
        int furthest_vertex = -1;

        // Find the vertex that is furthest away from all selected vertices
        for (int i = 0; i < distance_matrix.rows(); ++i) {
            double min_distance = std::numeric_limits<double>::infinity();

            // Compute the minimum distance to the selected vertices
            for (int j : selected_vertices) {
                double distance = distance_matrix(i, j);

                if (distance < min_distance) {
                    min_distance = distance;
                }
            }

            // If this vertex is further away than the current furthest vertex, update it
            if (min_distance > max_distance) {
                max_distance = min_distance;
                furthest_vertex = i;
            }
        }

        // Add the furthest vertex to the selected vertices
        selected_vertices.push_back(furthest_vertex);
    }

    return selected_vertices;
}





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

std::vector<VertexData> update_vertex_data(
    const Eigen::MatrixXd& old_r_3D_coord,
    const Eigen::MatrixXd& new_r_3D_coord,
    const std::vector<int>& inside_uv_ids,
    int start_id
){
    int num_r = old_r_3D_coord.rows();
    std::vector<VertexData> vertex_data(num_r);

    // Initialize the vertex data
    for (int i = 0; i < num_r; ++i) {
        VertexData& vd = vertex_data[i];

        vd.old_particle_pos = old_r_3D_coord.row(i);
        vd.next_particle_pos = old_r_3D_coord.row(i);
        vd.valid = false;
        vd.uv_mesh_id = start_id;
    }

    // Update the vertex data based on inside_uv_ids
    for (int i : inside_uv_ids) {

        if (!vertex_data[i].valid) {
            // Get the vertex data
            // ? VertexData& vd = vertex_data[inside_uv_ids[i]];
            VertexData& vd = vertex_data[i];

            vd.next_particle_pos = new_r_3D_coord.row(i);
            vd.uv_mesh_id = start_id;
            vd.valid = true;
        }
    }
    return vertex_data;
}


void update_if_valid(
    std::vector<VertexData>& vertex_data,
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_coord,
    const Eigen::MatrixXd& r_3D_coord,
    int start_id
){
    // Find out which particles are inside the mesh
    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_UV_coord);

    for (int i : inside_uv_ids) {
        if (!vertex_data[i].valid) {
            VertexData& vd = vertex_data[i];

            vd.next_particle_pos = r_3D_coord.row(i);
            vd.uv_mesh_id = start_id;
            vd.valid = true;
        }
    }
}

// Check if the given point r is inside the UV parametrization bounds
bool is_inside_uv(const Eigen::Vector2d& r) {
    return (0 <= r[0] && r[0] <= 1) && (0 <= r[1] && r[1] <= 1);
}

std::vector<int> find_inside_uv_vertices_id(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r) {
    int nrows = r.rows();
    std::vector<int> inside_id;

    for (int i = 0; i < nrows; ++i) {
        // Check if the point is inside the UV parametrization bounds
        Eigen::Vector2d first_two_columns = r.row(i).head<2>();
        if (is_inside_uv(first_two_columns)) {
            inside_id.push_back(i);
        }
    }

    return inside_id;
}


// class UVVertices {
// public:
//     UVVertices(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r)
//         : r_(r)
//     {
//         int nrows = r.rows();
//         for (int i = 0; i < nrows; ++i) {
//             Eigen::Vector2d first_two_columns = r.row(i).head<2>();
//             if (is_inside_uv(first_two_columns)) {
//                 inside_uv_ids_.insert(i);
//             }
//             else {
//                 outside_uv_ids_.insert(i);
//             }
//         }
//     }

//     const std::set<int>& get_inside_uv_ids() const {
//         return inside_uv_ids_;
//     }

//     const std::set<int>& get_outside_uv_ids() const {
//         return outside_uv_ids_;
//     }

// private:
//     // Check if the given point r is inside the UV parametrization bounds
//     static bool is_inside_uv(const Eigen::Vector2d& r) {
//         return (0 <= r[0] && r[0] <= 1) && (0 <= r[1] && r[1] <= 1);
//     }

//     Eigen::Matrix<double, Eigen::Dynamic, 2> r_;
//     std::set<int> inside_uv_ids_;
//     std::set<int> outside_uv_ids_;
// };
