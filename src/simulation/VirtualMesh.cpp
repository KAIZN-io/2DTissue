// author: @Jan-Piotraschke
// date: 2023-07-31
// license: Apache License 2.0
// version: 0.1.1

#include <VirtualMesh.h>

VirtualMesh::VirtualMesh(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::MatrixXd& r_3D,
    Eigen::MatrixXd& halfedge_UV,
    Eigen::MatrixXi& face_UV,
    Eigen::MatrixXd& vertice_UV,
    std::vector<int64_t>& h_v_mapping,
    int particle_count,
    Eigen::VectorXd& n,
    Eigen::MatrixXi& face_3D,
    Eigen::MatrixXd& vertice_3D,
    Eigen::MatrixXd& distance_matrix,
    std::string mesh_path,
    int map_cache_count,
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    std::unique_ptr<GeometryProcessing> geometry_ptr,
    std::unique_ptr<Validation> validation_ptr
)
    : r_UV(r_UV),
      r_UV_old(r_UV_old),
      r_3D(r_3D),
      halfedge_UV(halfedge_UV),
      face_UV(face_UV),
      vertice_UV(vertice_UV),
      h_v_mapping(h_v_mapping),
      particle_count(particle_count),
      n(n),
      face_3D(face_3D),
      vertice_3D(vertice_3D),
      distance_matrix(distance_matrix),
      mesh_path(mesh_path),
      map_cache_count(map_cache_count),
      vertices_2DTissue_map(vertices_2DTissue_map),
      geometry_ptr(std::move(geometry_ptr)),
      validation_ptr(std::move(validation_ptr)),
      cell(particle_count, halfedge_UV, face_UV, face_3D, vertice_UV, vertice_3D, h_v_mapping, r_UV, r_3D, n),
      compass(northPole)
{
}

Eigen::Vector2d VirtualMesh::init_north_pole(){
    northPole_3D.resize(1, 3);
    northPole_3D << 1.7466, -0.220152, -2.40228;
    // Get the old 3D coordinates
    auto [new_r_3D, new_vertices_3D_active] = cell.get_r3d();
    r_3D = new_r_3D;
    r_3D.conservativeResize(r_3D.rows() + 1, 3);
    r_3D.row(r_3D.rows() - 1) = northPole_3D;

    // Map it to the 2D coordinates of the new map
    auto new_r_UV = cell.get_r2d();

    northPole = new_r_UV.row(new_r_UV.rows() - 1);  // The north pole (= v1 of the basic UV mesh)
    r_3D.conservativeResize(r_3D.rows() - 1, 3);

    return northPole;
}

std::vector<int> VirtualMesh::get_3D_splay_vertices(){
    std::vector<int> selected_vertices;

    // Start at a random vertex
    int start_vertex = rand() % distance_matrix.rows();
    selected_vertices.push_back(start_vertex);

    while (selected_vertices.size() < map_cache_count) {
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


void VirtualMesh::generate_virtual_mesh()
{
    auto splay_state_vertices_id = get_3D_splay_vertices();

    for (int i = 0; i < splay_state_vertices_id.size(); ++i) {
        int splay_state_v = splay_state_vertices_id[i];
        Eigen::MatrixXi face_UV_virtual;
        Eigen::MatrixXd halfedge_UV_virtual;

        auto [h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual] = geometry_ptr->create_uv_surface(mesh_path, splay_state_v);
        loadMeshVertices(mesh_file_path_virtual, halfedge_UV_virtual);
        loadMeshFaces(mesh_file_path_virtual, face_UV_virtual);

        // Store the virtual meshes
        vertices_2DTissue_map[splay_state_v] = Mesh_UV_Struct{splay_state_v, halfedge_UV_virtual, h_v_mapping_virtual, face_UV_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual};
    }
}


void VirtualMesh::prepare_virtual_mesh(int old_id) {
    // get the old UV coordinates
    r_UV = r_UV_old;

    // Get the old 3D coordinates
    auto [new_r_3D, new_vertices_3D_active] = cell.get_r3d();
    r_3D = new_r_3D;

    // Get the relative orientation
    auto n_pole = get_relative_orientation();

    // Find the next UV map
    change_UV_map(old_id);

    // Add the particles to the new map
    assign_particle_position();
    assign_particle_orientation(n_pole, northPole_virtual);
}


double VirtualMesh::get_n_orientation(Eigen::Vector2d position_, Eigen::Vector2d northPole_, double n_pole_) {
    return compass.calculate_n_orientation(position_, northPole_, n_pole_) ;
}


Eigen::VectorXd VirtualMesh::get_relative_orientation(){
    return compass.calculateRelativeAngle(r_UV, n);
}


void VirtualMesh::change_UV_map(int target_vertex) {
    // Get all the availabe 2D maps
    std::vector<int> vertices_2DTissue_map_keys;
    for (auto const& [key, val] : vertices_2DTissue_map) {
        vertices_2DTissue_map_keys.push_back(key);
    }

    double max_distance = std::numeric_limits<double>::lowest();
    int furthest_vertex = -1;

    for (int vertex : vertices_2DTissue_map_keys) {
        double distance = distance_matrix(target_vertex, vertex);
        if (distance > max_distance) {
            max_distance = distance;
            furthest_vertex = vertex;
        }
    }
    load_UV_map(furthest_vertex);
}


void VirtualMesh::assign_particle_position(){
    // Add the north pole to the 3D coordinates
    r_3D.conservativeResize(r_3D.rows() + 1, 3);
    r_3D.row(r_3D.rows() - 1) = northPole_3D;

    // Map it to the 2D coordinates of the new map
    r_UV = cell.get_r2d();

    // Remove the north pole from the coordinates
    r_3D.conservativeResize(r_3D.rows() - 1, 3);

    // Get the UV coordinates of the north pole
    northPole_virtual = r_UV.row(r_UV.rows() - 1);  // The north pole (= v1 of the UV mesh)

    // Remove the North Pole from the coordinates
    r_UV.conservativeResize(r_UV.rows() - 1, 3);
}


void VirtualMesh::assign_particle_orientation(Eigen::VectorXd n_pole, Eigen::Vector2d northPole_virtual_test){
    n = compass.assignOrientation(r_UV, n_pole, northPole_virtual_test);
}


void VirtualMesh::load_UV_map(int target_vertex){
    std::string mesh_file_path;

    auto it = vertices_2DTissue_map.find(target_vertex);
    if (it != vertices_2DTissue_map.end()) {
        // Load the mesh
        halfedge_UV = it->second.mesh;
        h_v_mapping = it->second.h_v_mapping;
        face_UV = it->second.face_UV;
        vertice_UV = it->second.vertices_UV;
    }
}



// void process_invalid_particle(
//     std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
//     int old_id,
//     std::vector<int> old_ids,
//     std::vector<VertexData>& vertex_struct,
//     const VertexData& particle,
//     int num_part,
//     const Eigen::MatrixXd& distance_matrix,
//     Eigen::MatrixXd& n,
//     double v0,
//     double k,
//     double k_next,
//     double v0_next,
//     double σ,
//     double μ,
//     double r_adh,
//     double k_adh,
//     double dt,
//     double current_step
// ) {
//     // Get the nearest vertice map
//     auto [halfedges_uv, h_v_mapping, vertices_UV, vertices_3D, mesh_file_path] = change_UV_map(old_id, distance_matrix, vertices_2DTissue_map);

//     // Find the new row indices of the used vertices
//     auto row_indices = find_vertice_rows_index(h_v_mapping, old_ids);

//     // Get the coordinates of the vertices based on the row indices
//     auto r_active = get_coordinates(row_indices, vertices_UV);

//     // Simulate the flight of the particle
//     auto [r_UV_virtual, r_dot, dist_length] = simulate_flight(r_active, n, old_ids, distance_matrix, v0, k, σ, μ, r_adh, k_adh, dt);

//     // Get the new vertice id
//     Eigen::MatrixXi faces_uv;
//     loadMeshFaces(mesh_file_path, faces_uv);

//     // Map them to the 3D coordinates
//     auto [r_3D_virtual, vertices_3D_active] = get_r3d(r_UV_virtual, halfedges_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

//     update_if_valid(vertex_struct, r_UV_virtual, r_3D_virtual, old_id);
// }


// void process_if_not_valid(
//     std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
//     std::vector<int> old_vertices_3D,
//     std::vector<VertexData>& vertex_struct,
//     int num_part,
//     Eigen::MatrixXd& distance_matrix_v,
//     Eigen::MatrixXd& n,
//     double v0,
//     double k,
//     double k_next,
//     double v0_next,
//     double σ,
//     double μ,
//     double r_adh,
//     double k_adh,
//     double step_size,
//     double current_step
// ) {
//     std::vector<int> invalid_ids;

//     for (int i = 0; i < vertex_struct.size(); ++i) {
//         if (!vertex_struct[i].valid) {
//             invalid_ids.push_back(i);
//         }
//     }

//     for (int invalid_id : invalid_ids) {
//         process_invalid_particle(vertices_2DTissue_map, invalid_id, old_vertices_3D, vertex_struct, vertex_struct[invalid_id], num_part, distance_matrix_v, n, v0, k, v0_next, k_next, σ, μ, r_adh, k_adh, step_size, current_step);

//         if (are_all_valid(vertex_struct)) {
//             break;
//         }
//     }

//     std::vector<int> still_invalid_ids;

//     for (int i = 0; i < vertex_struct.size(); ++i) {
//         if (!vertex_struct[i].valid) {
//             still_invalid_ids.push_back(i);
//         }
//     }
//     if (still_invalid_ids.size() > 0) {

//         // Create new 2D surfaces for the still invalid ids
//         for (int i = 0; i < still_invalid_ids.size(); ++i) {

//             int invalid_particle = still_invalid_ids[i];

//             Eigen::VectorXd particle_distance = distance_matrix_v.row(invalid_particle);  // 0-based indexing, so the "fifth" row is at index 4

//             Eigen::VectorXd::Index maxIndex;
//             double max_distance = particle_distance.maxCoeff(&maxIndex);

//             // transfrom the maxIndex to int
//             int maxIndex_int = static_cast<int>(maxIndex);

//             std::cout << "Creating new 2D surface for particle " << invalid_particle << " with the distance " << max_distance << " in the row " << maxIndex_int << std::endl;

//             // because it is on the Seam Edge line of its own mesh !!
//             auto [h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path] = create_uv_surface("Ellipsoid", maxIndex_int);
//             Eigen::MatrixXd halfedge_uv;
//             loadMeshVertices(mesh_file_path, halfedge_uv);

//             // Store the new meshes
//             vertices_2DTissue_map[maxIndex_int] = Mesh_UV_Struct{maxIndex_int, halfedge_uv, h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path};
//             process_invalid_particle(vertices_2DTissue_map, maxIndex_int, old_vertices_3D, vertex_struct, vertex_struct[maxIndex_int], num_part, distance_matrix_v, n, v0, k, v0_next, k_next, σ, μ, r_adh, k_adh, step_size, current_step);

//             if (are_all_valid(vertex_struct)) {
//                 break;
//             }
//         }
//     }
// }


