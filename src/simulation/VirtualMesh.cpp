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


Eigen::VectorXd VirtualMesh::get_n_orientation(Eigen::Matrix<double, Eigen::Dynamic, 2> position_, Eigen::Vector2d northPole_, Eigen::VectorXd n_pole_) {
    return compass.assign_n(position_, northPole_, n_pole_);
}

void VirtualMesh::assign_particle_orientation(Eigen::Vector2d northPole_, Eigen::VectorXd n_pole_){
    n = compass.assign_n(r_UV, northPole_, n_pole_);
}


Eigen::VectorXd VirtualMesh::get_relative_orientation(){
    return compass.calculate_n_pole(r_UV, n);
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

void VirtualMesh::change_UV_map(int target_vertex){
    auto it = vertices_2DTissue_map.find(target_vertex);
    if (it != vertices_2DTissue_map.end()) {
        // Load the mesh
        halfedge_UV = it->second.mesh;
        h_v_mapping = it->second.h_v_mapping;
        face_UV = it->second.face_UV;
        vertice_UV = it->second.vertices_UV;
    }
}
