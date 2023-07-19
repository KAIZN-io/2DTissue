// // VirtualMesh.h

// #pragma once

// #include <tuple>
// #include <vector>
// #include <Eigen/Dense>
// #include <cstdint>
// #include <map>

// using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;

// std::tuple<Eigen::MatrixXd, std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> find_nearest_vertice_map(
//     int target_vertex,
//     const Eigen::MatrixXd distance_matrix,
//     std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map
// );

// void process_if_not_valid(
//     std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
//     std::vector<int> old_vertices_3D,
//     std::vector<VertexData>& vertex_data,
//     int num_part,
//     Eigen::MatrixXd& distance_matrix_v,
//     Eigen::VectorXd& n,
//     double v0,
//     double k,
//     double k_next,
//     double v0_next,
//     double σ,
//     double μ,
//     double r_adh,
//     double k_adh,
//     double dt,
//     double tt
// );

// std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 2>, std::vector<int>> get_splay_state_vertices(
//     const Matrix3Xi& mesh_loaded_uv,
//     const Eigen::MatrixXd& halfedges_uv,
//     int modula_mode = 10
// );

// std::vector<int> get_3D_splay_vertices(
//     Eigen::MatrixXd distance_matrix,
//     int modula_mode
// );

// std::vector<VertexData> update_vertex_data(
//     const Eigen::MatrixXd& old_r_3D_coord,
//     const Eigen::MatrixXd& new_r_3D_coord,
//     const std::vector<int>& inside_uv_ids,
//     int start_id
// );

// void update_if_valid(
//     std::vector<VertexData>& vertex_data,
//     const Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_coord,
//     const Eigen::MatrixXd& r_3D_coord,
//     int start_id
// );

// std::vector<int> find_inside_uv_vertices_id(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);
