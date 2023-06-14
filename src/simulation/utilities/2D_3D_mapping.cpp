// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <Eigen/Dense>
#include <unordered_set>

#include <utilities/2D_3D_mapping.h>
#include <utilities/barycentric_coord.h>
#include <io/mesh_loader.h>
#include <io/csv.h>


// (2D Coordinates -> 3D Coordinates and Their Nearest 3D Vertice id (for the distance calculation on resimulations)) mapping
std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd vertices_uv,
    const Eigen::MatrixXd vertices_3D,
    std::vector<int64_t> h_v_mapping
){
    int num_r = r.rows();
    Eigen::MatrixXd new_3D_points(num_r, 3);
    std::vector<int> nearest_vertices_ids(num_r);

    for (int i = 0; i < num_r; ++i) {
        auto [barycentric_coord, nearest_vertex_id] = calculate_barycentric_3D_coord(r, halfedges_uv, faces_uv, vertices_uv, vertices_3D, h_v_mapping, i);
        new_3D_points.row(i) = barycentric_coord;
        nearest_vertices_ids[i] = nearest_vertex_id;
    }

    return std::make_pair(new_3D_points, nearest_vertices_ids);
}


// (3D Coordinates -> 2D Coordinates and Their Nearest 2D Vertice id) mapping
Eigen::MatrixXd get_r2d(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const Eigen::MatrixXi& faces_uv,
    const Eigen::MatrixXd& vertices_uv,
    const Eigen::MatrixXd& vertices_3D,
    std::vector<int64_t>& h_v_mapping
){
    int num_r = r.rows();
    Eigen::MatrixXd new_2D_points(num_r, 2);

    for (int i = 0; i < num_r; ++i) {
        auto [uv_coord, nearest_vertex_id] = calculate_barycentric_2D_coord(r, halfedges_uv, faces_uv, vertices_uv, vertices_3D, h_v_mapping, i);
        new_2D_points.row(i) = uv_coord;
    }

    return new_2D_points;
}


// (3D Vertice id -> 3D Vertice row position of the h-v map) mapping
std::vector<int> find_vertice_rows_index(
    std::vector<int64_t> h_v_mapping_vector,
    std::vector<int> r3d_vertices
){
    std::unordered_set<int> found_ids;
    std::vector<int> indices;

    for (int i = 0; i < h_v_mapping_vector.size(); ++i) {
        if (std::find(r3d_vertices.begin(), r3d_vertices.end(), h_v_mapping_vector[i]) != r3d_vertices.end()) {
            if (found_ids.find(h_v_mapping_vector[i]) == found_ids.end()) {
                indices.push_back(i);
                found_ids.insert(h_v_mapping_vector[i]);
            }
        }
    }
    return indices;
}


// (3D Vertice row position -> nD Vertice coordinates) mapping
Eigen::MatrixXd get_coordinates(
    std::vector<int> indices,
    Eigen::MatrixXd coord
){
    Eigen::MatrixXd found_coord(indices.size(), coord.cols());

    for (int i = 0; i < indices.size(); ++i) {
        found_coord.row(i) = coord.row(indices[i]);
    }

    return found_coord;
}
