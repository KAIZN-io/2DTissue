// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <Eigen/Dense>

#include <utilities/2D_3D_mapping.h>
#include <utilities/barycentric_coord.h>
#include <io/mesh_loader.h>
#include <io/csv.h>


// (3D Vertice id -> Halfedge id) mapping
std::vector<int64_t> get_first_uv_halfedge_from_3D_vertice_id(
    const std::vector<int64_t>& _vertice_3D_id,
    const std::vector<int64_t>& _halfedge_vertices_mapping
) {
    std::vector<int64_t> halfedge_id;
    halfedge_id.reserve(_vertice_3D_id.size());

    for (const auto& vertice_3D_id : _vertice_3D_id) {
        auto it = std::find(_halfedge_vertices_mapping.begin(), _halfedge_vertices_mapping.end(), vertice_3D_id);
        halfedge_id.push_back(static_cast<int64_t>(std::distance(_halfedge_vertices_mapping.begin(), it)));
    }

    return halfedge_id;
}


// (Halfedge id -> 2D Coordinates) mapping
Eigen::MatrixXd get_r_from_halfedge_id(
    const std::vector<int64_t>& halfedge_id,
    const Eigen::MatrixXd& halfedges_uv
){
    int num_halfedges = static_cast<int>(halfedge_id.size());
    Eigen::MatrixXd halfedge_uv_coord(num_halfedges, halfedges_uv.cols());

    for (int i = 0; i < num_halfedges; i++) {
        halfedge_uv_coord.row(i) = halfedges_uv.row(halfedge_id[i]);
    }

    return halfedge_uv_coord;
}


// (2D Coordinates -> 3D Coordinates and Their Nearest 3D Vertice id) mapping
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
