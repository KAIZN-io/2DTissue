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


// (2D Coordinates -> 3D Vertice id) mapping
Eigen::VectorXd get_vertice_id(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const std::vector<int64_t>& halfedge_vertices_mapping
){
    int num_r = r.rows();
    Eigen::VectorXd vertice_3D_id(num_r);

    for (int i = 0; i < num_r; ++i) {
        double min_distance = std::numeric_limits<double>::max();
        int64_t min_idx = -1;

        for (int j = 0; j < halfedges_uv.rows(); ++j) {
            Eigen::VectorXd diff = halfedges_uv.row(j) - r.row(i);
            double distance = diff.norm();

            if (distance < min_distance) {
                min_distance = distance;
                min_idx = j;
            }
        }

        vertice_3D_id(i) = halfedge_vertices_mapping[min_idx];
    }

    return vertice_3D_id;
}


Eigen::MatrixXd get_r3d(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const std::vector<int64_t>& halfedge_vertices_mapping,
    const Eigen::MatrixXi faces_uv
){
    int num_r = r.rows();
    Eigen::MatrixXd vertices_3D = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/vertice_3D_data.csv");
    Eigen::MatrixXd vertices_uv = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/vertice_uv_data.csv");

    Eigen::MatrixXd new_3D_points(num_r, 3);

    for (int i = 0; i < num_r; ++i) {

        new_3D_points.row(i) = calculate_barycentric_3D_coord(r, halfedges_uv, faces_uv, vertices_uv, vertices_3D, i);

    }

    return new_3D_points;
}


// // (3D Coordinates -> 3D Vertice id) mapping
// Eigen::VectorXd get_vertice3D_id(
//     const Eigen::MatrixXd r3d,
//     const Eigen::MatrixXd vertices_3D
// ){
//     int num_r = r3d.rows();
//     Eigen::VectorXd vertice_3D_id(num_r);

//     for (int i = 0; i < num_r; ++i) {
//         double min_distance = std::numeric_limits<double>::max();
//         int64_t min_idx = -1;

//         for (int j = 0; j < vertices_3D.rows(); ++j) {
//             Eigen::Vector3d diff = vertices_3D.row(j) - r3d.row(i).transpose();
//             double distance = diff.norm();

//             if (distance < min_distance) {
//                 min_distance = distance;
//                 min_idx = j;
//             }
//         }

//         vertice_3D_id(i) = min_idx;
//     }

//     return vertice_3D_id;
// }
