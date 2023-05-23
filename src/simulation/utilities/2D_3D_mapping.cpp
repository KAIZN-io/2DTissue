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
#include <io/mesh_loader.h>


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


// (2D Coordinates -> 3D Coordinates) mapping
Eigen::MatrixXd get_r3d(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const std::vector<int64_t>& halfedge_vertices_mapping
){
    int num_r = r.rows();
    Eigen::MatrixXd vertices_3D = loadMeshVertices("/Users/jan-piotraschke/git_repos/2DTissue/meshes/ellipsoid_x4.off");
    Eigen::MatrixXd vertices_3D_ids(num_r, 3);
    Eigen::MatrixXd distances_3D(num_r, 3);

    for (int i = 0; i < num_r; ++i) {
        // Use a vector of pairs: distance and index
        std::vector<std::pair<double, int>> distances(halfedges_uv.rows());

        // Compute all distances
        for (int j = 0; j < halfedges_uv.rows(); ++j) {
            Eigen::VectorXd diff = halfedges_uv.row(j) - r.row(i);
            distances[j] = {diff.norm(), j};
        }

        // Partially sort the distances to find the three smallest
        std::partial_sort(distances.begin(), distances.begin() + 3, distances.end());

        // Save the vertex IDs and distances of the three closest half-edges
        for (int j = 0; j < 3; ++j) {
            vertices_3D_ids(i, j) = halfedge_vertices_mapping[distances[j].second];
            distances_3D(i, j) = distances[j].first;
        }
    }

    // init empty Eigen::MatrixXd for the new points
    Eigen::MatrixXd new_3D_points(num_r, 3);

    // iterate over the rows of vertices_3D_active_eigen_test and access all 3 column values by asigning them to new variables
    for (int i = 0; i < vertices_3D_ids.rows(); ++i) {
        // assign the distances to the weights
        double w1 = distances_3D(i, 0);
        double w2 = distances_3D(i, 1);
        double w3 = distances_3D(i, 2);

        int v1 = vertices_3D_ids(i, 0);
        int v2 = vertices_3D_ids(i, 1);
        int v3 = vertices_3D_ids(i, 2);

        Eigen::Vector3d v1_3D = vertices_3D.row(v1);
        Eigen::Vector3d v2_3D = vertices_3D.row(v2);
        Eigen::Vector3d v3_3D = vertices_3D.row(v3);

        // Calculate the new point
        Eigen::Vector3d newPoint = (w1 * v1_3D + w2 * v2_3D + w3 * v3_3D) / (w1 + w2 + w3);

        // Assign the new point to the new_points matrix
        new_3D_points.row(i) = newPoint;
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
