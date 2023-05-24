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


double pointTriangleDistance(const Eigen::Vector3d& p, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) {
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;
    Eigen::Vector3d ap = p - a;

    double d1 = ab.dot(ap);
    double d2 = ac.dot(ap);

    if (d1 <= 0.0 && d2 <= 0.0)
        return ap.norm();

    Eigen::Vector3d bp = p - b;
    double d3 = ab.dot(bp);
    double d4 = ac.dot(bp);

    if (d3 >= 0.0 && d4 <= d3)
        return bp.norm();

    double vc = d1*d4 - d3*d2;

    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        return (a + v * ab - p).norm();
    }

    Eigen::Vector3d cp = p - c;
    double d5 = ab.dot(cp);
    double d6 = ac.dot(cp);

    if (d6 >= 0.0 && d5 <= d6)
        return cp.norm();

    double vb = d5*d2 - d1*d6;

    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        return (a + w * ac - p).norm();
    }

    double va = d3*d6 - d5*d4;

    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (b + w * (c - b) - p).norm();
    }

    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;

    return (a + ab * v + ac * w - p).norm();
}


int closestRow(const Eigen::MatrixXd& vertices_uv_test, const Eigen::Vector2d& halfedge_coord) {
    Eigen::VectorXd dists(vertices_uv_test.rows());
    for (int i = 0; i < vertices_uv_test.rows(); ++i) {
        dists[i] = (vertices_uv_test.row(i) - halfedge_coord.transpose()).squaredNorm();
    }

    Eigen::VectorXd::Index minRow;
    dists.minCoeff(&minRow);

    return minRow;
}


Eigen::MatrixXd get_r3d(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const std::vector<int64_t>& halfedge_vertices_mapping,
    const Eigen::MatrixXi faces_uv)
{
    int num_r = r.rows();
    Eigen::MatrixXd vertices_3D = loadMeshVertices("/Users/jan-piotraschke/git_repos/2DTissue/meshes/ellipsoid_x4.off");
    Eigen::MatrixXd vertices_3D_test = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/vertice_3D_data.csv");
    Eigen::MatrixXd vertices_uv_test = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/vertice_uv_data.csv");

    Eigen::MatrixXd new_3D_points(num_r, 3);

    for (int i = 0; i < num_r; ++i) {
        std::vector<std::pair<double, int>> distances(faces_uv.rows());

        for (int j = 0; j < faces_uv.rows(); ++j) {
            Eigen::Vector3d uv_a = halfedges_uv.row(faces_uv(j, 0));
            Eigen::Vector3d uv_b = halfedges_uv.row(faces_uv(j, 1));
            Eigen::Vector3d uv_c = halfedges_uv.row(faces_uv(j, 2));

            distances[j] = {pointTriangleDistance(r.row(i).head(2), uv_a, uv_b, uv_c), j};
        }

        std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

        int halfedge_a = faces_uv(min_distance.second, 0);
        int halfedge_b = faces_uv(min_distance.second, 1);
        int halfedge_c = faces_uv(min_distance.second, 2);

        Eigen::Vector2d halfedge_a_coord = halfedges_uv.row(halfedge_a);
        Eigen::Vector2d halfedge_b_coord = halfedges_uv.row(halfedge_b);
        Eigen::Vector2d halfedge_c_coord = halfedges_uv.row(halfedge_c);

        // Inside your loop...
        int closest_a = closestRow(vertices_uv_test, halfedge_a_coord);
        int closest_b = closestRow(vertices_uv_test, halfedge_b_coord);
        int closest_c = closestRow(vertices_uv_test, halfedge_c_coord);

        // Get the 3D coordinates of the 3 halfedges
        Eigen::Vector3d a = vertices_3D_test.row(closest_a);
        Eigen::Vector3d b = vertices_3D_test.row(closest_b);
        Eigen::Vector3d c = vertices_3D_test.row(closest_c);

        // Compute the weights (distances in UV space)
        double w_a = (r.row(i).head(2).transpose() - halfedge_a_coord).norm();
        double w_b = (r.row(i).head(2).transpose() - halfedge_b_coord).norm();
        double w_c = (r.row(i).head(2).transpose() - halfedge_c_coord).norm();

        // Compute the barycentric coordinates
        double sum_weights = w_a + w_b + w_c;
        w_a /= sum_weights;
        w_b /= sum_weights;
        w_c /= sum_weights;

        // Compute the new 3D point using the barycentric coordinates
        Eigen::Vector3d newPoint = w_a * a + w_b * b + w_c * c;

        // Store the new 3D point
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
