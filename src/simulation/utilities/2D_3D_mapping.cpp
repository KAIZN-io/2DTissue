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

Eigen::Vector3d cartesianToBarycentric(
    const Eigen::VectorXd& P,
    const Eigen::Vector3d& A,
    const Eigen::Vector3d& B,
    const Eigen::Vector3d& C)
{
    Eigen::Vector3d v0 = B - A;
    Eigen::Vector3d v1 = C - A;
    Eigen::Vector3d v2 = P.cast<double>() - A;

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);

    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    return Eigen::Vector3d(u, v, w);
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

        // // Get the 3 halfedges of the UV triangle
        // std::cout << "halfedge_a coordinates: " << halfedge_a_coord.transpose() << std::endl;
        // std::cout << "halfedge_b coordinates: " << halfedge_b_coord.transpose() << std::endl;
        // std::cout << "halfedge_c coordinates: " << halfedge_c_coord.transpose() << std::endl;

        // Inside your loop...
        int closest_a = closestRow(vertices_uv_test, halfedge_a_coord);
        int closest_b = closestRow(vertices_uv_test, halfedge_b_coord);
        int closest_c = closestRow(vertices_uv_test, halfedge_c_coord);

        // Get the 3D coordinates of the 3 halfedges
        Eigen::Vector3d a = vertices_3D_test.row(closest_a);
        Eigen::Vector3d b = vertices_3D_test.row(closest_b);
        Eigen::Vector3d c = vertices_3D_test.row(closest_c);

        std::cout << "a: " << a.transpose() << std::endl;
        std::cout << "b: " << b.transpose() << std::endl;
        std::cout << "c: " << c.transpose() << std::endl;

        // std::cout << "r.row(i): " << r.row(i) << std::endl;

        // Eigen::Vector3d bary_coords = cartesianToBarycentric(r.row(i), a, b, c);
        // Eigen::Vector3d newPoint = bary_coords.x() * a + bary_coords.y() * b + bary_coords.z() * c;
        // new_3D_points.row(i) = newPoint;


    }

    return new_3D_points;
}



// Eigen::MatrixXd get_r3d(
//     const Eigen::MatrixXd& r,
//     const Eigen::MatrixXd& halfedges_uv,
//     const std::vector<int64_t>& halfedge_vertices_mapping,
//     Eigen::MatrixXi faces_uv
// ){
//     int num_r = r.rows();
//     Eigen::MatrixXd vertices_3D = loadMeshVertices("/Users/jan-piotraschke/git_repos/2DTissue/meshes/ellipsoid_x4.off");
//     Eigen::MatrixXd vertices_3D_ids(num_r, 3);
//     Eigen::MatrixXd distances_3D(num_r, 3);

//     for (int i = 0; i < num_r; ++i) {
//         std::vector<std::pair<double, int>> distances(faces_uv.rows());

//         for (int j = 0; j < faces_uv.rows(); ++j) {
//             Eigen::Vector3d a = vertices_3D.row(halfedge_vertices_mapping[faces_uv(j, 0)]);
//             Eigen::Vector3d b = vertices_3D.row(halfedge_vertices_mapping[faces_uv(j, 1)]);
//             Eigen::Vector3d c = vertices_3D.row(halfedge_vertices_mapping[faces_uv(j, 2)]);

//             distances[j] = {pointTriangleDistance(r.row(i), a, b, c), j};
//         }

//         std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

//         for (int j = 0; j < 3; ++j) {
//             vertices_3D_ids(i, j) = halfedge_vertices_mapping[faces_uv(min_distance.second, j)];
//             Eigen::Vector3d vertex = vertices_3D.row(vertices_3D_ids(i, j));
//             distances_3D(i, j) = (vertex - Eigen::Vector3d(r.row(i))).norm();
//         }
//     }

//     Eigen::MatrixXd new_3D_points(num_r, 3);

//     for (int i = 0; i < vertices_3D_ids.rows(); ++i) {
//         double w1 = distances_3D(i, 0);
//         double w2 = distances_3D(i, 1);
//         double w3 = distances_3D(i, 2);

//         int v1 = vertices_3D_ids(i, 0);
//         int v2 = vertices_3D_ids(i, 1);
//         int v3 = vertices_3D_ids(i, 2);

//         Eigen::Vector3d v1_3D = vertices_3D.row(v1);
//         Eigen::Vector3d v2_3D = vertices_3D.row(v2);
//         Eigen::Vector3d v3_3D = vertices_3D.row(v3);

//         Eigen::Vector3d newPoint = (w1 * v1_3D + w2 * v2_3D + w3 * v3_3D) / (w1 + w2 + w3);
//         new_3D_points.row(i) = newPoint;
//     }

//     return new_3D_points;
// }


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
