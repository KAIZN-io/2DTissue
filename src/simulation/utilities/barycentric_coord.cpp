// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <Eigen/Dense>

#include <utilities/barycentric_coord.h>


double pointTriangleDistance(
    const Eigen::Vector3d p,
    const Eigen::Vector3d a,
    const Eigen::Vector3d b,
    const Eigen::Vector3d c
) {
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


int closestRow(const Eigen::MatrixXd& vertices_uv, const Eigen::Vector2d& halfedge_coord) {
    Eigen::VectorXd dists(vertices_uv.rows());
    for (int i = 0; i < vertices_uv.rows(); ++i) {
        dists[i] = (vertices_uv.row(i) - halfedge_coord.transpose()).squaredNorm();
    }

    Eigen::VectorXd::Index minRow;
    dists.minCoeff(&minRow);

    return minRow;
}


std::pair<Eigen::Vector3d, int>calculate_barycentric_3D_coord(
    const Eigen::MatrixXd r,
    const Eigen::MatrixXd halfedges_uv,
    const Eigen::MatrixXi faces_uv,
    Eigen::MatrixXd vertices_uv,
    Eigen::MatrixXd vertices_3D,
    std::vector<int64_t> h_v_mapping,
    int interator
){
    std::vector<std::pair<double, int>> distances(faces_uv.rows());

    for (int j = 0; j < faces_uv.rows(); ++j) {
        Eigen::Vector3d uv_a = halfedges_uv.row(faces_uv(j, 0));
        Eigen::Vector3d uv_b = halfedges_uv.row(faces_uv(j, 1));
        Eigen::Vector3d uv_c = halfedges_uv.row(faces_uv(j, 2));

        distances[j] = {pointTriangleDistance(r.row(interator).head(2), uv_a, uv_b, uv_c), j};
    }

    std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

    int halfedge_a = faces_uv(min_distance.second, 0);
    int halfedge_b = faces_uv(min_distance.second, 1);
    int halfedge_c = faces_uv(min_distance.second, 2);

    Eigen::Vector2d halfedge_a_coord = halfedges_uv.row(halfedge_a);
    Eigen::Vector2d halfedge_b_coord = halfedges_uv.row(halfedge_b);
    Eigen::Vector2d halfedge_c_coord = halfedges_uv.row(halfedge_c);

    // Inside your loop...
    int closest_a = closestRow(vertices_uv, halfedge_a_coord);
    int closest_b = closestRow(vertices_uv, halfedge_b_coord);
    int closest_c = closestRow(vertices_uv, halfedge_c_coord);

    // Get the 3D coordinates of the 3 halfedges
    Eigen::Vector3d a = vertices_3D.row(closest_a);
    Eigen::Vector3d b = vertices_3D.row(closest_b);
    Eigen::Vector3d c = vertices_3D.row(closest_c);

    // Compute the weights (distances in UV space)
    double w_a = (r.row(interator).head(2).transpose() - halfedge_a_coord).norm();
    double w_b = (r.row(interator).head(2).transpose() - halfedge_b_coord).norm();
    double w_c = (r.row(interator).head(2).transpose() - halfedge_c_coord).norm();

    // Compute the barycentric coordinates
    double sum_weights = w_a + w_b + w_c;
    w_a /= sum_weights;
    w_b /= sum_weights;
    w_c /= sum_weights;

    // // Closest UV vertice
    // double min_weight = std::min({w_a, w_b, w_c});
    // int closest_UV_vertice = (min_weight == w_a) ? closest_a : (min_weight == w_b) ? closest_b : closest_c;

    // Compute the new 3D point using the barycentric coordinates
    Eigen::Vector3d newPoint = w_a * a + w_b * b + w_c * c;

    // Compute distances from newPoint to a, b, c
    double dist_a = (newPoint - a).norm();
    double dist_b = (newPoint - b).norm();
    double dist_c = (newPoint - c).norm();

    // Find the vertex with the minimum distance to newPoint
    double min_dist = std::min({dist_a, dist_b, dist_c});

    int closest_row_id;
    if (min_dist == dist_a) closest_row_id = closest_a;
    else if (min_dist == dist_b) closest_row_id = closest_b;
    else closest_row_id = closest_c;

    // Get the vertice of h_v_mapping
    int closest_vertice_id = h_v_mapping[closest_row_id];

    return std::make_pair(newPoint, closest_vertice_id);
}
