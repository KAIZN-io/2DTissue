// author: @Jan-Piotraschke
// date: 2023-07-18
// license: Apache License 2.0
// version: 0.1.0

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <random>
#include <Eigen/Dense>
#include <unordered_set>
#include <boost/filesystem.hpp>

#include <Cell.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;


Eigen::Vector2d get_face_gravity_center_coord(
    const Eigen::MatrixXd& vertices,
    const Eigen::Vector3i& r_face
) {
    Eigen::Vector3d center_face_test(0, 0, 0);

    for (int j = 0; j < 3; ++j) {
        center_face_test += vertices.row(r_face[j]);
    }
    Eigen::Vector2d center_face = center_face_test.head(2);

    return center_face / 3.0;
}


void init_particle_position(
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd halfedges_uv,
    int num_part,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r,
    Eigen::VectorXd& n
) {
    int faces_length = faces_uv.rows();
    std::vector<int> faces_list(faces_length);
    std::iota(faces_list.begin(), faces_list.end(), 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    std::uniform_int_distribution<> dis_face(0, faces_length - 1);
    std::uniform_int_distribution<> dis_angle(0, 359);

    for (int i = 0; i < num_part; ++i) {
        int random_face = dis_face(gen);
        auto it = std::find(faces_list.begin(), faces_list.end(), random_face);
        if (it != faces_list.end()) {
            faces_list.erase(it);
        }

        Eigen::Vector3i r_face_uv = faces_uv.row(random_face);
        r.row(i) = get_face_gravity_center_coord(halfedges_uv, r_face_uv);

        n.row(i) << dis_angle(gen);
    }
}


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
    const Eigen::Matrix<double, Eigen::Dynamic, 2> r,
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


Eigen::Vector3d calculate_barycentric_2D_coord(
    const Eigen::MatrixXd& start_3D_points,
    const Eigen::MatrixXi& faces_3D_static,
    const Eigen::MatrixXd& vertices_UV,
    const Eigen::MatrixXd& vertices_3D,
    std::vector<int64_t>& h_v_mapping,
    int iterator
){
    std::vector<std::pair<double, int>> distances(faces_3D_static.rows());

    for (int j = 0; j < faces_3D_static.rows(); ++j) {

        // find first occurence of vertice_a in h_v_mapping_vector and get the row index
        int halfedge_a = std::find(h_v_mapping.begin(), h_v_mapping.end(), faces_3D_static(j, 0)) - h_v_mapping.begin();
        int halfedge_b = std::find(h_v_mapping.begin(), h_v_mapping.end(), faces_3D_static(j, 1)) - h_v_mapping.begin();
        int halfedge_c = std::find(h_v_mapping.begin(), h_v_mapping.end(), faces_3D_static(j, 2)) - h_v_mapping.begin();

        Eigen::Vector3d uv_a = vertices_3D.row(halfedge_a);
        Eigen::Vector3d uv_b = vertices_3D.row(halfedge_b);
        Eigen::Vector3d uv_c = vertices_3D.row(halfedge_c);

        distances[j] = {pointTriangleDistance(start_3D_points.row(iterator).head(2), uv_a, uv_b, uv_c), j};
    }

    std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

    int vertice_a = faces_3D_static(min_distance.second, 0);
    int vertice_b = faces_3D_static(min_distance.second, 1);
    int vertice_c = faces_3D_static(min_distance.second, 2);

    int nearest_halfedge_a = std::find(h_v_mapping.begin(), h_v_mapping.end(), vertice_a) - h_v_mapping.begin();
    int nearest_halfedge_b = std::find(h_v_mapping.begin(), h_v_mapping.end(), vertice_b) - h_v_mapping.begin();
    int nearest_halfedge_c = std::find(h_v_mapping.begin(), h_v_mapping.end(), vertice_c) - h_v_mapping.begin();

    // Get the 3D coordinates of the 3 halfedges
    Eigen::Vector3d a = vertices_3D.row(nearest_halfedge_a);
    Eigen::Vector3d b = vertices_3D.row(nearest_halfedge_b);
    Eigen::Vector3d c = vertices_3D.row(nearest_halfedge_c);

    // Get the 2D coordinates of the 3 halfedges
    Eigen::Vector3d uv_a = vertices_UV.row(nearest_halfedge_a);
    Eigen::Vector3d uv_b = vertices_UV.row(nearest_halfedge_b);
    Eigen::Vector3d uv_c = vertices_UV.row(nearest_halfedge_c);

    // Compute the weights (distances in 3D space)
    double w_a = (start_3D_points.row(iterator).head(2).transpose() - a).norm();
    double w_b = (start_3D_points.row(iterator).head(2).transpose() - b).norm();
    double w_c = (start_3D_points.row(iterator).head(2).transpose() - c).norm();

    // Compute the barycentric coordinates
    double sum_weights = w_a + w_b + w_c;
    w_a /= sum_weights;
    w_b /= sum_weights;
    w_c /= sum_weights;

    // Compute the new 3D point using the barycentric coordinates
    Eigen::Vector3d newPoint = w_a * uv_a + w_b * uv_b + w_c * uv_c;

    return newPoint;
}


// (2D Coordinates -> 3D Coordinates and Their Nearest 3D Vertice id (for the distance calculation on resimulations)) mapping
std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> r,
    const Eigen::MatrixXd halfedges_uv,
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
Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d(
    const Eigen::MatrixXd r,
    const Eigen::MatrixXd vertices_uv,
    const Eigen::MatrixXd vertices_3D,
    std::vector<int64_t> h_v_mapping
){
    // ! TODO: This is a temporary solution. The mesh file path should be passed as an argument.
    std::string mesh_3D_file_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    Eigen::MatrixXi faces_3D_static;
    loadMeshFaces(mesh_3D_file_path, faces_3D_static);

    int num_r = r.rows();
    Eigen::Matrix<double, Eigen::Dynamic, 2> new_2D_points(num_r, 2);

    for (int i = 0; i < num_r; ++i) {
        Eigen::Vector3d uv_coord_test = calculate_barycentric_2D_coord(r, faces_3D_static, vertices_uv, vertices_3D, h_v_mapping, i);
        Eigen::Vector2d uv_coord(uv_coord_test[0], uv_coord_test[1]);

        new_2D_points.row(i) = uv_coord;
    }

    return new_2D_points;
}


// // (3D Vertice row position -> nD Vertice coordinates) mapping
// Eigen::MatrixXd get_coordinates(
//     std::vector<int> indices,
//     Eigen::MatrixXd coord
// ){
//     Eigen::MatrixXd found_coord(indices.size(), coord.cols());

//     for (int i = 0; i < indices.size(); ++i) {
//         found_coord.row(i) = coord.row(indices[i]);
//     }

//     return found_coord;
// }

// (3D Vertice id -> 3D Vertice row position of the h-v map) mapping
// std::vector<int> find_vertice_rows_index(
//     std::vector<int64_t> h_v_mapping_vector,
//     std::vector<int> r3d_vertices
// ){
//     std::unordered_set<int> found_ids;
//     std::vector<int> indices;

//     for (int i = 0; i < h_v_mapping_vector.size(); ++i) {
//         if (std::find(r3d_vertices.begin(), r3d_vertices.end(), h_v_mapping_vector[i]) != r3d_vertices.end()) {
//             if (found_ids.find(h_v_mapping_vector[i]) == found_ids.end()) {
//                 indices.push_back(i);
//                 found_ids.insert(h_v_mapping_vector[i]);
//             }
//         }
//     }
//     return indices;
// }