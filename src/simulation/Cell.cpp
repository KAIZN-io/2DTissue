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

Cell::Cell(
    int particle_count,
    Eigen::MatrixXd& halfedge_UV,
    Eigen::MatrixXi& face_UV,
    Eigen::MatrixXi& face_3D,
    Eigen::MatrixXd& vertice_UV,
    Eigen::MatrixXd& vertice_3D,
    std::vector<int64_t>& h_v_mapping,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::MatrixXd& r_3D,
    Eigen::VectorXd& n
)
    : particle_count(particle_count),
    halfedge_UV(halfedge_UV),
    face_UV(face_UV),
    face_3D(face_3D),
    vertice_UV(vertice_UV),
    vertice_3D(vertice_3D),
    h_v_mapping(h_v_mapping),
    r_UV(r_UV),
    r_3D(r_3D),
    n(n)
{

}


Eigen::Vector2d Cell::get_face_gravity_center_coord(
    const Eigen::Vector3i r_face
) {
    Eigen::Vector3d center_face_test(0, 0, 0);

    for (int j = 0; j < 3; ++j) {
        center_face_test += halfedge_UV.row(r_face[j]);
    }
    Eigen::Vector2d center_face = center_face_test.head(2);

    return center_face / 3.0;
}


void Cell::init_particle_position() {
    int face_length = face_UV.rows();
    std::vector<int> face_list(face_length);
    std::iota(face_list.begin(), face_list.end(), 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    std::uniform_int_distribution<> dis_face(0, face_length - 1);
    std::uniform_int_distribution<> dis_angle(0, 359);

    for (int i = 0; i < particle_count; ++i) {
        int random_face = dis_face(gen);
        auto it = std::find(face_list.begin(), face_list.end(), random_face);
        if (it != face_list.end()) {
            face_list.erase(it);
        }

        Eigen::Vector3i r_face_uv = face_UV.row(random_face);
        r_UV.row(i) = get_face_gravity_center_coord(r_face_uv);
        n.row(i) << dis_angle(gen);
    }
}


double Cell::pointTriangleDistance(
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


int Cell::closestRow(const Eigen::Vector2d& halfedge_coord) {
    Eigen::VectorXd dists(vertice_UV.rows());
    for (int i = 0; i < vertice_UV.rows(); ++i) {
        dists[i] = (vertice_UV.row(i) - halfedge_coord.transpose()).squaredNorm();
    }

    Eigen::VectorXd::Index minRow;
    dists.minCoeff(&minRow);

    return minRow;
}


std::pair<Eigen::Vector3d, int> Cell::calculate_barycentric_3D_coord(int interator){
    std::vector<std::pair<double, int>> distances(face_UV.rows());

    for (int j = 0; j < face_UV.rows(); ++j) {
        Eigen::Vector3d uv_a = halfedge_UV.row(face_UV(j, 0));
        Eigen::Vector3d uv_b = halfedge_UV.row(face_UV(j, 1));
        Eigen::Vector3d uv_c = halfedge_UV.row(face_UV(j, 2));

        distances[j] = {pointTriangleDistance(r_UV.row(interator).head(2), uv_a, uv_b, uv_c), j};
    }

    std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

    int halfedge_a = face_UV(min_distance.second, 0);
    int halfedge_b = face_UV(min_distance.second, 1);
    int halfedge_c = face_UV(min_distance.second, 2);

    Eigen::Vector2d halfedge_a_coord = halfedge_UV.row(halfedge_a);
    Eigen::Vector2d halfedge_b_coord = halfedge_UV.row(halfedge_b);
    Eigen::Vector2d halfedge_c_coord = halfedge_UV.row(halfedge_c);

    // Inside your loop...
    int closest_a = closestRow(halfedge_a_coord);
    int closest_b = closestRow(halfedge_b_coord);
    int closest_c = closestRow(halfedge_c_coord);

    // Get the 3D coordinates of the 3 halfedges
    Eigen::Vector3d a = vertice_3D.row(closest_a);
    Eigen::Vector3d b = vertice_3D.row(closest_b);
    Eigen::Vector3d c = vertice_3D.row(closest_c);

    // Compute the weights (distances in UV space)
    double w_a = (r_UV.row(interator).head(2).transpose() - halfedge_a_coord).norm();
    double w_b = (r_UV.row(interator).head(2).transpose() - halfedge_b_coord).norm();
    double w_c = (r_UV.row(interator).head(2).transpose() - halfedge_c_coord).norm();

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


Eigen::Vector3d Cell::calculate_barycentric_2D_coord(int iterator){
    std::vector<std::pair<double, int>> distances(face_3D.rows());

    for (int j = 0; j < face_3D.rows(); ++j) {

        // find first occurence of vertice_a in h_v_mapping_vector and get the row index
        int halfedge_a = std::find(h_v_mapping.begin(), h_v_mapping.end(), face_3D(j, 0)) - h_v_mapping.begin();
        int halfedge_b = std::find(h_v_mapping.begin(), h_v_mapping.end(), face_3D(j, 1)) - h_v_mapping.begin();
        int halfedge_c = std::find(h_v_mapping.begin(), h_v_mapping.end(), face_3D(j, 2)) - h_v_mapping.begin();

        Eigen::Vector3d uv_a = vertice_3D.row(halfedge_a);
        Eigen::Vector3d uv_b = vertice_3D.row(halfedge_b);
        Eigen::Vector3d uv_c = vertice_3D.row(halfedge_c);

        distances[j] = {pointTriangleDistance(r_3D.row(iterator).head(2), uv_a, uv_b, uv_c), j};
    }

    std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

    int vertice_a = face_3D(min_distance.second, 0);
    int vertice_b = face_3D(min_distance.second, 1);
    int vertice_c = face_3D(min_distance.second, 2);

    int nearest_halfedge_a = std::find(h_v_mapping.begin(), h_v_mapping.end(), vertice_a) - h_v_mapping.begin();
    int nearest_halfedge_b = std::find(h_v_mapping.begin(), h_v_mapping.end(), vertice_b) - h_v_mapping.begin();
    int nearest_halfedge_c = std::find(h_v_mapping.begin(), h_v_mapping.end(), vertice_c) - h_v_mapping.begin();

    // Get the 3D coordinates of the 3 halfedges
    Eigen::Vector3d a = vertice_3D.row(nearest_halfedge_a);
    Eigen::Vector3d b = vertice_3D.row(nearest_halfedge_b);
    Eigen::Vector3d c = vertice_3D.row(nearest_halfedge_c);

    // Get the 2D coordinates of the 3 halfedges
    Eigen::Vector3d uv_a = vertice_UV.row(nearest_halfedge_a);
    Eigen::Vector3d uv_b = vertice_UV.row(nearest_halfedge_b);
    Eigen::Vector3d uv_c = vertice_UV.row(nearest_halfedge_c);

    // Compute the weights (distances in 3D space)
    double w_a = (r_3D.row(iterator).head(2).transpose() - a).norm();
    double w_b = (r_3D.row(iterator).head(2).transpose() - b).norm();
    double w_c = (r_3D.row(iterator).head(2).transpose() - c).norm();

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
std::pair<Eigen::MatrixXd, std::vector<int>> Cell::get_r3d(){
    int num_r = r_UV.rows();
    Eigen::MatrixXd new_3D_points(num_r, 3);
    std::vector<int> nearest_vertices_ids(num_r);

    for (int i = 0; i < num_r; ++i) {
        auto [barycentric_coord, nearest_vertex_id] = calculate_barycentric_3D_coord(i);
        new_3D_points.row(i) = barycentric_coord;
        nearest_vertices_ids[i] = nearest_vertex_id;
    }

    return std::make_pair(new_3D_points, nearest_vertices_ids);
}


// (3D Coordinates -> 2D Coordinates and Their Nearest 2D Vertice id) mapping
Eigen::Matrix<double, Eigen::Dynamic, 2> Cell::get_r2d(){
    // ! TODO: This is a temporary solution. The mesh file path should be passed as an argument.
    std::string mesh_3D_file_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    Eigen::MatrixXi face_3D;
    loadMeshFaces(mesh_3D_file_path, face_3D);

    int num_r = r_3D.rows();
    Eigen::Matrix<double, Eigen::Dynamic, 2> new_2D_points(num_r, 2);

    for (int i = 0; i < num_r; ++i) {
        Eigen::Vector3d uv_coord_test = calculate_barycentric_2D_coord(i);
        Eigen::Vector2d uv_coord(uv_coord_test[0], uv_coord_test[1]);

        new_2D_points.row(i) = uv_coord;
    }

    return new_2D_points;
}
