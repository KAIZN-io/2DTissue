/**
 * @file        CellHelper.cpp
 * @brief       Mapping 2D coordinates to 3D coordinates and vice versa; initializing the particle positions and orientations
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jul-18
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        mesh file path should be passed as an argument
 */

#include <CellHelper.h>

const std::filesystem::path MESH_CARTOGRAPHY = MeshCartographyLib_SOURCE_DIR;

CellHelper::CellHelper(
    int particle_count,
    Eigen::MatrixXi& face_UV,
    Eigen::MatrixXi& face_3D,
    Eigen::MatrixXd& vertice_UV,
    Eigen::MatrixXd& vertice_3D,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::MatrixXd& r_3D,
    Eigen::VectorXi& n
)
    : particle_count(particle_count),
    face_UV(face_UV),
    face_3D(face_3D),
    vertice_UV(vertice_UV),
    vertice_3D(vertice_3D),
    r_UV(r_UV),
    r_3D(r_3D),
    n(n)
{

}


// ========================================
// Public Functions
// ========================================

void CellHelper::init_particle_position() {
    std::vector<int> face_list(face_UV.rows());
    std::iota(face_list.begin(), face_list.end(), 0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    std::uniform_int_distribution<> dis_face(0, face_list.size() - 1);
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


// (2D Coordinates -> 3D Coordinates and Their Nearest 3D Vertice id (for the distance calculation on resimulations)) mapping
std::pair<Eigen::MatrixXd, std::vector<int>> CellHelper::get_r3d(){
    int num_r = r_UV.rows();
    Eigen::MatrixXd new_3D_points(num_r, 3);
    std::vector<int> nearest_vertices_ids(num_r);

    for (int i = 0; i < num_r; ++i) {
        auto [barycentric_coord, nearest_vertex_id] = calculate_barycentric_3D_coord(i);
        new_3D_points.row(i) = barycentric_coord;
        nearest_vertices_ids[i] = nearest_vertex_id;
    }

    return {new_3D_points, nearest_vertices_ids};
}



// ========================================
// Private Functions
// ========================================

Eigen::Vector2d CellHelper::get_face_gravity_center_coord(
    const Eigen::Vector3i r_face
) {
    Eigen::Vector3d center_face_test(0, 0, 0);

    for (int j = 0; j < 3; ++j) {
        center_face_test += vertice_UV.row(r_face[j]);
    }
    Eigen::Vector2d center_face = center_face_test.head(2);

    return center_face / 3.0;
}


std::pair<Eigen::Vector3d, int> CellHelper::calculate_barycentric_3D_coord(int iterator){
    std::vector<std::pair<double, int>> distances(face_UV.rows());

    for (int j = 0; j < face_UV.rows(); ++j) {
        Eigen::Vector3d uv_a = vertice_UV.row(face_UV(j, 0));
        Eigen::Vector3d uv_b = vertice_UV.row(face_UV(j, 1));
        Eigen::Vector3d uv_c = vertice_UV.row(face_UV(j, 2));

        distances[j] = {pointTriangleDistance(r_UV.row(iterator).head(2), uv_a, uv_b, uv_c), j};
    }

    std::pair<double, int> min_distance = *std::min_element(distances.begin(), distances.end());

    int triangle_vertice_a = face_UV(min_distance.second, 0);
    int triangle_vertice_b = face_UV(min_distance.second, 1);
    int triangle_vertice_c = face_UV(min_distance.second, 2);

    Eigen::Vector2d triangle_vertice_a_coord = vertice_UV.row(triangle_vertice_a);
    Eigen::Vector2d triangle_vertice_b_coord = vertice_UV.row(triangle_vertice_b);
    Eigen::Vector2d triangle_vertice_c_coord = vertice_UV.row(triangle_vertice_c);

    // Get the 3D coordinates of the 3 triangle_vertices
    Eigen::Vector3d a = vertice_3D.row(triangle_vertice_a);
    Eigen::Vector3d b = vertice_3D.row(triangle_vertice_b);
    Eigen::Vector3d c = vertice_3D.row(triangle_vertice_c);

    // Compute the weights (distances in UV space)
    double w_a = (r_UV.row(iterator).head(2).transpose() - triangle_vertice_a_coord).norm();
    double w_b = (r_UV.row(iterator).head(2).transpose() - triangle_vertice_b_coord).norm();
    double w_c = (r_UV.row(iterator).head(2).transpose() - triangle_vertice_c_coord).norm();

    // Normalize the weights
    normalize_weights(w_a, w_b, w_c);

    // Compute the new 3D point using the barycentric coordinates
    Eigen::Vector3d newPoint = w_a * a + w_b * b + w_c * c;

    // Compute distances from newPoint to a, b, c
    double dist_a = (newPoint - a).norm();
    double dist_b = (newPoint - b).norm();
    double dist_c = (newPoint - c).norm();

    // Find the vertex with the minimum distance to newPoint
    double min_dist = std::min({dist_a, dist_b, dist_c});

    int closest_vertice_id;
    if (min_dist == dist_a) closest_vertice_id = triangle_vertice_a;
    else if (min_dist == dist_b) closest_vertice_id = triangle_vertice_b;
    else closest_vertice_id = triangle_vertice_c;

    return std::make_pair(newPoint, closest_vertice_id);
}


double CellHelper::pointTriangleDistance(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c
) {
    // Edges of the triangle
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;
    // Line segments connecting the vertices with the point
    Eigen::Vector3d ap = p - a;
    Eigen::Vector3d bp = p - b;
    Eigen::Vector3d cp = p - c;

    double d_ab_ap = ab.dot(ap);
    double d_ac_ap = ac.dot(ap);
    double d_ab_bp = ab.dot(bp);
    double d_ac_bp = ac.dot(bp);
    double d_ab_cp = ab.dot(cp);
    double d_ac_cp = ac.dot(cp);

    // Vertex region outside A
    if (d_ab_ap <= 0.0 && d_ac_ap <= 0.0)
        return ap.norm();

    // Vertex region outside B
    if (d_ab_bp >= 0.0 && d_ac_bp <= d_ab_bp)
        return bp.norm();

    // Vertex region outside C
    if (d_ac_cp >= 0.0 && d_ab_cp <= d_ac_cp)
        return cp.norm();

    // Edge region of AB, and the projection of P onto AB.
    double vc = d_ab_ap * d_ac_bp - d_ab_bp * d_ac_ap;
    if (vc <= 0.0 && d_ab_ap >= 0.0 && d_ab_bp <= 0.0)
        return pointSegmentDistance(p, a, b);

    // Edge region of AC, and the projection of P onto AC.
    double vb = d_ab_cp * d_ac_ap - d_ab_ap * d_ac_cp;
    if (vb <= 0.0 && d_ac_ap >= 0.0 && d_ac_cp <= 0.0)
        return pointSegmentDistance(p, a, c);

    // Edge region of BC, and the projection of P onto BC.
    double va = d_ab_bp * d_ac_cp - d_ab_cp * d_ac_bp;
    if (va <= 0.0 && (d_ac_bp - d_ab_bp) >= 0.0 && (d_ab_cp - d_ac_cp) >= 0.0)
        return pointSegmentDistance(p, b, c);

    // Calculate the point's barycentric coordinates
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;

    // The projection of P onto the triangle's plane
    return (a + ab * v + ac * w - p).norm();
}


double CellHelper::pointSegmentDistance(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b
) {
    Eigen::Vector3d ab = b - a;
    double t = ab.dot(p - a) / ab.dot(ab);
    t = std::clamp(t, 0.0, 1.0);  // ensure t stays between 0 and 1
    return (a + ab * t - p).norm();
}


void CellHelper::normalize_weights(
    double& a,
    double& b,
    double& c
) {
    double sum_weights = a + b + c;
    a /= sum_weights;
    b /= sum_weights;
    c /= sum_weights;
}
