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
 * @todo        mesh file path should be passed as an argument; fix the isPointInsideTriangle function
 */

#include <CellHelper.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

CellHelper::CellHelper(
    int particle_count,
    Eigen::MatrixXd& halfedge_UV,
    Eigen::MatrixXi& face_UV,
    Eigen::MatrixXi& face_3D,
    Eigen::MatrixXd& vertice_UV,
    Eigen::MatrixXd& vertice_3D,
    std::vector<int64_t>& h_v_mapping,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::MatrixXd& r_3D,
    Eigen::VectorXi& n
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


// (3D Coordinates -> 2D Coordinates and Their Nearest 2D Vertice id) mapping
Eigen::Matrix<double, Eigen::Dynamic, 2> CellHelper::get_r2d(){
    const std::string mesh_3D_file_path = (PROJECT_PATH / "meshes/ellipsoid_x4.off").string();
    Eigen::MatrixXi face_3D;
    loadMeshFaces(mesh_3D_file_path, face_3D);

    int num_r = r_3D.rows();
    Eigen::Matrix<double, Eigen::Dynamic, 2> new_2D_points(num_r, 2);

    for (int i = 0; i < num_r; ++i) {
        Eigen::Vector3d uv_coord_test = calculate_barycentric_2D_coord(i);
        new_2D_points.row(i) = uv_coord_test.head<2>();
    }

    return new_2D_points;
}



// ========================================
// Private Functions
// ========================================

std::pair<Eigen::Vector3d, int> CellHelper::calculate_barycentric_3D_coord(int iterator){
    std::vector<std::pair<double, int>> distances(face_UV.rows());

    for (int j = 0; j < face_UV.rows(); ++j) {
        Eigen::Vector3d uv_a = halfedge_UV.row(face_UV(j, 0));
        Eigen::Vector3d uv_b = halfedge_UV.row(face_UV(j, 1));
        Eigen::Vector3d uv_c = halfedge_UV.row(face_UV(j, 2));

        distances[j] = {pointTriangleDistance(r_UV.row(iterator).head(2), uv_a, uv_b, uv_c), j};
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
    double w_a = (r_UV.row(iterator).head(2).transpose() - halfedge_a_coord).norm();
    double w_b = (r_UV.row(iterator).head(2).transpose() - halfedge_b_coord).norm();
    double w_c = (r_UV.row(iterator).head(2).transpose() - halfedge_c_coord).norm();

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

    int closest_row_id;
    if (min_dist == dist_a) closest_row_id = closest_a;
    else if (min_dist == dist_b) closest_row_id = closest_b;
    else closest_row_id = closest_c;

    // Get the vertice of h_v_mapping
    int closest_vertice_id = h_v_mapping[closest_row_id];

    return std::make_pair(newPoint, closest_vertice_id);
}


Eigen::Vector3d CellHelper::calculate_barycentric_2D_coord(int iterator){
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

    // Normalize the weights
    normalize_weights(w_a, w_b, w_c);

    // Compute the new 3D point using the barycentric coordinates
    Eigen::Vector3d newPoint = w_a * uv_a + w_b * uv_b + w_c * uv_c;

    return newPoint;
}


bool CellHelper::isPointInsideTriangle(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c
) {
    Eigen::Vector3d v0 = b - a, v1 = c - a, v2 = p - a;
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;

    double beta = (d11 * d20 - d01 * d21) / denom;
    double gamma = (d00 * d21 - d01 * d20) / denom;
    double alpha = 1.0 - beta - gamma;

    return alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1;
}


Eigen::Vector2d CellHelper::get_face_gravity_center_coord(
    const Eigen::Vector3i r_face
) {
    Eigen::Vector3d center_face_test(0, 0, 0);

    for (int j = 0; j < 3; ++j) {
        center_face_test += halfedge_UV.row(r_face[j]);
    }
    Eigen::Vector2d center_face = center_face_test.head(2);

    return center_face / 3.0;
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


int CellHelper::closestRow(
    const Eigen::Vector2d& halfedge_coord
) {
    Eigen::VectorXd dists(vertice_UV.rows());
    for (int i = 0; i < vertice_UV.rows(); ++i) {
        dists[i] = (vertice_UV.row(i) - halfedge_coord.transpose()).squaredNorm();
    }

    Eigen::VectorXd::Index minRow;
    dists.minCoeff(&minRow);

    return minRow;
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
