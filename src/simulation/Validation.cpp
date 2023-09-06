/**
 * @file        Validation.cpp
 * @brief       Error checking
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jul-19
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        -
 */

#include <Validation.h>

Validation::Validation(
    SurfaceParametrization& surface_parametrization,
    bool& original_mesh
)
    : surface_parametrization(surface_parametrization),
        original_mesh(original_mesh)
{

}


// ========================================
// Public Functions
// ========================================

bool Validation::checkForInvalidValues(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> matrix
) {
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            if (std::isnan(matrix(i, j)) || std::isinf(matrix(i, j))) {
                return true;
            }
        }
    }
    return false;
}

bool Validation::are_all_valid(std::vector<VertexData> particle_change) {
    for (const VertexData& data : particle_change) {
        if (!data.valid) {
            return false;
        }
    }
    return true;
}

void Validation::error_invalid_3D_values(std::vector<VertexData> particle_change){
    if (!are_all_valid(particle_change)) {
        throw std::runtime_error("There are still particles outside the mesh");
    }
}

void Validation::error_invalid_values(
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new
){
    if (checkForInvalidValues(r_UV_new)) {
        std::exit(1);  // stop script execution
    }
}

std::vector<int> Validation::find_inside_uv_vertices_id(
    const Eigen::Matrix<double,
    Eigen::Dynamic, 2>& r
){
    int nrows = r.rows();
    std::vector<int> inside_id;

    for (int i = 0; i < nrows; ++i) {
        // Check if the point is inside the UV parametrization bounds
        Eigen::Vector2d first_two_columns = r.row(i).head<2>();
        if (is_inside_uv(first_two_columns)) {
            inside_id.push_back(i);
        }
    }

    return inside_id;
}

void Validation::error_lost_particles(
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new,
    int num_part
){
    if (find_inside_uv_vertices_id(r_UV_new).size() != num_part) {
        throw std::runtime_error("We lost particles after getting the original UV mesh coord");
    }
}



// ========================================
// Private Functions
// ========================================

bool Validation::is_inside_uv(const Eigen::Vector2d& r) {
    return surface_parametrization.check_point_in_polygon(r, original_mesh);
}
