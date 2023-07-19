// author: @Jan-Piotraschke
// date: 2023-06-28
// license: Apache License 2.0
// version: 0.1.0

#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

#include <utilities/check_validity.h>
#include <utilities/sim_structs.h>

#include <utilities/error_checking.h>


void error_unvalid_vertices(
    std::vector<VertexData> vertex_data
){
    if (!are_all_valid(vertex_data)) {
        throw std::runtime_error("There are still particles outside the mesh");
    }
}

void error_invalid_values(
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new
){
    if (checkForInvalidValues(r_UV_new)) {
        std::exit(1);  // stop script execution
    }
}

bool is_inside_uv(const Eigen::Vector2d& r) {
    return (0 <= r[0] && r[0] <= 1) && (0 <= r[1] && r[1] <= 1);
}

std::vector<int> find_inside_uv_vertices_id(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r) {
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

void error_lost_particles(
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new,
    int num_part
){
    if (find_inside_uv_vertices_id(r_UV_new).size() != num_part) {
        throw std::runtime_error("We lost particles after getting the original UV mesh coord");
    }
}