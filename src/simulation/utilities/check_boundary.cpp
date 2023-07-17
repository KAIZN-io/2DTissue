// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0
// description: Contains UV mesh related functions

#include <vector>
#include <Eigen/Dense>

#include <utilities/check_boundary.h>


// Check if the given point r is inside the UV parametrization bounds
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

std::vector<int> set_difference(int num_part, const std::vector<int>& inside_uv_ids) {
    std::vector<int> outside_uv_ids;

    // Create a copy of inside_uv_ids to sort without modifying the input
    std::vector<int> sorted_inside_uv_ids = inside_uv_ids;
    std::sort(sorted_inside_uv_ids.begin(), sorted_inside_uv_ids.end()); // Sort sorted_inside_uv_ids for efficient lookup

    for (int i = 0; i < num_part; ++i) {
        // If i is not found in sorted_inside_uv_ids, add it to outside_uv_ids
        if (std::binary_search(sorted_inside_uv_ids.begin(), sorted_inside_uv_ids.end(), i) == false) {
            outside_uv_ids.push_back(i);
        }
    }

    return outside_uv_ids;
}
