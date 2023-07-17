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



// #include <vector>
// #include <Eigen/Dense>
// #include <algorithm>
// #include <set>

// #include <utilities/check_boundary.h>

// class UVVertices {
// public:
//     UVVertices(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r)
//         : r_(r)
//     {
//         int nrows = r.rows();
//         for (int i = 0; i < nrows; ++i) {
//             Eigen::Vector2d first_two_columns = r.row(i).head<2>();
//             if (is_inside_uv(first_two_columns)) {
//                 inside_uv_ids_.insert(i);
//             }
//             else {
//                 outside_uv_ids_.insert(i);
//             }
//         }
//     }

//     const std::set<int>& get_inside_uv_ids() const {
//         return inside_uv_ids_;
//     }

//     const std::set<int>& get_outside_uv_ids() const {
//         return outside_uv_ids_;
//     }

// private:
//     // Check if the given point r is inside the UV parametrization bounds
//     static bool is_inside_uv(const Eigen::Vector2d& r) {
//         return (0 <= r[0] && r[0] <= 1) && (0 <= r[1] && r[1] <= 1);
//     }

//     Eigen::Matrix<double, Eigen::Dynamic, 2> r_;
//     std::set<int> inside_uv_ids_;
//     std::set<int> outside_uv_ids_;
// };
