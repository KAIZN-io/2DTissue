// Validation.h

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstdint>

class Validation {
public:
    bool checkForInvalidValues(const Eigen::Matrix<double, Eigen::Dynamic, 2> matrix);
    void error_invalid_values(Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new);
    void error_lost_particles(
        Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new,
        int num_part
    );

private:
    std::vector<int> find_inside_uv_vertices_id(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);
    bool is_inside_uv(const Eigen::Vector2d& r);
};