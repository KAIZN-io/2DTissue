// Validation.h

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstdint>
#include <Struct.h>
#include <SurfaceParametrization.h>

class Validation {
public:
    Validation(
        SurfaceParametrization& surface_parametrization,
        bool& original_mesh
    );
    bool checkForInvalidValues(const Eigen::Matrix<double, Eigen::Dynamic, 2> matrix);
    void error_invalid_3D_values(std::vector<VertexData> particle_change);
    void error_invalid_values(Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new);
    void error_lost_particles(
        Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_new,
        int num_part
    );
    bool are_all_valid(std::vector<VertexData> particle_change);
    std::vector<int> find_inside_uv_vertices_id(
        const Eigen::Matrix<double,
        Eigen::Dynamic, 2>& r
    );

private:
    SurfaceParametrization& surface_parametrization;
    bool& original_mesh;

    bool is_inside_uv(const Eigen::Vector2d& r);
};