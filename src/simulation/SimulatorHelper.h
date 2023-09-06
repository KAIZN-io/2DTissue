// SimulatorHelper.h
#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include <Struct.h>
#include <SurfaceParametrization.h>

class SimulatorHelper {
public:
    SimulatorHelper(
        std::vector<VertexData>& particle_change,
        std::vector<bool>& simulated_particles,
        int particle_count,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
        Eigen::MatrixXd& r_3D,
        Eigen::MatrixXd& r_3D_old,
        Eigen::VectorXi& n,
        Eigen::VectorXi& n_pole,
        Eigen::VectorXi& n_pole_old,
        SurfaceParametrization& surface_parametrization,
        bool& original_mesh
    );

    void set_new_particle_data();
    void update_if_valid(std::vector<int> inside_UV_id);
    std::vector<int> get_inside_UV_id();
    std::vector<int> get_outside_UV_id();

private:
    SurfaceParametrization& surface_parametrization;

    std::vector<VertexData>& particle_change;
    std::vector<bool>& simulated_particles;
    int particle_count;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot;
    Eigen::MatrixXd& r_3D;
    Eigen::MatrixXd& r_3D_old;
    Eigen::VectorXi& n;
    Eigen::VectorXi& n_pole;
    Eigen::VectorXi& n_pole_old;

    bool& original_mesh;
    std::vector<int> outside_UV_id;

    // Check if the given point r is inside the UV parametrization bounds
    bool is_inside_uv(const Eigen::Vector2d r_UV_row) {
        return surface_parametrization.check_point_in_polygon(r_UV_row, original_mesh);
    };
};
