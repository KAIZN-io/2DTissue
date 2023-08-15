// SimulatorHelper.h

#pragma once

#include <Eigen/Dense>

#include <vector>
#include <Struct.h>

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
        Eigen::VectorXi& n_pole_old
    );

    void set_new_particle_data();
    void update_if_valid(std::vector<int> inside_UV_id);
    std::vector<int> get_inside_UV_id();
    std::vector<int> get_outside_UV_id(std::vector<int> inside_UV_id);

private:
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

    // Check if the given point r is inside the UV parametrization bounds
    static bool is_inside_uv(const Eigen::Vector2d& r_UV) {
        return (0 <= r_UV[0] && r_UV[0] <= 1) && (0 <= r_UV[1] && r_UV[1] <= 1);
    };
};
