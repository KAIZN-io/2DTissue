// SimulatorHelper.h

#pragma once

#include <Eigen/Dense>

#include <vector>
#include <Struct.h>

class SimulatorHelper {
public:
    SimulatorHelper(
        std::vector<VertexData>& particle_change,
        int particle_count,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::MatrixXd& r_3D,
        Eigen::MatrixXd& r_3D_old,
        Eigen::VectorXi& n,
        Eigen::VectorXi& n_pole_old
    );

    void set_new_particle_data();

private:
    std::vector<VertexData>& particle_change;
    int particle_count;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::MatrixXd& r_3D;
    Eigen::MatrixXd& r_3D_old;
    Eigen::VectorXi& n;
    Eigen::VectorXi& n_pole_old;
};
