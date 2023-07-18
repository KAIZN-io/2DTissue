// Cell.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

#include <IO.h>

class Cell {
pulbic:
    Cell(
        int num_part,
        const Eigen::MatrixXd halfedge_UV,
        const Eigen::MatrixXi face_UV,
        const Eigen::MatrixXd vertice_UV,
        const Eigen::MatrixXd vertice_3D,
        std::vector<int64_t> h_v_mapping
    );

    void init_particle_position(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::VectorXd& n
    );
    std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV
    );
    Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d(
        const Eigen::MatrixXd r_3D
    );

private:
    const int num_part;
    const Eigen::MatrixXd halfedge_UV,
    const Eigen::MatrixXi face_UV,
    const Eigen::MatrixXd vertice_UV,
    const Eigen::MatrixXd vertice_3D,
    std::vector<int64_t> h_v_mapping

    std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV,
        int interator
    );

    Eigen::Vector3d calculate_barycentric_2D_coord(
        const Eigen::MatrixXd r_3D,
        int iterator
    );
};







