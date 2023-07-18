// Cell.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

#include <IO.h>

// Cell.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

#include <IO.h>

class Cell {
public:
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

    std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d();
    Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d();

    // getters
    Eigen::Matrix<double, Eigen::Dynamic, 2> get_r_UV() const { return r_UV; }
    Eigen::VectorXd get_n() const { return n; }
    Eigen::MatrixXd get_r_3D() const { return r_3D; }

private:
    const int num_part;
    const Eigen::MatrixXd halfedge_UV,
    const Eigen::MatrixXi face_UV,
    const Eigen::MatrixXd vertice_UV,
    const Eigen::MatrixXd vertice_3D,
    std::vector<int64_t> h_v_mapping;

    // additional private members
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV;
    Eigen::VectorXd n;
    Eigen::MatrixXd r_3D;

    std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(int interator);
    Eigen::Vector3d calculate_barycentric_2D_coord(int iterator);
};



