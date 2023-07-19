// Cell.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

#include "IO.h"

// Cell.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

class Cell {
public:
    Cell(
        int num_part,
        Eigen::MatrixXd halfedge_UV,
        Eigen::MatrixXi face_UV,
        Eigen::MatrixXi face_3D,
        Eigen::MatrixXd vertice_UV,
        Eigen::MatrixXd vertice_3D,
        std::vector<int64_t> h_v_mapping
    );

    void init_particle_position(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::VectorXd& n
    );

    std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d(Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV);
    Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d(Eigen::MatrixXd r_3D);

private:
    int num_part;
    Eigen::MatrixXd halfedge_UV;
    Eigen::MatrixXi face_UV;
    Eigen::MatrixXi face_3D;
    Eigen::MatrixXd vertice_UV;
    Eigen::MatrixXd vertice_3D;
    std::vector<int64_t> h_v_mapping;

    std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV, int interator);
    Eigen::Vector3d calculate_barycentric_2D_coord(Eigen::MatrixXd r_3D, int iterator);
    Eigen::Vector2d get_face_gravity_center_coord(
        const Eigen::Vector3i r_face
    );
    double pointTriangleDistance(
        const Eigen::Vector3d p,
        const Eigen::Vector3d a,
        const Eigen::Vector3d b,
        const Eigen::Vector3d c
    );
    int closestRow(const Eigen::MatrixXd& vertice_UV, const Eigen::Vector2d& halfedge_coord);
};



