// CellHelper.h
#pragma once

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <random>
#include <unordered_set>
#include <boost/filesystem.hpp>

#include "IO.h"

class CellHelper {
public:
    CellHelper(
        int particle_count,
        Eigen::MatrixXd& halfedge_UV,
        Eigen::MatrixXi& face_UV,
        Eigen::MatrixXi& face_3D,
        Eigen::MatrixXd& vertice_UV,
        Eigen::MatrixXd& vertice_3D,
        std::vector<int64_t>& h_v_mapping,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::MatrixXd& r_3D,
        Eigen::VectorXi& n
    );

    void init_particle_position();
    std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d();
    Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d();

private:
    int particle_count;
    Eigen::MatrixXd& halfedge_UV;
    Eigen::MatrixXi& face_UV;
    Eigen::MatrixXi& face_3D;
    Eigen::MatrixXd& vertice_UV;
    Eigen::MatrixXd& vertice_3D;
    std::vector<int64_t>& h_v_mapping;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::MatrixXd& r_3D;
    Eigen::VectorXi& n;

    std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(int interator);
    Eigen::Vector3d calculate_barycentric_2D_coord(int iterator);
    Eigen::Vector2d get_face_gravity_center_coord(
        const Eigen::Vector3i r_face
    );
    double pointTriangleDistance(
        const Eigen::Vector3d& p,
        const Eigen::Vector3d& a,
        const Eigen::Vector3d& b,
        const Eigen::Vector3d& c
    );
    bool isPointInsideTriangle(
        const Eigen::Vector3d& p,
        const Eigen::Vector3d& a,
        const Eigen::Vector3d& b,
        const Eigen::Vector3d& c
    );
    static double pointSegmentDistance(
        const Eigen::Vector3d& p,
        const Eigen::Vector3d& a,
        const Eigen::Vector3d& b
    );
    int closestRow(
        const Eigen::Vector2d& halfedge_coord
    );
    void normalize_weights(
        double& a,
        double& b,
        double& c
    );
};
