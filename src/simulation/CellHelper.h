// CellHelper.h
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <random>
#include <unordered_set>
#include <vector>

#include "IO.h"
#include "IO/IO.h"

class CellHelper
{
  public:
    CellHelper(
        int particle_count,
        Eigen::MatrixXi& face_UV,
        Eigen::MatrixXi& face_3D,
        Eigen::MatrixXd& vertice_UV,
        Eigen::MatrixXd& vertice_3D,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::MatrixXd& r_3D,
        Eigen::VectorXi& n);

    void init_particle_position();
    std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d();
    Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d();

  private:
    int particle_count;
    Eigen::MatrixXi& face_UV;
    Eigen::MatrixXi& face_3D;
    Eigen::MatrixXd& vertice_UV;
    Eigen::MatrixXd& vertice_3D;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::MatrixXd& r_3D;
    Eigen::VectorXi& n;

    std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(int interator);
    Eigen::Vector2d get_face_gravity_center_coord(const Eigen::Vector3i r_face);
    double pointTriangleDistance(
        const Eigen::Vector3d& p, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);
    static double pointSegmentDistance(const Eigen::Vector3d& p, const Eigen::Vector3d& a, const Eigen::Vector3d& b);
    void normalize_weights(double& a, double& b, double& c);
};
