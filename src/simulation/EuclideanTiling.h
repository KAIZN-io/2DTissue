// EuclideanTiling.h

#pragma once

#include <tuple>
#include <Eigen/Dense>

class EuclideanTiling {
public:
    EuclideanTiling(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::VectorXi& n
    );

    void opposite_seam_edges_square_border();
    void diagonal_seam_edges_square_border();

private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::VectorXi& n;

    std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> processPoints(
        const Eigen::Vector2d& pointA,
        const Eigen::Vector2d& point_outside,
        double n
    );
    int calculateSteepnessSwitch(double steepness);
    double interpolateX(
        const Eigen::Vector2d& pointA,
        const Eigen::Vector2d& pointB,
        double y
    );
    double interpolateY(
        const Eigen::Vector2d& pointA,
        const Eigen::Vector2d& pointB,
        double x
    );

    static constexpr double DEG_TO_RAD = M_PI / 180.0;
    static constexpr double RAD_TO_DEG = 180.0 / M_PI;
    static constexpr double FULL_CIRCLE = 360.0;
    static constexpr double QUARTER_CIRCLE = 90.0;
};
