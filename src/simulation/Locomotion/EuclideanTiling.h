#pragma once

#include <tuple>
#include <Eigen/Dense>
#include "SurfaceParametrization/SurfaceParametrization.h"
#include "SurfaceParametrization/TessellationHelper.h"

class EuclideanTiling {
public:
    EuclideanTiling(
        SurfaceParametrization& surface_parametrization,
        Tessellation& tessellation,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::VectorXi& n
    );

    void diagonal_seam_edges_square_border();
    std::pair<std::string, Point_2_eigen> check_border_crossings(
        const Point_2_eigen& start_eigen,
        const Point_2_eigen& end_eigen
    );

private:
    // Member variables
    SurfaceParametrization& surface_parametrization;
    Tessellation& tessellation;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::VectorXi& n;
    const bool original_mesh = true;
    std::vector<Point_2_eigen> left, right, up, down;

    // Helper functions
    std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> processPoints(
        const Eigen::Vector2d& pointA,
        const Eigen::Vector2d& point_outside,
        double n
    );

    bool is_point_on_segment(const Point_2_eigen& P, const Point_2_eigen& A, const Point_2_eigen& B);
    std::optional<Point_2_eigen> intersection_point(const Segment_2_eigen& line, const std::vector<Point_2_eigen>& border);

    // Constants
    static constexpr double KACHEL_ROTATION = 90.0;
    static constexpr double BORDER_THRESHOLD = 1e-3;
};
