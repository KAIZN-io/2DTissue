#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <random>
#include <vector>

#include "LocomotionHelperInterface.h"

class OrientationHelper : public OrientationHelperInterface
{
  public:
    OrientationHelper(
        const std::vector<Eigen::MatrixXd>& dist_vect, const Eigen::MatrixXd& dist_length, Eigen::VectorXi& n, double σ);

    void calculate_average_n_within_distance() override;

  private:
    const std::vector<Eigen::MatrixXd>& dist_vect;
    const Eigen::MatrixXd& dist_length;
    Eigen::VectorXi& n;
    double σ;

    static double mean_unit_circle_vector_angle_degrees(std::vector<double> angles);

    static constexpr double DEG_TO_RAD = M_PI / 180.0;
    static constexpr double RAD_TO_DEG = 180.0 / M_PI;
    static constexpr double FULL_CIRCLE = 360.0;
    static constexpr double QUARTER_CIRCLE = 90.0;
};
