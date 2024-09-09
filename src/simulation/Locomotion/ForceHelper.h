#pragma once

#include <Eigen/Dense>
#include <vector>

#include "LocomotionHelperInterface.h"

class ForceHelper : public LocomotionHelperInterface
{
  public:
    ForceHelper(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& F_track,
        double k,
        double σ,
        double r_adh,
        double k_adh,
        Eigen::MatrixXd& dist_length,
        const std::vector<Eigen::MatrixXd>& dist_vect);

    void calculate_forces_between_particles() override;

  private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& F_track;
    double k, σ, r_adh, k_adh;
    Eigen::MatrixXd& dist_length;
    const std::vector<Eigen::MatrixXd>& dist_vect;

    Eigen::Vector2d repulsive_adhesion_motion(
        double k, double σ, double dist, double r_adh, double k_adh, const Eigen::Vector2d dist_v);
};
