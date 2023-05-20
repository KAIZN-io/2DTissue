// cell_cell_interactions.h
#pragma once

#include <Eigen/Dense>

Eigen::Vector3d repulsive_adhesion_motion(
    double k,
    double σ,
    double dist,
    double r_adh,
    double k_adh,
    const Eigen::Vector3d& dist_v
);