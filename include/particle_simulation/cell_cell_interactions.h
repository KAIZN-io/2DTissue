// cell_cell_interactions.h
#pragma once

#include <Eigen/Dense>

Eigen::Vector2d repulsive_adhesion_motion(
    double k,
    double σ,
    double dist,
    double r_adh,
    double k_adh,
    const Eigen::Vector2d dist_v
);