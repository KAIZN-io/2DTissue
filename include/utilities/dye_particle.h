// dye_particle.h
#pragma once

#include <Eigen/Core>
#include <vector>

Eigen::VectorXi dye_particles(const Eigen::MatrixXd dist_length, double σ);