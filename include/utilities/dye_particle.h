// dye_particle.h
#pragma once

#include <Eigen/Core>
#include <vector>

std::vector<int> dye_particles(const Eigen::MatrixXd dist_length, double σ);