// dye_particle.h
#pragma once

#include <Eigen/Core>
#include <vector>

Eigen::VectorXd dye_particles(const Eigen::VectorXd& dist_length, double σ);