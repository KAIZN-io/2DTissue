// check_validity.h
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstdint>

#include "utilities/sim_structs.h"

bool are_all_valid(const std::vector<VertexData>& vertex_data);

bool checkForInvalidValues(const Eigen::Matrix<double, Eigen::Dynamic, 2> matrix);