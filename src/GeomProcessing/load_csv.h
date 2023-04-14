// load_csv.h
#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

template<typename M>
M load_csv(const std::string &path);