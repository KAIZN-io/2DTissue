// csv.h
#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <boost/filesystem.hpp>

const boost::filesystem::path PROJECT_PATH_IO = PROJECT_SOURCE_DIR;


// We need do define it in the header file or otherwise the template specialization will not be available at link time
template<typename M>
M load_csv(const std::string &path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}


template <typename MatrixType>
void save_matrix_to_csv(const MatrixType& matrix, const std::string& file_name, int num_particles) {
    std::string path = PROJECT_PATH_IO.string() + "/data/" + file_name;
    std::ofstream file(path, std::ios::app); // Open the file in append mode

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    // Write the number of particles first
    // file << num_particles << ",";

    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            file << std::setprecision(15) << matrix(i, j);
            if (j < matrix.cols() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}
