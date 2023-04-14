// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/jlcxx.hpp"
#include "Eigen/Dense"

#include "julia_handler.h"

// Jlcxx type aliases
using JuliaArray = jlcxx::ArrayRef<int64_t, 1>;
using JuliaArray2D = jlcxx::ArrayRef<double, 2>;


Eigen::MatrixXd reshape_vertices_array(
    const JuliaArray2D& vertices_stl,
    int num_rows,
    int num_cols
) {
    Eigen::MatrixXd vertices(num_rows, num_cols);

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            vertices(i, j) = vertices_stl[j * num_rows + i];
        }
    }

    return vertices;
}


Eigen::VectorXd jlcxxArrayRefToEigenVectorXd(
    const jlcxx::ArrayRef<int64_t, 1>& inputArray
){
    int arraySize = inputArray.size();
    Eigen::VectorXd outputVector(arraySize);

    for (int i = 0; i < arraySize; i++) {
        outputVector(i) = static_cast<double>(inputArray[i]);
    }

    return outputVector;
}