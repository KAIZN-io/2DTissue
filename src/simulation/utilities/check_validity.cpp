// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

#include <utilities/check_validity.h>

bool are_all_valid(const std::vector<VertexData>& vertex_data) {
    for (const VertexData& data : vertex_data) {
        if (!data.valid) {
            return false;
        }
    }
    return true;
}

bool checkForInvalidValues(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> matrix
) {
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            if (std::isnan(matrix(i, j)) || std::isinf(matrix(i, j))) {
                return true;
            }
        }
    }
    return false;
}
