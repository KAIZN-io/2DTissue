// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <vector>
#include <Eigen/Dense>

#include <utilities/matrix_algebra.h>

/*
normalize a 3D matrix A
*/
Eigen::MatrixXd normalize_3D_matrix(const Eigen::MatrixXd &A) {
    Eigen::VectorXd row_norms_t = A.rowwise().norm(); // Compute row-wise Euclidean norms
    Eigen::MatrixXd row_norms = row_norms_t.replicate(1, 3); // Repeat each norm for each column

    return row_norms; // Normalize each row
}


/**
 * Calculate the cross product of two 3D matrices A and B.
 *
 * @param A The first input matrix.
 * @param B The second input matrix.
 *
 * @return The cross product of A and B, computed for each row of the matrices.
 */
Eigen::MatrixXd calculate_3D_cross_product(
    const Eigen::MatrixXd &A,
    const Eigen::MatrixXd &B
){ 
    // Ensure that A and B have the correct size
    assert(A.cols() == 3 && B.cols() == 3 && A.rows() == B.rows());

    // Preallocate output matrix with the same size and type as A
    const int num_rows = A.rows();
    Eigen::MatrixXd new_A = Eigen::MatrixXd::Zero(num_rows, 3);

    // Compute cross product for each row and directly assign the result to the output matrix
    for (int i = 0; i < num_rows; ++i) {
        // Get the i-th row of matrices A and B
        Eigen::Vector3d A_row = A.row(i);
        Eigen::Vector3d B_row = B.row(i);

        // Compute the cross product of the i-th rows of A and B
        new_A.row(i) = A_row.cross(B_row);
    }

    return new_A;
}
