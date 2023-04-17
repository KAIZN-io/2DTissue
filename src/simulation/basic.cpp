// author: @Jan-Piotraschke
// date: 2023-03-24
// license: Apache License 2.0
// version: 0.1.0

// This file is examplary of how to pass a Julia Array to C++ and get a Vector back


#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"

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


Eigen::Vector3d compute_normal(
    const Eigen::MatrixXd& vertices,
    const JuliaArray& faces_stl
) {
    int a = faces_stl[0];
    int b = faces_stl[1];
    int c = faces_stl[2];

    Eigen::Vector3d A = vertices.row(a);
    Eigen::Vector3d B = vertices.row(b);
    Eigen::Vector3d C = vertices.row(c);

    Eigen::Vector3d BA = B - A;
    Eigen::Vector3d CA = C - A;

    return BA.cross(CA);
}


jlcxx::ArrayRef<double, 1> calculate_vertex_normals(
    JuliaArray faces_stl,
    JuliaArray2D vertices_stl
) {
    int num_entry = vertices_stl.size();
    int num_rows = num_entry / 3;

    Eigen::MatrixXd vertices = reshape_vertices_array(vertices_stl, num_rows, 3);
    Eigen::Vector3d normal = compute_normal(vertices, faces_stl);

    return jlcxx::ArrayRef<double, 1>(normal.data(), normal.size());
}


int main() {
    std::vector<double> data = {0, 3, 6, 9, 12, 15, 1, 4, 7, 10, 13, 16, 2, 5, 8, 11, 14, 17};

    int num_rows = 6;
    int num_cols = 3;

    Eigen::MatrixXd matrix(num_rows, num_cols);

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            matrix(i, j) = data[j * num_rows + i];
        }
    }

    std::cout << "Matrix:\n" << matrix << std::endl;
    return 0;
}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    // register a standard C++ function
    mod.method("calculate_vertex_normals", calculate_vertex_normals);
}
