// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <utilities/dye_particle.h>


std::vector<int> count_particle_neighbors(const Eigen::MatrixXd dist_length, double σ) {
    int num_rows = dist_length.rows();
    std::vector<int> num_partic(num_rows); // create an empty vector

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < dist_length.cols(); j++) {
            if (dist_length(i, j) != 0 && dist_length(i, j) <= 2.4 * σ) {
                num_partic[i] += 1;
            }
        }
    }

    return num_partic;
}


std::vector<int> dye_particles(const Eigen::MatrixXd dist_length, double σ) {
    // Count the number of neighbours for each particle
    std::vector<int> number_neighbors = count_particle_neighbors(dist_length, σ);

    return number_neighbors;
}