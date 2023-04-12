// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Core>
#include <vector>

#include "dye_particle.h"

Eigen::VectorXd count_particle_neighbours(const Eigen::VectorXd& dist_length, double σ) {
    Eigen::VectorXd num_partic(dist_length.size()); // create an empty vector
    num_partic.setZero(); // initialize to zero

    for (int i = 0; i < dist_length.size(); i++) {
        if (dist_length(i) == 0 || dist_length(i) > 2.4 * σ) {
            num_partic(i) = 0;
        }
    }
    Eigen::VectorXd num_neighbors = num_partic.colwise().sum();
    return num_neighbors;
}


Eigen::VectorXd dye_particles(const Eigen::VectorXd& dist_length, double σ) {
    // Count the number of neighbours for each particle
    Eigen::VectorXd number_neighbours = count_particle_neighbours(dist_length, σ);
    int num_part = number_neighbours.size();
    std::vector<int> N_color_temp;
    // #pragma omp parallel for
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < number_neighbours(i); j++) {
            N_color_temp.push_back(number_neighbours(i));
        }
    }

    // Convert std::vector<int> to Eigen::VectorXd
    Eigen::VectorXd N_color(N_color_temp.size());
    for (size_t i = 0; i < N_color_temp.size(); ++i) {
        N_color(i) = N_color_temp[i];
    }

    return N_color;
}