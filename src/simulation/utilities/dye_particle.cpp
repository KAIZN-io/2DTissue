// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Core>
#include <vector>

#include <utilities/dye_particle.h>


void count_particle_neighbors(std::vector<int>& num_partic, const Eigen::MatrixXd dist_length, double σ) {
    const int num_rows = dist_length.rows();

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < dist_length.cols(); j++) {
            if (dist_length(i, j) != 0 && dist_length(i, j) <= 2.4 * σ) {
                num_partic[i] += 1;
            }
        }
    }
}
