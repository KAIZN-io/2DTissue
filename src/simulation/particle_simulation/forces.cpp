// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <vector>
#include <Eigen/Dense>

#include <particle_simulation/forces.h>


Eigen::MatrixXd calculate_forces_between_particles(
    const std::vector<Eigen::MatrixXd>& dist_vect,
    const Eigen::MatrixXd& dist_length,
    double k,
    double σ,
    double r_adh,
    double k_adh
){
    // Get the number of particles
    int num_part = dist_vect[0].rows();

    // Initialize force matrix with zeros
    Eigen::MatrixXd F(num_part, 3);
    F.setZero();

    // Loop over all particle pairs
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < num_part; j++) {

            // Skip if particle itself
            if (i == j) continue;

            // Distance between particles A and B
            double dist = dist_length(i, j);

            // No force if particles too far from each other 
            if (dist >= 2 * σ) continue;

            // Add a small value if the distance is zero or you get nan values due to 'Fij * (dist_v / dist)' (division by zero)
            if (dist == 0) {
                dist += 0.0001;
            }

            // Eigen::Vector3d for the 3D distance vector
            Eigen::Vector3d dist_v(dist_vect[0](i, j), dist_vect[1](i, j), dist_vect[2](i, j));

            double Fij_rep = (-k * (2 * σ - dist)) / (2 * σ);
            double Fij_adh = (dist > r_adh) ? 0 : (k_adh * (2 * σ - dist)) / (2 * σ - r_adh);
            double Fij = Fij_rep + Fij_adh;

            F.row(i) += Fij * (dist_v / dist);
        }
    }

    // Actual force felt by each particle
    return F;
}
