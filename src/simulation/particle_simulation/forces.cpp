// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <vector>
#include <Eigen/Dense>

#include <particle_simulation/cell_cell_interactions.h>
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

            // Calculate the force between particles A and B
            F.row(i) += repulsive_adhesion_motion(k, σ, dist, r_adh, k_adh, dist_v);
        }
    }

    // Actual force felt by each particle
    return F;
}
