/**
 * @file        ForceHelper.cpp
 * @brief       Calculate the force acting between the particles
 *
 * @author      Jan-Piotraschke
 * @date        2023-Sep-20
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        -
 */

#include "ForceHelper.h"

ForceHelper::ForceHelper(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& F_track,
    double k,
    double σ,
    double r_adh,
    double k_adh,
    Eigen::MatrixXd& dist_length,
    const std::vector<Eigen::MatrixXd>& dist_vect
)
 : F_track(F_track),
    k(k),
    σ(σ),
    r_adh(r_adh),
    k_adh(k_adh),
    dist_length(dist_length),
    dist_vect(dist_vect) {
};


// ========================================
// Public Functions
// ========================================

/**
* @brief: Calculate the force that each particle feels due to all the other particles
*
* @info: Unittest implemented
*/
void ForceHelper::calculate_forces_between_particles(){
    // Get the number of particles
    int num_part = dist_vect[0].rows();
    F_track.setZero();

    // Loop over all particle pairs
    // #pragma omp parallel for
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
                dist += 0.001;
            }

            // Eigen::Vector3d for the 3D distance vector
            Eigen::Vector2d dist_v(dist_vect[0](i, j), dist_vect[1](i, j));

            // Calculate the force between particles A and B
            // #pragma omp critical
            F_track.row(i) += repulsive_adhesion_motion(k, σ, dist, r_adh, k_adh, dist_v);
        }
    }
}



// ========================================
// Private Functions
// ========================================

/**
* @brief: Locomotions can attract or repel each other as they move depending on their distance.
*
* @info: Unittest implemented
*/
Eigen::Vector2d ForceHelper::repulsive_adhesion_motion(
    double k,
    double σ,
    double dist,
    double r_adh,
    double k_adh,
    const Eigen::Vector2d dist_v
) {
    double Fij_rep = 0;
    double Fij_adh = 0;

    if (dist < 2*σ)
    {
        Fij_rep = (-k * (2 * σ - dist)) / (2 * σ);
    }

    if (dist >= 2*σ && dist <= r_adh)
    {
        Fij_adh = (k_adh * (2 * σ - dist)) / (2 * σ - r_adh);
    }

    // double Fij = 0.1 * Fij_rep + 0.1 * Fij_adh;
    double Fij = Fij_rep + Fij_adh;

    return Fij * (dist_v / dist);
}
