/**
 * @file        OrientationHelper.cpp
 * @brief       Calculate the orientation of the particles
 *
 * @author      Jan-Piotraschke
 * @date        2023-Sep-20
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        -
 */

#include "OrientationHelper.h"

OrientationHelper::OrientationHelper(
    const std::vector<Eigen::MatrixXd>& dist_vect,
    const Eigen::MatrixXd& dist_length,
    Eigen::VectorXi& n,
    double σ
)
 : dist_vect(dist_vect),
    dist_length(dist_length),
    n(n),
    σ(σ) {
};

// ========================================
// Public Functions
// ========================================

/**
* {\displaystyle \Theta _{i}(t+\Delta t)=\langle \Theta _{j}\rangle _{|r_{i}-r_{j}|<r}+\eta _{i}(t)}
*
* @brief At each time step, each particle aligns with its neighbours within a given distance with an uncertainity due to a noise.
*/
void OrientationHelper::calculate_average_n_within_distance(){
    // Get the number of particles
    int num_part = dist_vect[0].rows();

    // Initialize an Eigen::MatrixXd to store average n values
    Eigen::VectorXi avg_n(num_part);
    avg_n.setZero();

     // Set up the random number generation for noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_noise(0.00, 0.00);  // ! No noise for now

    // Loop over all particles
    for (int i = 0; i < num_part; i++) {
        // Initialize a vector to accumulate angle values and a counter for the number of valid pairs
        std::vector<double> angles;
        int count = 0;

        // Loop over all other particles
        for (int j = 0; j < num_part; j++) {

            // Distance between particles A and B
            double dist = dist_length(i, j);

            // Only consider particles within the specified distance
            if (dist < 2 * σ) {
                // Add the angle of particle j to the angles vector
                angles.push_back(n(j, 0));
                count++;
            }
        }

        // Calculate the average angle for particle i
        avg_n(i) = mean_unit_circle_vector_angle_degrees(angles);

        // Add noise to the average angle
        avg_n(i) += dist_noise(gen);
    }

    n = avg_n;
}



// ========================================
// Private Functions
// ========================================

/**
 * @brief Calculate the mean direction angle of a set of angles in degrees
*/

double OrientationHelper::mean_unit_circle_vector_angle_degrees(std::vector<double> angles) {
    if (angles.empty()) {
        throw std::invalid_argument("The input vector should not be empty.");
    }

    Eigen::Vector2d mean_vector(0.0, 0.0);

    for (const auto& angle_degrees : angles) {
        double angle_radians = angle_degrees * DEG_TO_RAD;

        // Convert the angle to a 2D unit vector
        Eigen::Vector2d vec(cos(angle_radians), sin(angle_radians));
        mean_vector += vec;
    }

    // Normalize the mean vector to keep it on the unit circle
    mean_vector.normalize();

    // Calculate the angle in radians using atan2 and convert it to degrees
    double angle_radians = std::atan2(mean_vector.y(), mean_vector.x());
    double angle_degrees = angle_radians * RAD_TO_DEG;

    // Make sure the angle is in the range [0, 360)
    if (angle_degrees < 0) {
        angle_degrees += FULL_CIRCLE;
    }

    return angle_degrees;
}
