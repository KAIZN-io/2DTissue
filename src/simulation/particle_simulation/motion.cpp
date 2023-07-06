// author: @Jan-Piotraschke
// date: 2023-04-12
// license: Apache License 2.0
// version: 0.1.0

#include <tuple>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>

#include <utilities/angles_to_unit_vectors.h>

#include <particle_simulation/forces.h>
#include <particle_simulation/motion.h>


/**
* @brief Transform a matrix into a symmetric matrix where the choosen value is the minimum of the two values
*
* @info: Unitests implemented
*/
void transform_into_symmetric_matrix(Eigen::MatrixXd &A) {
    int n = A.rows();

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            A(i, j) = A(j, i) = std::min(A(i, j), A(j, i));
        }
    }
}


/**
 * @info: Unitests implemented
*/
std::vector<Eigen::MatrixXd> get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r) {
    Eigen::VectorXd dist_x = r.col(0);
    Eigen::VectorXd dist_y = r.col(1);

    // Replicate each column into a square matrix
    Eigen::MatrixXd square_x = dist_x.replicate(1, r.rows());
    Eigen::MatrixXd square_y = dist_y.replicate(1, r.rows());

    // Compute the difference between each pair of rows
    Eigen::MatrixXd diff_x = square_x.array().rowwise() - dist_x.transpose().array();
    Eigen::MatrixXd diff_y = square_y.array().rowwise() - dist_y.transpose().array();

    // Store the difference matrices in a vector
    std::vector<Eigen::MatrixXd> dist_vect;
    dist_vect.push_back(diff_x);
    dist_vect.push_back(diff_y);

    return dist_vect;
}


/**
 * @brief Calculate the distance between each pair of particles
*/
Eigen::MatrixXd get_distances_between_particles(
    Eigen::Matrix<double, Eigen::Dynamic, 2> r,
    Eigen::MatrixXd distance_matrix,
    std::vector<int> vertice_3D_id
){
    int num_part = r.rows();

    // Get the distances from the distance matrix
    Eigen::MatrixXd dist_length = Eigen::MatrixXd::Zero(num_part, num_part);

    // Use the #pragma omp parallel for directive to parallelize the outer loop
    // The directive tells the compiler to create multiple threads to execute the loop in parallel, splitting the iterations among them
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < num_part; j++) {
            dist_length(i, j) = distance_matrix(vertice_3D_id[i], vertice_3D_id[j]);
        }
    }

    dist_length.diagonal().array() = 0.0;

    return dist_length;
}


/**
 * @brief Calculate the mean direction angle of a set of angles in degrees
 *
 * @info: Unitests implemented
*/
double mean_unit_circle_vector_angle_degrees(std::vector<double> angles) {
    if (angles.empty()) {
        throw std::invalid_argument("The input vector should not be empty.");
    }

    Eigen::Vector2d mean_vector(0.0, 0.0);

    for (const auto& angle_degrees : angles) {
        double angle_radians = angle_degrees * M_PI / 180.0;

        // Convert the angle to a 2D unit vector
        Eigen::Vector2d vec(cos(angle_radians), sin(angle_radians));
        mean_vector += vec;
    }

    // Normalize the mean vector to keep it on the unit circle
    mean_vector.normalize();

    // Calculate the angle in radians using atan2 and convert it to degrees
    double angle_radians = std::atan2(mean_vector.y(), mean_vector.x());
    double angle_degrees = angle_radians * 180.0 / M_PI;

    // Make sure the angle is in the range [0, 360)
    if (angle_degrees < 0) {
        angle_degrees += 360;
    }

    return angle_degrees;
}


/**
* {\displaystyle \Theta _{i}(t+\Delta t)=\langle \Theta _{j}\rangle _{|r_{i}-r_{j}|<r}+\eta _{i}(t)}
*
* @brief At each time step, each particle aligns with its neighbours within a given distance with an uncertainity due to a noise.
*/
void calculate_average_n_within_distance(
    const std::vector<Eigen::MatrixXd> dist_vect,
    const Eigen::MatrixXd dist_length,
    Eigen::VectorXd& n,
    double σ
){
    // Get the number of particles
    int num_part = dist_vect[0].rows();

    // Initialize an Eigen::MatrixXd to store average n values
    Eigen::VectorXd avg_n(num_part);
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




std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 2>, Eigen::MatrixXd, Eigen::MatrixXd> simulate_flight(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::VectorXd& n,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd distance_matrix_v,
    double v0,
    double k,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double step_size
){
    // Get distance vectors and calculate distances between particles
    auto dist_vect = get_dist_vect(r_UV);
    auto dist_length = get_distances_between_particles(r_UV, distance_matrix_v, vertices_3D_active);
    transform_into_symmetric_matrix(dist_length);

    // Calculate force between particles which pulls the particle in one direction within the 2D plane
    Eigen::MatrixXd F_track = calculate_forces_between_particles(dist_vect, dist_length, k, σ, r_adh, k_adh);
    Eigen::VectorXd abs_F = F_track.rowwise().norm();

    Eigen::Matrix<double, Eigen::Dynamic, 2> n_vec = angles_to_unit_vectors(n);

    // Velocity of each particle
    // 1. Every particle moves with a constant velocity v0 in the direction of the normal vector n
    // 2. Some particles are influenced by the force F_track
    abs_F = abs_F.array() + v0;

    // multiply elementwise the values of Eigen::VectorXd abs_F with the values of Eigen::MatrixXd n_vec to create a new matrix of size n_vec
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot = n_vec.array().colwise() * abs_F.array();

    // Calculate the new position of each particle
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_new = r_UV + r_dot * step_size;

    // Calculate the average for n for all particle pairs which are within dist < 2 * σ 
    calculate_average_n_within_distance(dist_vect, dist_length, n, σ);

    return std::make_tuple(r_new, r_dot, dist_length);
}
