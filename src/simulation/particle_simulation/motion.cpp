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

#include <particle_simulation/forces.h>
#include <particle_simulation/motion.h>


/*
Reminder: if you access the input variable with the '&' sign, you can change the variable in the function, without to return the new value.
The variable is Changed In the Memory and with that also in the main function
*/
void transform_into_symmetric_matrix(Eigen::MatrixXd &A) {
    int n = A.rows();

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (A(i, j) != 0 && A(j, i) != 0) {
                A(i, j) = A(j, i) = std::min(A(i, j), A(j, i));
            } else if (A(i, j) == 0) {
                A(i, j) = A(j, i);
            } else {
                A(j, i) = A(i, j);
            }
        }
    }
}


std::vector<Eigen::MatrixXd> get_dist_vect(const Eigen::MatrixXd& r) {
    Eigen::VectorXd dist_x = r.col(0);
    Eigen::VectorXd dist_y = r.col(1);
    Eigen::VectorXd dist_z = r.col(2);

    // Replicate each column into a square matrix
    Eigen::MatrixXd square_x = dist_x.replicate(1, r.rows());
    Eigen::MatrixXd square_y = dist_y.replicate(1, r.rows());
    Eigen::MatrixXd square_z = dist_z.replicate(1, r.rows());

    // Compute the difference between each pair of rows
    Eigen::MatrixXd diff_x = square_x.array().rowwise() - dist_x.transpose().array();
    Eigen::MatrixXd diff_y = square_y.array().rowwise() - dist_y.transpose().array();
    Eigen::MatrixXd diff_z = square_z.array().rowwise() - dist_z.transpose().array();

    // Store the difference matrices in a vector
    std::vector<Eigen::MatrixXd> dist_vect;
    dist_vect.push_back(diff_x);
    dist_vect.push_back(diff_y);
    dist_vect.push_back(diff_z);

    return dist_vect;
}


Eigen::MatrixXd get_distances_between_particles(
    Eigen::MatrixXd r,
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


double mean_unit_circle_vector_angle_degrees(std::vector<double> angles) {
    if (angles.empty()) {
        throw std::invalid_argument("The input vector should not be empty.");
    }

    Eigen::Vector2d sum_vector(0.0, 0.0);

    for (const auto& angle_degrees : angles) {
        double angle_radians = angle_degrees * M_PI / 180.0;

        // Convert the angle to a 2D unit vector
        Eigen::Vector2d vec(cos(angle_radians), sin(angle_radians));
        sum_vector += vec;
    }

    Eigen::Vector2d mean_vector = sum_vector / static_cast<double>(angles.size());

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


/*
{\displaystyle \Theta _{i}(t+\Delta t)=\langle \Theta _{j}\rangle _{|r_{i}-r_{j}|<r}+\eta _{i}(t)}

At each time step, each particle aligns with its neighbours within a given distance with an uncertainity due
to a noise.
*/
Eigen::MatrixXd calculate_average_n_within_distance(
    const std::vector<Eigen::MatrixXd> dist_vect,
    const Eigen::MatrixXd dist_length,
    Eigen::MatrixXd& n,
    double σ
){
    // Get the number of particles
    int num_part = dist_vect[0].rows();

    // Initialize an Eigen::MatrixXd to store average n values
    Eigen::MatrixXd avg_n(num_part, 1);
    avg_n.setZero();

    // Loop over all particles
    for (int i = 0; i < num_part; i++) {
        // Initialize a vector to accumulate angle values and a counter for the number of valid pairs
        std::vector<double> angles;
        // Eigen::Vector3d sum_n(0, 0, 0);
        int count = 0;

        // Loop over all other particles
        for (int j = 0; j < num_part; j++) {

            // Distance between particles A and B
            double dist = dist_length(i, j);

            // Only consider particles within the specified distance
            if (dist < 2 * σ) {
                // Add the angle of particle j to the angles vector
                angles.push_back(n(j, 0));

                // // Add the n vector of particle j to the sum
                // sum_n += n.row(j);
                count++;
            }
        }

        // Calculate the average angle for particle i
        avg_n(i, 0) = mean_unit_circle_vector_angle_degrees(angles);
        // avg_n.row(i) = sum_n / count;

        // // Add random noise to the average angle
        // double noise = std::abs(Eigen::MatrixXd::Random(1, 1)(0, 0)) * 0.099 + 0.001; // Range: [0.001, 0.1)
        // avg_n(i, 0) += noise;

        // // Add random noise to the average n
        // Eigen::MatrixXd noise = Eigen::MatrixXd::Random(1, 3).cwiseAbs();
        // Eigen::MatrixXd range(1, 3);
        // range << 0.099, 0.099, 0.099; // Set the range (0.1 - 0.001)
        // noise = (noise.cwiseProduct(range)).array() + 0.001;

        // avg_n.row(i) += noise.row(0);
    }

    n = avg_n;
    return avg_n;
}


Eigen::MatrixXd angles_to_unit_vectors(const Eigen::MatrixXd& avg_n) {
    if (avg_n.cols() != 1) {
        throw std::invalid_argument("The input matrix must have exactly 1 column.");
    }

    // Initialize an Eigen::MatrixXd to store the 2D unit vectors
    Eigen::MatrixXd n_vec(avg_n.rows(), 3);

    for (int i = 0; i < avg_n.rows(); ++i) {
        double angle_degrees = avg_n(i, 0);
        double angle_radians = angle_degrees * M_PI / 180.0;

        // Convert the angle to a 2D unit vector
        Eigen::Vector2d vec(cos(angle_radians), sin(angle_radians));
        n_vec.row(i) = vec;
    }
    n_vec.col(2).setZero();

    return n_vec;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> simulate_flight(
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& n,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd distance_matrix_v,
    double v0,
    double k,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt
){
    // Get distance vectors and calculate distances between particles
    auto dist_vect = get_dist_vect(r);
    auto dist_length = get_distances_between_particles(r, distance_matrix_v, vertices_3D_active);
    transform_into_symmetric_matrix(dist_length);

    // Calculate force between particles which pulls the particle in one direction within the 2D plane
    Eigen::MatrixXd F_track = calculate_forces_between_particles(dist_vect, dist_length, k, σ, r_adh, k_adh);
    Eigen::VectorXd abs_F = F_track.rowwise().norm();

    Eigen::MatrixXd n_vec = angles_to_unit_vectors(n);

    // Velocity of each particle
    // 1. Every particle moves with a constant velocity v0 in the direction of the normal vector n
    // 2. Some particles are influenced by the force F_track
    abs_F = abs_F.array() + v0;
    
    // ? Should I use the code in particle_vector.cpp for n_vec -> "a small correction of n"

    // multiply elementwise the values of Eigen::VectorXd abs_F with the values of Eigen::MatrixXd n_vec to create a new matrix of size n_vec
    Eigen::MatrixXd r_dot = n_vec.array().colwise() * abs_F.array();
    // Eigen::MatrixXd r_dot = v0 * n_vec + μ * F_track;
    // Eigen::MatrixXd r_dot = v0 * n_vec;
    r_dot.col(2).setZero();

    // Calculate the new position of each particle
    Eigen::MatrixXd r_new = r + r_dot * dt;
    r_new.col(2).setZero();

    // Calculate the average for n for all particle pairs which are within dist < 2 * σ 
    calculate_average_n_within_distance(dist_vect, dist_length, n, σ);

    return std::make_tuple(r_new, r_dot, dist_length);
}
