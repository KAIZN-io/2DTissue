// author: @Jan-Piotraschke
// date: 2023-07-17
// license: Apache License 2.0
// version: 0.1.0

#include <tuple>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>

#include <LinearAlgebra.h>
#include <Simulator.h>

Simulator::Simulator(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
    Eigen::VectorXd& n,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd& distance_matrix,
    Eigen::MatrixXd& dist_length,
    double v0,
    double k,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double step_size,
    std::unique_ptr<LinearAlgebra> linear_algebra_ptr
)
    : r_UV(r_UV),
      r_UV_old(r_UV_old),
      r_dot(r_dot),
      n(n),
      vertices_3D_active(vertices_3D_active),
      distance_matrix(distance_matrix),
      dist_length(dist_length),
      v0(v0),
      k(k),
      σ(σ),
      μ(μ),
      r_adh(r_adh),
      k_adh(k_adh),
      step_size(step_size),
      linear_algebra_ptr(std::move(linear_algebra_ptr)) {
}

void Simulator::resize_F_track() {
    F_track.resize(r_UV.rows(), 2);
}
    // linear_algebra_ptr(std::make_unique<LinearAlgebra>())

void Simulator::simulate_flight() {
    resize_F_track();

    // Get distance vectors and calculate distances between particles
    auto dist_vect = get_dist_vect(r_UV);

    get_distances_between_particles(dist_length, r_UV, distance_matrix, vertices_3D_active);

    // Calculate force between particles which pulls the particle in one direction within the 2D plane
    calculate_forces_between_particles(dist_vect);
    Eigen::VectorXd abs_F = F_track.rowwise().norm();
    Eigen::Matrix<double, Eigen::Dynamic, 2> n_vec = linear_algebra_ptr->angles_to_unit_vectors(n);

    // Velocity of each particle
    // 1. Every particle moves with a constant velocity v0 in the direction of the normal vector n
    // 2. Some particles are influenced by the force F_track
    abs_F = abs_F.array() + v0;

    // multiply elementwise the values of Eigen::VectorXd abs_F with the values of Eigen::MatrixXd n_vec to create a new matrix of size n_vec
    r_dot = n_vec.array().colwise() * abs_F.array();

    // Calculate the new position of each particle
    r_UV += r_dot * step_size;

    // Calculate the average for n for all particle pairs which are within dist < 2 * σ 
    calculate_average_n_within_distance(dist_vect, dist_length, n, σ);
}


/**
* @brief Transform a matrix into a symmetric matrix where the choosen value is the minimum of the two values
*
* @info: Unittest implemented
*/
void Simulator::transform_into_symmetric_matrix(Eigen::MatrixXd &A) {
    int n = A.rows();

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            A(i, j) = A(j, i) = std::min(A(i, j), A(j, i));
        }
    }
}


/**
 * @brief Calculate the mean direction angle of a set of angles in degrees
 *
 * @info: Unittest implemented
*/
double Simulator::mean_unit_circle_vector_angle_degrees(std::vector<double> angles) {
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
* @brief: Simulators can attract or repel each other as they move depending on their distance.
*
* @info: Unittest implemented
*/
Eigen::Vector2d Simulator::repulsive_adhesion_motion(
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

    double Fij = Fij_rep + Fij_adh;

    return Fij * (dist_v / dist);
}


/**
* @brief: Calculate the force that each particle feels due to all the other particles
*
* @info: Unittest implemented
*/
void Simulator::calculate_forces_between_particles(const std::vector<Eigen::MatrixXd> dist_vect){
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


/**
 * @info: Unittest implemented
*/
std::vector<Eigen::MatrixXd> Simulator::get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r) {
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
void Simulator::get_distances_between_particles(
    Eigen::MatrixXd& dist_length,
    Eigen::Matrix<double, Eigen::Dynamic, 2> r,
    Eigen::MatrixXd distance_matrix,
    std::vector<int> vertice_3D_id
){
    int num_part = r.rows();

    // Use the #pragma omp parallel for directive to parallelize the outer loop
    // The directive tells the compiler to create multiple threads to execute the loop in parallel, splitting the iterations among them
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < num_part; j++) {
            dist_length(i, j) = distance_matrix(vertice_3D_id[i], vertice_3D_id[j]);
        }
    }
    dist_length.diagonal().array() = 0.0;
    transform_into_symmetric_matrix(dist_length);
}


/**
* {\displaystyle \Theta _{i}(t+\Delta t)=\langle \Theta _{j}\rangle _{|r_{i}-r_{j}|<r}+\eta _{i}(t)}
*
* @brief At each time step, each particle aligns with its neighbours within a given distance with an uncertainity due to a noise.
*
* @info: Unittest implemented
*/
void Simulator::calculate_average_n_within_distance(
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


// function to calculate y given x
double Simulator::interpolateY(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double x) {
    double y = pointA[1] + ((x - pointA[0]) * (pointB[1] - pointA[1])) / (pointB[0] - pointA[0]);
    // For the case that the point is on the edge
    if (y == 1) {
        y -= 0.0001;
    }
    else if (y == 0) {
        y += 0.0001;
    }
    return y;
}

// function to calculate x given y
double Simulator::interpolateX(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double y) {
    double x = pointA[0] + ((y - pointA[1]) * (pointB[0] - pointA[0])) / (pointB[1] - pointA[1]);
    // For the case that the point is on the edge
    if (x == 1) {
        x -= 0.0001;
    }
    else if (x == 0) {
        x += 0.0001;
    }
    return x;
}

// Function to calculate steepness switch
int Simulator::calculateSteepnessSwitch(double steepness) {
    if (steepness > 0) return -1;
    else if (steepness < 0) return 1;
    return 0;
}

std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> Simulator::processPoints(const Eigen::Vector2d& pointA, const Eigen::Vector2d& point_outside, double n) {
    Eigen::Vector2d entry_angle(1, 1);
    Eigen::Vector2d entry_point(1, 1);
    Eigen::Vector2d new_point(2, 1);
    auto delta_x = point_outside[0] - pointA[0];
    auto delta_y = point_outside[1] - pointA[1];
    auto steepness = delta_y / delta_x;
    int steepness_switch = calculateSteepnessSwitch(steepness);

    if (point_outside[0] < 0 || point_outside[0] > 1 || point_outside[1] < 0 || point_outside[1] > 1) {
        // oben oder rechts
        if (delta_x >= 0 && delta_y >= 0){
            double x = 1;
            double y = interpolateY(pointA, point_outside, x);

            // rechte Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(1, y);
                Eigen::Vector2d entry_point(y, 1);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y, because we also switched the entry coordinates before
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(0) *= steepness_switch;  // has to be variable
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n -= 90;
            }
            // obere Grenze passiert
            else {
                double y_back = 1;
                double x_back = interpolateX(pointA, point_outside, y_back);

                Eigen::Vector2d exit_point(x_back, 1);
                Eigen::Vector2d entry_point(1, x_back);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(1) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n += 90;
            }
        }
        // unten oder rechts
        else if (delta_x > 0 && delta_y < 0){
            double x = 1;
            double y = interpolateY(pointA, point_outside, x);

            // rechte Grenze passiert
            if (y < 1 && y > 0){

                Eigen::Vector2d exit_point(1, y);
                Eigen::Vector2d entry_point(y, 1);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(0) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array()  * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n -= 90;
            }
            // unten Grenze passiert
            else {
                double y_back_neg = 0;
                double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

                Eigen::Vector2d exit_point(x_back_neg, 0);
                Eigen::Vector2d entry_point(0, x_back_neg);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n += 90;
            }
        }
        // oben oder links
        else if (delta_x < 0 && delta_y > 0){
            double x = 0;
            double y = interpolateY(pointA, point_outside, x);

            // linke Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(0, y);
                Eigen::Vector2d entry_point(y, 0);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(0) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n -= 90;
            }
            // obere Grenze passiert
            else {
                double y_back = 1;
                double x_back = interpolateX(pointA, point_outside, y_back);

                Eigen::Vector2d exit_point(x_back, 1);
                Eigen::Vector2d entry_point(1, x_back);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n += 90;
            }
        }
        // unten oder links
        else {
            double x = 0;
            double y = interpolateY(pointA, point_outside, x);

            // linke Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(0, y);
                Eigen::Vector2d entry_point(y, 0);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();

                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(0) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n -= 90;
            }
            // unten Grenze passiert
            else {
                double y_back_neg = 0;
                double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

                Eigen::Vector2d exit_point(x_back_neg, 0);
                Eigen::Vector2d entry_point(0, x_back_neg);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n += 90;
            }
        }
    }
    else {
        new_point = point_outside;
    }
    return std::tuple(new_point, n, entry_point);
}


/**
 * @brief Because we have a mod(2) seam edge cute line, pairing edges are on the exact same opposite position in the UV mesh with the same lenght
*/
void Simulator::opposite_seam_edges_square_border(){
    r_UV.col(0) = r_UV.col(0).array() - r_UV.col(0).array().floor();  // Wrap x values
    r_UV.col(1) = r_UV.col(1).array() - r_UV.col(1).array().floor();  // Wrap y values
}


/**
 * @brief By using the '&' we pass the reference of the variable to the function, so we can change the value of the variable inside the function
*/
void Simulator::diagonal_seam_edges_square_border(){
    bool valid;
    do {
        valid = true;
        for (int i = 0; i < r_UV_old.rows(); ++i) {
            Eigen::Vector2d pointA = r_UV_old.row(i).head<2>(); // only takes the first two columns for the ith row
            Eigen::Vector2d point_outside = r_UV.row(i).head<2>(); // only takes the first two columns for the ith row
            double n_double = n(i);

            auto results = processPoints(pointA, point_outside, n_double);
            auto new_point = std::get<0>(results);
            n(i) = std::get<1>(results);  // ? Nachschauen, ob diese Zeile nicht hinter dem 'else' sein sollte
            auto entry_point = std::get<2>(results);

            // ! TODO: this logic can be improved
            // As soon as one of the conditions is not met (i.e., a value is outside the [0,1] interval), it breaks the for loop and starts another iteration of the while loop.
            if (new_point[0] < 0 || new_point[0] > 1 || new_point[1] < 0 || new_point[1] > 1) {
                r_UV_old.row(i).head<2>() = entry_point;
                r_UV.row(i).head<2>().noalias() = new_point;
                valid = false;
                break;
            }
            else {
                r_UV.row(i).head<2>().noalias() = new_point;
            }
        }
    } while (!valid);
}
