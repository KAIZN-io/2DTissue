/**
 * @file        Locomotion.cpp
 * @brief       Simulating the movement ("locomotion") of the particles on the surface
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jul-17
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        nicht nur muss beim Kachelmuster die Partikel die Distanz der Partikeln auf den Kachelmuster wissen, sondern die Kraft und Orientierung von den Partikeln ihres Seam Kanten Zwilling, der neben ihr der Monolite ist.
 */

#include <Locomotion.h>
#include "Locomotion/ForceHelper.h"
#include "Locomotion/OrientationHelper.h"

Locomotion::Locomotion(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
    Eigen::VectorXi& n,
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


// ========================================
// Public Functions
// ========================================

void Locomotion::simulate_flight() {
    resize_F_track();

    // Get distance vectors and calculate distances between particles
    auto dist_vect = get_dist_vect(r_UV);

    get_distances_between_particles(dist_length, distance_matrix, vertices_3D_active);

    // Calculate force between particles which pulls the particle in one direction within the 2D plane
    ForceHelper helper = ForceHelper(F_track, k, σ, r_adh, k_adh, dist_length, dist_vect);
    LocomotionHelperInterface& locomotion_helper = helper;
    locomotion_helper.calculate_forces_between_particles();

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

    OrientationHelper helper_orientation = OrientationHelper(dist_vect, dist_length, n, σ);
    OrientationHelperInterface& orientation_helper = helper_orientation;
    orientation_helper.calculate_average_n_within_distance();
}


/**
 * @brief Calculate the distance between each pair of particles
*/
void Locomotion::get_distances_between_particles(
    Eigen::MatrixXd& dist_length,
    Eigen::MatrixXd distance_matrix,
    std::vector<int> vertice_3D_id
){
    int num_part = vertice_3D_id.size();

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



// ========================================
// Private Functions
// ========================================

void Locomotion::resize_F_track() {
    F_track.resize(r_UV.rows(), 2);
}


/**
* @brief Transform a matrix into a symmetric matrix where the choosen value is the minimum of the two values
*
* @info: Unittest implemented
*/
void Locomotion::transform_into_symmetric_matrix(Eigen::MatrixXd &A) {
    int n = A.rows();

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            A(i, j) = A(j, i) = std::min(A(i, j), A(j, i));
        }
    }
}


/**
 * @info: Unittest implemented
*/
std::vector<Eigen::MatrixXd> Locomotion::get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r) {
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
