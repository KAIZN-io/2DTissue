// author: @Jan-Piotraschke
// date: 2023-02-17
// license: Apache License 2.0
// version: 0.1.0

/*
Shortest paths on a terrain using one source point 
The heat method package returns an Approximation of the Geodesic Distance for all vertices of a triangle mesh to the closest vertex in a given set of source vertices.
As a rule of thumb, the method works well on triangle meshes, which are Delaunay.
*/

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstddef>
#include <vector>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"


using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;

using vertex_descriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using Vertex_distance_map = Triangle_mesh::Property_map<vertex_descriptor, double>;

//  The Intrinsic Delaunay Triangulation algorithm is switched off by the template parameter Heat_method_3::Direct.
using Heat_method_idt = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct>;
using Heat_method = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh>;


std::vector<double> geo_distance(
    int32_t start_node
){
    std::ifstream filename(CGAL::data_file_path("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/ellipsoid_x4.off"));
    Triangle_mesh tm;
    filename >> tm;

    //property map for the distance values to the source set
    Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

    //pass in the idt object and its vertex_distance_map
    Heat_method hm_idt(tm);

    //add the first vertex as the source set
    vertex_descriptor source = *(vertices(tm).first + start_node);
    hm_idt.add_source(source);
    hm_idt.estimate_geodesic_distances(vertex_distance);

    std::vector<double> distances_list;
    for(vertex_descriptor vd : vertices(tm)){
        distances_list.push_back(get(vertex_distance, vd));
    }

    return distances_list;
}


Eigen::MatrixXd get_distances_between_particles(Eigen::MatrixXd r, Eigen::MatrixXd distance_matrix, std::vector<int> vertice_3D_id) {
    int num_part = r.rows();

    // Get the distances from the distance matrix
    Eigen::MatrixXd dist_length = Eigen::MatrixXd::Zero(num_part, num_part);

    // ! BUG: it is interesting that using the heat method the distance from a -> b is not the same as b -> a
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < num_part; j++) {
            dist_length(i, j) = distance_matrix(vertice_3D_id[i], vertice_3D_id[j]);
        }
    }

    // Set the diagonal elements to 0.0
    dist_length.diagonal().array() = 0.0;

    return dist_length;
}


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


void fill_distance_matrix(Eigen::MatrixXd &distance_matrix, int closest_vertice){
    if (distance_matrix.row(closest_vertice).head(2).isZero()) {
        // get the distance of all vertices to all other vertices
        std::vector<double> vertices_3D_distance_map = geo_distance(closest_vertice);

        // fill the distance matrix
        for (int i = 0; i < vertices_3D_distance_map.size(); i++) {
            distance_matrix.coeffRef(closest_vertice, i) = vertices_3D_distance_map[i];
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


Eigen::MatrixXd calculate_forces_between_particles(
    std::vector<Eigen::MatrixXd>& dist_vect,
    Eigen::MatrixXd& dist_length,
    double k,
    double σ,
    double r_adh,
    double k_adh
){
    int num_part = dist_vect[0].rows();
    Eigen::MatrixXd F(num_part, 3);
    F.setZero();

    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < num_part; j++) {

            // No force if particle itself
            if (i != j) {
                // Distance between particles A and B
                double dist = dist_length(i, j);

                // Eigen::Vector3d for the 3D distance vector
                Eigen::Vector3d dist_v(dist_vect[0](i, j), dist_vect[1](i, j), dist_vect[2](i, j));

                // No force if particles too far from each other
                if (dist < 2 * σ) {
                    double Fij_rep = (-k * (2 * σ - dist)) / (2 * σ);
                    double Fij_adh = (dist > r_adh) ? 0 : (k_adh * (2 * σ - dist)) / (2 * σ - r_adh);
                    double Fij = Fij_rep + Fij_adh;

                    F.row(i) += Fij * (dist_v / dist);
                }
            }
        }
    }

    // Actual force felt by each particle
    return F;
}


Eigen::MatrixXd calculate_velocity(
    std::vector<Eigen::MatrixXd>& dist_vect,
    Eigen::MatrixXd& dist_length,
    Eigen::MatrixXd& n,
    double v0,
    double k,
    double σ,
    double μ,
    double r_adh,
    double k_adh
) {
    // Calculate force between particles
    Eigen::MatrixXd F_track = calculate_forces_between_particles(dist_vect, dist_length, k, σ, r_adh, k_adh);

    // Velocity of each particle
    Eigen::MatrixXd r_dot = v0 * n + μ * F_track;  // TODO: this isn't correct yet
    r_dot.col(2).setZero();

    return r_dot;
}


Eigen::MatrixXd calculate_next_position(
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& r_dot,
    double dt
){
    Eigen::MatrixXd r_new = r + r_dot * dt;
    r_new.col(2).setZero();
    return r_new;
}


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

std::vector<int> dye_particles(const Eigen::VectorXd& dist_length, int num_part, double σ) {
    // Count the number of neighbours for each particle
    Eigen::VectorXd number_neighbours = count_particle_neighbours(dist_length, σ);

    std::vector<int> N_color;
    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < number_neighbours(i); j++) {
            N_color.push_back(number_neighbours(i));
        }
    }

    return N_color;
}


// Eigen::MatrixXd calculate_3D_cross_product(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
//     Eigen::MatrixXd c(a.rows(), 3);
//     c.col(0) = a.col(1).array() * b.col(2).array() - a.col(2).array() * b.col(1).array();
//     c.col(1) = a.col(2).array() * b.col(0).array() - a.col(0).array() * b.col(2).array();
//     c.col(2) = a.col(0).array() * b.col(1).array() - a.col(1).array() * b.col(0).array();
//     return c;
// }


// Eigen::MatrixXd correct_n(
//     const Eigen::MatrixXd& r_dot,
//     const Eigen::MatrixXd& n,
//     double τ,
//     double dt
// ){
//     Eigen::MatrixXd ncross = calculate_3D_cross_product(n, r_dot).array().rowwise() / r_dot.rowwise().norm().array();
//     Eigen::MatrixXd n_cross_correction = (1.0 / τ) * ncross * dt;
//     Eigen::MatrixXd new_n = n - calculate_3D_cross_product(n, n_cross_correction);
//     return normalize_3D_matrix(new_n);
// }


Eigen::MatrixXd normalize_3D_matrix(const Eigen::MatrixXd &A) {
    Eigen::VectorXd row_norms_t = A.rowwise().norm(); // Compute row-wise Euclidean norms
    Eigen::MatrixXd row_norms = row_norms_t.replicate(1, 3); // Repeat each norm for each column

    return row_norms; // Normalize each row
}


int calculate_particle_vectors(Eigen::MatrixXd &r_dot, Eigen::MatrixXd &n) {

    // TODO:     n = correct_n(r_dot, n, τ, dt)

    n = normalize_3D_matrix(n);

    std::cout << "n: " << n << std::endl;
    Eigen::MatrixXd nr_dot = normalize_3D_matrix(r_dot);

    return 0;
}


int main()
{
    Eigen::MatrixXd distance_matrix(4670, 4670);
    std::vector<int> vertices_3D_active = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int num_part = 10;
    Eigen::MatrixXd n(5, 3);
    n <<  -0.999984,  -0.232996,   0.0,
            -0.736924,    0.0388327,  0.0,
            0.511211,  0.661931,   0.0,
            -0.0826997, -0.930856,   0.0,
            0.0655345, -0.893077,   0.0; 
    n.col(2).setZero();

    double v0 = 1.0;
    double k = 0.1;
    double σ = 0.5;
    double μ = 0.2;
    double r_adh = 0.4;
    double k_adh = 0.2;
    double dt = 0.1;

    // iterate over the active vertices and fill the distance matrix
    for (int i = 0; i < vertices_3D_active.size(); i++) {
        fill_distance_matrix(distance_matrix, vertices_3D_active[i]);
    }

    Eigen::MatrixXd r = Eigen::MatrixXd::Random(10, 3);  // create a 10x3 matrix with random values between -1 and 1
    r.block(0, 0, 10, 2) = (r.block(0, 0, 10, 2).array() + 1.0) / 2.0;  // rescale the values in the first 2 columns to be between 0 and 1
    r.block(0, 2, 10, 1) = Eigen::MatrixXd::Zero(10, 1);

    auto dist_vect = get_dist_vect(r);
    auto dist_length = get_distances_between_particles(r, distance_matrix, vertices_3D_active);
    transform_into_symmetric_matrix(dist_length);

    // calculate the next position and velocity of each particle based on the distances
    auto r_dot = calculate_velocity(dist_vect, dist_length, n, v0, k, σ, μ, r_adh, k_adh);
    auto r_new = calculate_next_position(r, r_dot, dt);
    // std::cout << r_new << std::endl;
    // std::cout << r_dot << std::endl;

    auto particles_color = dye_particles(dist_length, num_part, σ);

    calculate_particle_vectors(r_dot, n);
    return 0;
}






/*
normalize a 3D matrix A by dividing each row by its norm
*/
// Eigen::MatrixXd normalize_3D_matrix(
//     const Eigen::MatrixXd &A
// ){
//     return A.rowwise().normalized();
// }
