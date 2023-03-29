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

using JuliaArray = jlcxx::ArrayRef<int64_t, 1>;
using JuliaArray2D = jlcxx::ArrayRef<double, 2>;


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


/*
normalize a 3D matrix A
*/
Eigen::MatrixXd normalize_3D_matrix(const Eigen::MatrixXd &A) {
    Eigen::VectorXd row_norms_t = A.rowwise().norm(); // Compute row-wise Euclidean norms
    Eigen::MatrixXd row_norms = row_norms_t.replicate(1, 3); // Repeat each norm for each column

    return row_norms; // Normalize each row
}


/**
 * Calculate the cross product of two 3D matrices A and B.
 *
 * @param A The first input matrix.
 * @param B The second input matrix.
 *
 * @return The cross product of A and B, computed for each row of the matrices.
 */
Eigen::MatrixXd calculate_3D_cross_product(
    const Eigen::MatrixXd &A,
    const Eigen::MatrixXd &B
){ 
    // Ensure that A and B have the correct size
    assert(A.cols() == 3 && B.cols() == 3 && A.rows() == B.rows());

    // Preallocate output matrix with the same size and type as A
    const int num_rows = A.rows();
    Eigen::MatrixXd new_A = Eigen::MatrixXd::Zero(num_rows, 3);

    // Compute cross product for each row and directly assign the result to the output matrix
    for (int i = 0; i < num_rows; ++i) {
        // Get the i-th row of matrices A and B
        Eigen::Vector3d A_row = A.row(i);
        Eigen::Vector3d B_row = B.row(i);

        // Compute the cross product of the i-th rows of A and B
        new_A.row(i) = A_row.cross(B_row);
    }

    return new_A;
}


void calculate_order_parameter(
    Eigen::VectorXd& v_order, 
    Eigen::MatrixXd r, 
    Eigen::MatrixXd r_dot, 
    int num_part,
    double tt,
    double plotstep
) {
    // Define a vector normal to position vector and velocity vector
    Eigen::MatrixXd v_tp = calculate_3D_cross_product(r, r_dot);

    // Normalize v_tp
    Eigen::MatrixXd v_norm = v_tp.rowwise().normalized();

    // Sum v_tp vectors and divide by number of particle to obtain order parameter of collective motion for spheroids
    v_order((int)(tt / plotstep)) = (1.0 / num_part) * v_norm.colwise().sum().norm();
}


/*

Visceck-type n correction adapted from "Phys. Rev. E 74, 061908"
*/
Eigen::MatrixXd correct_n(
    const Eigen::MatrixXd& r_dot,
    const Eigen::MatrixXd& n,
    double τ,
    double dt
){
    // cross product of n and r_dot
    auto ncross = calculate_3D_cross_product(n, r_dot).cwiseQuotient(normalize_3D_matrix(r_dot));

    Eigen::MatrixXd n_cross_correction = (1.0 / τ) * ncross * dt;

    Eigen::MatrixXd new_n = n - calculate_3D_cross_product(n, n_cross_correction);

    return new_n.rowwise().normalized();
}


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> calculate_particle_vectors(
    Eigen::MatrixXd &r_dot,
    Eigen::MatrixXd &n
){
    // make a small correct for n according to Vicsek
    n = correct_n(r_dot, n, 1, 1);

    // Project the orientation of the corresponding faces using normal vectors
    n = n.rowwise().normalized();
    Eigen::MatrixXd nr_dot = r_dot.rowwise().normalized();

    return std::pair(n, nr_dot);
}


Eigen::MatrixXd reshape_vertices_array(
    const JuliaArray2D& vertices_stl,
    int num_rows,
    int num_cols
) {
    Eigen::MatrixXd vertices(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            vertices(i, j) = vertices_stl[j * num_rows + i];
        }
    }

    return vertices;
}


Eigen::VectorXd jlcxxArrayRefToEigenVectorXd(
    const jlcxx::ArrayRef<int64_t, 1>& inputArray
){
    int arraySize = inputArray.size();
    Eigen::VectorXd outputVector(arraySize);

    for (int i = 0; i < arraySize; i++) {
        outputVector(i) = static_cast<double>(inputArray[i]);
    }

    return outputVector;
}


void particle_simulation(
    jl_function_t* f,
    JuliaArray2D r_v,
    JuliaArray2D n_v,
    JuliaArray vertices_3D_active_id,
    JuliaArray2D distance_matrix_v,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
){

    int num_entry = r_v.size();
    int num_rows = num_entry / 3;
    int num_rows_dist = sqrt(distance_matrix_v.size());
    Eigen::MatrixXd r = reshape_vertices_array(r_v, num_rows, 3);
    Eigen::MatrixXd n = reshape_vertices_array(n_v, num_rows, 3);
    Eigen::MatrixXd distance_matrix = reshape_vertices_array(distance_matrix_v, num_rows_dist, num_rows_dist);
    Eigen::VectorXd vertices_3D_active_eigen = jlcxxArrayRefToEigenVectorXd(vertices_3D_active_id);

    // convert the active vertices to a vector -> only neccessary as long as the other functions of this script depend von std::vector
    std::vector<int> vertices_3D_active(vertices_3D_active_eigen.data(), vertices_3D_active_eigen.data() + vertices_3D_active_eigen.size());

    int num_part = 10;
    double plotstep = 0.1;

    // iterate over the active vertices and fill the distance matrix
    for (int i = 0; i < vertices_3D_active.size(); i++) {
        fill_distance_matrix(distance_matrix, vertices_3D_active[i]);
    }

    auto dist_vect = get_dist_vect(r);
    auto dist_length = get_distances_between_particles(r, distance_matrix, vertices_3D_active);
    transform_into_symmetric_matrix(dist_length);

    // calculate the next position and velocity of each particle based on the distances
    auto r_dot = calculate_velocity(dist_vect, dist_length, n, v0, k, σ, μ, r_adh, k_adh);
    auto r_new = calculate_next_position(r, r_dot, dt);

    auto particles_color = dye_particles(dist_length, num_part, σ);
    auto [ntest, nr_dot] = calculate_particle_vectors(r_dot, n);
    std::cout << ntest << std::endl;
    std::cout << nr_dot << std::endl;

    // Define the output vector v_order
    Eigen::VectorXd v_order((int)(tt / plotstep) + 1);

    // Calculate the order parameter
    calculate_order_parameter(v_order, r, r_dot, num_part, tt, plotstep);
    std::cout << v_order << std::endl;

    // transform the data into a Julia array
    auto r_julia = jlcxx::ArrayRef<double, 2>(r.data(), r.rows(), r.cols());
    auto r_dot_julia = jlcxx::ArrayRef<double, 2>(r_dot.data(), r_dot.rows(), r_dot.cols());
    auto dist_length_julia = jlcxx::ArrayRef<double, 2>(dist_length.data(), dist_length.rows(), dist_length.cols());
    auto distance_matrix_julia = jlcxx::ArrayRef<double, 2>(distance_matrix.data(), distance_matrix.rows(), distance_matrix.cols());


    auto v_order_julia = jlcxx::ArrayRef<double, 1>(v_order.data(), v_order.size());

    // Prepare to call the function defined in Julia
    jlcxx::JuliaFunction fnClb(f);

    // Fill the Julia Function with the inputs
    fnClb((jl_value_t*)r_julia.wrapped(), (jl_value_t*)r_dot_julia.wrapped(), (jl_value_t*)dist_length_julia.wrapped(), (jl_value_t*)distance_matrix_julia.wrapped());
}


int main()
{
    // Eigen::MatrixXd distance_matrix(4670, 4670);
    // std::vector<int> vertices_3D_active = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    // Eigen::MatrixXd n(5, 3);
    // n <<  -0.999984,  -0.232996,   0.0,
    //         -0.736924,    0.0388327,  0.0,
    //         0.511211,  0.661931,   0.0,
    //         -0.0826997, -0.930856,   0.0,
    //         0.0655345, -0.893077,   0.0; 
    // n.col(2).setZero();
    // Eigen::MatrixXd r = Eigen::MatrixXd::Random(5, 3);  // create a 10x3 matrix with random values between -1 and 1
    // r.block(0, 0, 10, 2) = (r.block(0, 0, 10, 2).array() + 1.0) / 2.0;  // rescale the values in the first 2 columns to be between 0 and 1
    // r.block(0, 2, 10, 1) = Eigen::MatrixXd::Zero(10, 1);
    return 0;
}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // Register a standard C++ function
    mod.method("particle_simulation", particle_simulation);
}
