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


int main()
{
    Eigen::MatrixXd distance_matrix(4670, 4670);
    std::vector<int> vertices_3D_active = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    // iterate over the active vertices and fill the distance matrix
    for (int i = 0; i < vertices_3D_active.size(); i++) {
        fill_distance_matrix(distance_matrix, vertices_3D_active[i]);
    }

    Eigen::MatrixXd r = Eigen::MatrixXd::Random(10, 3);  // create a 10x3 matrix with random values between -1 and 1
    r.block(0, 0, 10, 2) = (r.block(0, 0, 10, 2).array() + 1.0) / 2.0;  // rescale the values in the first 2 columns to be between 0 and 1
    r.block(0, 2, 10, 1) = Eigen::MatrixXd::Zero(10, 1);

    auto dist_vect = get_dist_vect(r);
    Eigen::MatrixXd dist_length = get_distances_between_particles(r, distance_matrix, vertices_3D_active);
    transform_into_symmetric_matrix(dist_length);
    std::cout << dist_vect << std::endl;

    return 0;
}
