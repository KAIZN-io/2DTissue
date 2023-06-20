// author: @Jan-Piotraschke
// date: 2023-02-17
// license: Apache License 2.0
// version: 0.1.0

/*
Shortest paths on a terrain using one source point 
The heat method package returns an Approximation of the Geodesic Distance for all vertices of a triangle mesh to the closest vertex in a given set of source vertices.
As a rule of thumb, the method works well on triangle meshes, which are Delaunay.

Disclaimer: The heat method solver is the bottle neck of the algorithm.
*/

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <utilities/paths_lib.h>
#include <utilities/distance.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;
using vertex_descriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using Vertex_distance_map = Triangle_mesh::Property_map<vertex_descriptor, double>;

//  The Intrinsic Delaunay Triangulation algorithm is switched off by the template parameter Heat_method_3::Direct.
using Heat_method_idt = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct>;
using Heat_method = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh>;


std::vector<double> geo_distance(const std::string mesh_path, int32_t start_node){
    std::ifstream filename(CGAL::data_file_path(mesh_path));
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


/*
atomic variable to keep track of the current index of the vector of distances, and each thread processes a
different index until all the distances have been added to the distance matrix.
*/
void fill_distance_matrix(
    const std::string mesh_path,
    Eigen::MatrixXd &distance_matrix,
    int closest_vertice
){
    if (distance_matrix.row(closest_vertice).head(2).isZero()) {
        // get the distance of all vertices to all other vertices
        std::vector<double> vertices_3D_distance_map = geo_distance(mesh_path, closest_vertice);
        distance_matrix.row(closest_vertice) = Eigen::Map<Eigen::VectorXd>(vertices_3D_distance_map.data(), vertices_3D_distance_map.size());
    }
}


int get_all_distances(){
    auto mesh_path = get_full_path("2DTissue/meshes/ellipsoid_x4.off");
    std::cout << mesh_path << std::endl;

    std::ifstream filename(CGAL::data_file_path(mesh_path));
    Triangle_mesh tm;
    filename >> tm;

    Eigen::MatrixXd distance_matrix_v(num_vertices(tm), num_vertices(tm));
    // ! dieser Schritt ist der Bottleneck der Simulation!
    // ! wir müssen nämlich n mal die geo distance ausrechnen und die kostet jeweils min 25ms pro Start Vertex
    // loop over all vertices and fill the distance matrix
    for (auto vi = vertices(tm).first; vi != vertices(tm).second; ++vi) {
        fill_distance_matrix(mesh_path, distance_matrix_v, *vi);
    }

    // save the distance matrix to a csv file using comma as delimiter
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    const std::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

    std::cout << "Saving distance matrix to file..." << std::endl;
    std::string distance_matrix_path = PROJECT_PATH.string() + "/meshes/data/ellipsoid_x4_distance_matrix_static.csv";
    std::ofstream file(distance_matrix_path);
    file << distance_matrix_v.format(CSVFormat);
    file.close();
    std::cout << "saved" << std::endl;

    return 0;
}
