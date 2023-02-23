// author: @Jan-Piotraschke
// date: 2023-02-17
// license: Apache License 2.0
// version: 0.1.0


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstddef>
#include <vector>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;

typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;

typedef jlcxx::ArrayRef<double, 1> JuliaArray;


JuliaArray geo_distance(int32_t start_node)
{
    std::ifstream filename(CGAL::data_file_path("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/ellipsoid_x4.off"));
    Triangle_mesh tm;
    filename >> tm;

    //property map for the distance values to the source set
    Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor,double>("v:distance",0).first;

    //pass in the idt object and its vertex_distance_map
    Heat_method_idt hm_idt(tm);

    //add the first vertex as the source set
    vertex_descriptor source = *(vertices(tm).first + start_node);

    hm_idt.add_source(source);
    hm_idt.estimate_geodesic_distances(vertex_distance);

    std::vector<double> distances_list;
    for(vertex_descriptor vd : vertices(tm)){
        distances_list.push_back(get(vertex_distance, vd));
    }

    // The ArrayRef type is provided to work conveniently with array data from Julia.
    JuliaArray distances(distances_list.data(), distances_list.size());

    return distances;
}


int main()
{
    int32_t start_node = 0;
    geo_distance(start_node);
    return 0;
}


// make this function visible to Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // register a standard C++ function
    mod.method("geo_distance", geo_distance);
}