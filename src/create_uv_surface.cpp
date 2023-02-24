// author: @Jan-Piotraschke
// date: 2023-02-13
// license: Apache License 2.0
// version: 0.1.0

// known Issue: https://github.com/CGAL/cgal/issues/2994


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

// borders
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>

// surface parameterization methods
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

// boost graph
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/property_map/property_map.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/format.hpp>

#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>
#include <format>

#include <algorithm>
#include <cstddef>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"


typedef CGAL::Simple_cartesian<double>            Kernel;
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Point_3                           Point_3;
typedef CGAL::Surface_mesh<Point_3>               SurfaceMesh;

namespace My {
  struct Mesh: public CGAL::Surface_mesh<Point_3> {
    typedef CGAL::Surface_mesh<Point_3> Base;
    std::string name;
  };
} // namespace My

#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My::Mesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

typedef boost::graph_traits<My::Mesh>::vertex_descriptor          my_vertex_descriptor;
typedef boost::graph_traits<My::Mesh>::edge_descriptor          my_edge_descriptor;
// We use a std::map to store the index
typedef std::map<my_vertex_descriptor, int>                              VertexIndexMap;
VertexIndexMap my_vertex_id_map;
// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;
VertexIdPropertyMap my_vertex_index_pmap(my_vertex_id_map);

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     SM_vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_iterator                      SM_edge_iterator;
typedef boost::graph_traits<SurfaceMesh>::halfedge_iterator                      SM_halfedge_iterator;

typedef SurfaceMesh::Property_map<SM_edge_descriptor, bool>           Seam_edge_pmap;
typedef SurfaceMesh::Property_map<SM_vertex_descriptor, bool>         Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap>  Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor                    vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor                  halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator                      vertex_iterator;

typedef SurfaceMesh::Property_map<SM_halfedge_descriptor, Point_2>      UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;


/*
Get the mesh name from its path
*/
std::string get_mesh_name(std::string mesh_3D)
{
    if (mesh_3D.find('.') < mesh_3D.length())
    {
        size_t pos = mesh_3D.find_last_of(".");
        size_t pos_slash = mesh_3D.find_last_of("/");
        std::string mesh_3D_name = mesh_3D.substr(pos_slash+1, pos-pos_slash-1);
        std::cout << "We extract the mesh name from the path string: " << mesh_3D_name << std::endl;
        return mesh_3D_name;
    }
    else
    {
        std::string mesh_3D_name = mesh_3D;
        return mesh_3D_name;
    }
}


/*
Get the mesh from the file.
Add new meshes with there path in this function
*/
std::ifstream get_mesh_obj(std::string mesh_3D)
{
    if (mesh_3D == "Ellipsoid")
    {
        std::ifstream in(CGAL::data_file_path("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/ellipsoid_x4.off"));
        return in;
    }
    else if (mesh_3D == "Sphere")
    {
        std::ifstream in(CGAL::data_file_path("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/sphere.off"));
        return in;
    }
    else if (mesh_3D == "Bear")
    {
        std::ifstream in(CGAL::data_file_path("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/bear.off"));
        return in;
    }
    else
    {
        std::ifstream in(CGAL::data_file_path(mesh_3D));
        return in;
    }
}


/*
Calculate the virtual border of the mesh

NOTE: We have this function not in a separate file because the C Language doesn't support returning a vector of our Edge data
*/
std::vector<my_edge_descriptor> calc_virtual_border(std::string mesh_3D)
{
    My::Mesh mesh;
    auto in = get_mesh_obj(mesh_3D);
    in >> mesh;

    std::cout << "number of vertices in the 3D mesh: " << mesh.number_of_vertices() << std::endl;

    typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;
    Point_property_map ppm = get(CGAL::vertex_point, mesh);

    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<my_vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
    auto indexmap = get(boost::vertex_index, mesh);
    std::vector<int> distance(num_vertices(mesh));  // record the distance from the root
    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);

    // to use two visitors, you need to put them in a pair (from https://theboostcpplibraries.com/boost.graph-algorithms)
    auto vis = boost::make_bfs_visitor(
        std::make_pair(
            boost::record_distances(dist_pmap, boost::on_tree_edge{}),
            boost::record_predecessors(&predecessor_pmap[0], boost::on_tree_edge{})
        )
    );


    /*
    Find the target node
    */
    my_vertex_descriptor start_node = *(vertices(mesh).first);
    boost::breadth_first_search(mesh, start_node, visitor(vis));

    // Traverse all vertices and show at what distance they are
    int max_distances = 0;
    decltype(start_node) target_node;
    for(my_vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
            // std::cout << vd << " at " << get(ppm, vd) << " is " << distance[vd] << " hops away" << std::endl;
            if (distance[vd] > max_distances) {
                max_distances = distance[vd];
                target_node = vd;
            }
        }
    }
    std::cout << "max distance: " << max_distances << std::endl;
    std::cout << "got the following start point: " << start_node << " at " << get(ppm, start_node) << std::endl;
    std::cout << "got the following target point: " << target_node << " at " << get(ppm, target_node) << std::endl;

    // get the edges of the path between the start and the target node
    std::vector<my_edge_descriptor> path_list;
    my_vertex_descriptor current = target_node;
    while (current != start_node) {
        my_vertex_descriptor predecessor = predecessor_pmap[current];
        // std::cout << predecessor << " -> " << current << std::endl;
        std::pair<my_edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        my_edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        current = predecessor;
    }

    // overgive the path_list to the vector b because handling vectors is easier for me
    std::vector<my_edge_descriptor> b(path_list.begin(), path_list.end());

    return b;
}


int create_uv_surface(std::string mesh_3D)
{
    // Start a timer
    CGAL::Timer task_timer;
    task_timer.start();

    // Load the 3D mesh
    SurfaceMesh sm;
    auto filename = get_mesh_obj(mesh_3D);
    filename >> sm;

    // save the STL file as an OFF file
    // std::ofstream out_OFF("git_repos/Confined_active_particles/meshes/ellipsoid_x4.off");
    // CGAL::IO::write_OFF(out_OFF, sm);
    // out_OFF.close();

    // Create property maps to store seam edges and vertices
    Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
    Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

    // Create the seam mesh
    Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);


    /*
    Calculate the virtual border
    1. mit der Funktion 'calc_virtual_border' wird die virtuelle Kante berechnet
    */
    auto calc_edges = calc_virtual_border(mesh_3D);
    for(SM_edge_descriptor e : calc_edges) {
    mesh.add_seam(source(e, sm), target(e, sm));  // Add the seams to the seam mesh
    }

    // Print the number of seam edges in the input
    std::cout << mesh.number_of_seam_edges() << " seam edges in input" << std::endl;

    // Create an index map for the seam mesh
    typedef std::unordered_map<vertex_descriptor, int> Indices;
    Indices indices;
    boost::associative_property_map<Indices> vimap(indices);
    int counter = 0;
    for(vertex_descriptor vd : vertices(mesh)) {
        put(vimap, vd, counter++);
        // std::cout << "vertex " << vd << " has index " << get(vimap, vd) << std::endl;
    }

    // The 2D points of the uv parametrisation will be written into this map
    // NOTE that this is a halfedge property map, and that uv values are only stored for the canonical Halfedges Representing a Vertex
    UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

    // Choose the border of the uv parametrisation
    typedef SMP::Circular_border_arc_length_parameterizer_3<Mesh> Border_parameterizer;
    // typedef SMP::Square_border_uniform_parameterizer_3<Mesh> Border_parameterizer;
    Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices

    // Iterative Authalic Parameterization:
    // from https://doi.org/10.1109/ICCVW.2019.00508
    // This parameterization is a fixed border parameterization and is part of the authalic parameterization family,
    // meaning that it aims to minimize area distortion between the input surface mesh and the parameterized output.
    typedef SMP::Iterative_authalic_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;
    Parameterizer parameterizer(border_parameterizer);

    // Other parameterization algorithms:
    // typedef SMP::Discrete_authalic_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;
    // typedef SMP::Mean_value_coordinates_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;

    // Choose a halfedge on the (possibly virtual) border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    // Perform the parameterization
    const unsigned int ITERATIONS = 9;

    /*
    computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
    The mapping is piecewise linear on the triangle mesh. The result is a pair (u,v) of parameter coordinates for each vertex of the input mesh.
    ! A one-to-one mapping may be guaranteed or not, depending on the chosen Parameterizer algorithm
    */
    SMP::Error_code err = parameterizer.parameterize(mesh, bhd, uvmap, ITERATIONS);
    // SMP::Error_code err = SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);

    if(err != SMP::OK){
        std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
    }

    // print the number of vertices of the uvmap
    std::cout << "Number of vertices in uvmap: " << indices.size() << std::endl;
    // get the mapping of vertices between the 3D mesh and the uvmap
    // for(vertex_descriptor vd : vertices(mesh)) {
    //     std::cout << "Input point: " << vd << " is mapped to " << get(uvmap, vd) << '\n';
    // }

    auto mesh_3D_name = get_mesh_name(mesh_3D);

    auto path_uv = str(boost::format("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/%s_uv.off") % mesh_3D_name);
    std::cout << "The UV mesh is saved to the following path: " << path_uv << "\n" << std::endl;
    std::ofstream out(path_uv);
    SMP::IO::output_uvmap_to_off(mesh, bhd, uvmap, out);

    std::cout << "Finished in " << task_timer.time() << " seconds" << std::endl;
    return 0;
}


int main()
{
    std::string mesh_3D = "Ellipsoid";
    create_uv_surface(mesh_3D);

    return 0;
}


// make this function visible to Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // register a standard C++ function
    mod.method("create_uv_surface", create_uv_surface);
}