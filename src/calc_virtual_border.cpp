// author: @janpiotraschke
// date: 2023-01-27
// license: Apache License 2.0
// version: 0.1.0


// g++ -std=c++14 -lpthread -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/calc_virtual_border.cpp -o src/calc_virtual_border
// g++ -std=c++14 -lpthread -I /usr/local/Cellar/CGAL/5.5.1/include -I /usr/local/Cellar/boost/1.80.0/include  -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3 src/calc_virtual_border.cpp -o src/calc_virtual_border

// ! Look at this again : https://theboostcpplibraries.com/boost.graph-algorithms 
// ! also look at this: https://doc.cgal.org/4.14.3/BGL/classCGAL_1_1Seam__mesh.html#a27d5997959e6a348e7206e359a378d85  

#include "calc_virtual_border.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <boost/graph/depth_first_search.hpp>
#include <vector>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <fstream>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <boost/property_map/property_map.hpp>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <iterator>
#include <string>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3                                          Point_3;
namespace My {
  struct Mesh: public CGAL::Surface_mesh<Point_3> {
    typedef CGAL::Surface_mesh<Point_3> Base;
    std::string name;
  };
} // namespace My
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My::Mesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

typedef boost::graph_traits<My::Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<My::Mesh>::edge_descriptor          edge_descriptor;

// We use a std::map to store the index
typedef std::map<vertex_descriptor, int>                              VertexIndexMap;
VertexIndexMap vertex_id_map;
// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;
VertexIdPropertyMap vertex_index_pmap(vertex_id_map);


/*
! Solve the graph using the node notation and not the Point_3 notation
*/
// int main(int argc, char** argv)
// {
//     My::Mesh mesh;
//     std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/bear.off"));
//     in >> mesh;

//     // typedef boost::graph_traits<My::Mesh>::vertex_descriptor vertex_descriptor;
//     typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;

//     Point_property_map ppm = get(CGAL::vertex_point, mesh);
//     // Create vectors to store the predecessors (p) and the distances from the root (d)
//     std::vector<vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
//     std::vector<int> distance(num_vertices(mesh));

//     auto indexmap = get(boost::vertex_index, mesh);
//     auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);
//     std::vector<vertex_descriptor> predecessor(num_vertices(mesh));  // We first declare a vector

//     boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmapa(predecessor.begin(), vertex_index_pmap);  // and then turn it into a property map

//     // to use two visitors, you need to put them in a pair (from https://theboostcpplibraries.com/boost.graph-algorithms)
//     auto vis = boost::make_bfs_visitor(std::make_pair(boost::record_distances(dist_pmap, boost::on_tree_edge{}),boost::record_predecessors(predecessor_pmapa, boost::on_tree_edge{})));

//     std::cout << "number of vertices: " << mesh.number_of_vertices() << std::endl;
//     std::cout << "number of edges in the mesh: " << mesh.number_of_edges() << std::endl;

//     vertex_descriptor start_node = *(vertices(mesh).first);


//     /*
//     Find the target node
//     */
//     boost::breadth_first_search(mesh, start_node, visitor(vis));

//     // Traverse all vertices and show at what distance they are
//     int max_distances = 0;
//     decltype(start_node) target_node;
//     for(vertex_descriptor vd : vertices(mesh)){
//         if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
//             // std::cout << vd << " at " << get(ppm, vd) << " is " << distance[vd] << " hops away" << std::endl;
//             if (distance[vd] > max_distances) {
//                 max_distances = distance[vd];
//                 target_node = vd;
//             }
//         }
//     }
//     std::cout << "max distance: " << max_distances << std::endl;
//     std::cout << "got the following start point: " << start_node << " at " << get(ppm, start_node) << std::endl;
//     std::cout << "got the following target point: " << target_node << " at " << get(ppm, target_node) << std::endl;


//     // ! From: https://www.technical-recipes.com/2015/getting-started-with-the-boost-graph-library/
//     // Evaluate Dijkstra on graph mesh with source start_node, predecessor_map predecessor_pmap and distance_map d
//     boost::dijkstra_shortest_paths(mesh, start_node, boost::predecessor_map(&predecessor_pmap[0]).distance_map(&distance[0]));

//     // get the edges of the path between the start and the target node
//     std::vector<edge_descriptor> path_list;
//     vertex_descriptor current = target_node;
//     while (current != start_node) {
//         vertex_descriptor predecessor = predecessor_pmap[current];
//         std::pair<edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
//         edge_descriptor edge = edge_pair.first;
//         path_list.push_back(edge);
//         // current = predecessor;
//         // std::cout << predecessor << " -> " << current << std::endl;
//         current = predecessor_pmap[current];
//     }

//     /*
//     The following code workes and is tested
//     -> get the vertices of the path between the start and the target node
//     */
//     // std::vector<boost::graph_traits<My::Mesh>::vertex_descriptor > path_list;
//     // boost::graph_traits<My::Mesh>::vertex_descriptor current = target_node;
//     // .selection.txt needs this weird format: https://stackoverflow.com/questions/62426624/cgal-seam-mesh-how-to-include-the-seam-in-the-halfedge-descriptor-of-the-bounda
//     // Note the duplicities (each vertex is there once as a starting point and once as an endpoint) and the two leading empty lines.
//     // "The edges to be marked as seams are described by the range [first, last) of vertices of the underlying mesh. Each edge to be marked is described by two consecutive iterators."
//     // ! TODO: remove the 'v' from the output nodes
//     // while(current!=start_node) 
//     // {
//     //     path_list.push_back(current);
//     //     current = predecessor_pmap[current];
//     //     path_list.push_back(current);
//     // }

//     // reverse the list
//     // std::reverse(path_list.begin(), path_list.end());

//     // Write the path to a file
//     // ! filename should be the name of a CGAL selection file with file extension "*.selection.txt": 
//     // ! TODO: edges are described by pairs of integers, on the THIRD line of the file.
//     std::ofstream output_files("git_repos/Confined_active_particles/meshes/bear2.selection.txt");
//     std::ostream_iterator<edge_descriptor> output_iterators(output_files, " ");
//     // write the list to the second line of the output_files
//     std::copy(path_list.begin(), path_list.end(), output_iterators);

//     return 0;
// }

int calc_virtual_border()
{
    My::Mesh mesh;
    std::ifstream in(CGAL::data_file_path("git_repos/Confined_active_particles/meshes/bear.off"));
    in >> mesh;

    // typedef boost::graph_traits<My::Mesh>::vertex_descriptor vertex_descriptor;
    typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;

    Point_property_map ppm = get(CGAL::vertex_point, mesh);
    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
    std::vector<int> distance(num_vertices(mesh));

    auto indexmap = get(boost::vertex_index, mesh);
    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);
    std::vector<vertex_descriptor> predecessor(num_vertices(mesh));  // We first declare a vector

    boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmapa(predecessor.begin(), vertex_index_pmap);  // and then turn it into a property map

    // to use two visitors, you need to put them in a pair (from https://theboostcpplibraries.com/boost.graph-algorithms)
    auto vis = boost::make_bfs_visitor(std::make_pair(boost::record_distances(dist_pmap, boost::on_tree_edge{}),boost::record_predecessors(predecessor_pmapa, boost::on_tree_edge{})));

    std::cout << "number of vertices: " << mesh.number_of_vertices() << std::endl;
    std::cout << "number of edges in the mesh: " << mesh.number_of_edges() << std::endl;

    vertex_descriptor start_node = *(vertices(mesh).first);


    /*
    Find the target node
    */
    boost::breadth_first_search(mesh, start_node, visitor(vis));

    // Traverse all vertices and show at what distance they are
    int max_distances = 0;
    decltype(start_node) target_node;
    for(vertex_descriptor vd : vertices(mesh)){
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


    // ! From: https://www.technical-recipes.com/2015/getting-started-with-the-boost-graph-library/
    // Evaluate Dijkstra on graph mesh with source start_node, predecessor_map predecessor_pmap and distance_map d
    boost::dijkstra_shortest_paths(mesh, start_node, boost::predecessor_map(&predecessor_pmap[0]).distance_map(&distance[0]));

    // get the edges of the path between the start and the target node
    std::vector<edge_descriptor> path_list;
    vertex_descriptor current = target_node;
    while (current != start_node) {
        vertex_descriptor predecessor = predecessor_pmap[current];
        std::pair<edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        // current = predecessor;
        // std::cout << predecessor << " -> " << current << std::endl;
        current = predecessor_pmap[current];
    }

    /*
    The following code workes and is tested
    -> get the vertices of the path between the start and the target node
    */
    // std::vector<boost::graph_traits<My::Mesh>::vertex_descriptor > path_list;
    // boost::graph_traits<My::Mesh>::vertex_descriptor current = target_node;
    // .selection.txt needs this weird format: https://stackoverflow.com/questions/62426624/cgal-seam-mesh-how-to-include-the-seam-in-the-halfedge-descriptor-of-the-bounda
    // Note the duplicities (each vertex is there once as a starting point and once as an endpoint) and the two leading empty lines.
    // "The edges to be marked as seams are described by the range [first, last) of vertices of the underlying mesh. Each edge to be marked is described by two consecutive iterators."
    // ! TODO: remove the 'v' from the output nodes
    // while(current!=start_node) 
    // {
    //     path_list.push_back(current);
    //     current = predecessor_pmap[current];
    //     path_list.push_back(current);
    // }

    // reverse the list
    // std::reverse(path_list.begin(), path_list.end());

    // Write the path to a file
    // ! filename should be the name of a CGAL selection file with file extension "*.selection.txt": 
    // ! TODO: edges are described by pairs of integers, on the THIRD line of the file.
    std::ofstream output_files("git_repos/Confined_active_particles/meshes/bear_now.selection.txt");
    std::ostream_iterator<edge_descriptor> output_iterators(output_files, " ");
    // write the list to the second line of the output_files
    std::copy(path_list.begin(), path_list.end(), output_iterators);

    // return path_list;
    return 0;
}
