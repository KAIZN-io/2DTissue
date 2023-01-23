// g++ -std=c++14 -lpthread -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/calc_path.cpp -o src/calc_path

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
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
namespace My {
  struct Mesh: public CGAL::Surface_mesh<Point_3> {
    typedef CGAL::Surface_mesh<Point_3> Base;
    std::string name;
  };
} // namespace My
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My::Mesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>


/*
! Solve the graph using the node notation and not the Point_3 notation
*/
int main(int argc, char** argv)
{
    My::Mesh mesh;
    std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
    in >> mesh;

    typedef boost::graph_traits<My::Mesh>::vertex_descriptor vertex_descriptor;
    typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;

    Point_property_map ppm = get(CGAL::vertex_point, mesh);
    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<vertex_descriptor> predecessor_pmap(num_vertices(mesh));
    std::vector<int> d(num_vertices(mesh));

    // Create descriptor for the source node
    // vertex_descriptor s = boost::vertex(vertices(mesh).first, mesh);
    //   vertex_descriptor goal = vertex(E, mesh);
    for(vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
            std::cout << vd << " at " << get(ppm, vd) << std::endl;
        }
    }

    auto indexmap = get(boost::vertex_index, mesh);
    std::vector<int> distance(mesh.number_of_vertices());
    std::cout << "P.size_of_vertices() = " << mesh.number_of_vertices() << std::endl;
    vertex_descriptor start_node = *(vertices(mesh).first);


    // the following code: DON'T ask me ... I will fix it in the future
    decltype(start_node) target_node;
    for(vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
            target_node = vd;
        }
    }

    std::cout << "starting node = " << target_node << std::endl;
    std::cout << "We compute the route from the following point: " << get(ppm, start_node) << " to " << get(ppm,target_node) <<  std::endl;

    // ! From: https://www.technical-recipes.com/2015/getting-started-with-the-boost-graph-library/
    // Evaluate Dijkstra on graph mesh with source start_node, predecessor_map predecessor_pmap and distance_map d
    boost::dijkstra_shortest_paths(mesh, start_node, boost::predecessor_map(&predecessor_pmap[0]).distance_map(&d[0]));



    // Polyhedron P;
    // std::ifstream input((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
    // input >> P ;

    // auto dist_pmap = boost::make_iterator_property_map(d.begin(), indexmap);

    // // record the predecessor of each vertex
    // std::vector<vertex_descriptor> predecessor(boost::num_vertices(P));  // We first declare a vector
    // boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmaps(predecessor.begin(), vertex_index_pmap);  // and then turn it into a property map
  
    // auto vis = boost::make_bfs_visitor(std::make_pair(boost::record_distances(&d[0], boost::on_tree_edge{}),boost::record_predecessors(&predecessor_pmaps[0], boost::on_tree_edge{})));
    // // run the search
    // boost::breadth_first_search(P, start_node, visitor(vis));




    //p[] is the predecessor map obtained through dijkstra
    //name[] is a vector with the names of the vertices
    //s and goal are vertex descriptors
    std::vector<boost::graph_traits<My::Mesh>::vertex_descriptor > path;
    boost::graph_traits<My::Mesh>::vertex_descriptor current = target_node;
 
    while(current!=start_node) 
    {
        path.push_back(current);
        current = predecessor_pmap[current];
    }
    path.push_back(start_node);

    // Prints the path obtained in reverse
    std::vector<boost::graph_traits<My::Mesh>::vertex_descriptor >::reverse_iterator it;
 
    int hope_number = 0;
    for (it = path.rbegin(); it != path.rend(); ++it) {
        std::cout << get(ppm, *it) << " \n";
        hope_number++;
    }
    std::cout << "\n " << hope_number << std::endl;
    std::cout << "number of edges in the mesh: " << mesh.number_of_edges() << '\n';

    return 0;
}
