// g++ -std=c++14 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/open-mesh/9.0/include -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/calc_distance.cpp -o src/calc_distance

// // ! anschauen: https://doc.cgal.org/latest/BGL/index.html#title30

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <boost/graph/depth_first_search.hpp>
#include <vector>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <fstream>
#include <iostream>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
// #include <CGAL/boost/graph/dijkstra_shortest_paths.h>


#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

// #include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
typedef std::map<vertex_descriptor, int>                              VertexIndexMap;
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;

int main(int argc, char** argv) {
    Polyhedron P;
    std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
    in >> P ;

    // associate indices to the vertices using the "id()" field of the vertex.
    vertex_iterator vertex_begin, vertex_end;
    int index = 0;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vit and ve
    for(boost::tie(vertex_begin,vertex_end)=boost::vertices(P); vertex_begin!=vertex_end; ++vertex_begin ){
        vertex_descriptor  vd = *vertex_begin;
        vd->id() = index++;
    }

    // This is the vector where the distance gets written to
    std::vector<int> distance(P.size_of_vertices());
    std::cout << "P.size_of_vertices() = " << P.size_of_vertices() << std::endl;

    // Here we start at an arbitrary vertex
    // Any other vertex could be the starting point
    boost::tie(vertex_begin,vertex_end)=boost::vertices(P);
    vertex_descriptor start = *vertex_begin;      // vertex_descriptor start = *(vertices(P).first);
    std::cout << "We compute the distances to the following point: " << start->point() << std::endl;


    // bfs = breadth first search explores the graph
    //      -> based on Moore, E. F. (1959). "The shortest path through a maze". Proceedings of an International Symposium on the Theory of Switching (Cambridge, Massachusetts, 2–5 April 1957). Cambridge: Harvard University Press. pp. 285–292.
    // Just as the distance_recorder there is a way to record the predecessor of a vertex
    auto indexmap = get(boost::vertex_index, P);

    std::vector<vertex_descriptor> path;
    auto record_path = [&path](const vertex_descriptor& v) {
        path.push_back(v);
    };

    auto dist_pmap = make_iterator_property_map(distance.begin(), indexmap);
    // auto vis = visitor(boost::make_bfs_visitor(record_path)).vertex_index_map(indexmap);
    auto vis = boost::make_bfs_visitor(boost::record_distances(dist_pmap, boost::on_tree_edge()));
    boost::breadth_first_search(P, start, visitor(vis));


    // Traverse all vertices and show at what distance they are
    int max_distances = 0;
    auto end_point = Kernel::Point_3(0,0,0);
    for(boost::tie(vertex_begin,vertex_end)=boost::vertices(P); vertex_begin!=vertex_end; ++vertex_begin ){
        start = *vertex_begin;
        auto v_distance = distance[start->id()];
        std::cout <<  start->point() << "  is " << v_distance << " hops away" << std::endl;
        if (v_distance > max_distances) {
            max_distances = v_distance;
            end_point = start->point();
        }
    }

    std::cout << "max_distances = " << max_distances << std::endl;
    Vertex_handle source = P.vertices_begin();
    std::cout << "source = " << source->point() << std::endl;
    std::cout << "target = " << end_point << std::endl;
    // TODO: gibt die Roote zum am weitesten entfernten Vertex aus   -> https://stackoverflow.com/questions/47518846/how-to-find-the-shortest-path-between-two-vertices-in-a-bgl-graph
    // TODO: nehme dann ein Nachbar Vertex und gehe zu diesem Vertex
    // TODO: verbinde die Vertices zu einer cut-line

    std::cout << "\n Start dijkstra_shortest_paths at " << source->point() <<"\n";

     // Create a vector to store the distances from the source vertex
    std::vector<double> distances(boost::num_vertices(P));
    std::vector<vertex_descriptor> predecessors(boost::num_vertices(P));     // std::vector<int> predecessors(boost::num_vertices(g));

    
    // boost::iterator_property_map<std::vector<vertex_descriptor>::iterator,boost::property_map<Polyhedron, boost::vertex_index_t>::const_type> predecessor_map(predecessors.begin(), get(boost::vertex_index, P));
    
    // boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), indexmap);
    // distance_pmap = boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, P));
    // Use Dijkstra's algorithm to find the shortest path
    // boost::dijkstra_shortest_paths(P, source, boost::distance_map(boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, P))).predecessor_map(predecessor_map));
    // boost::dijkstra_shortest_paths(P, source, distance_map(distance_pmap).predecessor_map(predecessor_map).vertex_index_map(vertex_index_pmap));
    // // Print the distance from the source to the target vertex
    // std::cout << "Shortest path distance: " << distances[target] << std::endl;

    // // Get the vertices of the shortest path
    // int vertex = target;
    // std::cout << "Shortest path: " << vertex;
    // while (predecessor_map[vertex] != vertex) {
    //     vertex = predecessor_map[vertex];
    //     std::cout << " <- " << vertex;
    // }
    // std::cout << std::endl;

    return 0;
}





// int main() {


//     // Print the distance from the source to the target vertex
//     std::cout << "Shortest path distance: " << distances[target] << std::endl;

//     // Get the vertices of the shortest path
//     vertex_descriptor vertex = target;
//     std::cout << "Shortest path: " << vertex;
//     while (predecessor_map[vertex] != vertex) {
//         vertex = predecessor_map[vertex];
//         std::cout << " <- " << vertex;
//     }
//     std::cout << std::endl;

//     return 0;
// }


// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
// #include <boost/graph/dijkstra_shortest_paths.hpp>
// #include <boost/property_map/property_map.hpp>

// typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
// typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
// typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;

// int main() {
//     // Create the polyhedron
//     Polyhedron polyhedron;

//     // Load the sphere mesh into the polyhedron
//     std::ifstream input("sphere.off");
//     input >> polyhedron;

//     // Define the source and target vertices
//     vertex_descriptor source = vertex(0, polyhedron);
//     vertex_descriptor target = vertex(7, polyhedron);

//     // Create a vector to store the distances from the source vertex
//     std::vector<double> distances(num_vertices(polyhedron));
//     std::vector<vertex_descriptor> predecessors(num_vertices(polyhedron));
//     boost::iterator_property_map<std::vector<vertex_descriptor>::iterator,
//                                  boost::property_map<Polyhedron, boost::vertex_index_t>::const_type>
//         predecessor_map(predecessors.begin(), get(boost::vertex_index, polyhedron));

//     // Use Dijkstra's algorithm to find the shortest path
//     boost::dijkstra_shortest_paths(polyhedron, source, boost::distance_map(boost::make_iterator_property_map(distances.begin(), get(boost::vertex_index, polyhedron))).predecessor_map(predecessor_map));

//     // Print the distance from the source to the target vertex
//     std::cout << "Shortest path distance: " << distances[target] << std::endl;

//     // Get the vertices of the shortest path
//     vertex_descriptor vertex = target;
//     std::cout << "Shortest path: " << vertex;
//     while (predecessor_map[vertex] != vertex) {
//         vertex = predecessor_map[vertex];
//         std::cout << " <- " << vertex;
//     }
//     std::cout << std::endl;

//     return 0;
// }

